!+ Data module for all data arrays, that are used by the radiation scheme
!-------------------------------------------------------------------------------

MODULE radiation_data

!-------------------------------------------------------------------------------
!
! Description:
!  This module declares and initializes all parametric data and data arrays 
!  that are used by the Ritter-Geleyn radiation scheme (i.e. the FESFT-routine). 
!  It has been derived from former data_radiation (which was used for the 
!  non-blocked scheme.
!
! Current Code Owner: DWD, ???
!  phone:  +49  69  8062 
!  fax:    +49  69  8062 3721
!  email:  @dwd.de
!
! History:
! Version      Date       Name
! ----------   ---------- ----
! V5_3         2015-10-09 Ulrich Schaettler
!  Initial release (based on data_radiation of V5_2)
!
! Code Description:
! Language: Fortran 90.
! Software Standards: "European Standards for Writing and
! Documenting Exchangeable Fortran 90 Code".
!=======================================================================
!
! Declarations:
!
! Modules used:
USE kind_parameters , ONLY :   &
           wp,       &! KIND-type parameters for real variables
           dp         ! KIND-type parameters for double precision variables

!-------------------------------------------------------------------------------

IMPLICIT NONE

! Local Declarations:

  INTEGER, PRIVATE ::    &
    i ! index for data initalization

! Global (i.e. public) Declarations:

! Global Parameters
! -----------------

  INTEGER, PARAMETER  ::    &

    ! Parameters for spectral resolution of radiative transfer code
    jpsol  = 3 ,     & ! Number of solar spectral intervals
    jpther = 5 ,     & ! Number of thermal spectral intervals
    jpspec = 8 ,     & ! (= jpsol+jpther) Total number of spectral intervals
 
    ! Parameters for gas absorbtion                                
    jpgas  = 3 ,     & ! Number of gases considered by radative transfer code
    jpabsc = 7         ! Maximum number of absorbtion coefficients in each   
                       ! spectral interval

! Global dimensions
! -----------------

  INTEGER                              ::    &
    idim_rad,        & ! ie-dimension of the coarser grid
    istartrad,       & ! start- and end-indices for computing the radiation
    iendrad,         & !   (when running on a coarser grid, the input values for
    jstartrad,       & !    fesft are computed on all grid points, to compute an
    jendrad,         & !    average input over several grid points)
    iendparrad,      & ! end-index just for preparations
    jendparrad         ! end-index just for preparations

  INTEGER                             ::  &
    istartradheat,   & !
    iendradheat,     & !
    jstartradheat,   & !
    jendradheat        !

! Global Arrays and Scalars
! -------------------------

! 1. absorption properties of atmospheric gases (..,1=h2o; ..,2 =CO2; ..,3=O3)
! --------------------------------------------

  REAL  (KIND=dp)         ::    &
    coai (jpabsc,jpspec,jpgas), & ! weigthing coefficients
    cobi (jpabsc,jpspec,jpgas), & ! absorption coefficients
    coali(jpabsc,jpspec,jpgas), & ! pressure correction coefficients
    cobti(jpabsc,jpspec,jpgas), & ! temperature correction coefficients
    pgas (       jpspec,jpgas), & ! reference pressure
    tgas (       jpspec,jpgas)    ! reference temperature

  INTEGER                 ::    &
    ncgas(jpspec,3),            & ! number of coefficients for each spectral  
                                  ! interval and gas (maximum=7)
    nfast(jpspec)                 ! control variable for choice between 
                                  ! ESFT/FESFT in each spectral interval

  ! Initialization of above arrays

  DATA nfast / 1, 1, 1, 1, 1, 1, 1, 1/
                                                                        
  ! coefficients for H2O in spectral interval 1
  DATA ncgas(1,1) /7/ ; DATA pgas(1,1) /101325.000_dp/ ; DATA tgas(1,1) /281.700_dp/
  DATA (coai (i,1,1),i=1,7)/ .114190E-01_dp,  .600200E-01_dp,  .111201E+00_dp,  .123340E+00_dp,  .902500E-01_dp,  .199632E+00_dp,  .404139E+00_dp/
  DATA (cobi (i,1,1),i=1,7)/ .209894E+02_dp,  .208930E+01_dp,  .184502E+00_dp,  .217771E-01_dp,  .279254E-02_dp,  .463447E-03_dp,  .000000E+00_dp/
  DATA (coali(i,1,1),i=1,7)/ .285370E-01_dp,  .688620E+00_dp,  .766031E+00_dp,  .833136E+00_dp,  .773491E+00_dp,  .768818E+00_dp,  .100000E+01_dp/
  DATA (cobti(i,1,1),i=1,7)/ .473006E+00_dp, -.468639E+00_dp, -.599601E+00_dp, -.162223E+01_dp, -.176002E+01_dp, -.153131E+01_dp,  .100000E+01_dp/

  ! coefficients for H2O in spectral interval 2
  DATA ncgas(2,1) /7/ ; DATA pgas(2,1) /101325.000_dp/ ; DATA tgas(2,1) /281.700_dp/
  DATA (coai (i,2,1),i=1,7)/ .201500E-02_dp,  .268530E-01_dp,  .598920E-01_dp,  .907740E-01_dp,  .102284E+00_dp,  .217298E+00_dp,  .500884E+00_dp/
  DATA (cobi (i,2,1),i=1,7)/ .508159E+01_dp,  .519996E+00_dp,  .465586E-01_dp,  .891251E-02_dp,  .159221E-02_dp,  .374973E-03_dp,  .000000E+00_dp/
  DATA (coali(i,2,1),i=1,7)/-.482300E-02_dp,  .529161E+00_dp,  .587751E+00_dp,  .756567E+00_dp,  .774607E+00_dp,  .733883E+00_dp,  .100000E+01_dp/
  DATA (cobti(i,2,1),i=1,7)/ .499755E+00_dp, -.529716E+00_dp, -.177970E-01_dp, -.746447E+00_dp, -.106191E+00_dp, -.727589E+00_dp,  .100000E+01_dp/
                                                                        
  ! coefficients for H2O in spectral interval 3
  DATA ncgas(3,1) /3/ ; DATA pgas(3,1) /101325.000_dp/ ; DATA tgas(3,1) /281.700_dp/
  DATA (coai (i,3,1),i=1,7)/ .566900E-02_dp,  .346720E-01_dp,  .959659E+00_dp,  .000000E+00_dp,  .000000E+00_dp,  .000000E+00_dp,  .000000E+00_dp/
  DATA (cobi (i,3,1),i=1,7)/ .716144E-03_dp,  .256449E-03_dp,  .000000E+00_dp,  .000000E+00_dp,  .000000E+00_dp,  .000000E+00_dp,  .000000E+00_dp/
  DATA (coali(i,3,1),i=1,7)/-.281669E+00_dp,  .611673E+00_dp,  .100000E+01_dp,  .000000E+00_dp,  .000000E+00_dp,  .000000E+00_dp,  .000000E+00_dp/
  DATA (cobti(i,3,1),i=1,7)/ .418657E+00_dp,  .405230E-01_dp,  .100000E+01_dp,  .000000E+00_dp,  .000000E+00_dp,  .000000E+00_dp,  .000000E+00_dp/
                                                                        
  ! coefficients for H2O in spectral interval 4
  DATA ncgas(4,1) /7/ ; DATA pgas(4,1) / 50662.500_dp/ ; DATA tgas(4,1) /255.800_dp/
  DATA (coai (i,4,1),i=1,7)/ .641200E-02_dp,  .362630E-01_dp,  .147064E+00_dp,  .285387E+00_dp,  .246376E+00_dp,  .226899E+00_dp,  .515980E-01_dp/
  DATA (cobi (i,4,1),i=1,7)/ .298538E+04_dp,  .139959E+03_dp,  .152405E+02_dp,  .144212E+01_dp,  .183654E+00_dp,  .283139E-01_dp,  .409261E-02_dp/
  DATA (coali(i,4,1),i=1,7)/ .183780E-01_dp,  .410557E+00_dp,  .808897E+00_dp,  .897332E+00_dp,  .932149E+00_dp,  .978389E+00_dp,  .100000E+01_dp/
  DATA (cobti(i,4,1),i=1,7)/ .413777E+00_dp, -.663704E+00_dp, -.953789E+00_dp, -.111883E+01_dp, -.156269E+01_dp, -.330557E+01_dp,  .100000E+01_dp/
                                                                        
  ! coefficients for H2O in spectral interval 5
  DATA ncgas(5,1) /7/ ; DATA pgas(5,1) / 86126.250_dp/ ; DATA tgas(5,1) /281.700_dp/
  DATA (coai (i,5,1),i=1,7)/ .147700E-02_dp,  .345020E-01_dp,  .865590E-01_dp,  .144237E+00_dp,  .218089E+00_dp,  .339440E+00_dp,  .175697E+00_dp/
  DATA (cobi (i,5,1),i=1,7)/ .126765E+02_dp,  .149624E+01_dp,  .147571E+00_dp,  .368129E-01_dp,  .792501E-02_dp,  .208930E-02_dp,  .000000E+00_dp/
  DATA (coali(i,5,1),i=1,7)/-.414300E-02_dp,  .504464E+00_dp,  .670985E+00_dp,  .920940E+00_dp,  .889089E+00_dp,  .966028E+00_dp,  .100000E+01_dp/
  DATA (cobti(i,5,1),i=1,7)/ .454691E+00_dp, -.423980E+01_dp, -.340869E+01_dp, -.410896E+01_dp, -.268068E+01_dp, -.250967E+01_dp,  .100000E+01_dp/

  ! coefficients for H2O in spectral interval 6
  DATA ncgas(6,1) /4/ ; DATA pgas(6,1) / 86126.250_dp/ ; DATA tgas(6,1) /281.700_dp/
  DATA (coai (i,6,1),i=1,7)/ .653200E-02_dp,  .700040E-01_dp,  .243768E+00_dp,  .679696E+00_dp,  .000000E+00_dp,  .000000E+00_dp,  .000000E+00_dp/
  DATA (cobi (i,6,1),i=1,7)/ .632412E+00_dp,  .473151E-02_dp,  .163305E-02_dp,  .000000E+00_dp,  .000000E+00_dp,  .000000E+00_dp,  .000000E+00_dp/
  DATA (coali(i,6,1),i=1,7)/ .794801E+00_dp,  .306898E+00_dp,  .100000E+01_dp,  .100000E+01_dp,  .000000E+00_dp,  .000000E+00_dp,  .000000E+00_dp/
  DATA (cobti(i,6,1),i=1,7)/-.100000E+02_dp, -.219711E+01_dp, -.369325E+01_dp,  .100000E+01_dp,  .000000E+00_dp,  .000000E+00_dp,  .000000E+00_dp/
                                                                        
  ! coefficients for H2O in spectral interval 7
  DATA ncgas(7,1) /3/ ; DATA pgas(7,1) / 86126.250_dp/ ; DATA tgas(7,1) /281.700_dp/
  DATA (coai (i,7,1),i=1,7)/ .138610E-01_dp,  .226595E+00_dp,  .759544E+00_dp, 0.000000E+00_dp,  .000000E+00_dp, 0.000000E+00_dp, 0.000000E+00_dp/
  DATA (cobi (i,7,1),i=1,7)/ .425598E-02_dp,  .155239E-02_dp, 0.000000E+00_dp, 0.000000E+00_dp,  .000000E+00_dp, 0.000000E+00_dp, 0.000000E+00_dp/
  DATA (coali(i,7,1),i=1,7)/-.736171E+00_dp,  .805828E+00_dp,  .100000E+01_dp, 0.000000E+00_dp,  .000000E+00_dp, 0.000000E+00_dp, 0.000000E+00_dp/
  DATA (cobti(i,7,1),i=1,7)/ .308301E+00_dp, -.267573E+01_dp,  .100000E+01_dp, 0.000000E+00_dp,  .000000E+00_dp, 0.000000E+00_dp, 0.000000E+00_dp/
  ! coefficients for H2O in spectral interval 8
  DATA ncgas(8,1) /7/ ; DATA pgas(8,1) / 75993.750_dp/ ; DATA tgas(8,1) /281.700_dp/
  DATA (coai (i,8,1),i=1,7)/ .181840E-01_dp,  .106586E+00_dp,  .237611E+00_dp,  .241085E+00_dp,  .157304E+00_dp,  .178767E+00_dp,  .604640E-01_dp/
  DATA (cobi (i,8,1),i=1,7)/ .822243E+02_dp,  .979490E+01_dp,  .905733E+00_dp,  .140281E+00_dp,  .193197E-01_dp,  .320627E-02_dp,  .000000E+00_dp/
  DATA (coali(i,8,1),i=1,7)/-.126888E+00_dp,  .701873E+00_dp,  .834941E+00_dp,  .920550E+00_dp,  .849506E+00_dp,  .931957E+00_dp,  .100000E+01_dp/
  DATA (cobti(i,8,1),i=1,7)/ .384580E+00_dp, -.187972E+01_dp, -.226834E+01_dp, -.475940E+01_dp, -.589531E+01_dp, -.395962E+01_dp,  .100000E+01_dp/
                    
  ! coefficients for CO2 in spectral interval 1
  DATA ncgas(1,2) /6/ ; DATA pgas(1,2) / 86126.250_dp/ ; DATA tgas(1,2) /255.800_dp/
  DATA (coai (i,1,2),i=1,7)/ .592000E-02_dp,  .667700E-02_dp,  .423020E-01_dp,  .732310E-01_dp,  .140143E+00_dp,  .731727E+00_dp,  .000000E+00_dp/
  DATA (cobi (i,1,2),i=1,7)/ .760326E+02_dp,  .480839E+01_dp,  .391742E+00_dp,  .133968E-01_dp,  .355631E-02_dp,  .000000E+00_dp,  .000000E+00_dp/
  DATA (coali(i,1,2),i=1,7)/ .659071E+00_dp,  .240858E+00_dp,  .694157E+00_dp,  .424843E+00_dp,  .694262E+00_dp,  .100000E+01_dp,  .000000E+00_dp/
  DATA (cobti(i,1,2),i=1,7)/ .467048E+00_dp,  .395422E+00_dp, -.902210E+00_dp, -.557526E+00_dp, -.473196E+00_dp,  .100000E+01_dp,  .000000E+00_dp/
                                                                        
  ! coefficients for CO2 in spectral interval 2
  DATA ncgas(2,2) /3/ ; DATA pgas(2,2)/ 86126.250_dp/ ;  DATA tgas(2,2) /255.800_dp/
  DATA (coai (i,2,2),i=1,7)/ .278000E-02_dp,  .197330E-01_dp,  .977487E+00_dp,  .000000E+00_dp,  .000000E+00_dp,  .000000E+00_dp,  .000000E+00_dp/
  DATA (cobi (i,2,2),i=1,7)/ .169434E+00_dp,  .103753E-01_dp,  .000000E+00_dp,  .000000E+00_dp,  .000000E+00_dp,  .000000E+00_dp,  .000000E+00_dp/
  DATA (coali(i,2,2),i=1,7)/ .138563E+00_dp,  .831359E+00_dp,  .100000E+01_dp,  .000000E+00_dp,  .000000E+00_dp,  .000000E+00_dp,  .000000E+00_dp/
  DATA (cobti(i,2,2),i=1,7)/ .475293E+00_dp, -.496213E+00_dp,  .100000E+01_dp,  .000000E+00_dp,  .000000E+00_dp,  .000000E+00_dp,  .000000E+00_dp/
                                                                        
  ! coefficients for CO2 in spectral interval 3
  DATA ncgas(3,2) /2/ ; DATA pgas(3,2) / 86126.250_dp/ ; DATA tgas(3,2) /255.800_dp/
  DATA (coai (i,3,2),i=1,7)/ .306100E-02_dp,  .996939E+00_dp,  .000000E+00_dp,  .000000E+00_dp,  .000000E+00_dp,  .000000E+00_dp,  .000000E+00_dp/
  DATA (cobi (i,3,2),i=1,7)/ .101625E+00_dp,  .000000E+00_dp,  .000000E+00_dp,  .000000E+00_dp,  .000000E+00_dp,  .000000E+00_dp,  .000000E+00_dp/
  DATA (coali(i,3,2),i=1,7)/ .100000E+01_dp,  .100000E+01_dp,  .000000E+00_dp,  .000000E+00_dp,  .000000E+00_dp,  .000000E+00_dp,  .000000E+00_dp/
  DATA (cobti(i,3,2),i=1,7)/-.100670E+01_dp,  .100000E+01_dp,  .000000E+00_dp,  .000000E+00_dp,  .000000E+00_dp,  .000000E+00_dp,  .000000E+00_dp/
                                                                        
  ! coefficients for CO2 in spectral interval 4
  DATA ncgas(4,2) /0/ ; DATA pgas(4,2) / 60795.000_dp/ ; DATA tgas(4,2) /255.800_dp/
  DATA (coai (i,4,2),i=1,7)/ .335800E-02_dp,  .996642E+00_dp,  .000000E+00_dp,  .000000E+00_dp,  .000000E+00_dp,  .000000E+00_dp,  .000000E+00_dp/ 
  DATA (cobi (i,4,2),i=1,7)/ .247172E+00_dp,  .000000E+00_dp,  .000000E+00_dp,  .000000E+00_dp,  .000000E+00_dp,  .000000E+00_dp,  .000000E+00_dp/
  DATA (coali(i,4,2),i=1,7)/ .100000E+01_dp,  .100000E+01_dp,  .000000E+00_dp,  .000000E+00_dp,  .000000E+00_dp,  .000000E+00_dp,  .000000E+00_dp/
  DATA (cobti(i,4,2),i=1,7)/-.807310E+00_dp,  .100000E+01_dp,  .000000E+00_dp,  .000000E+00_dp,  .000000E+00_dp,  .000000E+00_dp,  .000000E+00_dp/
                                                                        
  ! coefficients for CO2 in spectral interval 5
  DATA ncgas(5,2) /7/ ; DATA pgas(5,2) / 10132.500_dp/ ; DATA tgas(5,2) /229.900_dp/
  DATA (coai (i,5,2),i=1,7)/ .452500E-02_dp,  .321420E-01_dp,  .659180E-01_dp,  .101074E+00_dp,  .107224E+00_dp,  .186663E+00_dp,  .502454E+00_dp/
  DATA (cobi (i,5,2),i=1,7)/ .299226E+03_dp,  .364754E+02_dp,  .271644E+01_dp,  .570164E+00_dp,  .100231E+00_dp,  .224388E-01_dp,  .000000E+00_dp/
  DATA (coali(i,5,2),i=1,7)/ .466819E+00_dp,  .319510E+00_dp,  .596734E+00_dp,  .751216E+00_dp,  .708519E+00_dp,  .744381E+00_dp,  .100000E+01_dp/
  DATA (cobti(i,5,2),i=1,7)/ .358348E+00_dp, -.739332E+00_dp, -.183599E+01_dp, -.289470E+01_dp, -.214575E+01_dp, -.585028E+01_dp,  .100000E+01_dp/
                                                                        
  ! coefficients for CO2 in spectral interval 6
  DATA ncgas(6,2) /3/ ; DATA pgas(6,2) / 50662.500_dp/ ; DATA tgas(6,2) /255.800_dp/
  DATA (coai (i,6,2),i=1,7)/ .119551E+00_dp,  .899140E-01_dp,  .790535E+00_dp,  .000000E+00_dp,  .000000E+00_dp,  .000000E+00_dp,  .000000E+00_dp/
  DATA (cobi (i,6,2),i=1,7)/ .305492E-02_dp,  .148936E-02_dp,  .000000E+00_dp,  .000000E+00_dp,  .000000E+00_dp,  .000000E+00_dp,  .000000E+00_dp/
  DATA (coali(i,6,2),i=1,7)/ .783365E+00_dp, -.113116E+00_dp,  .100000E+01_dp,  .000000E+00_dp,  .000000E+00_dp,  .000000E+00_dp,  .000000E+00_dp/
  DATA (cobti(i,6,2),i=1,7)/-.447333E+01_dp,  .296352E+01_dp,  .100000E+01_dp,  .000000E+00_dp,  .000000E+00_dp,  .000000E+00_dp,  .000000E+00_dp/
                                                                        
  ! coefficients for CO2 in spectral interval 7
  DATA ncgas(7,2) /3/ ; DATA pgas(7,2) / 50662.500_dp/ ; DATA tgas(7,2) /255.800_dp/
  DATA (coai (i,7,2),i=1,7)/ .577890E-01_dp,  .321750E-01_dp,  .910036E+00_dp,  .000000E+00_dp,  .000000E+00_dp,  .000000E+00_dp,  .000000E+00_dp/
  DATA (cobi (i,7,2),i=1,7)/ .650130E-02_dp,  .309030E-02_dp,  .000000E+00_dp,  .000000E+00_dp,  .000000E+00_dp,  .000000E+00_dp,  .000000E+00_dp/
  DATA (coali(i,7,2),i=1,7)/ .295465E+00_dp,  .930860E-01_dp,  .100000E+01_dp,  .000000E+00_dp,  .000000E+00_dp,  .000000E+00_dp,  .000000E+00_dp/
  DATA (cobti(i,7,2),i=1,7)/-.562957E+01_dp, -.984577E+01_dp,  .100000E+01_dp,  .000000E+00_dp,  .000000E+00_dp,  .000000E+00_dp,  .000000E+00_dp/
                                                                        
  ! coefficients for CO2 in spectral interval 8
  DATA ncgas(8,2) /4/ ; DATA pgas(8,2) / 50662.500_dp/ ; DATA tgas(8,2) /255.800_dp/
  DATA (coai (i,8,2),i=1,7)/ .317000E-02_dp,  .127109E+00_dp,  .114118E+00_dp,  .755604E+00_dp,  .000000E+00_dp,  .000000E+00_dp,  .000000E+00_dp/
  DATA (cobi (i,8,2),i=1,7)/ .174181E+02_dp,  .495450E-01_dp,  .165196E-01_dp,  .000000E+00_dp,  .000000E+00_dp,  .000000E+00_dp,  .000000E+00_dp/
  DATA (coali(i,8,2),i=1,7)/ .511300E-02_dp,  .252848E+00_dp,  .851104E+00_dp,  .100000E+01_dp,  .000000E+00_dp,  .000000E+00_dp,  .000000E+00_dp/
  DATA (cobti(i,8,2),i=1,7)/ .495222E+00_dp,  .445084E+00_dp,  .117957E+01_dp,  .100000E+01_dp,  .000000E+00_dp,  .000000E+00_dp,  .000000E+00_dp/
                                                                        
  ! coefficients for O3  in spectral interval 1
  DATA ncgas(1,3) /0/ ; DATA pgas(1,3) /  3039.75_dp/ ; DATA tgas(1,3) /229.900_dp/
  DATA (coai (i,1,3),i=1,7)/ .306000E-03_dp,  .999694E+00_dp,  .000000E+00_dp,  .000000E+00_dp,  .000000E+00_dp,  .000000E+00_dp,  .000000E+00_dp/
  DATA (cobi (i,1,3),i=1,7)/ .409261E+02_dp,  .000000E+00_dp,  .000000E+00_dp,  .000000E+00_dp,  .000000E+00_dp,  .000000E+00_dp,  .000000E+00_dp/
  DATA (coali(i,1,3),i=1,7)/ .618332E+00_dp,  .100000E+01_dp,  .000000E+00_dp,  .000000E+00_dp,  .000000E+00_dp,  .000000E+00_dp,  .000000E+00_dp/
  DATA (cobti(i,1,3),i=1,7)/-.974847E+00_dp,  .100000E+01_dp,  .000000E+00_dp,  .000000E+00_dp,  .000000E+00_dp,  .000000E+00_dp,  .000000E+00_dp/
                                                                        
  ! coefficients for O3  in spectral interval 2
  DATA ncgas(2,3) /0/ ; DATA pgas(2,3) /  3039.75_dp/ ; DATA tgas(2,3) /229.900_dp/
  DATA (coai (i,2,3),i=1,7)/ .154800E-02_dp,  .998452E+00_dp,  .000000E+00_dp,  .000000E+00_dp,  .000000E+00_dp,  .000000E+00_dp,  .000000E+00_dp/
  DATA (cobi (i,2,3),i=1,7)/ .395367E+02_dp,  .000000E+00_dp,  .000000E+00_dp,  .000000E+00_dp,  .000000E+00_dp,  .000000E+00_dp,  .000000E+00_dp/
  DATA (coali(i,2,3),i=1,7)/ .592629E+00_dp,  .100000E+01_dp,  .000000E+00_dp,  .000000E+00_dp,  .000000E+00_dp,  .000000E+00_dp,  .000000E+00_dp/
  DATA (cobti(i,2,3),i=1,7)/-.106087E+01_dp,  .100000E+01_dp,  .000000E+00_dp,  .000000E+00_dp,  .000000E+00_dp,  .000000E+00_dp,  .000000E+00_dp/
                                                                        
  ! coefficients for O3  in spectral interval 3
  DATA ncgas(3,3) /5/ ; DATA pgas(3,3) /  3039.75_dp/ ; DATA tgas(3,3) /229.900_dp/
  DATA (coai (i,3,3),i=1,7)/ .564000E-03_dp,  .108690E-01_dp,  .124320E-01_dp,  .184417E+00_dp,  .791718E+00_dp,  .000000E+00_dp,  .000000E+00_dp/
  DATA (cobi (i,3,3),i=1,7)/ .191426E+05_dp,  .579429E+03_dp,  .717794E+02_dp,  .187068E+01_dp,  .000000E+00_dp,  .000000E+00_dp,  .000000E+00_dp/
  DATA (coali(i,3,3),i=1,7)/-.204400E-02_dp,  .776840E-01_dp, -.229667E+00_dp,  .994500E-01_dp,  .100000E+01_dp,  .000000E+00_dp,  .000000E+00_dp/
  DATA (cobti(i,3,3),i=1,7)/ .499912E+00_dp,  .485463E+00_dp,  .464581E+00_dp, -.254634E+00_dp,  .100000E+01_dp,  .000000E+00_dp,  .000000E+00_dp/
                                                                        
  ! coefficients for O3  in spectral interval 4
  DATA ncgas(4,3) /0/ ; DATA pgas(4,3) /  3039.75_dp/ ; DATA tgas(4,3) /229.900_dp/
  DATA (coai (i,4,3),i=1,7)/ .540000E-04_dp,  .999946E+00_dp,  .000000E+00_dp,  .000000E+00_dp,  .000000E+00_dp,  .000000E+00_dp,  .000000E+00_dp/
  DATA (cobi (i,4,3),i=1,7)/ .210378E+03_dp,  .000000E+00_dp,  .000000E+00_dp,  .000000E+00_dp,  .000000E+00_dp,  .000000E+00_dp,  .000000E+00_dp/
  DATA (coali(i,4,3),i=1,7)/ .490324E+00_dp,  .100000E+01_dp,  .000000E+00_dp,  .000000E+00_dp,  .000000E+00_dp,  .000000E+00_dp,  .000000E+00_dp/
  DATA (cobti(i,4,3),i=1,7)/ .500000E+00_dp,  .100000E+01_dp,  .000000E+00_dp,  .000000E+00_dp,  .000000E+00_dp,  .000000E+00_dp,  .000000E+00_dp/
                                                                        
  ! coefficients for O3  in spectral interval 5
  DATA ncgas(5,3) /2/ ; DATA pgas(5,3) /  3039.75_dp/ ; DATA tgas(5,3) /229.900_dp/
  DATA (coai (i,5,3),i=1,7)/ .587700E-02_dp,  .994123E+00_dp,  .000000E+00_dp,  .000000E+00_dp,  .000000E+00_dp,  .000000E+00_dp,  .000000E+00_dp/
  DATA (cobi (i,5,3),i=1,7)/ .223357E+03_dp,  .000000E+00_dp,  .000000E+00_dp,  .000000E+00_dp,  .000000E+00_dp,  .000000E+00_dp,  .000000E+00_dp/
  DATA (coali(i,5,3),i=1,7)/ .551312E+00_dp,  .100000E+01_dp,  .000000E+00_dp,  .000000E+00_dp,  .000000E+00_dp,  .000000E+00_dp,  .000000E+00_dp/
  DATA (cobti(i,5,3),i=1,7)/-.140025E+01_dp,  .100000E+01_dp,  .000000E+00_dp,  .000000E+00_dp,  .000000E+00_dp,  .000000E+00_dp,  .000000E+00_dp/
                                                                        
  ! coefficients for O3  in spectral interval 6
  DATA ncgas(6,3) /0/ ; DATA pgas(6,3) /  3039.75_dp/ ; DATA tgas(6,3) /229.900_dp/
  DATA (coai (i,6,3),i=1,7)/ .154100E-02_dp,  .998459E+00_dp,  .000000E+00_dp,  .000000E+00_dp,  .000000E+00_dp,  .000000E+00_dp,  .000000E+00_dp/
  DATA (cobi (i,6,3),i=1,7)/ .221820E+03_dp,  .000000E+00_dp,  .000000E+00_dp,  .000000E+00_dp,  .000000E+00_dp,  .000000E+00_dp,  .000000E+00_dp/
  DATA (coali(i,6,3),i=1,7)/ .546048E+00_dp,  .100000E+01_dp,  .000000E+00_dp,  .000000E+00_dp,  .000000E+00_dp,  .000000E+00_dp,  .000000E+00_dp/
  DATA (cobti(i,6,3),i=1,7)/-.273183E+01_dp,  .100000E+01_dp,  .000000E+00_dp,  .000000E+00_dp,  .000000E+00_dp,  .000000E+00_dp,  .000000E+00_dp/
                                                                        
  ! coefficients for O3  in spectral interval 7
  DATA ncgas(7,3) /7/ ; DATA pgas(7,3) / 10132.50_dp/ ; DATA tgas(7,3) /204.000_dp/
  DATA (coai (i,7,3),i=1,7)/ .220500E-02_dp,  .523500E-02_dp,  .951500E-02_dp,  .578800E-01_dp,  .277389E+00_dp,  .643850E-01_dp,  .583391E+00_dp/
  DATA (cobi (i,7,3),i=1,7)/ .434510E+03_dp,  .299916E+03_dp,  .121339E+03_dp,  .827942E+02_dp,  .157398E+02_dp,  .615177E+01_dp,  .000000E+00_dp/
  DATA (coali(i,7,3),i=1,7)/ .224000E-03_dp,  .100500E-02_dp,  .571600E-02_dp,  .508760E-01_dp,  .524641E+00_dp,  .896800E-01_dp,  .100000E+01_dp/
  DATA (cobti(i,7,3),i=1,7)/ .320370E+01_dp,  .130031E+01_dp, -.332851E+01_dp,  .105177E+01_dp, -.561714E+00_dp, -.357670E+01_dp,  .100000E+01_dp/
                                                                        
  ! coefficients for O3  in spectral interval 8
  DATA ncgas(8,3) /0/ ; DATA pgas(8,3) /  3039.75_dp/ ; DATA tgas(8,3) /229.900_dp/
  DATA (coai (I,8,3),I=1,7)/ .397000E-03_dp,  .999603E+00_dp,  .000000E+00_dp,  .000000E+00_dp,  .000000E+00_dp,  .000000E+00_dp,  .000000E+00_dp/
  DATA (cobi (I,8,3),I=1,7)/ .230675E+03_dp,  .000000E+00_dp,  .000000E+00_dp,  .000000E+00_dp,  .000000E+00_dp,  .000000E+00_dp,  .000000E+00_dp/
  DATA (coali(I,8,3),I=1,7)/ .564371E+00_dp,  .100000E+01_dp,  .000000E+00_dp,  .000000E+00_dp,  .000000E+00_dp,  .000000E+00_dp,  .000000E+00_dp/
  DATA (cobti(I,8,3),I=1,7)/-.479075E+01_dp,  .100000E+01_dp,  .000000E+00_dp,  .000000E+00_dp,  .000000E+00_dp,  .000000E+00_dp,  .000000E+00_dp/
 

! 2. Limits of spectral intervals in the radiation code (for information only)
! ------------------------------------------------------
 
  REAL  (KIND=dp)       ::   &
    grenze(2,2,jpspec)       ! limits of spectral intervals jpspec

  !              WMIN1    WMAX1     WMIN2    WMAX2
  DATA grenze/  1.5300_dp,   4.6420_dp    ,999._dp    ,  0._dp    , &
                0.7000_dp,   1.5300_dp    ,999._dp    ,  0._dp    , &
                0.2451_dp,   0.7000_dp    ,999._dp    ,  0._dp    , &
               20.0000_dp, 104.5150_dp    ,999._dp    ,  0._dp    , &
               12.5000_dp,  20.0000_dp    ,999._dp    ,  0._dp    , &
                8.3333_dp,   9.0090_dp    , 10.3093_dp, 12.5000_dp, &
                9.0090_dp,  10.3093_dp    ,999._dp    ,  0._dp    , &
                4.6420_dp,   8.3333_dp    ,999._dp    ,  0._dp    /


! 3. Rayleigh scattering coefficients in solar spectral intervals
! --------------------------------------------------------------

  REAL  (KIND=dp)       ::   &  ! US: could also be wp, if opt_so is run in wp
    zrsc(jpsol)              ! Coefficients for Rayleigh scattering in solar spectral intervals

  DATA zrsc / 0.59776370E-08_dp, 0.13266702E-06_dp, 0.20634412E-05_dp /

 
! 4. Fraction of solar energy at TOA contained in individual solar spectral intervals 
!    based on data from LABS and NECKEL (1970/1984):
! --------------------------------------------------------------

  REAL  (KIND=dp)       ::   &
    solant(jpsol)            ! Fraction of solar energy at TOA in individual spectral intervals

  DATA solant / 0.12888167_dp, 0.41683156_dp, 0.45428677_dp /

 
! 5. Coefficients for black body radiation and E-type coefficients                   
! --------------------------------------------------------------

  REAL  (KIND=dp)       ::   &
      planck(3,jpther), & !
      zketypr (jpther), & !
      ztetypr (jpther), & !
      zketypa (jpther), & !
      ztetypa (jpther), & !
      zteref
      ! planck: coefficients for the description of the fraction of the total 
      !         black body radiation contained in an thermal spectral interval 
      !         as a function (2.order polynomial) of temperature:
      !
      !         F(T) = PLANCK(1,ISPEC) + PLANCK(2,ISPEC)*T + PLANCK(3,ISPEC)*T**2
      !
      ! zketyp: e-type continuum-coefficient for all spectral intervals 
      !         (PA (H2O)**-2) at 296 K
      !         (r)  following ROBERTS ET AL. 1976
      !         (a)  implicitly derived from the AFGL spectral data
      ! ztetyp: constant for the temperature dependancy of e-type   
      !         absorption for all intervals
      ! zteref: reference temperaure 
 
  DATA planck / 0.157656E+01_dp, -0.711486E-02_dp,  0.908220E-05_dp, &
               -0.761337E-01_dp,  0.339014E-02_dp, -0.703246E-05_dp, &
               -0.353624E+00_dp,  0.321131E-02_dp, -0.472513E-05_dp, &
               -0.180726E+00_dp,  0.148131E-02_dp, -0.195189E-05_dp, &
                0.343422E-01_dp, -0.971936E-03_dp,  0.463714E-05_dp  /

  DATA zketypr / 0.0_dp       , 0.418E-06_dp , 0.859E-06_dp , 0.594E-06_dp , 0.767E-07_dp  /
  DATA ztetypr / 0.0_dp       , 0.0_dp       , 1800.0_dp    , 1800._dp     , 1800._dp      /
  DATA zketypa / 0.8426E-05_dp, 0.8982E-06_dp, 0.5489E-06_dp, 0.4743E-06_dp, 0.7040E-06_dp /
  DATA ztetypa / 1365.55_dp   , 1544.38_dp   , 1699.06_dp   , 1724.39_dp   , 1668.94_dp    /
  DATA zteref  / 296.0_dp     /
 
! 6. Aerosol optical properties for 8 spectral intervals
! ------------------------------------------------------
!US if opt_th/opt_so are run in wp, these variables have to be wp
 
  REAL  (KIND=dp)                  ::           &
  zaea  (jpspec,5),& ! ratio of optical thickness for the absorption in spectral
                     ! interval jpspec  and total optical thickness at 0.55m*1.E-06
                     ! for an aerosoltyp specified by second array index
  zaes  (jpspec,5),& ! analog for the optical thickness of scattering
                     !
  zaeg  (jpspec,5),& ! factor of asymetry for specified aerosoltyp in spectral
                     ! interval jspec

  zaef  (jpspec,5)   ! forward scatterd fraction from aerosols. This array is 
                     ! initialized with 0 but modified later in opt_th and opt_so.
 
                     ! the following aerosoltyps (second array index) are considered:
                     ! 1 : continental
                     ! 2 : maritim
                     ! 3 : urban
                     ! 4 : vulcano ashes
                     ! 5 : stratosphaeric background aerosol

  ! The DATA statements for zaea, zaes, zaeg, zaef have been removed. These 
  ! variables are now set in init_radiation

! 6a: parameters for the vertical distribution of aerosols
! --------------------------------------------------------

  REAL  (KIND=wp), ALLOCATABLE     ::           &
    pvdaes(:),      & ! normalized vertical distribution (sea)
    pvdael(:),      & ! normalized vertical distribution (land)
    pvdaeu(:),      & ! normalized vertical distribution (urban)
    pvdaed(:)         ! normalized vertical distrubution (desert)

  REAL  (KIND=wp)                  ::           &
    ptrbga        , & ! b. optical depths divided by pressure (tropospheric)
    pvobga        , & ! b. optical depths divided by pressure (volcanic)
    pstbga        , & ! b. optical depths divided by pressure (stratospheric)
    paeops        , & ! total optical depths for vertical varying aerosols (sea)
    paeopl        , & ! total optical depths for vertical varying aerosols (land)
    paeopu        , & ! total optical depths for vertical varying aerosols (urban)
    paeopd        , & ! total optical depths for vertical varying aerosols (desert)
    ptrpt         , & ! temperature exponent for the stratosperic definition
    paeadk(3)     , & ! constants for definition of the quantity of water 
    paeadm            ! vapour that will be adsorbed to the dry aerosols to
                      ! form moist aerosols

! 7. Optical properties of liquid water for all 8 spectral intervals (solar and thermal)
! --------------------------------------------------------------------------------------
  ! These data-arrays are used localy in routines opt_so and opt_th
!US if opt_th/opt_so are run in wp, these variables have to be wp

  ! For the calculation of the coefficients, spectral integration is performed 
  ! with NLA LSF and extinction = absorption + scattering is assumed.
  ! Single scattering albedo is used lateron.

  REAL  (KIND=dp)      ::  &
    zlwe(4,jpspec), &  ! 
    zlww(2,jpspec), &  !
    zlwg(2,jpspec), &  !
    zlwemn(jpspec), &  ! minimum values of the extinction coefficients in Pa(H2O)**-1
    zlwemx(jpspec)     ! maximum values of the extinction coefficients in Pa(H2O)**-1

  DATA zlwe / -23.014052_wp,     .173026_wp,     .811865_wp,     .000453_wp, &
              -28.122596_wp,     .172211_wp,     .705673_wp,     .000457_wp, &
              -28.162592_wp,     .198665_wp,     .810637_wp,     .000550_wp, &
             -170.687770_wp,     .498371_wp,     .356225_wp,     .001330_wp, &
              -68.573703_wp,     .263182_wp,     .568143_wp,     .000776_wp, &
             -122.833213_wp,     .297599_wp,     .306486_wp,     .000976_wp, &
             -192.594948_wp,     .440659_wp,     .317142_wp,     .001027_wp, &
              -62.018469_wp,     .281406_wp,     .732715_wp,     .000611_wp  /

  DATA zlww /    .989679_wp,  -22.291412_wp, &
                 .999529_wp,     .020875_wp, &
                 .999999_wp,     .000000_wp, &
                 .302657_wp,  102.711916_wp, &
                 .337398_wp,   80.596716_wp, &
                 .449308_wp,   52.823880_wp, &
                 .686930_wp,  -29.876242_wp, &
                 .804203_wp, -103.022685_wp  /
  DATA zlwg /    .804992_wp,   17.901033_wp, &
                 .814785_wp,   14.204375_wp, &
                 .843955_wp,    8.306586_wp, &
                 .279400_wp,  124.179115_wp, &
                 .499491_wp,  131.635343_wp, &
                 .696708_wp,   75.061613_wp, &
                 .704732_wp,   77.778408_wp, &
                 .784672_wp,   38.002913_wp  /

  DATA zlwemn / 5.000_wp,  5.000_wp,  4.930_wp,  5.800_wp,  5.400_wp,  5.200_wp,  5.500_wp,  5.5000_wp /
  DATA zlwemx /32.500_wp, 32.500_wp, 31.360_wp, 18.600_wp, 24.500_wp, 18.200_wp, 20.200_wp, 32.4000_wp /


! 8. Optical properties of ice clouds for all 8 spectral intervals (solar and thermal)
! --------------------------------------------------------------------------------------
  ! These data-arrays are used localy in routines opt_so and opt_th
!US if opt_th/opt_so are run in wp, these variables have to be wp

  ! The coefficients are derived using spectral averaging by weighted nonlinear LSF;
  ! Recombination of extinction, scattering and absorption after averaging of
  ! extinction = scattering + absorption

  REAL  (KIND=dp)      ::  &
    ziwe(4,jpspec), &  ! 
    ziww(2,jpspec), &  ! coefficients for logarithmic fit
    ziwg(2,jpspec), &  ! coefficients for logarithmic fit
    ziwemn(jpspec), &  ! minimum values of the extinction coefficients in Pa(H2O)**-1
    ziwemx(jpspec)     ! maximum values of the extinction coefficients in Pa(H2O)**-1
 
  DATA ziwe / 16.726535_wp,    0.007465_wp,    1.354626_wp,    0.000112_wp, &
              17.531261_wp,    0.003949_wp,    0.669605_wp,    0.000058_wp, &
              17.698999_wp,    0.003657_wp,    0.625067_wp,    0.000055_wp, &
              19.592746_wp,    0.008644_wp,    1.153213_wp,    0.000101_wp, &
              18.990998_wp,    0.006743_wp,    0.997361_wp,    0.000080_wp, &
              18.482156_wp,    0.004408_wp,    0.693883_wp,    0.000060_wp, &
              18.603168_wp,    0.005260_wp,    0.813026_wp,    0.000064_wp, &
              18.437818_wp,    0.004378_wp,    0.692678_wp,    0.000057_wp  /

  DATA ziww /  0.694631_wp,   -0.022160_wp, &
               0.998669_wp,   -0.000107_wp, &
               0.999993_wp,    0.000000_wp, &
               0.289966_wp,   -0.033855_wp, &
               0.555820_wp,    0.004491_wp, &
               0.554495_wp,    0.004904_wp, &
               0.375319_wp,   -0.017168_wp, &
               0.485290_wp,   -0.004358_wp  /

  DATA ziwg /  0.976960_wp,    0.007938_wp, &
               0.914842_wp,    0.003334_wp, &
               0.900536_wp,    0.001797_wp, &
               1.134025_wp,    0.032141_wp, &
               1.053136_wp,    0.015721_wp, &
               1.010632_wp,    0.006844_wp, &
               1.063545_wp,    0.011772_wp, &
               1.035725_wp,    0.008755_wp  /

  DATA ziwemn/ 2.000_wp, 2.000_wp, 2.000_wp, 2.000_wp, 2.000_wp, 2.000_wp, 2.000_wp, 2.0000_wp /
  DATA ziwemx/30.000_wp,30.000_wp,30.000_wp,30.000_wp,30.000_wp,30.000_wp,30.000_wp, 30.000_wp /

! 9. Optical properties of ice clouds for all 8 spectral intervals (solar and thermal)
! --------------------------------------------------------------------------------------

  !US: this variable is only used in organize_radiation, in working precision
  REAL  (KIND=wp)      ::  &
    rad_csalbw(10)     !  slope of solar albedo with respect to soil water content
                       ! as a function of depth of upper soil layer
                       ! (has been computed every radiation time step before,
                       !  and is now computed in init_radiation)

! 10. Organizational variables / types for a coarse radiation grid
! ----------------------------------------------------------------

! for a coarse radiation grid several grid points from the full COSMO grid are
! grouped together for one coarse grid point. The group can be characterized by
! the start- and end-indices in i- and j-direction. These indices, together with
! some more information, are stored in the following TYPE:

TYPE gp_radcoarse
  ! holds organizational variables for the coarser grid
  INTEGER        :: ilow, ihigh, jlow, jhigh  ! start and end indices of the 
                                              ! group for one coarse radiation
                                              ! grid point
  INTEGER        :: ntgp                      ! total number of grid points
                                              ! belonging to this group
  REAL (KIND=wp) :: rfacgp                    ! factor for averaging (1 / ntgp)
END TYPE gp_radcoarse

TYPE (gp_radcoarse), ALLOCATABLE :: gp_radcoarse_tot(:,:), gp_radcoarse_loc(:,:)

INTEGER  ::              &
  ietot_rc, jetot_rc,    & ! number of coarse radiation gropus in total domain
  ie_rc, je_rc             ! number of coarse radiation groups in a subdomain

INTEGER, ALLOCATABLE ::  &
  mind_ilon_rad(:,:,:),  & ! correspondance between blocked (coarse) radiation
  mind_jlat_rad(:,:,:)     ! grid and the usual COSMO grid (i,j,k)

!==============================================================================

CONTAINS

!==============================================================================
!+ Module procedure in "Radiation" 
!------------------------------------------------------------------------------

SUBROUTINE aerdis ( zetah,  zvdaes, zvdael, zvdaeu, zvdaed, klevp1, &
                    ztrbga, zvobga, zstbga, zaeops, zaeopl, zaeopu, &
                    zaeopd, ztrpt , zaeadk, zaeadm )

!------------------------------------------------------------------------------
!
! Description:
!
! The module procedure aerdis provides parameters for the vertical distribution
! of aerosols (based on the original code of J.F. Geleyn (ECMWF, 4.11.82).
!
! The routine computes the values PVDAE* (* = s, l, u or d for sea, land
! urban or desert) of a surface-normalised vertical distribution of aerosols'
! optical depth from the argument petah (vertical coordinate) at klevp1 levels.
! It also sets values for non-geograpically weighted total optical depths (at
! 55 micrometer wavelength) paeopn for the same four types and similar optical
! depths diveded by pressure for bachground well-mixed aerosols of three types
! p**bga (** = tr, vo or st for tropospheric, volcanic (stratosperic ashes) or
! stratosperic (sulfuric type)). It finally sets values for the power to be
! applied to a temperature ratio smaller than two in order to obtain an index
! one in the stratosphere and zero in the troposphere with a relatively smooth
! transistion (ptrpt), as well as for adsorption coefficients fo water to the
! three type of troposperic aerosols (paeadk) with a minimum value ( in the 
! whole atmosphere) for the sum of the products paeadk by the optical depths
! divided by pressure thickness: paeadm. 
!
! Method:
!
! Straightforward, equivalent heights are given in meters (8434 for the
! atmosphere) and tropospheric and stratospheric pressure boundary values
! are set at 101325 and 19330 Pascal. 
!
!------------------------------------------------------------------------------

! Subroutine arguments:
! --------------------

! Input data
! ----------
  INTEGER                 , INTENT (IN) ::  &
     klevp1           ! number of model layer interfaces      

  REAL    (KIND=wp)       , INTENT (IN) ::  &
     zetah(klevp1)    ! normalized vertical coordinate at half levels

! Output data
! -----------
  REAL    (KIND=wp)       , INTENT (OUT) ::  &
    zvdaes(klevp1), & ! normalized vertical distribution (sea)
    zvdael(klevp1), & ! normalized vertical distribution (land)
    zvdaeu(klevp1), & ! normalized vertical distribution (urban)
    zvdaed(klevp1), & ! normalized vertical distrubution (desert)
    ztrbga        , & ! b. optical depths div. by pressure (tropospheric)
    zvobga        , & ! b. optical depths div. by pressure (volcanic)
    zstbga        , & ! b. optical depths div. by pressure (stratospheric)
    zaeops        , & ! total opt. depths for ver. varying aerosols (sea)
    zaeopl        , & ! total opt. depths for ver. varying aerosols (land)
    zaeopu        , & ! total opt. depths for ver. varying aerosols (urban)
    zaeopd        , & ! total opt. depths for ver. varying aerosols (desert)
    ztrpt         , & ! temperature exponent for the stratosperic definition
    zaeadk(3)     , & ! constants for definition of the quantity of water 
    zaeadm            ! vapour that will be adsorbed to the dry aerosols to
                      ! form moist aerosols

! Local parameters:
! -------------
  REAL (KIND=wp),     PARAMETER  ::  &
    zhss = 8434.0_wp/1000.0_wp ,  & !
    zhsl = 8434.0_wp/1000.0_wp ,  & !
    zhsu = 8434.0_wp/1000.0_wp ,  & !
    zhsd = 8434.0_wp/3000.0_wp      !


!- End of header
!==============================================================================

!------------------------------------------------------------------------------
! Begin Subroutine aerdis              
!------------------------------------------------------------------------------

  zvdaes(1) = 0.0_wp
  zvdael(1) = 0.0_wp
  zvdaeu(1) = 0.0_wp
  zvdaed(1) = 0.0_wp

  IF(zetah(1) /= 0.0_wp) THEN
     zvdaes(1) = zetah(1)**zhss
     zvdael(1) = zetah(1)**zhsl
     zvdaeu(1) = zetah(1)**zhsu
     zvdaed(1) = zetah(1)**zhsd
  END IF

  zvdaes(2:klevp1) = zetah(2:klevp1)**zhss
  zvdael(2:klevp1) = zetah(2:klevp1)**zhsl
  zvdaeu(2:klevp1) = zetah(2:klevp1)**zhsu
  zvdaed(2:klevp1) = zetah(2:klevp1)**zhsd

  ztrbga = 0.03_wp  / (101325.0_wp - 19330.0_wp)
  zvobga = 0.007_wp /  19330.0_wp
  zstbga = 0.045_wp /  19330.0_wp

  zaeops = 0.05_wp
  zaeopl = 0.2_wp
  zaeopu = 0.1_wp
  zaeopd = 1.9_wp
  ztrpt  = 30.0_wp

  zaeadk(1) = 0.3876E-03_wp
  zaeadk(2) = 0.6693E-02_wp
  zaeadk(3) = 0.8563E-03_wp
  zaeadm    = 2.6000E-10_wp


!------------------------------------------------------------------------------
! End of the subroutine 
!------------------------------------------------------------------------------

END SUBROUTINE aerdis

!==============================================================================
!==============================================================================

END MODULE radiation_data
