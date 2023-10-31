module mkpftPIOMod

  use shr_kind_mod, only : r8 => shr_kind_r8
  use shr_sys_mod , only : shr_sys_flush
  use mkvarctl    , only : numpft
  use mkdomainMod , only : domain_checksame

  implicit none

  private           ! By default make data private
  !
  ! When pft_idx and pft_frc are set, they must be set together, and they will cause the
  ! entire area to be covered with vegetation and zero out other landunits.
  ! The sum of pft_frc must = 100%, and each pft_idx point in the array corresponds to
  ! the fraction in pft_frc. Only the first few points are used until pft_frc = 0.0.
  !
  integer            :: m                     ! index
  integer, parameter :: maxpft = 24           ! maximum # of PFT
  integer, public    :: num_natpft            ! number of PFTs on the natural vegetation landunit, NOT including bare ground (includes generic crops for runs with create_crop_landunit=false)
  integer, public    :: num_cft               ! number of CFTs on the crop landunit
  integer, public    :: natpft_lb             ! lower bound for natural pft arrays
  integer, public    :: natpft_ub             ! upper bound for natural pft arrays
  integer, public    :: cft_lb                ! lower bound for cft arrays
  integer, public    :: cft_ub                ! upper bound for cft arrays
  integer, public    :: pft_idx(0:maxpft) = & ! PFT vegetation index to override with
       (/ ( -1,  m = 0, maxpft )   /)
  real(r8), public   :: pft_frc(0:maxpft) = & ! PFT vegetation fraction to override with
       (/ ( 0.0, m = 0, maxpft ) /)
  integer, public :: baregroundindex = 0      ! index of bare ground in a natural pft array
  integer, public :: c3cropindex = 15
  integer, public :: c3irrcropindex = 16
  !
  ! !PRIVATE DATA MEMBERS:
  !
  logical, private :: zero_out      = .false. ! Flag to zero out PFT
  logical, private :: use_input_pft = .false. ! Flag to override PFT with input values
  integer, private :: nzero                   ! index of first zero fraction

  !
  ! !PUBLIC MEMBER FUNCTIONS:
  !
  public mkpftPIO

contains

  subroutine mkpftPIO(ldomain_pio, mapfname, fpft, ndiag, pctlnd_o, pctpft_o)

    use mkdomainPIOMod, only : domain_pio_type, domain_read_pio
    use mkvarpar      , only : numstdpft, numstdcft

    implicit none
    type(domain_pio_type) , intent(inout) :: ldomain_pio
    character(len=*)      , intent(in)    :: mapfname              ! input mapping file name
    character(len=*)      , intent(in)    :: fpft                  ! input pft dataset file name
    integer               , intent(in)    :: ndiag                 ! unit number for diag out
    real(r8)              , intent(out)   :: pctlnd_o(:)           ! output grid:%land/gridcell
    real(r8)              , pointer       :: pctpft_o(:,:)         ! PFT cover (% of vegetated area)

    type(domain_pio_type)    :: tdomain_pio     ! local domain
    !type(gridmap_type)    :: tgridmap           ! local gridmap
    real(r8), allocatable :: pctpft_i(:,:)      ! input grid: PFT percent
    integer  :: numpft_i                        ! num of plant types input data
    real(r8) :: sum_fldo                        ! global sum of dummy output fld
    real(r8) :: sum_fldi                        ! global sum of dummy input fld
    real(r8) :: wst(0:numpft)                   ! as pft_o at specific no
    real(r8) :: wst_sum                         ! sum of %pft
    real(r8) :: gpft_o(0:numpft)                ! output grid: global area pfts
    real(r8) :: garea_o                         ! output grid: global area
    real(r8) :: gpft_i(0:numpft)                ! input grid: global area pfts
    real(r8) :: garea_i                         ! input grid: global area
    integer  :: k,n,m,ni,no,ns_i,ns_o_loc       ! indices
    integer  :: ncid,dimid,varid                ! input netCDF id's
    integer  :: ier                             ! error status
    real(r8) :: relerr = 0.00001                ! max error: sum overlap wts ne 1

    character(len=35)  veg(0:maxpft)            ! vegetation types
    character(len=32) :: subname = 'mkpft'
    !-----------------------------------------------------------------------

    write (6,*)
    write (6,*) 'Attempting to make PFTs .....'
    call shr_sys_flush(6)

    ! -----------------------------------------------------------------
    ! Set the vegetation types
    ! -----------------------------------------------------------------
    if ( numpft >= numstdpft )then
       veg(0:maxpft) = (/                                   &
            'not vegetated                      ', &
            'needleleaf evergreen temperate tree', &
            'needleleaf evergreen boreal tree   ', &
            'needleleaf deciduous boreal tree   ', &
            'broadleaf evergreen tropical tree  ', &
            'broadleaf evergreen temperate tree ', &
            'broadleaf deciduous tropical tree  ', &
            'broadleaf deciduous temperate tree ', &
            'broadleaf deciduous boreal tree    ', &
            'broadleaf evergreen shrub          ', &
            'broadleaf deciduous temperate shrub', &
            'broadleaf deciduous boreal shrub   ', &
            'c3 arctic grass                    ', &
            'c3 non-arctic grass                ', &
            'c4 grass                           ', &
            'c3_crop                            ', &
            'c3_irrigated                       ', &
            'corn                               ', &
            'irrigated_corn                     ', &
            'spring_temperate_cereal            ', &
            'irrigated_spring_temperate_cereal  ', &
            'winter_temperate_cereal            ', &
            'irrigated_winter_temperate_cereal  ', &
            'soybean                            ', &
            'irrigated_soybean                  ' /)
    end if

    if (      numpft == numstdpft )then
       write(6,*)'Creating surface datasets with the standard # of PFTs =', numpft
    else if ( numpft > numstdpft )then
       write(6,*)'Creating surface datasets with extra types for crops; total pfts =', numpft
    else
       write(6,*) subname//': parameter numpft is NOT set to a known value (should be 16 or more) =',numpft
       call abort()
    end if

    ns_o_loc = ldomain_pio%ns_loc

    if (.not. use_input_pft) then
       call domain_read_pio(tdomain_pio, fpft)
    end if

  end subroutine mkpftPIO

end module mkpftPIOMod
