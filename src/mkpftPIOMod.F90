module mkpftPIOMod

  use shr_kind_mod, only : r8 => shr_kind_r8
  use shr_sys_mod , only : shr_sys_flush
  use mkvarctl    , only : numpft
  use mkdomainMod , only : domain_checksame
  use piofileutils
  use petsc

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
  public mkpft_pio

contains

  subroutine mkpft_pio(ldomain_pio, mapfname, fpft, ndiag, pctlnd_o, pctpft_o)

    use mkdomainPIOMod, only : domain_pio_type, domain_read_pio, domain_clean_pio
    use mkvarpar      , only : numstdpft, numstdcft
    use mkgridmapPIOMod
    use pio

    implicit none
    type(domain_pio_type) , intent(inout) :: ldomain_pio
    character(len=*)      , intent(in)    :: mapfname              ! input mapping file name
    character(len=*)      , intent(in)    :: fpft                  ! input pft dataset file name
    integer               , intent(in)    :: ndiag                 ! unit number for diag out
    real(r8)              , intent(out)   :: pctlnd_o(:)           ! output grid:%land/gridcell
    real(r8)              , pointer       :: pctpft_o(:,:)         ! PFT cover (% of vegetated area)

    type(domain_pio_type) :: tdomain_pio     ! local domain
    type(gridmap_pio_type):: tgridmap_pio
    real(r8) , pointer    :: pctpft3d_i(:,:,:)      ! input grid: PFT percent
    integer               :: numpft_i                        ! num of plant types input data
    real(r8)              :: sum_fldo                        ! global sum of dummy output fld
    real(r8)              :: sum_fldi                        ! global sum of dummy input fld
    real(r8)              :: wst(0:numpft)                   ! as pft_o at specific no
    real(r8)              :: wst_sum                         ! sum of %pft
    real(r8)              :: gpft_o(0:numpft)                ! output grid: global area pfts
    real(r8)              :: garea_o                         ! output grid: global area
    real(r8)              :: gpft_i(0:numpft)                ! input grid: global area pfts
    real(r8)              :: garea_i                         ! input grid: global area
    integer               :: k,n,m,ni,no,ns_loc_i,ns_loc_o   ! indices
    integer               :: dimid                           ! input netCDF id's
    type(var_desc_t)      :: varid
    integer               :: ier                             ! error status
    real(r8)              :: relerr = 0.00001                ! max error: sum overlap wts ne 1

    type(file_desc_t)     :: ncid
    type(iosystem_desc_t) :: pioIoSystem
    type(io_desc_t)       :: iodescNCells
    real     , pointer    :: dataBufferReal(:,:,:)
    real(r8) , pointer    :: pctpft1d_i(:)
    integer  , pointer    :: dataBuffer_int(:,:,:)
    integer  , pointer    :: compdof(:)
    integer  , pointer    :: var_dim_ids(:), dim_glb(:)
    integer               :: count, i, j, vartype, ierr
    integer               :: dim_idx(3,2)

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

    if (.not. use_input_pft) then

       ! Obtain input grid info
       call domain_read_pio(tdomain_pio, fpft)

       call OpenFilePIO(fpft, pioIoSystem, ncid, PIO_NOWRITE)

       call read_float_or_double_3d(tdomain_pio, pioIoSystem, ncid, 'PCT_PFT', 0, dim_idx, pctpft3d_i)

       call PIO_closefile(ncid)
       call PIO_finalize(pioIoSystem, ierr)

    end if


    if ( zero_out ) then

       pctpft_o(:,:) = 0._r8
       pctlnd_o(:)   = 100._r8

    else if ( use_input_pft ) then

       write(6,*) subname//': extend the code for the case use_input_pft = .true.'
       call abort()

    else

       ! Read the map
       call gridmap_mapread_pio(tgridmap_pio, mapfname)

       ns_loc_o = ldomain_pio%ns_loc

       do no = 1, ns_loc_o
          i = tgridmap_pio%dim_nb%begd + no - 1
          pctlnd_o(no) = tgridmap_pio%dst%frac(i) * 100._r8
       end do

       if (ldomain_pio%is_2d) then
          write(6,*)'Extend code to support 2D output grid'
          call abort()
       else
          do no = 1, ns_loc_o
             i = tgridmap_pio%dim_nb%begd + no - 1
             ldomain_pio%frac1d(ldomain_pio%begs + no - 1) = tgridmap_pio%dst%frac(i)
          end do
       end if

       ns_loc_i = (dim_idx(1,2) - dim_idx(1,1) + 1) * (dim_idx(2,2) - dim_idx(2,1) + 1)
       allocate(pctpft1d_i(ns_loc_i))

       do m = 0, numpft

          count = 0
          do j = dim_idx(2,1), dim_idx(2,2)
             do i = dim_idx(1,1), dim_idx(1,2)
                count = count + 1
                pctpft1d_i(count) = pctpft3d_i(i,j,m)
             end do
          end do

          call gridmap_areaave_pio(tgridmap_pio, pctpft1d_i(:), pctpft_o(:,m), nodata=0._r8)

          do no = 1, ns_loc_o
             if (pctlnd_o(no) < 1.0e-6) then
                if (m == 0) then
                   pctpft_o(no,m) = 100._r8
                else
                   pctpft_o(no,m) = 0._r8
                end if
             else
             end if
          end do
       end do

       deallocate(pctpft3d_i)
       deallocate(pctpft1d_i)

    end if

    ! Error check: percents should sum to 100 for land grid cells

    if ( .not. zero_out) then
       do no = 1,ns_loc_o
          wst_sum = 0.
          do m = 0,numpft
             wst_sum = wst_sum + pctpft_o(no,m)
          enddo
          if (abs(wst_sum-100._r8) > 0.00001_r8) then
             write (6,*) subname//'error: pft = ', &
                  (pctpft_o(no,m), m = 0, numpft), &
                  ' do not sum to 100. at no = ',no,' but to ', wst_sum
             stop
          end if
       end do
    end if

    ! Deallocate memory
    call domain_clean_pio(tdomain_pio)
    if ( .not. zero_out .and. .not. use_input_pft) then
       call gridmap_clean_pio(tgridmap_pio)
    end if

    write (6,*) 'Successfully made PFTs'
    write (6,*)

  end subroutine mkpft_pio

end module mkpftPIOMod
