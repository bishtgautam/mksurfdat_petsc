module mksoilMod
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: mksoilMod
!
! !DESCRIPTION:
! Make soil data (texture, color,order and organic)
!
! !REVISION HISTORY:
! Author: Erik Kluzek
!
!-----------------------------------------------------------------------
!!USES:
  use shr_kind_mod, only : r8 => shr_kind_r8, r4=>shr_kind_r4
  use shr_sys_mod , only : shr_sys_flush
  use mkdomainMod , only : domain_checksame
  implicit none

  SAVE
  private           ! By default make data private
!
! !PUBLIC MEMBER FUNCTIONS:
!
  public mksoilInit     ! Soil Initialization

  public mksoilAtt      ! Add attributes to output file
  public mksoilAttPIO   ! PIO-version to add attributes to output file

  public mksoiltex      ! Set soil texture
  public mkorganic      ! Set organic soil
  public mkorganic_pio  ! PIO-version of set organic soil
  public mksoilcol      ! Set soil color
  public mksoilcol_pio  ! PIO-version of set soil color
  public mksoilord      ! Set soil order
  public mksoilord_pio  ! PIO-version of set soil order
  public mkfmax         ! Make percent fmax
  public mkfmax_pio     ! PIO-version of make percent fmax

  public mksoiltex_pio  ! Set soil texture
!
! !PUBLIC DATA MEMBERS:
!
  real(r8), public, parameter :: unset = -999.99_r8 ! Flag to signify soil texture override not set
  real(r8), public    :: soil_sand = unset     ! soil texture sand % to override with
  real(r8), public    :: soil_clay = unset     ! soil texture clay % to override with
  real(r8), public    :: soil_fmax = unset     ! soil max saturation frac to override with
  integer , parameter :: unsetcol  = -999      ! flag to indicate soil color NOT set
  integer , public    :: soil_color= unsetcol  ! soil color to override with
  integer , parameter :: unsetord  = -999      ! flag to indicate soil order NOT set
  integer , public    :: soil_order= unsetord  ! soil order to override with
!
! !PRIVATE DATA MEMBERS:
!
! !PRIVATE MEMBER FUNCTIONS:
  private :: mkrank
  private :: mksoiltexInit  ! Soil texture Initialization
  private :: mksoilcolInit  ! Soil color Initialization
  private :: mksoilordInit  ! Soil order Initialization
  private :: mksoilfmaxInit ! Soil fmax Initialization

!EOP
!===============================================================
contains
!===============================================================

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: mksoilInit
!
! !INTERFACE:
subroutine mksoilInit( )
!
! !DESCRIPTION:
! Initialize the different soil types
! !USES:
!
! !ARGUMENTS:
  implicit none
!
! !REVISION HISTORY:
! Author: Erik Kluzek
!
!
! !LOCAL VARIABLES:
!EOP
!-----------------------------------------------------------------------
  call mksoiltexInit()
  call mksoilcolInit()
  call mksoilordInit()
  call mksoilfmaxInit()

end subroutine mksoilInit

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: mksoiltexInit
!
! !INTERFACE:
subroutine mksoiltexInit( )
!
! !DESCRIPTION:
! Initialize of make soil texture
! !USES:
!
! !ARGUMENTS:
  implicit none
!
! !REVISION HISTORY:
! Author: Erik Kluzek
!
!
! !LOCAL VARIABLES:
!EOP
  real(r8) :: sumtex
  character(len=32) :: subname = 'mksoiltexInit'
!-----------------------------------------------------------------------
    if ( soil_clay /= unset )then
       write(6,*) 'Replace soil clay % for all points with: ', soil_clay
       if ( soil_sand == unset )then
           write (6,*) subname//':error: soil_clay set, but NOT soil_sand'
           call abort()
       end if
    end if
    if ( soil_sand /= unset )then
       write(6,*) 'Replace soil sand % for all points with: ', soil_sand
       if ( soil_clay == unset )then
           write (6,*) subname//':error: soil_sand set, but NOT soil_clay'
           call abort()
       end if
       sumtex = soil_sand + soil_clay
       if ( sumtex < 0.0_r8 .or. sumtex > 100.0_r8 )then
           write (6,*) subname//':error: soil_sand and soil_clay out of bounds: sand, clay = ', &
                       soil_sand, soil_clay
           call abort()
       end if
    end if

end subroutine mksoiltexInit

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: mksoiltex
!
! !INTERFACE:
subroutine mksoiltex(ldomain, mapfname, datfname, ndiag, sand_o, clay_o)
!
! !DESCRIPTION:
! make %sand and %clay from IGBP soil data, which includes
! igbp soil 'mapunits' and their corresponding textures
!
! !USES:
  use mkdomainMod, only : domain_type, domain_clean, domain_read
  use mkgridmapMod
  use mkvarpar
  use mkvarctl    
  use mkncdio
!
! !ARGUMENTS:
  implicit none
  type(domain_type), intent(in) :: ldomain
  character(len=*)  , intent(in) :: mapfname      ! input mapping file name
  character(len=*)  , intent(in) :: datfname      ! input data file name
  integer           , intent(in) :: ndiag         ! unit number for diag out
  real(r8)          , intent(out):: sand_o(:,:)   ! % sand (output grid)
  real(r8)          , intent(out):: clay_o(:,:)   ! % clay (output grid)
!
! !CALLED FROM:
! subroutine mksrfdat in module mksrfdatMod
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
!
! !LOCAL VARIABLES:
!EOP
  type(gridmap_type)    :: tgridmap
  type(domain_type)    :: tdomain          ! local domain
  character(len=38)  :: typ                 ! soil texture based on ...
  integer  :: nlay                          ! number of soil layers
  integer  :: mapunitmax                    ! max value of igbp soil mapunits
  integer  :: mapunittemp                   ! temporary igbp soil mapunit
  integer  :: maxovr
  integer , allocatable :: novr(:)
  integer , allocatable :: kmap(:,:)
  real(r8), allocatable :: kwgt(:,:)
  integer , allocatable :: kmax(:)
  real(r8), allocatable :: wst(:)
  real(r8), allocatable :: sand_i(:,:)      ! input grid: percent sand
  real(r8), allocatable :: clay_i(:,:)      ! input grid: percent clay
  real(r8), allocatable :: mapunit_i(:)     ! input grid: igbp soil mapunits
  integer, parameter :: num=2               ! set soil mapunit number
  integer  :: wsti(num)                     ! index to 1st and 2nd largest wst
  integer, parameter :: nlsm=4              ! number of soil textures 
  character(len=38)  :: soil(0:nlsm)        ! name of each soil texture
  real(r8) :: gast_i(0:nlsm)                ! global area, by texture type
  real(r8) :: gast_o(0:nlsm)                ! global area, by texture type
  real(r8) :: wt                            ! map overlap weight
  real(r8) :: sum_fldi                      ! global sum of dummy input fld
  real(r8) :: sum_fldo                      ! global sum of dummy output fld
  integer  :: l,k,n,m,ni,no,ns_i,ns_o       ! indices
  integer  :: k1,k2                         ! indices
  integer  :: ncid,dimid,varid              ! input netCDF id's
  integer  :: ier                           ! error status
  integer  :: miss = 99999                  ! missing data indicator
  real(r8) :: relerr = 0.00001              ! max error: sum overlap wts ne 1
  logical  :: found                         ! temporary
  integer  :: kmap_max                      ! maximum overlap weights
  integer, parameter :: kmap_max_min   = 90 ! kmap_max mininum value
  integer, parameter :: km_mx_ns_prod = 160000 ! product of kmap_max*ns_o to keep constant
  character(len=32) :: subname = 'mksoiltex'
!-----------------------------------------------------------------------

  write (6,*) 'Attempting to make %sand and %clay .....'
  call shr_sys_flush(6)
  write(*,*)'mapfname:' ,trim(mapfname)
  write(*,*)'datfname:' ,trim(datfname)

  ! -----------------------------------------------------------------
  ! Define the model surface types: 0 to nlsm
  ! -----------------------------------------------------------------

  soil(0) = 'no soil: ocean, glacier, lake, no data'
  soil(1) = 'clays                                 '
  soil(2) = 'sands                                 '
  soil(3) = 'loams                                 '
  soil(4) = 'silts                                 '

  ! -----------------------------------------------------------------
  ! Read input file
  ! -----------------------------------------------------------------

  ! Obtain input grid info, read local fields

  call domain_read(tdomain,datfname)
  ns_i = tdomain%ns
  ns_o = ldomain%ns

  write (6,*) 'Open soil texture file: ', trim(datfname)
  call check_ret(nf_open(datfname, 0, ncid), subname)
  call check_ret(nf_inq_dimid  (ncid, 'number_of_layers', dimid), subname)
  call check_ret(nf_inq_dimlen (ncid, dimid, nlay), subname)

  call check_ret(nf_inq_dimid  (ncid, 'max_value_mapunit', dimid), subname)
  call check_ret(nf_inq_dimlen (ncid, dimid, mapunitmax), subname)

  allocate(sand_i(mapunitmax,nlay), &
           clay_i(mapunitmax,nlay), &
           mapunit_i(ns_i), stat=ier)
  if (ier/=0) call abort()

  call check_ret(nf_inq_varid (ncid, 'MAPUNITS', varid), subname)
  call check_ret(nf_get_var_double (ncid, varid, mapunit_i), subname)

  call check_ret(nf_inq_varid (ncid, 'PCT_SAND', varid), subname)
  call check_ret(nf_get_var_double (ncid, varid, sand_i), subname)

  call check_ret(nf_inq_varid (ncid, 'PCT_CLAY', varid), subname)
  call check_ret(nf_get_var_double (ncid, varid, clay_i), subname)

  call check_ret(nf_close(ncid), subname)

  ! Compute local fields _o
  if (soil_sand==unset .and. soil_clay==unset) then

     call gridmap_mapread(tgridmap, mapfname)

     ! Error checks for domain and map consistencies

     call domain_checksame( tdomain, ldomain, tgridmap )

     ! kmap_max are the maximum number of mapunits that will consider on
     ! any output gridcell - this is set currently above and can be changed
     ! kmap(:) are the mapunit values on the input grid
     ! kwgt(:) are the weights on the input grid

     allocate(novr(ns_o))
     novr(:) = 0
     do n = 1,tgridmap%ns
        ni = tgridmap%src_indx(n)
        no = tgridmap%dst_indx(n)
        wt = tgridmap%wovr(n)
        novr(no) = novr(no) + 1
     end do
     maxovr = maxval(novr(:))
     kmap_max = min(maxovr,max(kmap_max_min,km_mx_ns_prod/ns_o))
     deallocate(novr)

     write(6,*)'kmap_max= ',kmap_max,' maxovr= ',maxovr,' ns_o= ',ns_o,' size= ',(kmap_max+1)*ns_o

     allocate(kmap(0:kmap_max,ns_o), stat=ier)
     if (ier/=0) call abort()
     allocate(kwgt(0:kmap_max,ns_o), stat=ier)
     if (ier/=0) call abort()
     allocate(kmax(ns_o), stat=ier)
     if (ier/=0) call abort()
     allocate(wst(0:kmap_max), stat=ier)
     if (ier/=0) call abort()

     kwgt(:,:) = 0.
     kmap(:,:) = 0
     kmax(:) = 0

     do n = 1,tgridmap%ns
        ni = tgridmap%src_indx(n)
        no = tgridmap%dst_indx(n)
        wt = tgridmap%wovr(n)
        if (tgridmap%frac_src(ni) > 0) then
           k = int(mapunit_i(ni))
        else
           k = 0
        end if
        found = .false.
        do l = 0,kmax(no)
           if (k == kmap(l,no)) then
              kwgt(l,no) = kwgt(l,no) + wt
              kmap(l,no) = k
              found = .true.
              exit
           end if
        end do
        if (.not. found) then
           kmax(no) = kmax(no) + 1
           if (kmax(no) > kmap_max) then
              write(6,*)'kmax is > kmap_max= ',kmax(no), 'kmap_max = ', &
                         kmap_max,' for no = ',no
              write(6,*)'reset kmap_max in mksoilMod to a greater value'
              stop
           end if
           kmap(kmax(no),no) = k
           kwgt(kmax(no),no) = wt
        end if
     enddo

  end if

  do no = 1,ns_o

     if (soil_sand==unset .and. soil_clay==unset) then
        wst(:) = 0.
        wst(0:kmax(no)) = kwgt(0:kmax(no),no)

        ! Rank non-zero weights by soil mapunit.
        ! k1 is the most extensive mapunit.
        ! k2 is the second most extensive mapunit.

        if (maxval(wst(:)) > 0) then
           call mkrank (kmax(no)+1, wst(0:kmax(no)), miss, wsti, num)
           k1 = kmap(wsti(1),no)
           if (wsti(2) == miss) then
              k2 = miss
           else
              k2 = kmap(wsti(2),no)
           end if
        else
           k1 = 0
           k2 = 0
        end if

     end if

     ! Set soil texture as follows:
     !   a. Use dominant igbp soil mapunit based on area of overlap unless
     !     'no data' is dominant
     !   b. In this case use second most dominant mapunit if it has data
     !   c. If this has no data or if there isn't a second most dominant
     !      mapunit, use loam for soil texture

     if (soil_sand/=unset .and. soil_clay/=unset) then  !---soil texture is input
        do l = 1, nlay
           sand_o(no,l) = soil_sand
           clay_o(no,l) = soil_clay
        end do
     else if (k1 /= 0) then           !---not 'no data'
        do l = 1, nlay
           sand_o(no,l) = sand_i(k1,l)
           clay_o(no,l) = clay_i(k1,l)
        end do
     else                                  !---if (k1 == 0) then
        if (k2 == 0 .or. k2 == miss) then     !---no data
           do l = 1, nlay
              sand_o(no,l) = 43.           !---use loam
              clay_o(no,l) = 18.
           end do
        else                               !---if (k2 /= 0 and /= miss)
           do l = 1, nlay
              sand_o(no,l) = sand_i(k2,l)
              clay_o(no,l) = clay_i(k2,l)
           end do
        end if       !---end of k2 if-block
     end if          !---end of k1 if-block

  enddo

  if (soil_sand==unset .and. soil_clay==unset) then

     ! Global sum of output field 

     sum_fldi = 0.0_r8
     do ni = 1,ns_i
        sum_fldi = sum_fldi + tgridmap%area_src(ni)*tgridmap%frac_src(ni)*re**2
     enddo

     sum_fldo = 0.
     do no = 1,ns_o
        sum_fldo = sum_fldo + tgridmap%area_dst(no)*tgridmap%frac_dst(no)*re**2
     end do

     ! -----------------------------------------------------------------
     ! Error check1
     ! Compare global sum fld_o to global sum fld_i.
     ! -----------------------------------------------------------------

     if ( trim(mksrf_gridtype) == 'global') then
        if ( abs(sum_fldo/sum_fldi-1.) > relerr ) then
           write (6,*) 'MKSOILTEX error: input field not conserved'
           write (6,'(a30,e20.10)') 'global sum output field = ',sum_fldo
           write (6,'(a30,e20.10)') 'global sum input  field = ',sum_fldi
           stop
        end if
     end if

     ! -----------------------------------------------------------------
     ! Error check2
     ! Compare global area of each soil type on input and output grids
     ! -----------------------------------------------------------------

     ! input grid: global areas by texture class

     gast_i(:) = 0.
     do l = 1, nlay
        do ni = 1,ns_i
           mapunittemp = nint(mapunit_i(ni))
           if (mapunittemp==0) then
              typ = 'no soil: ocean, glacier, lake, no data'
           else if (clay_i(mapunittemp,l) >= 40.) then
              typ = 'clays'
           else if (sand_i(mapunittemp,l) >= 50.) then
              typ = 'sands'
           else if (clay_i(mapunittemp,l)+sand_i(mapunittemp,l) < 50.) then
              if (tdomain%mask(ni) /= 0.) then
                 typ = 'silts'
              else            !if (tdomain%mask(ni) == 0.) then no data
                 typ = 'no soil: ocean, glacier, lake, no data'
              end if
           else
              typ = 'loams'
           end if
           do m = 0, nlsm
              if (typ == soil(m)) go to 101
           end do
           write (6,*) 'MKSOILTEX error: sand = ',sand_i(mapunittemp,l), &
             ' clay = ',clay_i(mapunittemp,l), &
             ' not assigned to soil type for input grid lon,lat,layer = ',ni,l
           call abort()
101        continue
           gast_i(m) = gast_i(m) + tgridmap%area_src(ni)*tgridmap%frac_src(ni)*re**2
        end do
     end do

     ! output grid: global areas by texture class

     gast_o(:) = 0.
     do l = 1, nlay
        do no = 1,ns_o
           if (clay_o(no,l)==0. .and. sand_o(no,l)==0.) then
              typ = 'no soil: ocean, glacier, lake, no data'
           else if (clay_o(no,l) >= 40.) then
              typ = 'clays'
           else if (sand_o(no,l) >= 50.) then
              typ = 'sands'
           else if (clay_o(no,l)+sand_o(no,l) < 50.) then
              typ = 'silts'
           else
              typ = 'loams'
           end if
           do m = 0, nlsm
              if (typ == soil(m)) go to 102
           end do
           write (6,*) 'MKSOILTEX error: sand = ',sand_o(no,l), &
             ' clay = ',clay_o(no,l), &
             ' not assigned to soil type for output grid lon,lat,layer = ',no,l
           call abort()
102        continue
           gast_o(m) = gast_o(m) + tgridmap%area_dst(no)*tgridmap%frac_dst(no)*re**2
        end do
     end do

     ! Diagnostic output

     write (ndiag,*)
     write (ndiag,'(1x,70a1)') ('=',l=1,70)
     write (ndiag,*) 'Soil Texture Output'
     write (ndiag,'(1x,70a1)') ('=',l=1,70)
     write (ndiag,*)

     write (ndiag,*) 'The following table of soil texture classes is for comparison only.'
     write (ndiag,*) 'The actual data is continuous %sand, %silt and %clay not textural classes'
     write (ndiag,*)

     write (ndiag,*)
     write (ndiag,'(1x,70a1)') ('.',l=1,70)
     write (ndiag,1001)
1001 format (1x,'soil texture class',17x,' input grid area output grid area',/ &
             1x,33x,'     10**6 km**2','      10**6 km**2')
     write (ndiag,'(1x,70a1)') ('.',l=1,70)
     write (ndiag,*)

     do l = 0, nlsm
        write (ndiag,1002) soil(l),gast_i(l)*1.e-6,gast_o(l)*1.e-6
1002    format (1x,a38,f16.3,f17.3)
     end do

  end if

  ! Deallocate dynamic memory

  call domain_clean(tdomain)
  if (soil_sand==unset .and. soil_clay==unset) then
     call gridmap_clean(tgridmap)
     deallocate (kmap, kwgt, kmax, wst)
     deallocate (sand_i,clay_i,mapunit_i)
  end if


  write (6,*) 'Successfully made %sand and %clay'
  write (6,*)
  call shr_sys_flush(6)

end subroutine mksoiltex

!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: mksoiltex_pio
!
! !INTERFACE:
subroutine mksoiltex_pio(ldomain_pio, mapfname, datfname, ndiag, sand_o, clay_o)
  !
  ! !DESCRIPTION:
  ! make %sand and %clay from IGBP soil data, which includes
  ! igbp soil 'mapunits' and their corresponding textures
  !
  ! !USES:
  use mkdomainPIOMod, only : domain_pio_type, domain_read_pio, domain_clean_pio
  use mkgridmapMod
  use mkvarpar
  use mkvarctl    
  use mkncdio
  use pio
  use piofileutils
  !
  ! !ARGUMENTS:
  implicit none
  type(domain_pio_type) , intent(inout) :: ldomain_pio
  character(len=*)  , intent(in) :: mapfname      ! input mapping file name
  character(len=*)  , intent(in) :: datfname      ! input data file name
  integer           , intent(in) :: ndiag         ! unit number for diag out
  real(r8)          , intent(out):: sand_o(:,:)   ! % sand (output grid)
  real(r8)          , intent(out):: clay_o(:,:)   ! % clay (output grid)
  !
  ! !CALLED FROM:
  ! subroutine mksrfdat in module mksrfdatMod
  !
  ! !REVISION HISTORY:
  ! Author: Gautam Bisht
  !
  !
  ! !LOCAL VARIABLES:
  !EOP
  type(gridmap_type)    :: tgridmap
  type(domain_pio_type) :: tdomain_pio     ! local domain
  character(len=38)  :: typ                 ! soil texture based on ...
  integer  :: nlay                          ! number of soil layers
  integer  :: mapunittemp                   ! temporary igbp soil mapunit
  integer  :: maxovr
  integer , allocatable :: novr(:)
  integer , allocatable :: kmap(:,:)
  real(r8), allocatable :: kwgt(:,:)
  integer , allocatable :: kmax(:)
  real(r8), allocatable :: wst(:)
  real(r8), pointer :: sand3d_i(:,:,:)      ! input grid: percent sand
  real(r8), pointer :: clay3d_i(:,:,:)      ! input grid: percent clay
  real(r8), pointer :: sand2d_i(:,:)      ! input grid: percent sand
  real(r8), pointer :: clay2d_i(:,:)      ! input grid: percent clay
  real(r8), pointer :: tmp1d_i(:)
  real(r8), pointer :: mapunit2d_i(:,:)     ! input grid: igbp soil mapunits
  real(r8), pointer :: mapunit1d_i(:)       ! input grid: igbp soil mapunits
  integer, parameter :: num=2               ! set soil mapunit number
  integer  :: wsti(num)                     ! index to 1st and 2nd largest wst
  integer, parameter :: nlsm=4              ! number of soil textures 
  character(len=38)  :: soil(0:nlsm)        ! name of each soil texture
  real(r8) :: gast_i(0:nlsm)                ! global area, by texture type
  real(r8) :: gast_o(0:nlsm)                ! global area, by texture type
  real(r8) :: wt                            ! map overlap weight
  real(r8) :: sum_fldi                      ! global sum of dummy input fld
  real(r8) :: sum_fldo                      ! global sum of dummy output fld
  integer  :: l,k,n,m,ni,no,ns_loc_i,ns_loc_o       ! indices
  integer  :: k1,k2                         ! indices
  integer  :: ier                           ! error status
  integer  :: miss = 99999                  ! missing data indicator
  real(r8) :: relerr = 0.00001              ! max error: sum overlap wts ne 1
  logical  :: found                         ! temporary
  integer  :: kmap_max                      ! maximum overlap weights
  integer, parameter :: kmap_max_min   = 90 ! kmap_max mininum value
  integer, parameter :: km_mx_ns_prod = 160000 ! product of kmap_max*ns_o to keep constant
  integer, pointer    :: vec_row_indices(:)

  type(file_desc_t)     :: ncid
  type(iosystem_desc_t) :: pioIoSystem
  integer               :: count, i, j, ierr
  integer               :: dim_idx_3d(3,2)
  integer               :: dim_idx_2d(3,2)
  !-----------------------------------------------------------------------

  write (6,*) 'Attempting to make %sand and %clay .....'
  call shr_sys_flush(6)
  write(*,*)'mapfname:' ,trim(mapfname)
  write(*,*)'datfname:' ,trim(datfname)

  ! -----------------------------------------------------------------
  ! Define the model surface types: 0 to nlsm
  ! -----------------------------------------------------------------

  soil(0) = 'no soil: ocean, glacier, lake, no data'
  soil(1) = 'clays                                 '
  soil(2) = 'sands                                 '
  soil(3) = 'loams                                 '
  soil(4) = 'silts                                 '

  ! -----------------------------------------------------------------
  ! Read input file
  ! -----------------------------------------------------------------

  ! Obtain input grid info, read local fields

  write (6,*) 'Open soil texture file: ', trim(datfname)

  call domain_read_pio(tdomain_pio, datfname)

  call OpenFilePIO(datfname, pioIoSystem, ncid, PIO_NOWRITE)

  call read_float_or_double_2d(tdomain_pio, pioIoSystem, ncid, 'MAPUNITS', dim_idx_2d, vec_row_indices, mapunit2d_i)
  call read_float_or_double_3d(tdomain_pio, pioIoSystem, ncid, 'PCT_SAND', 0, dim_idx_3d, vec_row_indices, sand3d_i)
  call read_float_or_double_3d(tdomain_pio, pioIoSystem, ncid, 'PCT_CLAY', 0, dim_idx_3d, vec_row_indices, clay3d_i)

  call PIO_closefile(ncid)
  call PIO_finalize(pioIoSystem, ierr)

  call gridmap_mapread(tgridmap, mapfname)

  ! Convert 2D vector to 1D vector
  ns_loc_i = (dim_idx_2d(1,2) - dim_idx_2d(1,1) + 1) * (dim_idx_2d(2,2) - dim_idx_2d(2,1) + 1)

  nlay = dim_idx_3d(3,2)
  allocate(sand2d_i(ns_loc_i,nlay))
  allocate(clay2d_i(ns_loc_i,nlay))
  allocate(tmp1d_i(ns_loc_i))

  do l = 1, nlay
     count = 0
     do j = dim_idx_3d(2,1), dim_idx_3d(2,2)
        do i = dim_idx_3d(1,1), dim_idx_3d(1,2)
           count = count + 1
           sand2d_i(count,l) = sand3d_i(i,j,l)
           clay2d_i(count,l) = clay3d_i(i,j,l)
        end do
     end do
  end do
  
  ! Compute local fields _o
  if (soil_sand==unset .and. soil_clay==unset) then

     allocate(mapunit1d_i(ns_loc_i))

     count = 0
     do j = dim_idx_2d(2,1), dim_idx_2d(2,2)
        do i = dim_idx_2d(1,1), dim_idx_2d(1,2)
           count = count + 1
           mapunit1d_i(count) = mapunit2d_i(i,j)
        end do
     end do

     ! kmap_max are the maximum number of mapunits that will consider on
     ! any output gridcell - this is set currently above and can be changed
     ! kmap(:) are the mapunit values on the input grid
     ! kwgt(:) are the weights on the input grid

     ns_loc_o = ldomain_pio%ns_loc

     allocate(novr(ns_loc_o))
     novr(:) = 0._r8
     do n = 1,tgridmap%ns
        ni = tgridmap%src_indx(n)
        no = tgridmap%dst_indx(n)
        wt = tgridmap%wovr(n)
        novr(no) = novr(no) + 1
     end do
     maxovr = maxval(novr(:))
     kmap_max = min(maxovr,max(kmap_max_min,km_mx_ns_prod/ns_loc_o))
     deallocate(novr)

     write(6,*)'kmap_max= ',kmap_max,' maxovr= ',maxovr,' ns_loc_o= ',ns_loc_o,' size= ',(kmap_max+1)*ns_loc_o

     allocate(kmap(0:kmap_max,ns_loc_o), stat=ier)
     if (ier/=0) call abort()
     allocate(kwgt(0:kmap_max,ns_loc_o), stat=ier)
     if (ier/=0) call abort()
     allocate(kmax(ns_loc_o), stat=ier)
     if (ier/=0) call abort()
     allocate(wst(0:kmap_max), stat=ier)
     if (ier/=0) call abort()

     kwgt(:,:) = 0.
     kmap(:,:) = 0
     kmax(:) = 0

     do n = 1,tgridmap%ns
        ni = tgridmap%src_indx(n)
        no = tgridmap%dst_indx(n)
        wt = tgridmap%wovr(n)
        if (tgridmap%frac_src(ni) > 0) then
           k = int(mapunit1d_i(ni))
        else
           k = 0
        end if
        found = .false.
        do l = 0,kmax(no)
           if (k == kmap(l,no)) then
              kwgt(l,no) = kwgt(l,no) + wt
              kmap(l,no) = k
              found = .true.
              exit
           end if
        end do
        if (.not. found) then
           kmax(no) = kmax(no) + 1
           if (kmax(no) > kmap_max) then
              write(6,*)'kmax is > kmap_max= ',kmax(no), 'kmap_max = ', &
                         kmap_max,' for no = ',no
              write(6,*)'reset kmap_max in mksoilMod to a greater value'
              stop
           end if
           kmap(kmax(no),no) = k
           kwgt(kmax(no),no) = wt
        end if
     enddo

  end if

#if 1
  do no = 1,ns_loc_o

     if (soil_sand==unset .and. soil_clay==unset) then
        wst(:) = 0.
        wst(0:kmax(no)) = kwgt(0:kmax(no),no)

        ! Rank non-zero weights by soil mapunit.
        ! k1 is the most extensive mapunit.
        ! k2 is the second most extensive mapunit.

        if (maxval(wst(:)) > 0) then
           call mkrank (kmax(no)+1, wst(0:kmax(no)), miss, wsti, num)
           k1 = kmap(wsti(1),no)
           if (wsti(2) == miss) then
              k2 = miss
           else
              k2 = kmap(wsti(2),no)
           end if
        else
           k1 = 0
           k2 = 0
        end if

     end if

     ! Set soil texture as follows:
     !   a. Use dominant igbp soil mapunit based on area of overlap unless
     !     'no data' is dominant
     !   b. In this case use second most dominant mapunit if it has data
     !   c. If this has no data or if there isn't a second most dominant
     !      mapunit, use loam for soil texture

     if (soil_sand/=unset .and. soil_clay/=unset) then  !---soil texture is input
        do l = 1, nlay
           sand_o(no,l) = soil_sand
           clay_o(no,l) = soil_clay
        end do
     else if (k1 /= 0) then           !---not 'no data'
        do l = 1, nlay
           sand_o(no,l) = sand2d_i(k1,l)
           clay_o(no,l) = clay2d_i(k1,l)
        end do
     else                                  !---if (k1 == 0) then
        if (k2 == 0 .or. k2 == miss) then     !---no data
           do l = 1, nlay
              sand_o(no,l) = 43.           !---use loam
              clay_o(no,l) = 18.
           end do
        else                               !---if (k2 /= 0 and /= miss)
           do l = 1, nlay
              sand_o(no,l) = sand2d_i(k2,l)
              clay_o(no,l) = clay2d_i(k2,l)
           end do
        end if       !---end of k2 if-block
     end if          !---end of k1 if-block

  enddo

  if (soil_sand==unset .and. soil_clay==unset) then

     ! Global sum of output field 

     sum_fldi = 0.0_r8
     do ni = 1,ns_loc_i
        sum_fldi = sum_fldi + tgridmap%area_src(ni)*tgridmap%frac_src(ni)*re**2
     enddo

     sum_fldo = 0.
     do no = 1,ns_loc_o
        sum_fldo = sum_fldo + tgridmap%area_dst(no)*tgridmap%frac_dst(no)*re**2
     end do

     ! -----------------------------------------------------------------
     ! Error check1
     ! Compare global sum fld_o to global sum fld_i.
     ! -----------------------------------------------------------------

     if ( trim(mksrf_gridtype) == 'global') then
        if ( abs(sum_fldo/sum_fldi-1.) > relerr ) then
           write (6,*) 'MKSOILTEX error: input field not conserved'
           write (6,'(a30,e20.10)') 'global sum output field = ',sum_fldo
           write (6,'(a30,e20.10)') 'global sum input  field = ',sum_fldi
           stop
        end if
     end if

     ! -----------------------------------------------------------------
     ! Error check2
     ! Compare global area of each soil type on input and output grids
     ! -----------------------------------------------------------------

     ! input grid: global areas by texture class

     gast_i(:) = 0.
     do l = 1, nlay
        do ni = 1,ns_loc_i
           mapunittemp = nint(mapunit1d_i(ni))
           if (mapunittemp==0) then
              typ = 'no soil: ocean, glacier, lake, no data'
           else if (clay2d_i(mapunittemp,l) >= 40.) then
              typ = 'clays'
           else if (sand2d_i(mapunittemp,l) >= 50.) then
              typ = 'sands'
           else if (clay2d_i(mapunittemp,l)+sand2d_i(mapunittemp,l) < 50.) then
              !if (tdomain_pio%mask(ni) /= 0.) then
                 typ = 'silts'
              !else            !if (tdomain%mask(ni) == 0.) then no data
              !   typ = 'no soil: ocean, glacier, lake, no data'
              !end if
           else
              typ = 'loams'
           end if
           do m = 0, nlsm
              if (typ == soil(m)) go to 101
           end do
           write (6,*) 'MKSOILTEX error: sand = ',sand2d_i(mapunittemp,l), &
             ' clay = ',clay2d_i(mapunittemp,l), &
             ' not assigned to soil type for input grid lon,lat,layer = ',ni,l
           call abort()
101        continue
           gast_i(m) = gast_i(m) + tgridmap%area_src(ni)*tgridmap%frac_src(ni)*re**2
        end do
     end do

     ! output grid: global areas by texture class

     gast_o(:) = 0.
     do l = 1, nlay
        do no = 1,ns_loc_o
           if (clay_o(no,l)==0. .and. sand_o(no,l)==0.) then
              typ = 'no soil: ocean, glacier, lake, no data'
           else if (clay_o(no,l) >= 40.) then
              typ = 'clays'
           else if (sand_o(no,l) >= 50.) then
              typ = 'sands'
           else if (clay_o(no,l)+sand_o(no,l) < 50.) then
              typ = 'silts'
           else
              typ = 'loams'
           end if
           do m = 0, nlsm
              if (typ == soil(m)) go to 102
           end do
           write (6,*) 'MKSOILTEX error: sand = ',sand_o(no,l), &
             ' clay = ',clay_o(no,l), &
             ' not assigned to soil type for output grid lon,lat,layer = ',no,l
           call abort()
102        continue
           gast_o(m) = gast_o(m) + tgridmap%area_dst(no)*tgridmap%frac_dst(no)*re**2
        end do
     end do

     ! Diagnostic output

     write (ndiag,*)
     write (ndiag,'(1x,70a1)') ('=',l=1,70)
     write (ndiag,*) 'Soil Texture Output'
     write (ndiag,'(1x,70a1)') ('=',l=1,70)
     write (ndiag,*)

     write (ndiag,*) 'The following table of soil texture classes is for comparison only.'
     write (ndiag,*) 'The actual data is continuous %sand, %silt and %clay not textural classes'
     write (ndiag,*)

     write (ndiag,*)
     write (ndiag,'(1x,70a1)') ('.',l=1,70)
     write (ndiag,1001)
1001 format (1x,'soil texture class',17x,' input grid area output grid area',/ &
             1x,33x,'     10**6 km**2','      10**6 km**2')
     write (ndiag,'(1x,70a1)') ('.',l=1,70)
     write (ndiag,*)

     do l = 0, nlsm
        write (ndiag,1002) soil(l),gast_i(l)*1.e-6,gast_o(l)*1.e-6
1002    format (1x,a38,f16.3,f17.3)
     end do

  end if
#endif

  ! Deallocate dynamic memory

  call domain_clean_pio(tdomain_pio)
  if (soil_sand==unset .and. soil_clay==unset) then
     call gridmap_clean(tgridmap)
     deallocate (kmap, kwgt, kmax, wst)
     deallocate (sand2d_i,clay2d_i,mapunit2d_i)
     deallocate (sand3d_i,clay3d_i,mapunit1d_i,tmp1d_i)
     deallocate (vec_row_indices)
  end if


  write (6,*) 'Successfully made %sand and %clay'
  write (6,*)
  call shr_sys_flush(6)

end subroutine mksoiltex_pio

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: mksoilcolInit
!
! !INTERFACE:
subroutine mksoilcolInit( )
!
! !DESCRIPTION:
! Initialize of make soil color
! !USES:
!
! !ARGUMENTS:
  implicit none
!
! !REVISION HISTORY:
! Author: Erik Kluzek
!
!
! !LOCAL VARIABLES:
!EOP
!-----------------------------------------------------------------------

  ! Error check soil_color if it is set
  if ( soil_color /= unsetcol )then
     if ( soil_color < 0 .or. soil_color > 20 )then
        write(6,*)'soil_color is out of range = ', soil_color
        call abort()
     end if
     write(6,*) 'Replace soil color for all points with: ', soil_color
  end if
end subroutine mksoilcolInit


!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: mksoilcol
!
! !INTERFACE:
subroutine mksoilcol(ldomain, mapfname, datfname, ndiag, &
                    soil_color_o, nsoicol)
!
! !DESCRIPTION:
! make %sand and %clay from IGBP soil data, which includes
! igbp soil 'mapunits' and their corresponding textures
!
! !USES:
  use mkdomainMod, only : domain_type, domain_clean, domain_read
  use mkgridmapMod
  use mkvarpar
  use mkvarctl    
  use mkncdio
!
! !ARGUMENTS:
  implicit none
  type(domain_type), intent(in) :: ldomain
  character(len=*)  , intent(in) :: mapfname           ! input mapping file name
  character(len=*)  , intent(in) :: datfname           ! input data file name
  integer           , intent(in) :: ndiag              ! unit number for diag out
  integer           , intent(out):: soil_color_o(:)    ! soil color classes
  integer           , intent(out):: nsoicol            ! number of soil colors 
!
! !CALLED FROM:
! subroutine mksrfdat in module mksrfdatMod
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
!
! !LOCAL VARIABLES:
!EOP
  type(gridmap_type)    :: tgridmap
  type(domain_type)    :: tdomain          ! local domain
  integer, parameter :: num=2               ! set soil mapunit number
  integer  :: wsti(num)                     ! index to 1st and 2nd largest wst
  real(r8), allocatable :: wst(:,:)         ! overlap weights, by surface type
  real(r8), allocatable :: gast_i(:)        ! global area, by surface type
  real(r8), allocatable :: gast_o(:)        ! global area, by surface type
  integer , allocatable :: soil_color_i(:)  ! input grid: BATS soil color
  integer , allocatable :: color(:)         ! 0: none; 1: some
  real(r8) :: wt                            ! map overlap weight
  real(r8) :: sum_fldi                      ! global sum of dummy input fld
  real(r8) :: sum_fldo                      ! global sum of dummy output fld
  character(len=35), allocatable :: col(:)  ! name of each color
  integer  :: k,n,ni,no,ns_i,ns_o           ! indices
  integer  :: ncid,varid                    ! input netCDF id's
  integer  :: ier                           ! error status
  integer  :: miss = 99999                  ! missing data indicator
  real(r8) :: relerr = 0.00001              ! max error: sum overlap wts ne 1
  character(len=32) :: subname = 'mksoilcol'
!-----------------------------------------------------------------------

  write (6,*) 'Attempting to make soil color classes .....'
  call shr_sys_flush(6)
  write(*,*)'mapfname:' ,trim(mapfname)
  write(*,*)'datfname:' ,trim(datfname)

  ! -----------------------------------------------------------------
  ! Read input file
  ! -----------------------------------------------------------------

  ns_o = ldomain%ns

  ! Obtain input grid info, read local fields

  call domain_read(tdomain,datfname)
  ns_i = tdomain%ns
  allocate(soil_color_i(ns_i), stat=ier)
  if (ier/=0) call abort()

  write (6,*) 'Open soil color file: ', trim(datfname)
  call check_ret(nf_open(datfname, 0, ncid), subname)
  call check_ret(nf_inq_varid (ncid, 'SOIL_COLOR', varid), subname)
  call check_ret(nf_get_var_int (ncid, varid, soil_color_i), subname)
  call check_ret(nf_close(ncid), subname)

  nsoicol = maxval(soil_color_i)
  write(6,*)'nsoicol = ',nsoicol

  allocate(gast_i(0:nsoicol),gast_o(0:nsoicol),col(0:nsoicol))

  ! -----------------------------------------------------------------
  ! Define the model color classes: 0 to nsoicol
  ! -----------------------------------------------------------------

  if (nsoicol == 20) then
     col(0)  = 'no soil                            '
     col(1)  = 'class 1: light                     '
     col(2)  = 'class 2:                           '
     col(3)  = 'class 3:                           '
     col(4)  = 'class 4:                           '
     col(5)  = 'class 5:                           '
     col(6)  = 'class 6:                           '
     col(7)  = 'class 7:                           '
     col(8)  = 'class 8:                           '
     col(9)  = 'class 9:                           '
     col(10) = 'class 10:                          '
     col(11) = 'class 11:                          '
     col(12) = 'class 12:                          '
     col(13) = 'class 13:                          '
     col(14) = 'class 14:                          '
     col(15) = 'class 15:                          '
     col(16) = 'class 16:                          '
     col(17) = 'class 17:                          '
     col(18) = 'class 18:                          '
     col(19) = 'class 19:                          '
     col(20) = 'class 20: dark                     '
  else if (nsoicol == 8) then
     col(0) = 'no soil                            '
     col(1) = 'class 1: light                     '
     col(2) = 'class 2:                           '
     col(3) = 'class 3:                           '
     col(4) = 'class 4:                           '
     col(5) = 'class 5:                           '
     col(6) = 'class 6:                           '
     col(7) = 'class 7:                           '
     col(8) = 'class 8: dark                      '
  else
     write(6,*)'nsoicol value of ',nsoicol,' is not currently supported'
     call abort()
  end if

  ! Error check soil_color if it is set
  if ( soil_color /= unsetcol )then
     if ( soil_color > nsoicol )then
        write(6,*)'soil_color is out of range = ', soil_color
        call abort()
     end if

     do no = 1,ns_o
        soil_color_o(no) = soil_color
     end do

  else

     call gridmap_mapread(tgridmap, mapfname)

     ! Error checks for domain and map consistencies

     call domain_checksame( tdomain, ldomain, tgridmap )

     ! find area of overlap for each soil color for each no

     allocate(wst(0:nsoicol,ns_o))
     wst(0:nsoicol,:) = 0
     allocate(color(ns_o))
     color(:) = 0
     
     ! TODO: need to do a loop to determine
     ! the maximum number of over lap cells throughout the grid 
     ! first get an array that is novr(ns_o) and fill this in - then set
     ! maxovr - to max(novr) - then allocate the array wst to be size of
     ! maxovr,ns_o or 0:nsoilcol,ns_o

     do n = 1,tgridmap%ns
        ni = tgridmap%src_indx(n)
        no = tgridmap%dst_indx(n)
        wt = tgridmap%wovr(n)
        k  = soil_color_i(ni) * tdomain%mask(ni)
        wst(k,no) = wst(k,no) + wt
        if (k>0 .and. wst(k,no)>0.) then
           color(no) = 1
           wst(0,no) = 0.0
        end if
     enddo

     soil_color_o(:) = 0
     do no = 1,ns_o

        ! Rank non-zero weights by color type. wsti(1) is the most extensive
        ! color type. 

        if (color(no) == 1) then
           call mkrank (nsoicol, wst(0:nsoicol,no), miss, wsti, num)
           soil_color_o(no) = wsti(1)
        end if

        ! If land but no color, set color to 15 (in older dataset generic 
        ! soil color 4)
        
        if (nsoicol == 8) then
           if (soil_color_o(no)==0) soil_color_o(no) = 4
        else if (nsoicol == 20) then
           if (soil_color_o(no)==0) soil_color_o(no) = 15
        end if
        
        ! Error checks

        if (soil_color_o(no) < 0 .or. soil_color_o(no) > nsoicol) then
           write (6,*) 'MKSOILCOL error: land model soil color = ', &
                soil_color_o(no),' is not valid for lon,lat = ',no
           call abort()
        end if

     enddo
     deallocate (wst)
     deallocate (color)

     ! Global sum of output field 

     sum_fldi = 0.0_r8
     do ni = 1,ns_i
       sum_fldi = sum_fldi + tgridmap%area_src(ni) * tgridmap%frac_src(ni)
     enddo

     sum_fldo = 0.
     do no = 1,ns_o
        sum_fldo = sum_fldo + tgridmap%area_dst(no) * tgridmap%frac_dst(no)
     end do

     ! -----------------------------------------------------------------
     ! Error check1
     ! Compare global sum fld_o to global sum fld_i.
     ! -----------------------------------------------------------------

     if ( trim(mksrf_gridtype) == 'global') then
        if ( abs(sum_fldo/sum_fldi-1.) > relerr ) then
           write (6,*) 'MKSOILCOL error: input field not conserved'
           write (6,'(a30,e20.10)') 'global sum output field = ',sum_fldo
           write (6,'(a30,e20.10)') 'global sum input  field = ',sum_fldi
           stop
        end if
     end if

     ! -----------------------------------------------------------------
     ! Error check2
     ! Compare global area of each soil color on input and output grids
     ! -----------------------------------------------------------------

     gast_i(:) = 0.
     do ni = 1,ns_i
        k = soil_color_i(ni)
        gast_i(k) = gast_i(k) + tgridmap%area_src(ni)*tgridmap%frac_src(ni)*re**2
     end do

     gast_o(:) = 0.
     do no = 1,ns_o
        k = soil_color_o(no)
        gast_o(k) = gast_o(k) + tgridmap%area_dst(no)*tgridmap%frac_dst(no)*re**2
     end do

     ! area comparison

     write (ndiag,*)
     write (ndiag,'(1x,70a1)') ('=',k=1,70)
     write (ndiag,*) 'Soil Color Output'
     write (ndiag,'(1x,70a1)') ('=',k=1,70)

     write (ndiag,*)
     write (ndiag,'(1x,70a1)') ('.',k=1,70)
     write (ndiag,1001)
1001 format (1x,'soil color type',20x,' input grid area output grid area',/ &
             1x,33x,'     10**6 km**2','      10**6 km**2')
     write (ndiag,'(1x,70a1)') ('.',k=1,70)
     write (ndiag,*)

     do k = 0, nsoicol
        write (ndiag,1002) col(k),gast_i(k)*1.e-6,gast_o(k)*1.e-6
1002    format (1x,a35,f16.3,f17.3)
     end do

  end if

  ! Deallocate dynamic memory

  call domain_clean(tdomain)
  if ( soil_color == unsetcol )then
     call gridmap_clean(tgridmap)
  end if
  deallocate (soil_color_i,gast_i,gast_o,col)

  write (6,*) 'Successfully made soil color classes'
  write (6,*)
  call shr_sys_flush(6)

end subroutine mksoilcol

!-----------------------------------------------------------------------

subroutine mksoilcol_pio(ldomain_pio, mapfname, datfname, ndiag, &
     soil_color_o, nsoicol)
  !
  ! !DESCRIPTION:
  ! make %sand and %clay from IGBP soil data, which includes
  ! igbp soil 'mapunits' and their corresponding textures
  !
  ! !USES:
  use mkdomainPIOMod, only : domain_pio_type
  use mkdataPIOMod
  use mkgridmapMod
  use mkvarpar
  use mkvarctl    
  use mkncdio
  !
  ! !ARGUMENTS:
  implicit none
  type(domain_pio_type) , intent(in) :: ldomain_pio
  character(len=*)      , intent(in) :: mapfname           ! input mapping file name
  character(len=*)      , intent(in) :: datfname           ! input data file name
  integer               , intent(in) :: ndiag              ! unit number for diag out
  integer               , intent(out):: soil_color_o(:)    ! soil color classes
  integer               , intent(out):: nsoicol            ! number of soil colors 
  !
  ! !CALLED FROM:
  ! subroutine mksrfdat in module mksrfdatMod
  !
  ! !REVISION HISTORY:
  ! Author: Mariana Vertenstein
  !
  !
  ! !LOCAL VARIABLES:
  integer :: nodata_value
  !-----------------------------------------------------------------------

  write (6,*) 'Attempting to make soil color classes .....'
  call shr_sys_flush(6)

  ! -----------------------------------------------------------------
  ! Define the model color classes: 0 to nsoicol
  ! -----------------------------------------------------------------
  nsoicol = 20

  if (nsoicol == 8) then
     nodata_value = 4
  else if (nsoicol == 20) then
     nodata_value = 15
  else
     write(*,*)'Unsupported nsoilcol = ', nsoicol
     call exit(0)
  end if

  call mkdata_dominant_int_2d_pio(ldomain_pio, mapfname=mapfname, datfname=datfname, varname='SOIL_COLOR', &
       data_descrip='soil_color', ndiag=ndiag, zero_out=.false., nodata_value=nodata_value, &
       max_value = nsoicol, data_o=soil_color_o)


end subroutine mksoilcol_pio

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: mksoilordInit
!
! !INTERFACE:
subroutine mksoilordInit( )
!
! !DESCRIPTION:
! Initialize of make soil order 
! !USES:
!
! !ARGUMENTS:
  implicit none
!
! !REVISION HISTORY:
! Author: Xiaoying Shi 
!
!
! !LOCAL VARIABLES:
!EOP
!-----------------------------------------------------------------------

  ! Error check soil_order if it is set
  if ( soil_order /= unsetord )then
     if ( soil_order < 0 .or. soil_order > 16 )then
        write(6,*)'soil_order is out of range = ', soil_order
        call abort()
     end if
     write(6,*) 'Replace soil order for all points with: ', soil_order
  end if
end subroutine mksoilordInit
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: mksoilord
!
! !INTERFACE:
subroutine mksoilord(ldomain, mapfname, datfname, ndiag, &
                    pctglac_o, soil_order_o, nsoiord)
!
! !DESCRIPTION:
! make soil order from soil order data
!
! !USES:
  use mkdomainMod, only : domain_type, domain_clean, domain_read
  use mkgridmapMod
  use mkvarpar
  use mkvarctl
  use mkncdio
!
! !ARGUMENTS:
  implicit none
  type(domain_type), intent(in) :: ldomain
    character(len=*)  , intent(in) :: mapfname           ! input mapping file name
  character(len=*)  , intent(in) :: datfname           ! input data file name
  integer           , intent(in) :: ndiag              ! unit number for diagout
  real(r8)          , intent(in) :: pctglac_o(:)       ! % glac (output grid)
  integer           , intent(out):: soil_order_o(:)    ! soil order classes
  integer           , intent(out):: nsoiord            ! number of soil order 

!
! !CALLED FROM:
! subroutine mksrfdat in module mksrfdatMod
!
! !REVISION HISTORY:
! Author: Xiaoying Shi 
!
!
! !LOCAL VARIABLES:
!EOP
  type(gridmap_type)    :: tgridmap
  type(domain_type)    :: tdomain          ! local domain
  integer, parameter :: num=2               ! set soil mapunit number
  integer  :: wsti(num)                     ! index to 1st and 2nd largest wst
  real(r8), allocatable :: wst(:,:)         ! overlap weights, by surface type
  real(r8), allocatable :: gast_i(:)        ! global area, by surface type
  real(r8), allocatable :: gast_o(:)        ! global area, by surface type
  integer , allocatable :: soil_order_i(:)  ! input grid:  soil order 
  integer , allocatable :: order(:)         ! 0: none; 1: some
  real(r8) :: wt                            ! map overlap weight
  real(r8) :: sum_fldi                      ! global sum of dummy input fld
  real(r8) :: sum_fldo                      ! global sum of dummy output fld
  character(len=35), allocatable :: ord(:)  ! name of each ord 
  integer  :: k,n,ni,no,ns_i,ns_o           ! indices
  integer  :: ncid,varid                    ! input netCDF id's
  integer  :: ier                           ! error status
  integer  :: miss = 99999                  ! missing data indicator
  real(r8) :: relerr = 0.00001              ! max error: sum overlap wts ne 1
  character(len=32) :: subname = 'mksoilord'
!-----------------------------------------------------------------------

  write (6,*) 'Attempting to make soil order classes .....'
  call shr_sys_flush(6)
  write(*,*)'mapfname:' ,trim(mapfname)
  write(*,*)'datfname:' ,trim(datfname)

  ! -----------------------------------------------------------------
  ! Read input file
  ! -----------------------------------------------------------------

  ns_o = ldomain%ns
  ! Obtain input grid info, read local fields

  call domain_read(tdomain,datfname)
  ns_i = tdomain%ns
  allocate(soil_order_i(ns_i), stat=ier)
  if (ier/=0) call abort()

  write (6,*) 'Open soil order file: ', trim(datfname)
  call check_ret(nf_open(datfname, 0, ncid), subname)
  call check_ret(nf_inq_varid (ncid, 'SOILORDER', varid), subname)
  call check_ret(nf_get_var_int (ncid, varid, soil_order_i), subname)
  call check_ret(nf_close(ncid), subname)

  nsoiord = maxval(soil_order_i)
  write(6,*)'nsoiord = ',nsoiord

  allocate(gast_i(0:nsoiord),gast_o(0:nsoiord),ord(0:nsoiord))
  ! -----------------------------------------------------------------
  ! Define the model order classes: 0 to nsoiord
  ! -----------------------------------------------------------------

  if (nsoiord == 16) then
     ord(0)  = 'no soil                            '
     ord(1)  = 'class 1: light                     '
     ord(2)  = 'class 2:                           '
     ord(3)  = 'class 3:                           '
     ord(4)  = 'class 4:                           '
     ord(5)  = 'class 5:                           '
     ord(6)  = 'class 6:                           '
     ord(7)  = 'class 7:                           '
     ord(8)  = 'class 8:                           '
     ord(9)  = 'class 9:                           '
     ord(10) = 'class 10:                          '
     ord(11) = 'class 11:                          '
     ord(12) = 'class 12: dark                     '
  else
     write(6,*)'nsoiord value of ',nsoiord,' is not currently supported'
     call abort()
  end if

   ! Error check soil_order if it is set
   if ( soil_order /= unsetord )then
     if ( soil_order > nsoiord )then
        write(6,*)'soil_order is out of range = ', soil_order
        call abort()
     end if
        do no = 1,ns_o
        soil_order_o(no) = soil_order
     end do

  else
      call gridmap_mapread(tgridmap, mapfname)
     ! Error checks for domain and map consistencies

     call domain_checksame( tdomain, ldomain, tgridmap )

     ! find area of overlap for each soil order for each no

     allocate(wst(0:nsoiord,ns_o))
     wst(0:nsoiord,:) = 0
     allocate(order(ns_o))
     order(:) = 0

     ! TODO: need to do a loop to determine
     ! the maximum number of over lap cells throughout the grid 
     ! first get an array that is novr(ns_o) and fill this in - then set
     ! maxovr - to max(novr) - then allocate the array wst to be size of
     ! maxovr,ns_o or 0:nsoilord,ns_o
     
     do n = 1,tgridmap%ns
        ni = tgridmap%src_indx(n)
        no = tgridmap%dst_indx(n)
        wt = tgridmap%wovr(n)
        k  = soil_order_i(ni) * tdomain%mask(ni)
        wst(k,no) = wst(k,no) + wt
        if (k>0 .and. wst(k,no)>0.) then
           order(no) = 1
           wst(0,no) = 0.0
        end if
     enddo

     soil_order_o(:) = 0
     do no = 1,ns_o

        ! Rank non-zero weights by order type. wsti(1) is the most extensive
        ! order type. 

        if (order(no) == 1) then
           call mkrank (nsoiord, wst(0:nsoiord,no), miss, wsti, num)
           soil_order_o(no) = wsti(1)
        end if

        ! If land but no order, set order to  6?

           if (soil_order_o(no)==0) soil_order_o(no) = 15
        ! If ocean, change the soil order to a valid value (the valid value for soil order is 1-15,
        ! but the values of ocean gridcells are 16 in the raw data) 
           if (soil_order_o(no)==16) soil_order_o(no) = 15 
 

        ! Set order for grid cells that are 100% glacier to zero. Otherwise,
        ! must have a soil order for the non-glacier portion of grid cell.

        if (abs(pctglac_o(no)-100.)<1.e-06) soil_order_o(no)= 15
         ! Error checks

        if (soil_order_o(no) < 0 .or. soil_order_o(no) > nsoiord) then
           write (6,*) 'MKSOILCOL error: land model soil order = ', &
                soil_order_o(no),' is not valid for lon,lat = ',no
           call abort()
        end if

     enddo
     deallocate (wst)
     deallocate (order)

     ! Global sum of output field 

     sum_fldi = 0.0_r8
     do ni = 1,ns_i
       sum_fldi = sum_fldi + tgridmap%area_src(ni) * tgridmap%frac_src(ni)
     enddo

     sum_fldo = 0.
     do no = 1,ns_o
        sum_fldo = sum_fldo + tgridmap%area_dst(no) * tgridmap%frac_dst(no)
     end do
     ! -----------------------------------------------------------------
     ! Error check1
     ! Compare global sum fld_o to global sum fld_i.
     ! -----------------------------------------------------------------

     if ( trim(mksrf_gridtype) == 'global') then
        if ( abs(sum_fldo/sum_fldi-1.) > relerr ) then
           write (6,*) 'MKSOILORD error: input field not conserved'
           write (6,'(a30,e20.10)') 'global sum output field = ',sum_fldo
           write (6,'(a30,e20.10)') 'global sum input  field = ',sum_fldi
           stop
        end if
     end if

       ! -----------------------------------------------------------------
     ! Error check2
     ! Compare global area of each soil order on input and output grids
     ! -----------------------------------------------------------------

     gast_i(:) = 0.
     do ni = 1,ns_i
        k = soil_order_i(ni)
        gast_i(k) = gast_i(k) + tgridmap%area_src(ni)*tgridmap%frac_src(ni)*re**2
     end do

     gast_o(:) = 0.
     do no = 1,ns_o
        k = soil_order_o(no)
        gast_o(k) = gast_o(k) + tgridmap%area_dst(no)*tgridmap%frac_dst(no)*re**2
     end do
       ! area comparison

     write (ndiag,*)
     write (ndiag,'(1x,70a1)') ('=',k=1,70)
     write (ndiag,*) 'Soil order Output'
     write (ndiag,'(1x,70a1)') ('=',k=1,70)

     write (ndiag,*)
     write (ndiag,'(1x,70a1)') ('.',k=1,70)
     write (ndiag,1001)
1001 format (1x,'soil order type',20x,' input grid area output grid area',/ &
             1x,33x,'     10**6 km**2','      10**6 km**2')
     write (ndiag,'(1x,70a1)') ('.',k=1,70)
     write (ndiag,*)

     do k = 0, nsoiord
        write (ndiag,1002) ord(k),gast_i(k)*1.e-6,gast_o(k)*1.e-6
1002    format (1x,a35,f16.3,f17.3)
     end do

  end if
        ! Deallocate dynamic memory

  call domain_clean(tdomain)
  if ( soil_order == unsetord )then
     call gridmap_clean(tgridmap)
  end if
  deallocate (soil_order_i,gast_i,gast_o,ord)

  write (6,*) 'Successfully made soil order classes'

  write (6,*)
  call shr_sys_flush(6)

end subroutine mksoilord

!-----------------------------------------------------------------------
subroutine mksoilord_pio(ldomain_pio, mapfname, datfname, ndiag, &
     pctglac_o, soil_order_o, nsoiord)
  !
  use mkdomainPIOMod, only : domain_pio_type
  use mkdataPIOMod
  use mkgridmapMod
  use mkvarpar
  use mkvarctl
  use mkncdio
  !
  ! !ARGUMENTS:
  implicit none
  type(domain_pio_type) , intent(in) :: ldomain_pio
  character(len=*)  , intent(in) :: mapfname           ! input mapping file name
  character(len=*)  , intent(in) :: datfname           ! input data file name
  integer           , intent(in) :: ndiag              ! unit number for diagout
  real(r8)          , intent(in) :: pctglac_o(:)       ! % glac (output grid)
  integer           , intent(out):: soil_order_o(:)    ! soil order classes
  integer           , intent(out):: nsoiord            ! number of soil order 
  ! !LOCAL VARIABLES:
  integer, parameter :: nodata_value = 0

  write (6,*) 'Attempting to make soil order classes .....'
  call shr_sys_flush(6)

  nsoiord = 16
  call mkdata_dominant_int_2d_pio(ldomain_pio, mapfname=mapfname, datfname=datfname, varname='SOILORDER', &
       data_descrip='soil_color', ndiag=ndiag, zero_out=.false., nodata_value=nodata_value, &
       max_value = nsoiord, data_o=soil_order_o)

  write (6,*) 'Successfully made soil order classes'

  write (6,*)
  call shr_sys_flush(6)

end subroutine mksoilord_pio
!-----------------------------------------------------------------------

!BOP
!
! !IROUTINE: mkorganic
!
! !INTERFACE:
subroutine mkorganic(ldomain, mapfname, datfname, ndiag, organic_o)
!
! !DESCRIPTION:
! make organic matter dataset
!
! !USES:
  use mkdomainMod, only : domain_type, domain_clean, domain_read
  use mkgridmapMod
  use mkvarpar
  use mkvarctl    
  use mkncdio
!
! !ARGUMENTS:
  implicit none
  type(domain_type), intent(in) :: ldomain
  character(len=*)  , intent(in) :: mapfname       ! input mapping file name
  character(len=*)  , intent(in) :: datfname       ! input data file name
  integer           , intent(in) :: ndiag          ! unit number for diag out
  real(r8)          , intent(out):: organic_o(:,:) ! output grid:
!
! !CALLED FROM:
! subroutine mksrfdat in module mksrfdatMod
!
! !REVISION HISTORY:
! 
! Author: David Lawrence
!
!
! !LOCAL VARIABLES:
!EOP
  type(gridmap_type)    :: tgridmap
  type(domain_type)    :: tdomain         ! local domain
  real(r8), allocatable :: organic_i(:,:)  ! input grid: total column organic matter
  integer  :: no,ns_i                      ! indices
  integer  :: lev                          ! level index
  integer  :: nlay                         ! number of soil layers
  integer  :: ncid,dimid,varid             ! input netCDF id's
  integer  :: ier                          ! error status
  character(len=32) :: subname = 'mkorganic'
!-----------------------------------------------------------------------

  write (6,*) 'Attempting to make organic matter dataset .....'
  call shr_sys_flush(6)

  ! -----------------------------------------------------------------
  ! Read input file
  ! -----------------------------------------------------------------

  ! Obtain input grid info, read local fields

  call domain_read(tdomain,datfname)
  ns_i = tdomain%ns

  write (6,*) 'Open soil organic file: ', trim(datfname)
  call check_ret(nf_open(datfname, 0, ncid), subname)

  call check_ret(nf_inq_dimid  (ncid, 'number_of_layers', dimid), subname)
  call check_ret(nf_inq_dimlen (ncid, dimid, nlay), subname)

  allocate(organic_i(ns_i,nlay),stat=ier)
  if (ier/=0) call abort()
  if (nlay /= nlevsoi) then
     write(6,*)'nlay, nlevsoi= ',nlay,nlevsoi,' do not match'
     stop
  end if

  call check_ret(nf_inq_varid (ncid, 'ORGANIC', varid), subname)
  call check_ret(nf_get_var_double (ncid, varid, organic_i), subname)

  call check_ret(nf_close(ncid), subname)

  ! Area-average percent cover on input grid to output grid 
  ! and correct according to land landmask
  ! Note that percent cover is in terms of total grid area.

  call gridmap_mapread(tgridmap, mapfname )

  call domain_checksame( tdomain, ldomain, tgridmap )

  do lev = 1,nlay
     call gridmap_areaave(tgridmap, organic_i(:,lev), organic_o(:,lev), nodata=0._r8)
  end do

  do lev = 1,nlevsoi

     ! Check for conservation

     do no = 1,ldomain%ns
        if ((organic_o(no,lev)) > 130.000001_r8) then
           write (6,*) 'MKORGANIC error: organic = ',organic_o(no,lev), &
                ' greater than 130.000001 for column, row = ',no
           call shr_sys_flush(6)
           stop
        end if
     enddo

!    ! Diagnostic output

     ! TODO: there is nothing being written out here currently - all zeroes
     ! So for now these are commented out
!!$     write (ndiag,*)
!!$     write (ndiag,'(1x,70a1)') ('.',k=1,70)
!!$     write (ndiag,2001)
!!$2001 format (1x,'surface type   input grid area  output grid area'/ &
!!$             1x,'                 10**6 km**2      10**6 km**2   ')
!!$     write (ndiag,'(1x,70a1)') ('.',k=1,70)
!!$     write (ndiag,*)
!!$     write (ndiag,2002) gomlev_i*1.e-06,gomlev_o*1.e-06
!!$     write (ndiag,2004) garea_i*1.e-06,garea_o*1.e-06
!!$2002 format (1x,'organic    ',f14.3,f17.3)
!!$2004 format (1x,'all surface ',f14.3,f17.3)
!!$
     call shr_sys_flush(ndiag)

     write (6,*) 'Successfully made organic matter, level = ', lev
     call shr_sys_flush(6)

  end do   ! lev

  ! Deallocate dynamic memory

  call domain_clean(tdomain)
  call gridmap_clean(tgridmap)
  deallocate (organic_i)

  write (6,*) 'Successfully made organic matter'
  call shr_sys_flush(6)
  write(6,*)

end subroutine mkorganic

!-----------------------------------------------------------------------
subroutine mkorganic_pio(ldomain_pio, mapfname, datfname, ndiag, organic_o)
!
! !DESCRIPTION:
! make organic matter dataset
!
! !USES:
  use mkdomainPIOMod, only : domain_pio_type
  use mkdataPIOMod
!
! !ARGUMENTS:
  implicit none
  type(domain_pio_type), intent(in) :: ldomain_pio
  character(len=*)  , intent(in) :: mapfname       ! input mapping file name
  character(len=*)  , intent(in) :: datfname       ! input data file name
  integer           , intent(in) :: ndiag          ! unit number for diag out
  real(r8)          , intent(out):: organic_o(:,:) ! output grid:
  !
  real(r8), parameter :: nodata_value = 0._r8
  real(r8), parameter :: max_valid    = 130.000001_r8

  write (6,*) 'Attempting to make organic matter dataset .....'
  call shr_sys_flush(6)

  call mkdata_double_3d_pio(ldomain_pio, mapfname=mapfname, datfname=datfname, varname='ORGANIC', &
       data_descrip='%organic', ndiag=ndiag, zero_out=.false., nodata_value=nodata_value, start_id_for_dim3=1, &
       data_o=organic_o, max_valid_value=max_valid)

  write (6,*) 'Successfully made organic matter'
  call shr_sys_flush(6)
  write(6,*)

end subroutine mkorganic_pio

!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: mkrank
!
! !INTERFACE:
subroutine mkrank (n, a, miss, iv, num)
!
! !DESCRIPTION:
! Return indices of largest [num] values in array [a]. Private method
! only used for soil color and soil texture.
!
! !USES:
!
! !ARGUMENTS:
  implicit none
  integer , intent(in) :: n        !array length
  real(r8), intent(in) :: a(0:n)   !array to be ranked
  integer , intent(in) :: miss     !missing data value
  integer , intent(in) :: num      !number of largest values requested
  integer , intent(out):: iv(num)  !index to [num] largest values in array [a]
!
! !CALLED FROM:
! subroutine mksoilcol 
! subroutine mksoiltex
!
! !REVISION HISTORY:
! Author: Gordon Bonan
!
!
! !LOCAL VARIABLES:
!EOP
  real(r8) a_max       !maximum value in array
  integer i            !array index
  real(r8) delmax      !tolerance for finding if larger value
  integer m            !do loop index
  integer k            !do loop index
  logical exclude      !true if data value has already been chosen
!-----------------------------------------------------------------------

  delmax = 1.e-06

  ! Find index of largest non-zero number

  iv(1) = miss
  a_max = -9999.

  do i = 0, n
     if (a(i)>0. .and. (a(i)-a_max)>delmax) then
        a_max = a(i)
        iv(1)  = i
     end if
  end do

  ! iv(1) = miss indicates no values > 0. this is an error

  if (iv(1) == miss) then
     write (6,*) 'MKRANK error: iv(1) = missing'
     call abort()
  end if

  ! Find indices of the next [num]-1 largest non-zero number.
  ! iv(m) = miss if there are no more values > 0

  do m = 2, num
     iv(m) = miss
     a_max = -9999.
     do i = 0, n

        ! exclude if data value has already been chosen

        exclude = .false.
        do k = 1, m-1
           if (i == iv(k)) exclude = .true.
        end do

        ! if not already chosen, see if it is the largest of
        ! the remaining values

        if (.not. exclude) then
           if (a(i)>0. .and. (a(i)-a_max)>delmax) then
              a_max = a(i)
              iv(m)  = i
           end if
        end if
     end do
  end do

end subroutine mkrank

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: mksoilfmaxInit
!
! !INTERFACE:
subroutine mksoilfmaxInit( )
!
! !DESCRIPTION:
! Initialize of make soil fmax
! !USES:
!
! !ARGUMENTS:
  implicit none
!
! !REVISION HISTORY:
! Author: Erik Kluzek
!
!
! !LOCAL VARIABLES:
!EOP
!-----------------------------------------------------------------------

  ! Error check soil_fmax if it is set
  if ( soil_fmax /= unset )then
     if ( soil_fmax < 0.0 .or. soil_fmax > 1.0 )then
        write(6,*)'soil_fmax is out of range = ', soil_fmax
        stop
     end if
     write(6,*) 'Replace soil fmax for all points with: ', soil_fmax
  end if

end subroutine mksoilfmaxInit

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: mkfmax
!
! !INTERFACE:
subroutine mkfmax(ldomain, mapfname, datfname, ndiag, fmax_o)
!
! !DESCRIPTION:
! make percent fmax
!
! !USES:
  use mkdomainMod, only : domain_type, domain_clean, domain_read
  use mkgridmapMod
  use mkvarpar
  use mkvarctl
  use mkncdio
!
! !ARGUMENTS:
  implicit none
  type(domain_type), intent(in) :: ldomain
  character(len=*)  , intent(in) :: mapfname  ! input mapping file name
  character(len=*)  , intent(in) :: datfname  ! input data file name
  integer           , intent(in) :: ndiag     ! unit number for diag out
  real(r8)          , intent(out):: fmax_o(:) ! output grid: %fmax
!
! !CALLED FROM:
! subroutine mksrfdat in module mksrfdatMod
!
! !REVISION HISTORY:
! Revised: Nan Rosenbloom - used mkglacier.F90 as template.
! Original Author: Mariana Vertenstein
!
!
! !LOCAL VARIABLES:
!EOP
  type(gridmap_type)    :: tgridmap
  type(domain_type)    :: tdomain         ! local domain
  real(r8), allocatable :: fmax_i(:)       ! input grid: percent fmax
  real(r8) :: sum_fldi                     ! global sum of dummy input fld
  real(r8) :: sum_fldo                     ! global sum of dummy output fld
  real(r8) :: gfmax_i                      ! input  grid: global fmax
  real(r8) :: garea_i                      ! input  grid: global area
  real(r8) :: gfmax_o                      ! output grid: global fmax
  real(r8) :: garea_o                      ! output grid: global area
  integer  :: k,ni,no,ns_i,ns_o            ! indices
  integer  :: ncid,varid                   ! input netCDF id's
  integer  :: ier                          ! error status
  real(r8) :: relerr = 0.00001             ! max error: sum overlap wts ne 1
  character(len=32) :: subname = 'mkfmax'
!-----------------------------------------------------------------------

  write (6,*) 'Attempting to make %fmax .....'
  call shr_sys_flush(6)
  write(*,*)'mapfname:' ,trim(mapfname)
  write(*,*)'datfname:' ,trim(datfname)

  ! -----------------------------------------------------------------
  ! Read input file
  ! -----------------------------------------------------------------

  ! Obtain input grid info, read local fields

  call domain_read(tdomain,datfname)
  ns_i = tdomain%ns
  allocate(fmax_i(ns_i), stat=ier)
  if (ier/=0) call abort()
  ns_o = ldomain%ns

  write (6,*) 'Open soil fmax file: ', trim(datfname)
  call check_ret(nf_open(datfname, 0, ncid), subname)
  call check_ret(nf_inq_varid (ncid, 'FMAX', varid), subname)
  call check_ret(nf_get_var_double (ncid, varid, fmax_i), subname)
  call check_ret(nf_close(ncid), subname)

  ! Area-average percent cover on input grid to output grid 
  ! and correct according to land landmask
  ! Note that percent cover is in terms of total grid area.

  call gridmap_mapread(tgridmap, mapfname )

  ! Error checks for domain and map consistencies

  call domain_checksame( tdomain, ldomain, tgridmap )

  ! Determine fmax_o on output grid
 
  ! In points with no data, use globalAvg
  ! (WJS (3-11-13): use real(.365783,r8) rather than .365783_r8 to maintain bfb results
  ! with old code)
  call gridmap_areaave(tgridmap, fmax_i, fmax_o, nodata=real(.365783,r8))

  ! Check for conservation

  do no = 1, ns_o
     if ((fmax_o(no)) > 1.000001_r8) then
        write (6,*) 'MKFMAX error: fmax = ',fmax_o(no), &
                ' greater than 1.000001 for column, row = ',no
        call shr_sys_flush(6)
        stop
     end if
  enddo

  ! Global sum of output field -- must multiply by fraction of
  ! output grid that is land as determined by input grid

  sum_fldi = 0.0_r8
  do ni = 1,ns_i
    sum_fldi = sum_fldi + tgridmap%area_src(ni) * tgridmap%frac_src(ni)
  enddo

  sum_fldo = 0.
  do no = 1,ns_o
     sum_fldo = sum_fldo + tgridmap%area_dst(no) * tgridmap%frac_dst(no)
  end do

  ! -----------------------------------------------------------------
  ! Error check1
  ! Compare global sum fld_o to global sum fld_i.
  ! -----------------------------------------------------------------

  if ( trim(mksrf_gridtype) == 'global') then
     if ( abs(sum_fldo/sum_fldi-1.) > relerr ) then
        write (6,*) 'MKFMAX error: input field not conserved'
        write (6,'(a30,e20.10)') 'global sum output field = ',sum_fldo
        write (6,'(a30,e20.10)') 'global sum input  field = ',sum_fldi
        stop
     end if
  end if

  ! -----------------------------------------------------------------
  ! Error check2
  ! Compare global areas on input and output grids
  ! -----------------------------------------------------------------

  gfmax_i = 0.
  garea_i = 0.
  do ni = 1,ns_i
     garea_i = garea_i + tgridmap%area_src(ni)*re**2
     gfmax_i = gfmax_i + fmax_i(ni)*(tgridmap%area_src(ni)/100.)* &
                                     tgridmap%frac_src(ni)*re**2
  end do

  gfmax_o = 0.
  garea_o = 0.
  do no = 1,ns_o
     garea_o = garea_o + tgridmap%area_dst(no)*re**2
     gfmax_o = gfmax_o + fmax_o(no)*(tgridmap%area_dst(no)/100.) * &
                                     tgridmap%frac_dst(no)*re**2
     if ((tgridmap%mask_dst(no) > 0)) then 
        if ((tgridmap%frac_dst(no) < 0.0) .or. (tgridmap%frac_dst(no) > 1.0001)) then
           write(6,*) "ERROR:: frac out of range: ", tgridmap%frac_dst(no),no
           stop 
        end if
     end if
  end do

  ! Diagnostic output

  write (ndiag,*)
  write (ndiag,'(1x,70a1)') ('=',k=1,70)
  write (ndiag,*) 'Maximum Fractional Saturated Area Output'
  write (ndiag,'(1x,70a1)') ('=',k=1,70)

  write (ndiag,*)
  write (ndiag,'(1x,70a1)') ('.',k=1,70)
  write (ndiag,2001)
2001 format (1x,'surface type   input grid area  output grid area'/ &
             1x,'                 10**6 km**2      10**6 km**2   ')
  write (ndiag,'(1x,70a1)') ('.',k=1,70)
  write (ndiag,*)
  write (ndiag,2002) gfmax_i*1.e-06,gfmax_o*1.e-06
  write (ndiag,2004) garea_i*1.e-06,garea_o*1.e-06
2002 format (1x,'fmax    ',f14.3,f17.3)
2004 format (1x,'all surface ',f14.3,f17.3)

  write (6,*) 'Successfully made %fmax'
  write (6,*)
  call shr_sys_flush(6)

  ! Deallocate dynamic memory

  call domain_clean(tdomain)
  call gridmap_clean(tgridmap)
  deallocate (fmax_i)

end subroutine mkfmax

!-----------------------------------------------------------------------
subroutine mkfmax_pio(ldomain_pio, mapfname, datfname, ndiag, fmax_o)
  !
  ! !DESCRIPTION:
  ! make percent fmax
  !
  ! !USES:
  use mkdomainPIOMod, only : domain_pio_type
  use mkdataPIOMod
  !
  ! !ARGUMENTS:
  implicit none
  type(domain_pio_type) , intent(in) :: ldomain_pio
  character(len=*)      , intent(in) :: mapfname  ! input mapping file name
  character(len=*)      , intent(in) :: datfname  ! input data file name
  integer               , intent(in) :: ndiag     ! unit number for diag out
  real(r8)              , intent(out):: fmax_o(:) ! output grid: %fmax

  real(r8), parameter :: nodata_value = 0.365783_r8
  real(r8), parameter :: min_valid    = 0._r8
  real(r8), parameter :: max_valid    = 1.000001_r8
  
  write (6,*) 'Attempting to make %fmax .....'
  call shr_sys_flush(6)

  call mkdata_double_2d_pio(ldomain_pio, mapfname=mapfname, datfname=datfname, varname='FMAX', &
       data_descrip='%fmax', ndiag=ndiag, zero_out=.false., nodata_value=nodata_value, data_o=fmax_o, &
       threshold_o=0._r8)

  write (6,*) 'Successfully made %fmax'
  write (6,*)
  call shr_sys_flush(6)

end subroutine mkfmax_pio

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: mksoilAtt
!
! !INTERFACE:
subroutine mksoilAtt( ncid, dynlanduse, xtype )
!
! !DESCRIPTION:
! add atttributes to output file regarding the soil module
!
! !USES:
  use fileutils  , only : get_filename
  use mkncdio    , only : check_ret, ncd_defvar
  use mkvarpar   
  use mkvarctl   

! !ARGUMENTS:
  implicit none
  include 'netcdf.inc'
  integer, intent(in) :: ncid         ! NetCDF file ID to write out to
  logical, intent(in) :: dynlanduse   ! if dynamic land-use file
  integer, intent(in) :: xtype        ! external type to output real data as
!
! !CALLED FROM:
! subroutine mkfile in module mkfileMod
!
! !REVISION HISTORY:
! Original Author: Erik Kluzek
!
!
! !LOCAL VARIABLES:
!EOP
  integer :: dimid                ! temporary
  character(len=512) :: str       ! global attribute string
  character(len=32) :: subname = 'mksoilAtt'
!-----------------------------------------------------------------------

  if (.not. dynlanduse) then

     ! Define dimensions unique to soil

     call check_ret(nf_def_dim (ncid, 'nlevsoi',  &
                                       nlevsoi    , dimid), subname)

     ! Add global attributes to file

     if ( soil_clay /= unset .and. soil_sand /= unset )then
        str = 'TRUE'
        call check_ret(nf_put_att_text (ncid, NF_GLOBAL, &
             'soil_clay_override', len_trim(str), trim(str)), subname)
        str = 'TRUE'
        call check_ret(nf_put_att_text (ncid, NF_GLOBAL, &
             'soil_sand_override', len_trim(str), trim(str)), subname)
     else
        str = get_filename(mksrf_fsoitex)
        call check_ret(nf_put_att_text(ncid, NF_GLOBAL, &
             'Soil_texture_raw_data_file_name', len_trim(str), trim(str)), subname)
     end if
     if ( soil_color /= unsetcol )then
        str = 'TRUE'
        call check_ret(nf_put_att_text (ncid, NF_GLOBAL, &
             'soil_color_override', len_trim(str), trim(str)), subname)
     else
        str = get_filename(mksrf_fsoicol)
        call check_ret(nf_put_att_text(ncid, NF_GLOBAL, &
             'Soil_color_raw_data_file_name', len_trim(str), trim(str)), subname)
     end if
     if ( soil_order /= unsetord )then
        str = 'TRUE'
        call check_ret(nf_put_att_text (ncid, NF_GLOBAL, &
             'soil_order_override', len_trim(str), trim(str)), subname)
     else
        str = get_filename(mksrf_fsoiord)
        call check_ret(nf_put_att_text(ncid, NF_GLOBAL, &
             'Soil_order_raw_data_file_name', len_trim(str), trim(str)),subname)
     end if

     if ( soil_fmax /= unset )then
        str = 'TRUE'
        call check_ret(nf_put_att_text (ncid, NF_GLOBAL, &
             'soil_fmax_override', len_trim(str), trim(str)), subname)
     else
        str = get_filename(mksrf_fmax)
        call check_ret(nf_put_att_text(ncid, NF_GLOBAL, &
             'Fmax_raw_data_file_name', len_trim(str), trim(str)), subname)
     end if
     str = get_filename(mksrf_forganic)
     call check_ret(nf_put_att_text(ncid, NF_GLOBAL, &
          'Organic_matter_raw_data_file_name', len_trim(str), trim(str)), subname)
     
     ! Define variables

     call ncd_defvar(ncid=ncid, varname='mxsoil_color', xtype=nf_int, &
          long_name='maximum numbers of soil colors', units='unitless')
     call ncd_defvar(ncid=ncid, varname='mxsoil_order', xtype=nf_int, &
          long_name='maximum numbers of soil order', units='unitless')
     
     if (outnc_1d) then
        call ncd_defvar(ncid=ncid, varname='SOIL_COLOR', xtype=nf_int, &
             dim1name='gridcell',&
             long_name='soil color', units='unitless')
     else
        call ncd_defvar(ncid=ncid, varname='SOIL_COLOR', xtype=nf_int, &
             dim1name='lsmlon', dim2name='lsmlat', &
             long_name='soil color', units='unitless')
     end if
     
     if (outnc_1d) then
        call ncd_defvar(ncid=ncid, varname='SOIL_ORDER', xtype=nf_int, &
             dim1name='gridcell',&
             long_name='soil order', units='unitless')
     else
        call ncd_defvar(ncid=ncid, varname='SOIL_ORDER', xtype=nf_int, &
             dim1name='lsmlon', dim2name='lsmlat', &
             long_name='soil order', units='unitless')
     end if

     if (outnc_1d) then
        call ncd_defvar(ncid=ncid, varname='PCT_SAND', xtype=xtype, &
             dim1name='gridcell', dim2name='nlevsoi', &
             long_name='percent sand', units='unitless')
     else
        call ncd_defvar(ncid=ncid, varname='PCT_SAND', xtype=xtype, &
             dim1name='lsmlon', dim2name='lsmlat', dim3name='nlevsoi', &
             long_name='percent sand', units='unitless')
     end if

     if (outnc_1d) then
        call ncd_defvar(ncid=ncid, varname='PCT_CLAY', xtype=xtype, &
             dim1name='gridcell', dim2name='nlevsoi', &
             long_name='percent clay', units='unitless')
     else
        call ncd_defvar(ncid=ncid, varname='PCT_CLAY', xtype=xtype, &
             dim1name='lsmlon', dim2name='lsmlat', dim3name='nlevsoi', &
             long_name='percent clay', units='unitless')
     end if

     if (outnc_1d) then
        call ncd_defvar(ncid=ncid, varname='ORGANIC', xtype=xtype, &
             dim1name='gridcell', dim2name='nlevsoi', &
             long_name='organic matter density at soil levels', &
             units='kg/m3 (assumed carbon content 0.58 gC per gOM)')
     else
        call ncd_defvar(ncid=ncid, varname='ORGANIC', xtype=xtype, &
             dim1name='lsmlon', dim2name='lsmlat', dim3name='nlevsoi', &
             long_name='organic matter density at soil levels', &
             units='kg/m3 (assumed carbon content 0.58 gC per gOM)')
     end if

     if (outnc_1d) then
        call ncd_defvar(ncid=ncid, varname='FMAX', xtype=xtype, &
             dim1name='gridcell', &
             long_name='maximum fractional saturated area', units='unitless')
     else
        call ncd_defvar(ncid=ncid, varname='FMAX', xtype=xtype, &
             dim1name='lsmlon', dim2name='lsmlat', &
             long_name='maximum fractional saturated area', units='unitless')
     end if

  end if

end subroutine mksoilAtt

!-----------------------------------------------------------------------
subroutine mksoilAttPIO( ncid, dynlanduse, xtype, dim_id_gridcell, dim_id_lsmlon, dim_id_lsmlat, dim_id_nlevsoi)

  use fileutils  , only : get_filename
  use pio
  use piofileutils
  use mkvarpar
  use mkvarctl

  implicit none

  type(file_desc_t) , intent(in)    :: ncid
  logical           , intent(in)    :: dynlanduse
  integer           , intent(in)    :: xtype          ! external type to output real data as
  integer           , intent(in)    :: dim_id_gridcell
  integer           , intent(in)    :: dim_id_lsmlon
  integer           , intent(in)    :: dim_id_lsmlat
  integer           , intent(inout) :: dim_id_nlevsoi

  character(len=512)                :: str            ! global attribute string
  character(len=512)                :: att_name, att_value
  character(len=32)                 :: subname = 'mksoilAttPIO'
  type(var_desc_t)                  :: pioVar
  integer                           :: dim1d(1), dim2d(2), dim3d(3)

  if (.not.dynlanduse) then

     ! Define dimensions unique to soil

     call DefineDimPIO(ncid, 'nlevsoi' , nlevsoi, dim_id_nlevsoi)

     ! Add global attributes to file

     if ( soil_clay /= unset .and. soil_sand /= unset )then
        str = 'TRUE'
        call check_ret(PIO_put_att (ncid, PIO_GLOBAL, 'soil_clay_override', trim(str)), subname)
        str = 'TRUE'
        call check_ret(PIO_put_att (ncid, PIO_GLOBAL, 'soil_sand_override', trim(str)), subname)
     else
        str = get_filename(mksrf_fsoitex)
        call check_ret(PIO_put_att(ncid, PIO_GLOBAL, 'Soil_texture_raw_data_file_name', trim(str)), subname)
     end if
     if ( soil_color /= unsetcol )then
        str = 'TRUE'
        call check_ret(PIO_put_att (ncid, PIO_GLOBAL, 'soil_color_override', trim(str)), subname)
     else
        str = get_filename(mksrf_fsoicol)
        call check_ret(PIO_put_att(ncid, PIO_GLOBAL, 'Soil_color_raw_data_file_name', trim(str)), subname)
     end if

     if ( soil_order /= unsetord )then
        str = 'TRUE'
        call check_ret(PIO_put_att (ncid, PIO_GLOBAL, 'soil_order_override', trim(str)), subname)
     else
        str = get_filename(mksrf_fsoiord)
        call check_ret(PIO_put_att(ncid, PIO_GLOBAL, 'Soil_order_raw_data_file_name', trim(str)),subname)
     end if

     if ( soil_fmax /= unset )then
        str = 'TRUE'
        call check_ret(PIO_put_att (ncid, PIO_GLOBAL, 'soil_fmax_override', trim(str)), subname)
     else
        str = get_filename(mksrf_fmax)
        call check_ret(PIO_put_att(ncid, PIO_GLOBAL, 'Fmax_raw_data_file_name', trim(str)), subname)
     end if
     str = get_filename(mksrf_forganic)
     call check_ret(PIO_put_att(ncid, PIO_GLOBAL, 'Organic_matter_raw_data_file_name', trim(str)), subname)

     ! Define variables

     call check_ret(PIO_def_var(ncid, 'mxsoil_color', PIO_INT, pioVar), subname)
     att_name = 'long_name'; att_value = 'maximum numbers of soil colors';
     call check_ret(PIO_put_att(ncid, pioVar, att_name, trim(att_value)), subname)
     att_name = 'units'; att_value = 'unitless';
     call check_ret(PIO_put_att(ncid, pioVar, att_name, trim(att_value)), subname)

     call check_ret(PIO_def_var(ncid, 'mxsoil_order', PIO_INT, pioVar), subname)
     att_name = 'long_name'; att_value = 'maximum numbers of soil orders';
     call check_ret(PIO_put_att(ncid, pioVar, att_name, trim(att_value)), subname)
     att_name = 'units'; att_value = 'unitless';
     call check_ret(PIO_put_att(ncid, pioVar, att_name, trim(att_value)), subname)

     if (outnc_1d) then
        dim1d(1) = dim_id_gridcell
        call DefineVarPIO_1d(ncid, 'SOIL_COLOR', PIO_INT, dim1d, longName='soil color', units='unitless')
        call DefineVarPIO_1d(ncid, 'SOIL_ORDER', PIO_INT, dim1d, longName='soil order', units='unitless')

        dim2d(1) = dim_id_gridcell; dim2d(2) = dim_id_nlevsoi
        call DefineVarPIO_2d(ncid, 'PCT_CLAY', xtype, dim2d, longName='percent sand', units='unitless')
        call DefineVarPIO_2d(ncid, 'PCT_SAND', xtype, dim2d, longName='percent clay', units='unitless')
        call DefineVarPIO_2d(ncid, 'ORGANIC', xtype, dim2d, longName='organic matter density at soil levels', &
             units='kg/m3 (assumed carbon content 0.58 gC per gOM)')

        dim1d(1) = dim_id_gridcell
        call DefineVarPIO_1d(ncid, 'FMAX', xtype, dim1d, longName='maximum fractional saturated area', units='unitless')

     else
        dim2d(1) = dim_id_lsmlon; dim2d(2) = dim_id_lsmlat;
        call DefineVarPIO_2d(ncid, 'SOIL_COLOR', PIO_INT, dim2d, longName='soil color', units='unitless')
        call DefineVarPIO_2d(ncid, 'SOIL_ORDER', PIO_INT, dim2d, longName='soil order', units='unitless')

        dim3d(1) = dim_id_lsmlon; dim3d(2) = dim_id_lsmlat; dim3d(3) = dim_id_nlevsoi
        call DefineVarPIO_3d(ncid, 'PCT_CLAY', xtype, dim3d, longName='percent sand', units='unitless')
        call DefineVarPIO_3d(ncid, 'PCT_SAND', xtype, dim3d, longName='percent clay', units='unitless')
        call DefineVarPIO_3d(ncid, 'ORGANIC', xtype, dim3d, longName='organic matter density at soil levels', &
             units='kg/m3 (assumed carbon content 0.58 gC per gOM)')

        dim2d(1) = dim_id_lsmlon; dim2d(2) = dim_id_lsmlat;
        call DefineVarPIO_2d(ncid, 'FMAX', xtype, dim2d, longName='maximum fractional saturated area', units='unitless')
     end if

  end if

end subroutine mksoilAttPIO

end module mksoilMod
