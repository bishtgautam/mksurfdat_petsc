module mklaiMod
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: mklai
!
! !DESCRIPTION:
! Make LAI/SAI/height data
!
! !REVISION HISTORY:
! Author: Sam Levis
!
!EOP
!-----------------------------------------------------------------------
  use shr_kind_mod, only : r8 => shr_kind_r8
  use shr_sys_mod , only : shr_sys_flush
  use mkdomainMod , only : domain_checksame
  use spmdMod     , only : masterproc
  use mkvarctl    

  implicit none

  private

  public  :: mklai
  public  :: mklai_pio
  private :: pft_laicheck

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: mklai
!
! !INTERFACE:
subroutine mklai(ldomain, mapfname, datfname, ndiag, ncido)
!
! !DESCRIPTION:
! Make LAI/SAI/height data
! Portions of this code could be moved out of the month loop
! for improved efficiency
!
! !USES:
  use mkdomainMod, only : domain_type, domain_clean, domain_read
  use mkgridmapMod
  use mkvarpar    , only : re
  use mkvarctl    
  use mkncdio
  use mkpftMod    , only : c3cropindex, c3irrcropindex
!
! !ARGUMENTS:
  implicit none
  type(domain_type), intent(in) :: ldomain
  character(len=*)  , intent(in) :: mapfname     ! input mapping file name
  character(len=*)  , intent(in) :: datfname     ! input data file name
  integer           , intent(in) :: ndiag        ! unit number for diag out
  integer           , intent(in) :: ncido        ! output netcdf file id
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
  integer  :: numpft_i                      ! number of plant types on input
  real(r8) :: glai_o(0:numpft)              ! output grid: global area pfts
  real(r8) :: gsai_o(0:numpft)              ! output grid: global area pfts
  real(r8) :: ghgtt_o(0:numpft)             ! output grid: global area pfts
  real(r8) :: ghgtb_o(0:numpft)             ! output grid: global area pfts
  real(r8) :: glai_i(0:numpft)              ! input grid: global area pfts
  real(r8) :: gsai_i(0:numpft)              ! input grid: global area pfts
  real(r8) :: ghgtt_i(0:numpft)             ! input grid: global area pfts
  real(r8) :: ghgtb_i(0:numpft)             ! input grid: global area pfts

  real(r8), allocatable :: mlai_o(:,:)      ! monthly lai
  real(r8), allocatable :: msai_o(:,:)      ! monthly sai
  real(r8), allocatable :: mhgtt_o(:,:)     ! monthly height (top)
  real(r8), allocatable :: mhgtb_o(:,:)     ! monthly height (bottom)
  real(r8), allocatable :: mlai_i(:,:)      ! monthly lai in
  real(r8), allocatable :: msai_i(:,:)      ! monthly sai in
  real(r8), allocatable :: mhgtt_i(:,:)     ! monthly height (top) in
  real(r8), allocatable :: mhgtb_i(:,:)     ! monthly height (bottom) in
  real(r8), allocatable :: mask_src(:)      ! input grid: mask (0, 1)
  integer,  pointer     :: laimask(:,:)     ! lai+sai output mask for each plant function type
  integer,  allocatable :: pft_index(:)     ! PFT index for each grid cell
  real(r8), allocatable :: mlai_i_tmp(:)      ! monthly lai in
  real(r8), allocatable :: msai_i_tmp(:)      ! monthly sai in
  real(r8), allocatable :: mhgtt_i_tmp(:)     ! monthly height (top) in
  real(r8), allocatable :: mhgtb_i_tmp(:)     ! monthly height (bottom) in
  real(r8) :: garea_i                       ! input  grid: global area
  real(r8) :: garea_o                       ! output grid: global area
  integer  :: ni,no,ns_i,ns_o               ! indices
  integer  :: i,k,l,m                         ! indices
  integer  :: ncidi,dimid,varid             ! input netCDF id's
  integer  :: ndimsi,ndimso                 ! netCDF dimension sizes 
  integer  :: dimids(4)                     ! netCDF dimension ids
  integer  :: bego(4),leno(4)               ! netCDF bounds
  integer  :: begi(4),leni(4)               ! netCDF bounds 
  integer  :: ntim                          ! number of input time samples
  logical  :: use_pft_index                 ! whether to use PFT index
  integer  :: ier                           ! error status
  character(len= 32) :: subname = 'mklai'
!-----------------------------------------------------------------------

  if (masterproc) write (6,*) 'Attempting to make LAIs/SAIs/heights .....'
  call shr_sys_flush(6)

  ! -----------------------------------------------------------------
  ! Read input file
  ! -----------------------------------------------------------------

  ! Obtain input grid info, read local fields

  ns_o = ldomain%ns

  call domain_read(tdomain,datfname)
  ns_i = tdomain%ns

  if (masterproc) write (6,*) 'Open LAI file: ', trim(datfname)
  call check_ret(nf_open(datfname, 0, ncidi), subname)
  call check_ret(nf_inq_dimid(ncidi, 'pft', dimid), subname)
  call check_ret(nf_inq_dimlen(ncidi, dimid, numpft_i), subname)
  call check_ret(nf_inq_dimid(ncidi, 'time', dimid), subname)
  call check_ret(nf_inq_dimlen(ncidi, dimid, ntim), subname)

  if (numpft_i == 1) then
     use_pft_index = .true.
  else
     use_pft_index = .false.
     if (numpft_i /= numpft+1) then
        write(6,*)'MKLAI: parameter numpft+1= ',numpft+1, &
             'does not equal input dataset numpft= ',numpft_i
        stop
     endif
  endif

  if (ntim /= 12) then
     write(6,*)'MKLAI: must have 12 time samples on input data'
     call abort()
  endif

  ! NOTE - close data set at bottom of routine

  ! Dynamic allocation of variables

  allocate(mlai_i(ns_i,0:numpft),  &
           msai_i(ns_i,0:numpft),  &
           mhgtt_i(ns_i,0:numpft), &
           mhgtb_i(ns_i,0:numpft), &
           mask_src(ns_i),         &
           mlai_o(ns_o,0:numpft),  &
           msai_o(ns_o,0:numpft),  &
           mhgtt_o(ns_o,0:numpft), &
           mhgtb_o(ns_o,0:numpft), &
           laimask(ns_i,0:numpft), stat=ier )
  if (ier /= 0) then
     write(6,*)'mklai allocation error'; call abort()
  end if

  ! Determine mapping weights and map

  call gridmap_mapread(tgridmap, mapfname)

  ! Error checks for domain and map consistencies

  call domain_checksame( tdomain, ldomain, tgridmap )

  ! Determine number of dimensions in input by querying MONTHLY_LAI

  call check_ret(nf_inq_varid(ncidi, 'MONTHLY_LAI', varid), subname)
  call check_ret(nf_inq_vardimid(ncidi, varid, dimids), subname)
  call check_ret(nf_inq_varndims(ncidi, varid, ndimsi), subname)
  if (ndimsi ==4) then
     begi(1) = 1
     begi(2) = 1
     begi(3) = 1
     leni(4) = 1
     call check_ret(nf_inq_dimlen(ncidi, dimids(1), leni(1)), subname)
     call check_ret(nf_inq_dimlen(ncidi, dimids(2), leni(2)), subname)
     call check_ret(nf_inq_dimlen(ncidi, dimids(3), leni(3)), subname)
  else if (ndimsi== 3) then
     begi(1) = 1
     begi(2) = 1
     leni(3) = 1
     call check_ret(nf_inq_dimlen(ncidi, dimids(1), leni(1)), subname)
     call check_ret(nf_inq_dimlen(ncidi, dimids(2), leni(2)), subname)
  end if

  ! Determine number of dimensions in output by querying MONTHLY_LAI

  call check_ret(nf_inq_varid(ncido, 'MONTHLY_LAI', varid), subname)
  call check_ret(nf_inq_varndims(ncido, varid, ndimso), subname)
  call check_ret(nf_inq_vardimid(ncido, varid, dimids), subname)
  if (ndimso ==4) then
     bego(1) = 1
     bego(2) = 1
     bego(3) = 1
     leno(4) = 1
     call check_ret(nf_inq_dimlen(ncido, dimids(1), leno(1)), subname)
     call check_ret(nf_inq_dimlen(ncido, dimids(2), leno(2)), subname)
     call check_ret(nf_inq_dimlen(ncido, dimids(3), leno(3)), subname)
  else if (ndimso== 3) then
     bego(1) = 1
     bego(2) = 1
     leno(3) = 1
     call check_ret(nf_inq_dimlen(ncido, dimids(1), leno(1)), subname)
     call check_ret(nf_inq_dimlen(ncido, dimids(2), leno(2)), subname)
  end if

  if (use_pft_index) then
     allocate(pft_index(ns_i) , &
          mlai_i_tmp(ns_i)    , &
          msai_i_tmp(ns_i)    , &
          mhgtt_i_tmp(ns_i)   , &
          mhgtb_i_tmp(ns_i)   , stat=ier)
     call check_ret(nf_inq_varid (ncidi, 'NAT_PFT_INDEX', varid), subname)
     call check_ret(nf_get_vara_int (ncidi, varid, begi(1:ndimsi-2), leni(1:ndimsi-2), &
          pft_index), subname)
  endif

  ! Loop over months 

  do m = 1, ntim

     if (ndimsi == 4) begi(4)=m
     if (ndimsi == 3) begi(3)=m
     
     call check_ret(nf_inq_varid (ncidi, 'MONTHLY_LAI', varid), subname)
     call check_ret(nf_get_vara_double (ncidi, varid, begi(1:ndimsi), leni(1:ndimsi), &
          mlai_i), subname)

     call check_ret(nf_inq_varid (ncidi, 'MONTHLY_SAI', varid), subname)
     call check_ret(nf_get_vara_double (ncidi, varid, begi(1:ndimsi), leni(1:ndimsi), &
          msai_i(:,0:numpft)), subname)

     call check_ret(nf_inq_varid (ncidi, 'MONTHLY_HEIGHT_TOP', varid), subname)
     call check_ret(nf_get_vara_double (ncidi, varid, begi(1:ndimsi), leni(1:ndimsi), &
          mhgtt_i), subname)

     call check_ret(nf_inq_varid (ncidi, 'MONTHLY_HEIGHT_BOT', varid), subname)
     call check_ret(nf_get_vara_double (ncidi, varid, begi(1:ndimsi), leni(1:ndimsi), &
          mhgtb_i), subname)

     mlai_o(:,:)  = 0.
     msai_o(:,:)  = 0.
     mhgtt_o(:,:) = 0.
     mhgtb_o(:,:) = 0.
  
     ! Loop over pft types to do mapping
    
     if (.not.use_pft_index) then
        do l = 0,numpft
           mask_src(:) = 1._r8 
           call gridmap_areaave(tgridmap, mlai_i(:,l) , mlai_o(:,l) , nodata=0._r8, mask_src=mask_src)
           call gridmap_areaave(tgridmap, msai_i(:,l) , msai_o(:,l) , nodata=0._r8, mask_src=mask_src)
           call gridmap_areaave(tgridmap, mhgtt_i(:,l), mhgtt_o(:,l), nodata=0._r8, mask_src=mask_src)
           call gridmap_areaave(tgridmap, mhgtb_i(:,l), mhgtb_o(:,l), nodata=0._r8, mask_src=mask_src)
        enddo
     else
        mask_src(:) = 1._r8 
        do l = 0, numpft
           do i = 1, ns_i
              if (l == pft_index(i)) then
                 mlai_i_tmp(i)  = mlai_i(i,0)
                 msai_i_tmp(i)  = msai_i(i,0)
                 mhgtt_i_tmp(i) = mhgtt_i(i,0)
                 mhgtb_i_tmp(i) = mhgtb_i(i,0)
              else
                 mlai_i_tmp(i)  = 0._r8
                 msai_i_tmp(i)  = 0._r8
                 mhgtt_i_tmp(i) = 0._r8
                 mhgtb_i_tmp(i) = 0._r8
              end if
           enddo
           call gridmap_areaave(tgridmap, mlai_i_tmp(:) , mlai_o(:,l) , nodata=0._r8, mask_src=mask_src)
           call gridmap_areaave(tgridmap, msai_i_tmp(:) , msai_o(:,l) , nodata=0._r8, mask_src=mask_src)
           call gridmap_areaave(tgridmap, mhgtt_i_tmp(:), mhgtt_o(:,l), nodata=0._r8, mask_src=mask_src)
           call gridmap_areaave(tgridmap, mhgtb_i_tmp(:), mhgtb_o(:,l), nodata=0._r8, mask_src=mask_src)
        enddo
     endif

     ! Determine laimask
     
     laimask(:,:) = 0
     
     ! copy LAI, SAI, & heights from the C3 crop (pft15)
     ! to the irrigated (pft16) whether crop is on or off
     mlai_o(:,c3irrcropindex)  = mlai_o(:,c3cropindex)
     msai_o(:,c3irrcropindex)  = msai_o(:,c3cropindex)
     mhgtt_o(:,c3irrcropindex) = mhgtt_o(:,c3cropindex)
     mhgtb_o(:,c3irrcropindex) = mhgtb_o(:,c3cropindex)

     ! -----------------------------------------------------------------
     ! Output model resolution LAI/SAI/HEIGHT data
     ! -----------------------------------------------------------------

     ! Now write out all variables

     if (ndimso == 4) bego(4)=m
     if (ndimso == 3) bego(3)=m

     call check_ret(nf_inq_varid(ncido, 'MONTHLY_LAI', varid), subname)
     call check_ret(nf_put_vara_double(ncido, varid, bego, leno, mlai_o), subname)
     
     call check_ret(nf_inq_varid(ncido, 'MONTHLY_SAI', varid), subname)
     call check_ret(nf_put_vara_double(ncido, varid, bego, leno, msai_o), subname)
     
     call check_ret(nf_inq_varid(ncido, 'MONTHLY_HEIGHT_TOP', varid), subname)
     call check_ret(nf_put_vara_double(ncido, varid, bego, leno, mhgtt_o), subname)
     
     call check_ret(nf_inq_varid(ncido, 'MONTHLY_HEIGHT_BOT', varid), subname)
     call check_ret(nf_put_vara_double(ncido, varid, bego, leno, mhgtb_o), subname)

     call check_ret(nf_inq_varid(ncido, 'time', varid), subname)
     call check_ret(nf_put_vara_int(ncido, varid, bego(ndimso), leno(ndimso), m), subname)

     call check_ret(nf_sync(ncido), subname)


     ! -----------------------------------------------------------------
     ! Error check2
     ! Compare global areas on input and output grids
     ! -----------------------------------------------------------------

     ! Input grid global area

     garea_i    = 0.
     do ni = 1,ns_i
        garea_i = garea_i + tgridmap%area_src(ni)
     end do

     glai_i(:)  = 0.
     gsai_i(:)  = 0.
     ghgtt_i(:) = 0.
     ghgtb_i(:) = 0.
     do l  = 0, numpft
     do ni = 1, ns_i
        glai_i(l)  = glai_i(l) + mlai_i(ni,l) *tgridmap%area_src(ni)*&
             tgridmap%frac_src(ni)*re**2
        gsai_i(l)  = gsai_i(l) + msai_i(ni,l) *tgridmap%area_src(ni)*&
             tgridmap%frac_src(ni)*re**2
        ghgtt_i(l) = ghgtt_i(l)+ mhgtt_i(ni,l)*tgridmap%area_src(ni)*&
             tgridmap%frac_src(ni)*re**2
        ghgtb_i(l) = ghgtb_i(l)+ mhgtb_i(ni,l)*tgridmap%area_src(ni)*&
             tgridmap%frac_src(ni)*re**2
     end do
     end do

     ! Output grid global area

     garea_o    = 0.
     do no = 1,ns_o
        garea_o = garea_o + tgridmap%area_dst(no)
     end do

     glai_o(:)  = 0.
     gsai_o(:)  = 0.
     ghgtt_o(:) = 0.
     ghgtb_o(:) = 0.
     do l  = 0, numpft
     do no = 1,ns_o
        glai_o(l)  = glai_o(l) + mlai_o(no,l)*tgridmap%area_dst(no)* &
             tgridmap%frac_dst(no)*re**2
        gsai_o(l)  = gsai_o(l) + msai_o(no,l)*tgridmap%area_dst(no)* &
             tgridmap%frac_dst(no)*re**2
        ghgtt_o(l) = ghgtt_o(l)+ mhgtt_o(no,l)*tgridmap%area_dst(no)* &
             tgridmap%frac_dst(no)*re**2
        ghgtb_o(l) = ghgtb_o(l)+ mhgtb_o(no,l)*tgridmap%area_dst(no)* &
             tgridmap%frac_dst(no)*re**2
     end do
     end do

     ! Comparison

     write (ndiag,*)
     write (ndiag,'(1x,70a1)') ('=',k=1,70)
     write (ndiag,*) 'LAI Output for month ',m
     write (ndiag,'(1x,70a1)') ('=',k=1,70)

     write (ndiag,*)
     write (ndiag,'(1x,70a1)') ('.',k=1,70)
     write (ndiag,1001)
1001 format (1x,'PFT input grid area output grid area',/ &
             1x,3x,'     10**6 km**2','      10**6 km**2')
     write (ndiag,'(1x,70a1)') ('.',k=1,70)
     write (ndiag,*)
     do l = 0, numpft
        write (ndiag,1002) l, glai_i(l)*1.e-06*1.e-02,glai_o(l)*1.e-06*1.e-02
1002    format (1x,i3,f16.3,f17.3)
     end do

     if (masterproc) then
        write (6,*) 'Successfully made LAIs/SAIs/heights for month ', m
        write (6,*) ' '
     end if
     call shr_sys_flush(6)

  enddo

  ! Close input file
  call check_ret(nf_close(ncidi), subname)

  ! consistency check that PFT and LAI+SAI make sense
  !call pft_laicheck( ni_s, pft_i, laimask )

  ! Deallocate dynamic memory

  deallocate(mlai_i,msai_i,mhgtt_i,mhgtb_i,&
             mask_src,mlai_o,msai_o,mhgtt_o,mhgtb_o,laimask)
  if (use_pft_index) then
     deallocate(pft_index, mlai_i_tmp, msai_i_tmp, mhgtt_i_tmp, mhgtb_i_tmp)
  end if

  call gridmap_clean(tgridmap)
  call domain_clean(tdomain) 

end subroutine mklai

!-----------------------------------------------------------------------
!BOP
!
! !INTERFACE:
subroutine pft_laicheck( ni_s, pctpft_i, laimask )

! !USES:
!
! !DESCRIPTION:
!
! consistency check that PFT and LAI+SAI make sense
!
! !ARGUMENTS:
  implicit none
  integer , intent(in) :: ni_s          ! input PFT grid resolution
  real(r8), pointer    :: pctpft_i(:,:)  ! % plant function types
  integer,  pointer    :: laimask(:,:)   ! mask where LAI+SAI > 0
!EOP

  character(len=*), parameter :: subName="pft_laicheck"
  integer :: ni,l,n,nc      ! Indices
!-----------------------------------------------------------------------

  do l  = 0, numpft
     n  = 0
     nc = 0
     do ni = 1,ni_s
        if ( pctpft_i(ni,l) > 0.0_r8 ) nc = nc + 1
        if ( (pctpft_i(ni,l) > 0.0_r8) .and. (laimask(ni,l) /= 1) )then
           write (6,*) subName//' :: warning: pft and LAI+SAI mask not consistent!'
           write (6,*) 'ni,l   = ', ni, l
           write (6,*) 'pctpft_i  = ',pctpft_i(ni,l)
           write (6,*) 'laimask   = ', laimask(ni,l)
           n = n + 1
        end if
     end do
     if ( n > max(4,nc/4) ) then
        write (6,*) subName//' :: pft/LAI+SAI inconsistency over more than 25% land-cover'
        write (6,*) '# inconsistent points, total PFT pts, total LAI+SAI pts = ', &
                     n, nc, sum(laimask(:,l))
        stop
     end if
  end do

end subroutine pft_laicheck

!-----------------------------------------------------------------------
subroutine mklai_pio(ldomain_pio, mapfname, datfname, ndiag, pioIoSystem_o, ncid_o)

  use mkdomainPIOMod, only : domain_pio_type, domain_read_pio, domain_clean_pio
  use mkvarpar      , only : numstdpft, numstdcft
  use mkvarctl      , only : outnc_1d
  use mkpftMod      , only : c3cropindex, c3irrcropindex
  use mkgridmapPIOMod
  use piofileutils
  use pio

  implicit none

  type(domain_pio_type) , intent(inout) :: ldomain_pio
  character(len=*)      , intent(in)    :: mapfname              ! input mapping file name
  character(len=*)      , intent(in)    :: datfname              ! input dataset file name
  integer               , intent(in)    :: ndiag                 ! unit number for diag out
  type(iosystem_desc_t) , intent(in)    :: pioIoSystem_o
  type(file_desc_t)                     :: ncid_o                ! file descriptor of the output file

  type(domain_pio_type)                 :: tdomain_pio           ! local domain
  type(gridmap_pio_type)                :: tgridmap_pio
  type(file_desc_t)                     :: ncid_i
  type(iosystem_desc_t)                 :: pioIoSystem_i
  type(io_desc_t)                       :: iodesc
  integer               , pointer       :: compdof(:)
  integer               , pointer       :: pft_index(:,:)      ! PFT index for each grid cell
  real(r8)              , pointer       :: mlai4d_i(:,:,:,:)  , mlai1d_i(:)  , mlai3d_o(:,:,:)  ! monthly lai
  real(r8)              , pointer       :: msai4d_i(:,:,:,:)  , msai1d_i(:)  , msai3d_o(:,:,:)  ! monthly sai
  real(r8)              , pointer       :: mhgtt4d_i(:,:,:,:) , mhgtt1d_i(:) , mhgtt3d_o(:,:,:) ! monthly vegetation height top
  real(r8)              , pointer       :: mhgtb4d_i(:,:,:,:) , mhgtb1d_i(:) , mhgtb3d_o(:,:,:) ! monthly vegetation height bottom
  integer               , pointer       :: vec_row_indices(:), pft_vec_row_indices(:)
  integer                               :: dim_idx(4,2), dim2d(2), pft_dim_idx(2,2)
  integer                               :: i,j,p,t
  integer                               :: ns_loc_i,np,nt,count
  logical                               :: use_pft_index
  integer                               :: ierr

  !
  ! Obtain input grid info
  !

  call domain_read_pio(tdomain_pio, datfname)

  !
  ! Read the data from source file
  !

  call OpenFilePIO(datfname, pioIoSystem_i, ncid_i, PIO_NOWRITE)
  if (masterproc) write(6,*) 'Open LAI file: ', trim(datfname)

  call read_float_or_double_4d(tdomain_pio, pioIoSystem_i, ncid_i, 'MONTHLY_LAI', dim_idx, vec_row_indices, mlai4d_i)
  deallocate(vec_row_indices)

  call read_float_or_double_4d(tdomain_pio, pioIoSystem_i, ncid_i, 'MONTHLY_SAI', dim_idx, vec_row_indices, msai4d_i)
  deallocate(vec_row_indices)

  call read_float_or_double_4d(tdomain_pio, pioIoSystem_i, ncid_i, 'MONTHLY_HEIGHT_TOP', dim_idx, vec_row_indices, mhgtt4d_i)
  deallocate(vec_row_indices)

  call read_float_or_double_4d(tdomain_pio, pioIoSystem_i, ncid_i, 'MONTHLY_HEIGHT_BOT', dim_idx, vec_row_indices, mhgtb4d_i)

  np = dim_idx(3,2) - dim_idx(3,1) + 1
  nt = dim_idx(4,2) - dim_idx(4,1) + 1
  if (np == 1) then
     use_pft_index = .true.
     call read_integer_2d(tdomain_pio, pioIoSystem_i, ncid_i, 'NAT_PFT_INDEX', pft_dim_idx, pft_vec_row_indices, pft_index)
  else
     use_pft_index = .false.
  end if

  call PIO_closefile(ncid_i)
  call PIO_finalize(pioIoSystem_i, ierr)

  !
  ! Read the map
  !
  call gridmap_mapread_pio(tgridmap_pio, mapfname)

  ns_loc_i = (dim_idx(1,2) - dim_idx(1,1) + 1) * (dim_idx(2,2) - dim_idx(2,1) + 1)

  allocate(mlai1d_i (ns_loc_i))
  allocate(msai1d_i (ns_loc_i))
  allocate(mhgtt1d_i(ns_loc_i))
  allocate(mhgtb1d_i(ns_loc_i))


  allocate(mlai3d_o  (ldomain_pio%ns_loc, numpft + 1, nt))
  allocate(msai3d_o  (ldomain_pio%ns_loc, numpft + 1, nt))
  allocate(mhgtt3d_o (ldomain_pio%ns_loc, numpft + 1, nt))
  allocate(mhgtb3d_o (ldomain_pio%ns_loc, numpft + 1, nt))

  do t = 1, nt
     do p = 1, numpft + 1

        ! Convert data into 1D array

        if (.not.use_pft_index) then

           ! The input data contains information for ALL PFTs, so copy the values for the
           ! p-th PFT in the 1D array
           count = 0
           do j = dim_idx(2,1), dim_idx(2,2)
              do i = dim_idx(1,1), dim_idx(1,2)
                 count = count + 1
                 mlai1d_i  (count) = mlai4d_i  (i,j,p,t)
                 msai1d_i  (count) = msai4d_i  (i,j,p,t)
                 mhgtt1d_i (count) = mhgtt4d_i (i,j,p,t)
                 mhgtb1d_i (count) = mhgtb4d_i (i,j,p,t)
              end do
           end do

        else

           ! The input data contains information for only the dominant PFTs, so copy the values
           ! in the 1D array when the 'p' corresponds to the dominant PFT

           count = 0
           do j = dim_idx(2,1), dim_idx(2,2)
              do i = dim_idx(1,1), dim_idx(1,2)
                 count = count + 1
                 if (p == pft_index(i,j) + 1) then
                    mlai1d_i  (count) = mlai4d_i  (i,j,1,t)
                    msai1d_i  (count) = msai4d_i  (i,j,1,t)
                    mhgtt1d_i (count) = mhgtt4d_i (i,j,1,t)
                    mhgtb1d_i (count) = mhgtb4d_i (i,j,1,t)
                 else
                    mlai1d_i  (count) = 0._r8
                    msai1d_i  (count) = 0._r8
                    mhgtt1d_i (count) = 0._r8
                    mhgtb1d_i (count) = 0._r8
                 endif
              end do
           end do

        endif

        !
        ! Map the data onto output grid
        !
        call gridmap_areaave_pio(tgridmap_pio, ns_loc_i, vec_row_indices, mlai1d_i(:) , mlai3d_o(:,p,t) , nodata=0._r8)
        call gridmap_areaave_pio(tgridmap_pio, ns_loc_i, vec_row_indices, msai1d_i(:) , msai3d_o(:,p,t) , nodata=0._r8)
        call gridmap_areaave_pio(tgridmap_pio, ns_loc_i, vec_row_indices, mhgtt1d_i(:), mhgtt3d_o(:,p,t), nodata=0._r8)
        call gridmap_areaave_pio(tgridmap_pio, ns_loc_i ,vec_row_indices, mhgtb1d_i(:), mhgtb3d_o(:,p,t), nodata=0._r8)
     end do

     ! copy LAI, SAI, & heights from the C3 crop (pft15)
     ! to the irrigated (pft16) whether crop is on or off
     mlai3d_o (:,c3irrcropindex + 1, t) = mlai3d_o (:,c3cropindex + 1, t)
     msai3d_o (:,c3irrcropindex + 1, t) = msai3d_o (:,c3cropindex + 1, t)
     mhgtt3d_o(:,c3irrcropindex + 1, t) = mhgtt3d_o(:,c3cropindex + 1, t)
     mhgtb3d_o(:,c3irrcropindex + 1, t) = mhgtb3d_o(:,c3cropindex + 1, t)

   end do

   allocate(compdof(ldomain_pio%ns_loc * (numpft + 1) ))
   count = 0
   do p = 1, numpft + 1
      do i = 1, ldomain_pio%ns_loc
         count = count + 1
         compdof(count) = i + (p-1)*ldomain_pio%ns_glb + (ldomain_pio%begs - 1)
      end do
   end do

   if (outnc_1d) then
      dim2d(1) = ldomain_pio%ns_glb; dim2d(2) = numpft + 1;
      call PIO_initdecomp(pioIoSystem_o, PIO_DOUBLE, dim2d, compdof, iodesc)
   else
      write(6,*)'Extend mklai_pio to support the case when output mesh is not 1D'
      call abort()
   end if

   !
   ! Write data with unlimited time dimension one frame at a time in the file
   !
   do t = 1, nt
      call write_a_frame_of_double_3d(ncid_o, iodesc, 'MONTHLY_LAI'        , t, mlai3d_o  )
      call write_a_frame_of_double_3d(ncid_o, iodesc, 'MONTHLY_SAI'        , t, msai3d_o  )
      call write_a_frame_of_double_3d(ncid_o, iodesc, 'MONTHLY_HEIGHT_TOP' , t, mhgtt3d_o )
      call write_a_frame_of_double_3d(ncid_o, iodesc, 'MONTHLY_HEIGHT_BOT' , t, mhgtb3d_o )
   enddo

   call PIO_freedecomp(pioIoSystem_o, iodesc)
   deallocate(compdof)

 end subroutine mklai_pio

end module mklaiMod
