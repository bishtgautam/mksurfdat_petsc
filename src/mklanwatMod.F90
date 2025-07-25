module mklanwatMod

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: mklanwatMod
!
! !DESCRIPTION:
! make %lake and %wetland from input lake / wetland data
! also make lake parameters
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
!-----------------------------------------------------------------------
!
! !USES:
  use shr_kind_mod, only : r8 => shr_kind_r8
  use shr_sys_mod , only : shr_sys_flush
  use mkdomainMod , only : domain_checksame
  use spmdMod     , only : masterproc

  implicit none

  private

! !PUBLIC MEMBER FUNCTIONS:
  public mklakwat           ! make % lake
  public mklakwat_pio       ! make % lake PIO version
  public mkwetlnd           ! make % wetland
  public mkwetlnd_pio       ! make % wetland PIO version
  public mklakparams        ! make lake parameters
  public mklakparams_pio     ! make lake parameters PIO version

!EOP
!===============================================================
contains
!===============================================================

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: mklakwat
!
! !INTERFACE:
subroutine mklakwat(ldomain, mapfname, datfname, ndiag, zero_out, lake_o)
!
! !DESCRIPTION:
! make %lake
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
  logical           , intent(in) :: zero_out  ! if should zero glacier out
  real(r8)          , intent(out):: lake_o(:) ! output grid: %lake
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
  type(domain_type)    :: tdomain            ! local domain
  real(r8), allocatable :: lake_i(:)          ! input grid: percent lake
  real(r8) :: sum_fldi                        ! global sum of dummy input fld
  real(r8) :: sum_fldo                        ! global sum of dummy output fld
  real(r8) :: glake_i                         ! input  grid: global lake
  real(r8) :: garea_i                         ! input  grid: global area
  real(r8) :: glake_o                         ! output grid: global lake
  real(r8) :: garea_o                         ! output grid: global area
  integer  :: ni,no,k,ns_i,ns_o               ! indices
  integer  :: ncid,varid                      ! input netCDF id's
  integer  :: ier                             ! error status
  real(r8) :: relerr = 0.00001                ! max error: sum overlap wts ne 1
  character(len=32) :: subname = 'mklakwat'
!-----------------------------------------------------------------------

  if (masterproc) write (6,*) 'Attempting to make %lake and %wetland .....'
  call shr_sys_flush(6)

  ! -----------------------------------------------------------------
  ! Read input file
  ! -----------------------------------------------------------------

  ! Obtain input grid info, read local fields

  ns_o = ldomain%ns

  call domain_read(tdomain,datfname)
  ns_i = tdomain%ns

  if ( .not. zero_out )then
     allocate(lake_i(ns_i), stat=ier)
     if (ier/=0) call abort()

     if (masterproc) write(6,*)'Open lake file: ', trim(datfname)
     call check_ret(nf_open(datfname, 0, ncid), subname)
     call check_ret(nf_inq_varid (ncid, 'PCT_LAKE', varid), subname)
     call check_ret(nf_get_var_double (ncid, varid, lake_i), subname)
     call check_ret(nf_close(ncid), subname)

     ! Area-average percent cover on input grid to output grid 
     ! and correct according to land landmask
     ! Note that percent cover is in terms of total grid area.

     call gridmap_mapread(tgridmap, mapfname )

     ! Error checks for domain and map consistencies

     call domain_checksame( tdomain, ldomain, tgridmap )

     ! Determine lake_o on output grid

     call gridmap_areaave(tgridmap, lake_i,lake_o, nodata=0._r8)

     do no = 1,ns_o
        if (lake_o(no) < 1.) lake_o(no) = 0.
     enddo

     ! -----------------------------------------------------------------
     ! Error check prep
     ! Global sum of output field -- must multiply by fraction of
     ! output grid that is land as determined by input grid
     ! -----------------------------------------------------------------

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

     if ( .not. zero_out .and. (trim(mksrf_gridtype) == 'global') ) then
        if ( abs(sum_fldo/sum_fldi-1.) > relerr ) then
           write (6,*) 'MKLANWAT error: input field not conserved'
           write (6,'(a30,e20.10)') 'global sum output field = ',sum_fldo
           write (6,'(a30,e20.10)') 'global sum input  field = ',sum_fldi
           stop
        end if
     end if

     ! -----------------------------------------------------------------
     ! Error check2
     ! Compare global areas on input and output grids
     ! -----------------------------------------------------------------

     ! Input grid

     glake_i = 0.
     garea_i = 0.
     do ni = 1,ns_i
        garea_i = garea_i + tgridmap%area_src(ni)*re**2
        glake_i = glake_i + lake_i(ni)*tgridmap%area_src(ni)/100.*re**2
     end do

     ! Output grid

     glake_o = 0.
     garea_o = 0.
     do no = 1,ns_o
        garea_o = garea_o + tgridmap%area_dst(no)*re**2
        glake_o = glake_o + lake_o(no)*tgridmap%area_dst(no)/100.*re**2
     end do

     ! Diagnostic output

     write (ndiag,*)
     write (ndiag,'(1x,70a1)') ('=',k=1,70)
     write (ndiag,*) 'Inland Water Output'
     write (ndiag,'(1x,70a1)') ('=',k=1,70)

     write (ndiag,*)
     write (ndiag,'(1x,70a1)') ('.',k=1,70)
     write (ndiag,2001)
2001 format (1x,'surface type   input grid area  output grid area'/ &
             1x,'                 10**6 km**2      10**6 km**2   ')
     write (ndiag,'(1x,70a1)') ('.',k=1,70)
     write (ndiag,*)
     write (ndiag,2002) glake_i*1.e-06,glake_o*1.e-06
     write (ndiag,2004) garea_i*1.e-06,garea_o*1.e-06
2002 format (1x,'lakes       ',f14.3,f17.3)
2004 format (1x,'all surface ',f14.3,f17.3)
  else
     do no = 1,ns_o
        lake_o(no) = 0.
     enddo
  end if

  ! Deallocate dynamic memory

  call domain_clean(tdomain) 
  if ( .not. zero_out )then
     call gridmap_clean(tgridmap)
     deallocate (lake_i)
  end if

  if (masterproc) then
     write (6,*) 'Successfully made %lake'
     write (6,*)
  end if
  call shr_sys_flush(6)

end subroutine mklakwat

!-----------------------------------------------------------------------
subroutine mklakwat_pio(ldomain_pio, mapfname, datfname, ndiag, zero_out, lake_o)
!
! !DESCRIPTION:
! make %lake
!
! !USES:
  use mkdomainPIOMod, only : domain_pio_type, domain_clean_pio, domain_read_pio
  use mkgridmapMod
  use mkgridmapPIOMod
  use mkvarpar
  use mkvarctl
  use mkncdio
  use pio
  use piofileutils
!
! !ARGUMENTS:

  implicit none
  type(domain_pio_type), intent(in) :: ldomain_pio
  character(len=*)  , intent(in) :: mapfname  ! input mapping file name
  character(len=*)  , intent(in) :: datfname  ! input data file name
  integer           , intent(in) :: ndiag     ! unit number for diag out
  logical           , intent(in) :: zero_out  ! if should zero glacier out
  real(r8)          , intent(out):: lake_o(:) ! output grid: %lake
  !
  type(gridmap_pio_type) :: tgridmap_pio
  type(domain_pio_type)    :: tdomain_pio            ! local domain
  integer  :: no
  integer  :: ns_loc_i,ns_loc_o               !  indices

  type(file_desc_t)     :: ncid
  type(iosystem_desc_t) :: pioIoSystem
  real(r8) , pointer    :: lake2d_i(:,:)
  real(r8) , pointer    :: lake1d_i(:)
  integer               :: ierr
  integer               :: dim_idx(2,2)
  integer               :: i, j, count
  integer  , pointer    :: vec_row_indices(:)
  !-----------------------------------------------------------------------

  if (masterproc) write (6,*) 'Attempting to make %lake and %wetland .....'
  call shr_sys_flush(6)

  ! -----------------------------------------------------------------
  ! Read input file
  ! -----------------------------------------------------------------

  ! Obtain input grid info, read local fields

  call domain_read_pio(tdomain_pio, datfname)

  if ( .not. zero_out )then

     if (masterproc) write(6,*)'Open lake file: ', trim(datfname)

     ! Open the netcdf file
     call OpenFilePIO(datfname, pioIoSystem, ncid, PIO_NOWRITE)

     ! Read the variable
     call read_float_or_double_2d(tdomain_pio, pioIoSystem, ncid, 'PCT_LAKE', dim_idx, vec_row_indices, lake2d_i)

     call PIO_closefile(ncid)
     call PIO_finalize(pioIoSystem, ierr)

     ! Area-average percent cover on input grid to output grid
     ! and correct according to land landmask
     ! Note that percent cover is in terms of total grid area.

     call gridmap_mapread_pio(tgridmap_pio, mapfname )

     ! Convert 2D vector to 1D vector
     ns_loc_i = (dim_idx(1,2) - dim_idx(1,1) + 1) * (dim_idx(2,2) - dim_idx(2,1) + 1)

     allocate(lake1d_i(ns_loc_i))

     count = 0
     do j = dim_idx(2,1), dim_idx(2,2)
        do i = dim_idx(1,1), dim_idx(1,2)
           count = count + 1
           lake1d_i(count) = lake2d_i(i,j)
        end do
     end do

     ! Determine lake_o on output grid

     call gridmap_areaave_pio(tgridmap_pio, ns_loc_i, vec_row_indices, lake1d_i(:), lake_o, nodata=0._r8)

     ns_loc_o = ldomain_pio%ns_loc
     do no = 1, ns_loc_o
        if (lake_o(no) < 1._r8) lake_o(no) = 0._r8
     end do

  end if

  ! Deallocate dynamic memory

  call domain_clean_pio(tdomain_pio)
  if ( .not. zero_out )then
     call gridmap_clean_pio(tgridmap_pio)
     deallocate (lake2d_i)
     deallocate (lake1d_i)
     deallocate (vec_row_indices)
  end if

  if (masterproc) then
     write (6,*) 'Successfully made %lake'
     write (6,*)
  end if
  call shr_sys_flush(6)

end subroutine mklakwat_pio

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: mkwetlnd
!
! !INTERFACE:
subroutine mkwetlnd(ldomain, mapfname, datfname, ndiag, zero_out, swmp_o)
!
! !DESCRIPTION:
! make %wetland
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
  logical           , intent(in) :: zero_out  ! if should zero glacier out
  real(r8)          , intent(out):: swmp_o(:) ! output grid: %wetland
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
  type(domain_type)    :: tdomain            ! local domain
  real(r8), allocatable :: swmp_i(:)          ! input grid: percent swamp
  real(r8) :: sum_fldi                        ! global sum of dummy input fld
  real(r8) :: sum_fldo                        ! global sum of dummy output fld
  real(r8) :: gswmp_i                         ! input  grid: global swamp
  real(r8) :: garea_i                         ! input  grid: global area
  real(r8) :: gswmp_o                         ! output grid: global swamp
  real(r8) :: garea_o                         ! output grid: global area
  integer  :: ni,no,k,ns_i,ns_o               ! indices
  integer  :: ncid,varid                      ! input netCDF id's
  integer  :: ier                             ! error status
  real(r8) :: relerr = 0.00001                ! max error: sum overlap wts ne 1
  character(len=32) :: subname = 'mkwetlnd'
!-----------------------------------------------------------------------

  if (masterproc) write (6,*) 'Attempting to make %wetland .....'
  call shr_sys_flush(6)

  ! -----------------------------------------------------------------
  ! Read input file
  ! -----------------------------------------------------------------

  ! Obtain input grid info, read local fields

  ns_o = ldomain%ns

  call domain_read(tdomain,datfname)
  ns_i = tdomain%ns

  if ( .not. zero_out )then
     allocate(swmp_i(ns_i), stat=ier)
     if (ier/=0) call abort()

     if (masterproc) write(6,*)'Open wetland file: ', trim(datfname)
     call check_ret(nf_open(datfname, 0, ncid), subname)
     call check_ret(nf_inq_varid (ncid, 'PCT_WETLAND', varid), subname)
     call check_ret(nf_get_var_double (ncid, varid, swmp_i), subname)
     call check_ret(nf_close(ncid), subname)

     ! Area-average percent cover on input grid to output grid 
     ! and correct according to land landmask
     ! Note that percent cover is in terms of total grid area.

     call gridmap_mapread(tgridmap, mapfname )

     ! Error checks for domain and map consistencies

     call domain_checksame( tdomain, ldomain, tgridmap )
     ! Determine swmp_o on output grid

     call gridmap_areaave(tgridmap, swmp_i,swmp_o, nodata=0._r8)

     do no = 1,ns_o
        if (swmp_o(no) < 1.) swmp_o(no) = 0.
     enddo

     ! -----------------------------------------------------------------
     ! Error check prep
     ! Global sum of output field -- must multiply by fraction of
     ! output grid that is land as determined by input grid
     ! -----------------------------------------------------------------

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

     if ( .not. zero_out .and. (trim(mksrf_gridtype) == 'global') ) then
        if ( abs(sum_fldo/sum_fldi-1.) > relerr ) then
           write (6,*) 'MKLANWAT error: input field not conserved'
           write (6,'(a30,e20.10)') 'global sum output field = ',sum_fldo
           write (6,'(a30,e20.10)') 'global sum input  field = ',sum_fldi
           stop
        end if
     end if

     ! -----------------------------------------------------------------
     ! Error check2
     ! Compare global areas on input and output grids
     ! -----------------------------------------------------------------

     ! Input grid

     gswmp_i = 0.
     garea_i = 0.
     do ni = 1,ns_i
        garea_i = garea_i + tgridmap%area_src(ni)*re**2
        gswmp_i = gswmp_i + swmp_i(ni)*tgridmap%area_src(ni)/100.*re**2
     end do

     ! Output grid

     gswmp_o = 0.
     garea_o = 0.
     do no = 1,ns_o
        garea_o = garea_o + tgridmap%area_dst(no)*re**2
        gswmp_o = gswmp_o + swmp_o(no)*tgridmap%area_dst(no)/100.*re**2
     end do

     ! Diagnostic output

     write (ndiag,*)
     write (ndiag,'(1x,70a1)') ('=',k=1,70)
     write (ndiag,*) 'Inland Water Output'
     write (ndiag,'(1x,70a1)') ('=',k=1,70)

     write (ndiag,*)
     write (ndiag,'(1x,70a1)') ('.',k=1,70)
     write (ndiag,2001)
2001 format (1x,'surface type   input grid area  output grid area'/ &
             1x,'                 10**6 km**2      10**6 km**2   ')
     write (ndiag,'(1x,70a1)') ('.',k=1,70)
     write (ndiag,*)
     write (ndiag,2003) gswmp_i*1.e-06,gswmp_o*1.e-06
     write (ndiag,2004) garea_i*1.e-06,garea_o*1.e-06
2003 format (1x,'wetlands    ',f14.3,f17.3)
2004 format (1x,'all surface ',f14.3,f17.3)
  else
     do no = 1,ns_o
        swmp_o(no) = 0.
     enddo
  end if

  ! Deallocate dynamic memory

  call domain_clean(tdomain) 
  if ( .not. zero_out )then
     call gridmap_clean(tgridmap)
     deallocate (swmp_i)
  end if

  if (masterproc) write (6,*) 'Successfully made %wetland'
  call shr_sys_flush(6)

end subroutine mkwetlnd

!-----------------------------------------------------------------------
subroutine mkwetlnd_pio(ldomain_pio, mapfname, datfname, ndiag, zero_out, swmp_o)
  !
  ! !DESCRIPTION:
  ! make %wetland
  !
  ! !USES:
  use mkdomainPIOMod, only : domain_pio_type, domain_clean_pio, domain_read_pio
  use mkgridmapMod
  use mkgridmapPIOMod
  use mkvarpar
  use mkvarctl
  use mkncdio
  use pio
  use piofileutils
  !
  ! !ARGUMENTS:

  implicit none
  type(domain_pio_type), intent(in) :: ldomain_pio
  character(len=*)  , intent(in) :: mapfname  ! input mapping file name
  character(len=*)  , intent(in) :: datfname  ! input data file name
  integer           , intent(in) :: ndiag     ! unit number for diag out
  logical           , intent(in) :: zero_out  ! if should zero glacier out
  real(r8)          , intent(out):: swmp_o(:) ! output grid: %wetland
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
  type(gridmap_pio_type)    :: tgridmap_pio
  type(domain_pio_type)    :: tdomain_pio            ! local domain
  integer  :: no
  integer  :: ns_loc_i,ns_loc_o               ! indices

  type(file_desc_t)     :: ncid
  type(iosystem_desc_t) :: pioIoSystem
  real(r8) , pointer    :: swmp2d_i(:,:)
  real(r8) , pointer    :: swmp1d_i(:)
  integer               :: ierr
  integer               :: dim_idx(2,2)
  integer               :: i, j, count
  integer  , pointer    :: vec_row_indices(:)
  !-----------------------------------------------------------------------

  if (masterproc) write (6,*) 'Attempting to make %wetland .....'
  call shr_sys_flush(6)

  ! -----------------------------------------------------------------
  ! Read input file
  ! -----------------------------------------------------------------

  ! Obtain input grid info, read local fields

  if ( .not. zero_out )then

     if (masterproc) write(6,*)'Open wetland file: ', trim(datfname)

     ! Read the input domain
     call domain_read_pio(tdomain_pio, datfname)

     ! Open the netcdf file
     call OpenFilePIO(datfname, pioIoSystem, ncid, PIO_NOWRITE)

     ! Read the variable
     call read_float_or_double_2d(tdomain_pio, pioIoSystem, ncid, 'PCT_WETLAND', dim_idx, vec_row_indices, swmp2d_i)

     call PIO_closefile(ncid)
     call PIO_finalize(pioIoSystem, ierr)

     ! Area-average percent cover on input grid to output grid
     ! and correct according to land landmask
     ! Note that percent cover is in terms of total grid area.

     call gridmap_mapread_pio(tgridmap_pio, mapfname )

     ! Convert 2D vector to 1D vector
     ns_loc_i = (dim_idx(1,2) - dim_idx(1,1) + 1) * (dim_idx(2,2) - dim_idx(2,1) + 1)

     allocate(swmp1d_i(ns_loc_i))

     count = 0
     do j = dim_idx(2,1), dim_idx(2,2)
        do i = dim_idx(1,1), dim_idx(1,2)
           count = count + 1
           swmp1d_i(count) = swmp2d_i(i,j)
        end do
     end do

     ! Determine swmp_o on output grid

     call gridmap_areaave_pio(tgridmap_pio, ns_loc_i, vec_row_indices, swmp1d_i(:), swmp_o, nodata=0._r8)

     ns_loc_o = ldomain_pio%ns_loc
     do no = 1,ns_loc_o
        if (swmp_o(no) < 1._r8) swmp_o(no) = 0._r8
     enddo

  else

     ! set zero
     ns_loc_o = ldomain_pio%ns_loc
     do no = 1,ns_loc_o
        swmp_o(no) = 0._r8
     enddo

  end if

  ! Deallocate dynamic memory

  if ( .not. zero_out )then
     call domain_clean_pio(tdomain_pio)
     call gridmap_clean_pio(tgridmap_pio)
     deallocate (swmp2d_i)
     deallocate (swmp1d_i)
     deallocate (vec_row_indices)
  end if

  if (masterproc) write (6,*) 'Successfully made %wetland'
  call shr_sys_flush(6)

end subroutine mkwetlnd_pio

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: mklakparams
!
! !INTERFACE:
subroutine mklakparams(ldomain, mapfname, datfname, ndiag, &
                       lakedepth_o)
!
! !DESCRIPTION:
! make lake parameters (currently just lake depth)
!
! !USES:
  use mkdomainMod, only : domain_type, domain_clean, domain_read
  use mkgridmapMod
  use mkncdio
  use mkdiagnosticsMod, only : output_diagnostics_continuous
  use mkchecksMod, only : min_bad
!
! !ARGUMENTS:
  
  implicit none
  type(domain_type) , intent(in) :: ldomain
  character(len=*)  , intent(in) :: mapfname          ! input mapping file name
  character(len=*)  , intent(in) :: datfname          ! input data file name
  integer           , intent(in) :: ndiag             ! unit number for diag out
  real(r8)          , intent(out):: lakedepth_o(:)    ! output grid: lake depth (m)
!
! !CALLED FROM:
! subroutine mksrfdat in module mksrfdatMod
!
! !REVISION HISTORY:
! Author: Bill Sacks
!
!
! !LOCAL VARIABLES:
!EOP
  type(gridmap_type)    :: tgridmap
  type(domain_type)     :: tdomain            ! local domain
  real(r8), allocatable :: data_i(:)          ! data on input grid
  integer  :: ncid,varid                      ! input netCDF id's
  integer  :: ier                             ! error status
  
  real(r8), parameter :: min_valid_lakedepth = 0._r8

  character(len=32) :: subname = 'mklakparams'
!-----------------------------------------------------------------------

  write (6,*) 'Attempting to make lake parameters.....'
  call shr_sys_flush(6)
  write(*,*)'mapfname:' ,trim(mapfname)
  write(*,*)'datfname:' ,trim(datfname)

  ! -----------------------------------------------------------------
  ! Read domain and mapping information, check for consistency
  ! -----------------------------------------------------------------

  call domain_read(tdomain,datfname)
  
  call gridmap_mapread(tgridmap, mapfname )
  call gridmap_check( tgridmap, subname )

  call domain_checksame( tdomain, ldomain, tgridmap )

  ! -----------------------------------------------------------------
  ! Open input file, allocate memory for input data
  ! -----------------------------------------------------------------

  write(6,*)'Open lake parameter file: ', trim(datfname)
  call check_ret(nf_open(datfname, 0, ncid), subname)

  allocate(data_i(tdomain%ns), stat=ier)
  if (ier/=0) call abort()

  ! -----------------------------------------------------------------
  ! Regrid lake depth
  ! -----------------------------------------------------------------

  call check_ret(nf_inq_varid (ncid, 'LAKEDEPTH', varid), subname)
  call check_ret(nf_get_var_double (ncid, varid, data_i), subname)
  call gridmap_areaave(tgridmap, data_i, lakedepth_o, nodata=10._r8)

  ! Check validity of output data
  if (min_bad(lakedepth_o, min_valid_lakedepth, 'lakedepth')) then
     stop
  end if

  call output_diagnostics_continuous(data_i, lakedepth_o, tgridmap, "Lake Depth", "m", ndiag)

  ! -----------------------------------------------------------------
  ! Close files and deallocate dynamic memory
  ! -----------------------------------------------------------------

  call check_ret(nf_close(ncid), subname)
  call domain_clean(tdomain) 
  call gridmap_clean(tgridmap)
  deallocate (data_i)

  write (6,*) 'Successfully made lake parameters'
  write (6,*)
  call shr_sys_flush(6)

end subroutine mklakparams

!-----------------------------------------------------------------------
subroutine mklakparams_pio(ldomain_pio, mapfname, datfname, ndiag, lakedepth_o)
  !
  ! !USES:
  use mkdomainPIOMod, only : domain_pio_type
  use mkdataPIOMod
  !
  implicit none
  !
  type(domain_pio_type) , intent(in)  :: ldomain_pio
  character(len=*)      , intent(in)  :: mapfname       ! input mapping file name
  character(len=*)      , intent(in)  :: datfname       ! input data file name
  integer               , intent(in)  :: ndiag          ! unit number for diag out
  real(r8)              , intent(out) :: lakedepth_o(:) ! output grid: lake depth (m)

  real(r8), parameter :: nodata_value = 0._r8
  real(r8), parameter :: min_valid    = 0._r8
  !-----------------------------------------------------------------------

  if (masterproc) write (6,*) 'Attempting to make lake parameters.....'
  call shr_sys_flush(6)

  call mkdata_double_2d_pio(ldomain_pio, mapfname=mapfname, datfname=datfname, varname='LAKEDEPTH', &
       data_descrip='LAKEDEPTH', ndiag=ndiag, zero_out=.false., nodata_value=nodata_value, data_o=lakedepth_o, &
       min_valid_value=min_valid)

  if (masterproc) write (6,*) 'Successfully made lake parameters'
  call shr_sys_flush(6)

end subroutine mklakparams_pio

end module mklanwatMod
