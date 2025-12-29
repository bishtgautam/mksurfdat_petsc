module mktopradMod
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: mktopradMod
!
! !DESCRIPTION:
! Make topography data for TOP solar radiation parameterization
!
! !REVISION HISTORY:
! Author: Dalei Hao
!
!-----------------------------------------------------------------------
!!USES:
  use shr_kind_mod, only : r8 => shr_kind_r8
  use shr_sys_mod , only : shr_sys_flush
  use spmdMod     , only : masterproc
  implicit none

  SAVE
  private           ! By default make data private
!
! !PUBLIC MEMBER FUNCTIONS:
!
  public mktopradAtt      ! Add attributes to output file
  public mktopradAttPIO   ! PIO-version to add attributes to output file

  public mktoprad         ! Set topography
  public mktoprad_pio ! PIO-version of set topography data for TOP solar rad param
!
! !PUBLIC DATA MEMBERS:
!
!
! !PRIVATE DATA MEMBERS:
!
! !PRIVATE MEMBER FUNCTIONS:

!EOP
!===============================================================
contains
!===============================================================

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: mktoprad
!
! !INTERFACE:
subroutine mktoprad(ldomain, mapfname, datfname, sinsl_sinas_o, sinsl_cosas_o, sky_view_o, terrain_config_o)
!
! !DESCRIPTION:
! Make topography data for TOP solar radiation parameterization
!
! !USES:
  use mkdomainMod  , only : domain_type, domain_clean, domain_read, domain_checksame
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
  real(r8)          , intent(out):: sinsl_sinas_o(:)  ! output topography data
  real(r8)          , intent(out):: sinsl_cosas_o(:)  ! output topography data
  real(r8)          , intent(out):: sky_view_o(:)  ! output topography data
  real(r8)          , intent(out):: terrain_config_o(:)  ! output topography data
!
!
! !CALLED FROM:
! subroutine mksrfdat in module mksrfdatMod
!
! !REVISION HISTORY:
! Author: Dalei Hao
!
!
! !LOCAL VARIABLES:
!EOP
  type(domain_type)     :: tdomain            ! local domain
  type(gridmap_type)    :: tgridmap           ! local gridmap

  real(r8), allocatable :: sinsl_sinas_i(:)
  real(r8), allocatable :: sinsl_cosas_i(:)
  real(r8), allocatable :: sky_view_i(:)
  real(r8), allocatable :: terrain_config_i(:)
  real(r8), allocatable :: mask_i(:)          ! input grid: mask (0, 1)
  integer  :: ns_i,ns_o                       ! indices
  integer  :: ni                              ! indices
  integer  :: ncidi,varid                     ! input netCDF id's
  integer  :: ier                             ! error status
  character(len= 32) :: subname = 'mktop'
!-----------------------------------------------------------------------

  write (6,*) 'Attempting to make topography .....'
  call shr_sys_flush(6)
  write(*,*)'mapfname:' ,trim(mapfname)
  write(*,*)'datfname:' ,trim(datfname)

  ns_o = ldomain%ns

  ! -----------------------------------------------------------------
  ! Read input file
  ! -----------------------------------------------------------------

  ! Obtain input grid info, read local fields

  call domain_read(tdomain,datfname)

  ns_i = tdomain%ns
  allocate(sinsl_sinas_i(ns_i), sinsl_cosas_i(ns_i), sky_view_i(ns_i), terrain_config_i(ns_i), stat=ier)
  if (ier /= 0) then
     write(6,*)'mktoprad allocation error'; call abort()
  end if

  write (6,*) 'Open topography file: ', trim(datfname)
  call check_ret(nf_open(datfname, 0, ncidi), subname)

  call check_ret(nf_inq_varid (ncidi, 'SINSL_SINAS', varid), subname)
  call check_ret(nf_get_var_double (ncidi, varid, sinsl_sinas_i), subname)

  call check_ret(nf_inq_varid (ncidi, 'SINSL_COSAS', varid), subname)
  call check_ret(nf_get_var_double (ncidi, varid, sinsl_cosas_i), subname)

  call check_ret(nf_inq_varid (ncidi, 'SKY_VIEW', varid), subname)
  call check_ret(nf_get_var_double (ncidi, varid, sky_view_i), subname)

  call check_ret(nf_inq_varid (ncidi, 'TERRAIN_CONFIG', varid), subname)
  call check_ret(nf_get_var_double (ncidi, varid, terrain_config_i), subname)


  call check_ret(nf_close(ncidi), subname)

  ! set mask as 0 when topo data is filled value: -9999
  allocate(mask_i(ns_i), stat=ier)
  if (ier /= 0) then
     write(6,*)'mktoprad allocation error'; call abort()
  end if
  
  mask_i(:) = 1._r8
  do ni = 1,ns_i
      if (sinsl_sinas_i(ni) < -1000._r8 .or. sinsl_cosas_i(ni) < -1000._r8 .or. sky_view_i(ni) < -1000._r8 .or. terrain_config_i(ni) < -1000._r8) then
         mask_i(ni) = 0._r8
     end if
  enddo

  ! Read mapping file
  call gridmap_mapread(tgridmap, mapfname)

  ! Error checks for domain and map consistencies
  call domain_checksame( tdomain, ldomain, tgridmap )

  ! Determine top_o on output grid
  sinsl_sinas_o(:)    = 0._r8
  sinsl_cosas_o(:)    = 0._r8
  sky_view_o(:)       = 1._r8
  terrain_config_o(:) = 0._r8

  call gridmap_areaave(tgridmap, sinsl_sinas_i   , sinsl_sinas_o   , nodata=0._r8, mask_src=mask_i)
  call gridmap_areaave(tgridmap, sinsl_cosas_i   , sinsl_cosas_o   , nodata=0._r8, mask_src=mask_i)
  call gridmap_areaave(tgridmap, sky_view_i      , sky_view_o      , nodata=1._r8, mask_src=mask_i)
  call gridmap_areaave(tgridmap, terrain_config_i, terrain_config_o, nodata=0._r8, mask_src=mask_i)

  ! Deallocate dynamic memory

  call domain_clean(tdomain)
  call gridmap_clean(tgridmap)
  deallocate (sinsl_sinas_i, sinsl_cosas_i, sky_view_i, terrain_config_i)
  deallocate (mask_i)

  write (6,*) 'Successfully made topography parameters'
  write (6,*)
  call shr_sys_flush(6)

end subroutine mktoprad

!-----------------------------------------------------------------------

subroutine mktoprad_pio(ldomain_pio, mapfname, datfname, ndiag, sinsl_sinas_o, sinsl_cosas_o, sky_view_o, terrain_config_o)
!
! !DESCRIPTION:
! Make topography data for TOP solar radiation parameterization
!
! !USES:
  use mkdomainPIOMod, only : domain_pio_type
  use mkdataPIOMod
!
  type(domain_pio_type), intent(in) :: ldomain_pio
  character(len=*)  , intent(in) :: mapfname          ! input mapping file name
  character(len=*)  , intent(in) :: datfname          ! input data file name
  integer           , intent(in) :: ndiag             ! unit number for diag out
  real(r8)          , intent(out):: sinsl_sinas_o(:)  ! output topography data
  real(r8)          , intent(out):: sinsl_cosas_o(:)  ! output topography data
  real(r8)          , intent(out):: sky_view_o(:)     ! output topography data
  real(r8)          , intent(out):: terrain_config_o(:)  ! output topography data
!
  real(r8), parameter   :: nodata_value = 0._r8
  real(r8), parameter   :: min_valid    = 0._r8
!-----------------------------------------------------------------------

  if (masterproc) write (6,*) 'Attempting to make topography .....'
  call shr_sys_flush(6)
  
  call mkdata_double_2d_pio(ldomain_pio, mapfname=mapfname, datfname=datfname, varname='SINSL_SINAS', &
     data_descrip='SINSL_SINAS', ndiag=ndiag, zero_out=.false., nodata_value=nodata_value, data_o=sinsl_sinas_o, &
     min_valid_value=min_valid)

  call mkdata_double_2d_pio(ldomain_pio, mapfname=mapfname, datfname=datfname, varname='SINSL_COSAS', &
     data_descrip='SINSL_COSAS', ndiag=ndiag, zero_out=.false., nodata_value=nodata_value, data_o=sinsl_cosas_o, &
     min_valid_value=min_valid)

  call mkdata_double_2d_pio(ldomain_pio, mapfname=mapfname, datfname=datfname, varname='SKY_VIEW', &
     data_descrip='SKY_VIEW', ndiag=ndiag, zero_out=.false., nodata_value=nodata_value, data_o=sky_view_o, &
     min_valid_value=min_valid)

  call mkdata_double_2d_pio(ldomain_pio, mapfname=mapfname, datfname=datfname, varname='TERRAIN_CONFIG', &
     data_descrip='TERRAIN_CONFIG', ndiag=ndiag, zero_out=.false., nodata_value=nodata_value, data_o=terrain_config_o, &
     min_valid_value=min_valid)

    if (masterproc) then
       write (6,*) 'Successfully made topography parameters'
       write (6,*)
    end if
    call shr_sys_flush(6)

end subroutine mktoprad_pio

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: mktopradAtt
!
! !INTERFACE:
subroutine mktopradAtt( ncid, dynlanduse, xtype )
!
! !DESCRIPTION:
! add atttributes to output file regarding the topography module
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
! Original Author: Dalei Hao
!
!
! !LOCAL VARIABLES:
!EOP
  character(len=256) :: str       ! global attribute string
  character(len=32) :: subname = 'mktopAtt'
!-----------------------------------------------------------------------

  if (.not. dynlanduse) then

     ! Add global attributes to file

     str = get_filename(mksrf_fgrvl)
     call check_ret(nf_put_att_text(ncid, NF_GLOBAL, &
          'top_raw_data_file_name', len_trim(str), trim(str)), subname)
     
     ! Define variables

     if (outnc_1d) then
        call ncd_defvar(ncid=ncid, varname='SINSL_COSAS', xtype=xtype, &
             dim1name='gridcell',&
             long_name='sin(slope) * cos(aspect)', units='unitless')
     else
        call ncd_defvar(ncid=ncid, varname='SINSL_COSAS', xtype=xtype, &
             dim1name='lsmlon', dim2name='lsmlat', &
             long_name='sin(slope) * cos(aspect)', units='unitless')
     end if

     if (outnc_1d) then
        call ncd_defvar(ncid=ncid, varname='SINSL_SINAS', xtype=xtype, &
             dim1name='gridcell',&
             long_name='sin(slope) * sin(aspect)', units='unitless')
     else
        call ncd_defvar(ncid=ncid, varname='SINSL_SINAS', xtype=xtype, &
             dim1name='lsmlon', dim2name='lsmlat', &
             long_name='sin(slope) * sin(aspect)', units='unitless')
     end if

    if (outnc_1d) then
        call ncd_defvar(ncid=ncid, varname='SKY_VIEW', xtype=xtype, &
            dim1name='gridcell',&
            long_name='sky view factor', units='unitless')
    else
        call ncd_defvar(ncid=ncid, varname='SKY_VIEW', xtype=xtype, &
            dim1name='lsmlon', dim2name='lsmlat', &
            long_name='sky view factor', units='unitless')
    end if

    if (outnc_1d) then
        call ncd_defvar(ncid=ncid, varname='TERRAIN_CONFIG', xtype=xtype, &
            dim1name='gridcell',&
            long_name='terrain configuration factor', units='unitless')
    else
        call ncd_defvar(ncid=ncid, varname='TERRAIN_CONFIG', xtype=xtype, &
            dim1name='lsmlon', dim2name='lsmlat', &
            long_name='terrain configuration factor', units='unitless')
    end if

  end if

end subroutine mktopradAtt

!-----------------------------------------------------------------------

subroutine mktopradAttPIO( ncid, dynlanduse, xtype, dim_id_gridcell, dim_id_lsmlon, dim_id_lsmlat)

  use fileutils  , only : get_filename
  use pio
  use piofileutils
  use mkvarpar   
  use mkvarctl   

  implicit none

  type(file_desc_t) , intent(in)    :: ncid
  logical, intent(in) :: dynlanduse   ! if dynamic land-use file
  integer, intent(in) :: xtype        ! external type to output real data as
  integer, intent(in) :: dim_id_gridcell
  integer, intent(in) :: dim_id_lsmlon
  integer, intent(in) :: dim_id_lsmlat

  character(len=256) :: str       ! global attribute string
  character(len=32)  :: subname = 'mktopradAttPIO'
  integer            :: dim1d(1), dim2d(2)

  if (.not. dynlanduse) then

     ! Add global attributes to file

     str = get_filename(mksrf_fgrvl)
     call check_ret(PIO_put_att(ncid, PIO_GLOBAL, 'top_raw_data_file_name', trim(str)), subname)

     ! Define variables

     if (outnc_1d) then
        dim1d(1) = dim_id_gridcell
        call DefineVarPIO_1d(ncid, 'SINSL_COSAS', PIO_INT, dim1d, longName='sin(slope) * cos(aspect)', units='unitless')
        call DefineVarPIO_1d(ncid, 'SINSL_SINAS', PIO_INT, dim1d, longName='sin(slope) * sin(aspect)', units='unitless')
        call DefineVarPIO_1d(ncid, 'SKY_VIEW', PIO_INT, dim1d, longName='sky view factor', units='unitless')
        call DefineVarPIO_1d(ncid, 'TERRAIN_CONFIG', PIO_INT, dim1d, longName='terrain configuration factor', units='unitless')

     else
        dim2d(1) = dim_id_lsmlon; dim2d(2) = dim_id_lsmlat;
        call DefineVarPIO_2d(ncid, 'SINSL_COSAS', PIO_INT, dim2d, longName='sin(slope) * cos(aspect)', units='unitless')
        call DefineVarPIO_2d(ncid, 'SINSL_SINAS', PIO_INT, dim2d, longName='sin(slope) * sin(aspect)', units='unitless')
        call DefineVarPIO_2d(ncid, 'SKY_VIEW', PIO_INT, dim2d, longName='sky view factor', units='unitless')
        call DefineVarPIO_2d(ncid, 'TERRAIN_CONFIG', PIO_INT, dim2d, longName='terrain configuration factor', units='unitless')

     end if
  end if

end subroutine mktopradAttPIO

end module mktopradMod
