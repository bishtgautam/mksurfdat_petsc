module mksoilphosphorusMod
  !-----------------------------------------------------------------------
  !BOP
  !
  ! !MODULE: mksoilphosphorusMod
  !
  ! !DESCRIPTION:
  ! Make soil phosphorus data
  !
  ! !REVISION HISTORY:
  ! Author: Gautam Bisht
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
  public mksoilphosphorus     ! Set soil phosphorus
  public mksoilphosphorus_pio ! PIO-version of set soil phosphorus

  !EOP
  !===============================================================
contains

  !-----------------------------------------------------------------------
  !BOP
  !
  ! !IROUTINE: mksoilphosphorus
  !
  ! !INTERFACE:
  subroutine mksoilphosphorus(ldomain, mapfname, datfname, ndiag, apatiteP_o, &
       labileP_o, occludedP_o, secondaryP_o)
    !
    ! !DESCRIPTION:
    ! make soil phosphorus datasets
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
    !
    type(domain_type) , intent(in)  :: ldomain
    character(len=*)  , intent(in)  :: mapfname        ! input mapping file name
    character(len=*)  , intent(in)  :: datfname        ! input data file name
    integer           , intent(in)  :: ndiag           ! unit number for diag out
    real(r8)          , intent(out) :: apatiteP_o(:)   ! apatite phosphorus (output grid)
    real(r8)          , intent(out) :: labileP_o(:)    ! labile phosphorus (output grid)
    real(r8)          , intent(out) :: occludedP_o(:)  ! occluded phosphorus (output grid)
    real(r8)          , intent(out) :: secondaryP_o(:) ! secondaryP phosphorus (output grid)
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
    type(domain_type)     :: tdomain                      ! local domain
    integer               :: ncid,varid                   ! input netCDF id's
    integer               :: ier                          ! error status
    real(r8), allocatable :: data_i(:)                    ! data on input grid
    real(r8), parameter   :: apatiteP_nodata   = 0._r8
    real(r8), parameter   :: labileP_nodata    = 0._r8
    real(r8), parameter   :: occludedP_nodata  = 0._r8
    real(r8), parameter   :: secondaryP_nodata = 0._r8
    character(len=32)     :: subname = 'mksoilphosphorus' !
    !-----------------------------------------------------------------------

    write (6,*) 'Attempting to make soil phosphorus .....'
    call shr_sys_flush(6)
    write(*,*)'mapfname:' ,trim(mapfname)
    write(*,*)'datfname:' ,trim(datfname)
  
    ! -----------------------------------------------------------------
    ! Read input file
    ! -----------------------------------------------------------------

    ! Obtain input grid info, read local fields

    call domain_read(tdomain,datfname)

    call gridmap_mapread(tgridmap, mapfname )
    call gridmap_check(  tgridmap, subname )

    call domain_checksame( tdomain, ldomain, tgridmap )

    ! -----------------------------------------------------------------
    ! Open input file, allocate memory for input data
    ! -----------------------------------------------------------------
    write (6,*) 'Open soil phosphorus file: ', trim(datfname)
    call check_ret(nf_open(datfname, 0, ncid), subname)
    
    allocate(data_i(tdomain%ns), stat=ier)
    if (ier/=0) call abort()

    ! -----------------------------------------------------------------
    ! Regrid APATITE_P
    ! -----------------------------------------------------------------

    call check_ret(nf_inq_varid (ncid, 'APATITE_P', varid), subname)
    call check_ret(nf_get_var_double (ncid, varid, data_i), subname)
    call gridmap_areaave(tgridmap, data_i, apatiteP_o, nodata=apatiteP_nodata)

    ! -----------------------------------------------------------------
    ! Regrid LABILE_P
    ! -----------------------------------------------------------------

    call check_ret(nf_inq_varid (ncid, 'LABILE_P', varid), subname)
    call check_ret(nf_get_var_double (ncid, varid, data_i), subname)
    call gridmap_areaave(tgridmap, data_i, labileP_o, nodata=labileP_nodata)

    ! -----------------------------------------------------------------
    ! Regrid OCCLUDED_P
    ! -----------------------------------------------------------------

    call check_ret(nf_inq_varid (ncid, 'OCCLUDED_P', varid), subname)
    call check_ret(nf_get_var_double (ncid, varid, data_i), subname)
    call gridmap_areaave(tgridmap, data_i, occludedP_o, nodata=occludedP_nodata)

    ! -----------------------------------------------------------------
    ! Regrid SECONDARY_P
    ! -----------------------------------------------------------------

    call check_ret(nf_inq_varid (ncid, 'SECONDARY_P', varid), subname)
    call check_ret(nf_get_var_double (ncid, varid, data_i), subname)
    call gridmap_areaave(tgridmap, data_i, secondaryP_o, nodata=secondaryP_nodata)

    ! -----------------------------------------------------------------
    ! Close files and deallocate dynamic memory
    ! -----------------------------------------------------------------
    call check_ret(nf_close(ncid), subname)
    call domain_clean(tdomain)
    call gridmap_clean(tgridmap)
    deallocate (data_i)

    write (6,*) 'Successfully made phosphorus'
    write (6,*)
    call shr_sys_flush(6)

  end subroutine mksoilphosphorus

  !-----------------------------------------------------------------------
  subroutine mksoilphosphorus_pio(ldomain_pio, mapfname, datfname, ndiag, apatiteP_o, &
       labileP_o, occludedP_o, secondaryP_o)
    !
    ! !DESCRIPTION:
    ! make soil phosphorus datasets
    !
    ! !USES:
    use mkdomainPIOMod, only : domain_pio_type
    use mkdataPIOMod
    !
    type(domain_pio_type) , intent(in)  :: ldomain_pio
    character(len=*)  , intent(in)  :: mapfname        ! input mapping file name
    character(len=*)  , intent(in)  :: datfname        ! input data file name
    integer           , intent(in)  :: ndiag           ! unit number for diag out
    real(r8)          , intent(out) :: apatiteP_o(:)   ! apatite phosphorus (output grid)
    real(r8)          , intent(out) :: labileP_o(:)    ! labile phosphorus (output grid)
    real(r8)          , intent(out) :: occludedP_o(:)  ! occluded phosphorus (output grid)
    real(r8)          , intent(out) :: secondaryP_o(:) ! secondaryP phosphorus (output grid)
    !
    real(r8), parameter   :: nodata_value = 0._r8
    real(r8), parameter   :: min_valid    = 0._r8
    !-----------------------------------------------------------------------

    write (6,*) 'Attempting to make soil phosphorus .....'
    call shr_sys_flush(6)
    write(*,*)'mapfname:' ,trim(mapfname)
    write(*,*)'datfname:' ,trim(datfname)
  
  call mkdata_double_2d_pio(ldomain_pio, mapfname=mapfname, datfname=datfname, varname='APATITE_P', &
       data_descrip='APATITE_P', ndiag=ndiag, zero_out=.false., nodata_value=nodata_value, data_o=apatiteP_o, &
       min_valid_value=min_valid)

  call mkdata_double_2d_pio(ldomain_pio, mapfname=mapfname, datfname=datfname, varname='LABILE_P', &
       data_descrip='LABILE_P', ndiag=ndiag, zero_out=.false., nodata_value=nodata_value, data_o=labileP_o, &
       min_valid_value=min_valid)

  call mkdata_double_2d_pio(ldomain_pio, mapfname=mapfname, datfname=datfname, varname='OCCLUDED_P', &
       data_descrip='OCCLUDED_P', ndiag=ndiag, zero_out=.false., nodata_value=nodata_value, data_o=occludedP_o, &
       min_valid_value=min_valid)

  call mkdata_double_2d_pio(ldomain_pio, mapfname=mapfname, datfname=datfname, varname='SECONDARY_P', &
       data_descrip='SECONDARY_P', ndiag=ndiag, zero_out=.false., nodata_value=nodata_value, data_o=secondaryP_o, &
       min_valid_value=min_valid)

    write (6,*) 'Successfully made phosphorus'
    write (6,*)
    call shr_sys_flush(6)

  end subroutine mksoilphosphorus_pio

end module mksoilphosphorusMod
