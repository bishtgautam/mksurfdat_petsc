module piofileutils

  use shr_kind_mod, only : r8 => shr_kind_r8
  use mkdomainPIOMod
  use pio

  implicit none

  public OpenFilePio
  public read_float_or_double_2d
  public read_float_or_double_3d
  

contains

  !-----------------------------------------------------------------------
  subroutine check_ret(ret, calling)
    !
    ! !DESCRIPTION:
    ! Check return status from netcdf call
    !
    ! !ARGUMENTS:
    use pio, only : PIO_NOERR, PIO_STRERROR
    implicit none
    integer, intent(in) :: ret
    character(len=*)    :: calling
    !
    integer             :: ier
    character(len=256)  :: errmsg

    if (ret /= PIO_NOERR) then
       ier = PIO_STRERROR(ret, errmsg)
       write(6,*)'netcdf error from ',trim(calling), ' rcode = ', ret, &
            ' error = ', trim(errmsg)
       call abort()
    end if

  end subroutine check_ret

  !-----------------------------------------------------------------------
  subroutine OpenFilePIO(fname, pioIoSystem, ncid)

    use petsc
    use spmdMod, only : iam, npes, masterproc, mpicom
    use pio

    implicit none

    character (len=*)     :: fname
    type(file_desc_t)     :: ncid
    type(iosystem_desc_t) :: pioIoSystem
    !
    integer               :: niotasks
    integer               :: numAggregator
    integer               :: stride
    integer               :: optBase
    integer               :: iotype
    integer               :: retVal, ier
    character(len= 32)    :: subname = 'OpenFilePIO'

    stride        = 1
    numAggregator = 0
    iotype        = PIO_iotype_pnetcdf
    niotasks      = npes

    if (npes == 1) then
       optBase = 0
    else
       optBase = 1
    end if

    call PIO_init(iam,     & ! MPI rank
         MPI_COMM_WORLD,   & ! MPI communicator
         niotasks,         & ! Number of iotasks (ntasks/stride)
         numAggregator,    & ! number of aggregators to use
         stride,           & ! stride
         PIO_rearr_subset, & ! do not use any form of rearrangement
         pioIoSystem,      & ! iosystem
         optBase)

    call check_ret(PIO_openfile(pioIoSystem, ncid, iotype, trim(fname), PIO_NOWRITE), subname)

  end subroutine OpenFilePIO

  !-----------------------------------------------------------------------
  subroutine read_float_or_double_2d(domain, pioIoSystem, ncid, varname, dim_idx, data2d)

    use mkdomainPIOMod
    use pio

    implicit none

    type(domain_pio_type) , intent(in)           :: domain
    type(iosystem_desc_t) , intent(in)           :: pioIoSystem
    type(file_desc_t)     , intent(in)           :: ncid
    character(len=*)      , intent(in)           :: varname
    integer               , intent(out)          :: dim_idx(2,2)
    real(r8)              , pointer, intent(out) :: data2d(:,:)
    !
    type(var_desc_t)                             :: varid
    type(io_desc_t)                              :: iodescNCells
    integer                                      :: ndims, idim
    integer                                      :: begi, endi, begj, endj
    integer                                      :: numi, numj
    integer                                      :: i, j, count
    integer                                      :: vartype
    integer                                      :: ierr
    integer               , pointer              :: var_dim_ids(:), dim_glb(:), compdof(:)
    real                  , pointer              :: dataReal2d(:,:)
    character(len=32)                            :: subname = 'read_float_or_double_2d'

    call check_ret(PIO_inq_varid(ncid, varname, varid), subname)

    call check_ret(PIO_inq_varndims(ncid, varid, ndims), subname)
    if (ndims /= 2) then
       write(6,*)trim(varname),' is not a 2D variable.'
       call abort()
    end if

    allocate(var_dim_ids(ndims))
    allocate(dim_glb(ndims))

    call check_ret(PIO_inq_vardimid(ncid, varid, var_dim_ids), subname)

    do idim = 1, ndims
       call check_ret(PIO_inq_dimlen(ncid, var_dim_ids(idim), dim_glb(idim)), subname)
    end do

    call check_ret(PIO_inq_vartype(ncid, varid, vartype), subname)

    begi = domain%begi
    endi = domain%endi
    begj = domain%begj
    endj = domain%endj

    dim_idx(1,1) = begi; dim_idx(1,2) = endi
    dim_idx(2,1) = begj; dim_idx(2,2) = endj

    numi = domain%endi - domain%begi + 1
    numj = domain%endj - domain%begj + 1

    allocate(compdof(numi * numj))

    count = 0;
    do j = 1, numj
       do i = 1, numi
          count = count + 1
          compdof(count) = (domain%begi - 1) + i + (j-1)*dim_glb(1)
       end do
    end do

    call PIO_initdecomp(pioIoSystem, vartype, dim_glb, compdof, iodescNCells)

    allocate(data2d(begi:endi, begj:endj))

    if (vartype == PIO_REAL) then
       allocate(dataReal2d(begi:endi, begj:endj))
       call PIO_read_darray(ncid, varid, iodescNCells, dataReal2d, ierr)
       data2d = dble(dataReal2d)
       deallocate(dataReal2d)
    else
       call PIO_read_darray(ncid, varid, iodescNCells, data2d, ierr)
    end if

    ! Free up memory
    deallocate(compdof)
    call PIO_freedecomp(pioIoSystem, iodescNCells)
    
  end subroutine read_float_or_double_2d

  !-----------------------------------------------------------------------
  subroutine read_float_or_double_3d(domain, pioIoSystem, ncid, varname, start_id_for_dim3, dim_idx, data3d)

    use mkdomainPIOMod
    use pio

    implicit none

    type(domain_pio_type) , intent(in)           :: domain
    type(iosystem_desc_t) , intent(in)           :: pioIoSystem
    type(file_desc_t)     , intent(in)           :: ncid
    character(len=*)      , intent(in)           :: varname
    integer                                      :: start_id_for_dim3
    integer               , intent(out)          :: dim_idx(3,2)
    real(r8)              , pointer, intent(out) :: data3d(:,:,:)
    !
    type(var_desc_t)                             :: varid
    type(io_desc_t)                              :: iodescNCells
    integer                                      :: ndims, idim
    integer                                      :: begi, endi, begj, endj, begk, endk
    integer                                      :: numi, numj, numk
    integer                                      :: i, j, k, count
    integer                                      :: vartype
    integer                                      :: ierr
    integer               , pointer              :: var_dim_ids(:), dim_glb(:), compdof(:)
    real                  , pointer              :: dataReal3d(:,:,:)
    character(len=32)                            :: subname = 'read_float_or_double_3d'

    call check_ret(PIO_inq_varid(ncid, varname, varid), subname)

    call check_ret(PIO_inq_varndims(ncid, varid, ndims), subname)
    if (ndims /= 3) then
       write(6,*)trim(varname),' is not a 3D variable.'
       call abort()
    end if

    allocate(var_dim_ids(ndims))
    allocate(dim_glb(ndims))

    call check_ret(PIO_inq_vardimid(ncid, varid, var_dim_ids), subname)

    do idim = 1, ndims
       call check_ret(PIO_inq_dimlen(ncid, var_dim_ids(idim), dim_glb(idim)), subname)
    end do

    call check_ret(PIO_inq_vartype(ncid, varid, vartype), subname)

    begi = domain%begi
    endi = domain%endi
    begj = domain%begj
    endj = domain%endj
    begk = start_id_for_dim3
    endk = start_id_for_dim3 + dim_glb(3) - 1

    dim_idx(1,1) = begi; dim_idx(1,2) = endi
    dim_idx(2,1) = begj; dim_idx(2,2) = endj
    dim_idx(3,1) = begk; dim_idx(3,2) = endk

    numi = domain%endi - domain%begi + 1
    numj = domain%endj - domain%begj + 1
    numk = dim_glb(3)

    allocate(compdof(numi * numj * numk))

    count = 0;
    do k = 1, numk
       do j = 1, numj
          do i = 1, numi
             count = count + 1
             compdof(count) = (domain%begi - 1) + i + (j-1)*dim_glb(1) + (k-1)*dim_glb(1)*dim_glb(2)
          end do
       end do
    end do

    call PIO_initdecomp(pioIoSystem, vartype, dim_glb, compdof, iodescNCells)

    allocate(data3d(begi:endi, begj:endj, begk:endk))

    if (vartype == PIO_REAL) then
       allocate(dataReal3d(begi:endi, begj:endj, begk:endk))
       call PIO_read_darray(ncid, varid, iodescNCells, dataReal3d, ierr)
       data3d = dble(dataReal3d)
       deallocate(dataReal3d)
    else
       call PIO_read_darray(ncid, varid, iodescNCells, data3d, ierr)
    end if

    ! Free up memory
    deallocate(compdof)
    call PIO_freedecomp(pioIoSystem, iodescNCells)

  end subroutine read_float_or_double_3d

end module piofileutils
