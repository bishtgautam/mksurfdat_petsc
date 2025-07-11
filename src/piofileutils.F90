module piofileutils

  use shr_kind_mod, only : r8 => shr_kind_r8
  use mkdomainPIOMod
  use pio
  use spmdMod      , only : iam, npes, masterproc, mpicom

  implicit none

  public OpenFilePIO
  public CreateFilePIO
  public DefineDimPIO
  public read_integer_2d
  public read_float_or_double_2d
  public read_float_or_double_3d
  public read_float_or_double_4d
  public write_integer_1d
  public write_double_1d
  public write_double_2d
  public write_double_3d
  public write_double_4d
  public write_a_frame_of_double_3d
  public check_ret
  public DefineVarPIO_1d
  public DefineVarPIO_2d
  public DefineVarPIO_3d
  public DefineVarPIO_4d
  public get_dimlen
  public get_float_or_double_2d

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
  subroutine OpenFilePIO(fname, pioIoSystem, ncid, mode)

    use petsc
    use spmdMod, only : iam, npes
    use pio

    implicit none

    character (len=*)     :: fname
    type(file_desc_t)     :: ncid
    type(iosystem_desc_t) :: pioIoSystem
    integer, intent(in)   :: mode
    !
    integer               :: niotasks
    integer               :: numAggregator
    integer               :: stride
    integer               :: optBase
    integer               :: iotype
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

    call check_ret(PIO_openfile(pioIoSystem, ncid, iotype, trim(fname), mode), subname)

  end subroutine OpenFilePIO

  !-----------------------------------------------------------------------
  subroutine CreateFilePIO(fname, pioIoSystem, ncid)
    use petsc
    use spmdMod, only : iam, npes
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
    character(len= 32)    :: subname = 'CreateFilePIO'

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

    call check_ret(PIO_createfile(pioIoSystem, ncid, iotype, trim(fname), PIO_CLOBBER), subname)

  end subroutine CreateFilePIO

  !-----------------------------------------------------------------------
  subroutine DefineDimPIO(ncid, dimName, dimLen, dimID)

    implicit none

    type(file_desc_t)     , intent (in)  :: ncid
    character(len=*)      , intent (in)  :: dimName
    integer               , intent (in)  :: dimLen
    integer               , intent (out) :: dimID

    character(len= 32)    :: subname = 'DefineDimPIO'

    call check_ret(PIO_def_dim(ncid, dimName, dimLen, dimID), subname)

  end subroutine DefineDimPIO

  !-----------------------------------------------------------------------
  subroutine DefineVarPIO_1d(ncid, varName, xtype, dimID_1d, longName, units)

    use pio

    implicit none

    type(file_desc_t) , intent (in) :: ncid
    character(len=*)  , intent (in) :: varName
    integer           , intent (in) :: xtype
    integer           , intent (in) :: dimID_1d(1)
    character(len=*)  , optional    :: longName
    character(len=*)  , optional    :: units

    type(var_desc_t)                :: pioVar
    character(len=32)               :: subname = 'DefineVarPIO_1d'

    call check_ret(PIO_def_var(ncid, varName, xtype, dimID_1d, pioVar), subname)

    if (present(longName)) then
       call check_ret(PIO_put_att(ncid, pioVar, 'long_name', trim(longName)), subname)
    end if

    if (present(units)) then
       call check_ret(PIO_put_att(ncid, pioVar, 'units', trim(units)), subname)
    end if

  end subroutine DefineVarPIO_1d

  !-----------------------------------------------------------------------
  subroutine DefineVarPIO_2d(ncid, varName, xtype, dimIDs, longName, units)

    use pio

    implicit none

    type(file_desc_t) , intent (in) :: ncid
    character(len=*)  , intent (in) :: varName
    integer           , intent (in) :: xtype
    integer           , intent (in) :: dimIDs(2)
    character(len=*)  , optional    :: longName
    character(len=*)  , optional    :: units

    type(var_desc_t)                :: pioVar
    character(len=32)               :: subname = 'DefineVarPIO_1d'

    call check_ret(PIO_def_var(ncid, varName, xtype, dimIDs, pioVar), subname)

    if (present(longName)) then
       call check_ret(PIO_put_att(ncid, pioVar, 'long_name', trim(longName)), subname)
    end if

    if (present(units)) then
       call check_ret(PIO_put_att(ncid, pioVar, 'units', trim(units)), subname)
    end if

  end subroutine DefineVarPIO_2d

  !-----------------------------------------------------------------------
  subroutine DefineVarPIO_3d(ncid, varName, xtype, dimIDs, longName, units)

    use pio

    implicit none

    type(file_desc_t) , intent (in) :: ncid
    character(len=*)  , intent (in) :: varName
    integer           , intent (in) :: xtype
    integer           , intent (in) :: dimIDs(3)
    character(len=*)  , optional    :: longName
    character(len=*)  , optional    :: units

    type(var_desc_t)                :: pioVar
    character(len=32)               :: subname = 'DefineVarPIO_1d'

    call check_ret(PIO_def_var(ncid, varName, xtype, dimIDs, pioVar), subname)

    if (present(longName)) then
       call check_ret(PIO_put_att(ncid, pioVar, 'long_name', trim(longName)), subname)
    end if

    if (present(units)) then
       call check_ret(PIO_put_att(ncid, pioVar, 'units', trim(units)), subname)
    end if

  end subroutine DefineVarPIO_3d

  !-----------------------------------------------------------------------
  subroutine DefineVarPIO_4d(ncid, varName, xtype, dimIDs, longName, units)

    use pio

    implicit none

    type(file_desc_t) , intent (in) :: ncid
    character(len=*)  , intent (in) :: varName
    integer           , intent (in) :: xtype
    integer           , intent (in) :: dimIDs(4)
    character(len=*)  , optional    :: longName
    character(len=*)  , optional    :: units

    type(var_desc_t)                :: pioVar
    character(len=32)               :: subname = 'DefineVarPIO_1d'

    call check_ret(PIO_def_var(ncid, varName, xtype, dimIDs, pioVar), subname)

    if (present(longName)) then
       call check_ret(PIO_put_att(ncid, pioVar, 'long_name', trim(longName)), subname)
    end if

    if (present(units)) then
       call check_ret(PIO_put_att(ncid, pioVar, 'units', trim(units)), subname)
    end if

  end subroutine DefineVarPIO_4d

  !-----------------------------------------------------------------------
  subroutine read_integer_2d(domain, pioIoSystem, ncid, varname, dim_idx, vec_row_indices, data2d)

    use mkdomainPIOMod
    use pio

    implicit none

    type(domain_pio_type) , intent(in)           :: domain
    type(iosystem_desc_t) , intent(in)           :: pioIoSystem
    type(file_desc_t)     , intent(in)           :: ncid
    character(len=*)      , intent(in)           :: varname
    integer               , intent(out)          :: dim_idx(2,2)
    integer               , pointer, intent(out) :: vec_row_indices(:)
    integer               , pointer, intent(out) :: data2d(:,:)
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
    allocate(vec_row_indices(numi * numj))

    count = 0;
    do j = 1, numj
       do i = 1, numi
          count = count + 1
          compdof(count) = (domain%begi - 1) + i + (j-1)*dim_glb(1)
          vec_row_indices(count) = compdof(count) - 1
       end do
    end do

    call PIO_initdecomp(pioIoSystem, vartype, dim_glb, compdof, iodescNCells)

    allocate(data2d(begi:endi, begj:endj))

    if (vartype == PIO_INT) then
      call PIO_read_darray(ncid, varid, iodescNCells, data2d, ierr)
    elseif (vartype == PIO_DOUBLE .or. vartype == PIO_REAL) then
      allocate(dataReal2d(begi:endi, begj:endj))
      call PIO_read_darray(ncid, varid, iodescNCells, dataReal2d, ierr)
      data2d = int(dataReal2d)
      deallocate(dataReal2d)
    else
      write(*,*) 'Unknown vartype = ', vartype
      call exit(0)
    endif

    ! Free up memory
    deallocate(compdof)
    call PIO_freedecomp(pioIoSystem, iodescNCells)
    
  end subroutine read_integer_2d

  !-----------------------------------------------------------------------
  subroutine read_float_or_double_2d(domain, pioIoSystem, ncid, varname, dim_idx, vec_row_indices, data2d)

   use mkdomainPIOMod
   use pio

   implicit none

   type(domain_pio_type) , intent(in)           :: domain
   type(iosystem_desc_t) , intent(in)           :: pioIoSystem
   type(file_desc_t)     , intent(in)           :: ncid
   character(len=*)      , intent(in)           :: varname
   integer               , intent(out)          :: dim_idx(2,2)
   integer               , pointer, intent(out) :: vec_row_indices(:)
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
   allocate(vec_row_indices(numi * numj))

   count = 0;
   do j = 1, numj
      do i = 1, numi
         count = count + 1
         compdof(count) = (domain%begi - 1) + i + (j-1)*dim_glb(1)
         vec_row_indices(count) = compdof(count) - 1
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
 subroutine get_float_or_double_2d(ncid, varname, dim_idx, data2d)

   use mkdomainPIOMod
   use pio
   use spmdMod

   implicit none

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
   character(len=32)                            :: subname = 'get_float_or_double_2d'

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

   begi = 1; endi = dim_glb(1)
   begj = 1; endj = dim_glb(2)

   dim_idx(1,1) = begi; dim_idx(1,2) = endi
   dim_idx(2,1) = begj; dim_idx(2,2) = endj

   numi = dim_glb(1)
   numj = dim_glb(2)

   allocate(data2d(begi:endi, begj:endj))

   if (vartype == PIO_REAL) then
      allocate(dataReal2d(begi:endi, begj:endj))
      ierr = PIO_get_var(ncid, varid, dataReal2d)
      data2d = dble(dataReal2d)
      deallocate(dataReal2d)
   else
      ierr = PIO_get_var(ncid, varid, data2d)
   end if


 end subroutine get_float_or_double_2d

 !-----------------------------------------------------------------------
  subroutine read_float_or_double_3d(domain, pioIoSystem, ncid, varname, start_id_for_dim3, dim_idx, vec_row_indices, data3d)

    use mkdomainPIOMod
    use pio

    implicit none

    type(domain_pio_type) , intent(in)           :: domain
    type(iosystem_desc_t) , intent(in)           :: pioIoSystem
    type(file_desc_t)     , intent(in)           :: ncid
    character(len=*)      , intent(in)           :: varname
    integer                                      :: start_id_for_dim3
    integer               , intent(out)          :: dim_idx(3,2)
    integer               , pointer, intent(out) :: vec_row_indices(:)
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

    allocate(vec_row_indices(numi * numj))
    count = 0;
    do k = 1, numk
       do j = 1, numj
          do i = 1, numi
             count = count + 1
             compdof(count) = (domain%begi - 1) + i + (j-1)*dim_glb(1) + (k-1)*dim_glb(1)*dim_glb(2)
             if (k == 1) then
                vec_row_indices(count) = compdof(count) - 1
             end if
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

 !-----------------------------------------------------------------------
  subroutine read_float_or_double_4d(domain, pioIoSystem, ncid, varname, dim_idx, vec_row_indices, data4d)

    use mkdomainPIOMod
    use pio

    implicit none

    type(domain_pio_type) , intent(in)           :: domain
    type(iosystem_desc_t) , intent(in)           :: pioIoSystem
    type(file_desc_t)     , intent(in)           :: ncid
    character(len=*)      , intent(in)           :: varname
    integer               , intent(out)          :: dim_idx(4,2)
    integer               , pointer, intent(out) :: vec_row_indices(:)
    real(r8)              , pointer, intent(out) :: data4d(:,:,:,:)
    !
    type(var_desc_t)                             :: varid
    type(io_desc_t)                              :: iodescNCells
    integer                                      :: ndims, idim
    integer                                      :: begi, endi, begj, endj, begk, endk, begl, endl
    integer                                      :: numi, numj, numk, numl
    integer                                      :: i, j, k, l, count
    integer                                      :: vartype
    integer                                      :: ierr
    integer               , pointer              :: var_dim_ids(:), dim_glb(:), compdof(:)
    real                  , pointer              :: dataReal4d(:,:,:,:)
    character(len=32)                            :: subname = 'read_float_or_double_4d'

    call check_ret(PIO_inq_varid(ncid, varname, varid), subname)

    call check_ret(PIO_inq_varndims(ncid, varid, ndims), subname)
    if (ndims /= 4) then
       write(6,*)trim(varname),' is not a 4D variable.'
       call abort()
    end if

    allocate(var_dim_ids(ndims))
    allocate(dim_glb(ndims))

    call check_ret(PIO_inq_vardimid(ncid, varid, var_dim_ids), subname)

    do idim = 1, ndims
       call check_ret(PIO_inq_dimlen(ncid, var_dim_ids(idim), dim_glb(idim)), subname)
    end do

    call check_ret(PIO_inq_vartype(ncid, varid, vartype), subname)

    begi = domain%begi; endi = domain%endi
    begj = domain%begj; endj = domain%endj
    begk = 1;           endk = dim_glb(3)
    begl = 1;           endl = dim_glb(4)

    dim_idx(1,1) = begi; dim_idx(1,2) = endi
    dim_idx(2,1) = begj; dim_idx(2,2) = endj
    dim_idx(3,1) = begk; dim_idx(3,2) = endk
    dim_idx(4,1) = begl; dim_idx(4,2) = endl

    numi = domain%endi - domain%begi + 1
    numj = domain%endj - domain%begj + 1
    numk = dim_glb(3)
    numl = dim_glb(4)

    allocate(compdof        (numi * numj * numk * numl))
    allocate(vec_row_indices(numi * numj))

    count = 0;
    do l = 1, numl
       do k = 1, numk
          do j = 1, numj
             do i = 1, numi
                count = count + 1
                compdof(count) = (domain%begi - 1) + i + (j-1)*dim_glb(1) + (k-1)*dim_glb(1)*dim_glb(2) + (l-1)*dim_glb(1)*dim_glb(2)*dim_glb(3)
                if (k == 1 .and. l == 1) then
                   vec_row_indices(count) = compdof(count) - 1
                end if
             end do
          end do
       end do
    end do

    call PIO_initdecomp(pioIoSystem, vartype, dim_glb, compdof, iodescNCells)

    allocate(data4d(begi:endi, begj:endj, begk:endk, begl:endl))

    if (vartype == PIO_REAL) then
       allocate(dataReal4d(begi:endi, begj:endj, begk:endk, begl:endl))
       call PIO_read_darray(ncid, varid, iodescNCells, dataReal4d, ierr)
       data4d = dble(dataReal4d)
       deallocate(dataReal4d)
    else
       call PIO_read_darray(ncid, varid, iodescNCells, data4d, ierr)
    end if

    ! Free up memory
    deallocate(compdof)
    call PIO_freedecomp(pioIoSystem, iodescNCells)

  end subroutine read_float_or_double_4d

  !-----------------------------------------------------------------------
  subroutine write_integer_1d(ncid, iodesc, varname, data1d)
    !
    use pio

    implicit none

    type(file_desc_t)     , intent(inout)        :: ncid
    type(io_desc_t)       , intent (in)          :: iodesc
    character(len=*)      , intent(in)           :: varname
    integer               , pointer, intent(in)  :: data1d(:)

    type(var_desc_t)      :: varid
    integer               :: ier
    character(len=32)                            :: subname = 'write_double_2d'

    call check_ret(PIO_inq_varid(ncid, varname, varid), subname)
    call PIO_write_darray(ncid, varid, iodesc, data1d, ier)
    call PIO_syncfile(ncid)

  end subroutine write_integer_1d

  !-----------------------------------------------------------------------
  subroutine write_double_1d(ncid, iodesc, varname, data1d)
    !
    use pio

    implicit none

    type(file_desc_t)     , intent(inout)        :: ncid
    type(io_desc_t)       , intent (in)          :: iodesc
    character(len=*)      , intent(in)           :: varname
    real(r8)              , pointer, intent(in)  :: data1d(:)

    type(var_desc_t)      :: varid
    integer               :: ier
    character(len=32)                            :: subname = 'write_double_2d'

    call check_ret(PIO_inq_varid(ncid, varname, varid), subname)
    call PIO_write_darray(ncid, varid, iodesc, data1d, ier)
    call PIO_syncfile(ncid)

  end subroutine write_double_1d

  !-----------------------------------------------------------------------
  subroutine write_double_2d(ncid, iodesc, varname, data2d)
    !
    use pio

    implicit none

    type(file_desc_t)     , intent(inout)        :: ncid
    type(io_desc_t)       , intent (in)          :: iodesc
    character(len=*)      , intent(in)           :: varname
    real(r8)              , pointer, intent(in)  :: data2d(:,:)

    type(var_desc_t)      :: varid
    integer               :: ier
    character(len=32)                            :: subname = 'write_double_2d'

    call check_ret(PIO_inq_varid(ncid, varname, varid), subname)
    call PIO_write_darray(ncid, varid, iodesc, data2d, ier)
    call PIO_syncfile(ncid)

  end subroutine write_double_2d

  !-----------------------------------------------------------------------
  subroutine write_double_3d(ncid, iodesc, varname, data3d)
    !
    use pio

    implicit none

    type(file_desc_t)     , intent(inout)        :: ncid
    type(io_desc_t)       , intent (in)          :: iodesc
    character(len=*)      , intent(in)           :: varname
    real(r8)              , pointer, intent(in)  :: data3d(:,:,:)

    type(var_desc_t)      :: varid
    integer               :: ier
    character(len=32)                            :: subname = 'write_double_2d'

    call check_ret(PIO_inq_varid(ncid, varname, varid), subname)
    call PIO_write_darray(ncid, varid, iodesc, data3d, ier)
    call PIO_syncfile(ncid)

  end subroutine write_double_3d

  !-----------------------------------------------------------------------
  subroutine write_a_frame_of_double_3d(ncid, iodesc, varname, frame_id, data3d)
   !
   use pio

   implicit none

   type(file_desc_t)     , intent(inout)        :: ncid
   type(io_desc_t)       , intent (in)          :: iodesc
   character(len=*)      , intent(in)           :: varname
   integer               , intent(in)           :: frame_id
   real(r8)              , pointer, intent(in)  :: data3d(:,:,:)

   type(var_desc_t)                             :: varid
   integer                                      :: ier
   character(len=32)                            :: subname = 'write_a_frame_of_double_3d'

   call check_ret(PIO_inq_varid(ncid, varname, varid), subname)
   call PIO_setframe(ncid, varid, int(frame_id,kind=PIO_OFFSET_KIND))
   call PIO_write_darray(ncid, varid, iodesc, data3d(:,:,frame_id), ier)
   call PIO_syncfile(ncid)

 end subroutine write_a_frame_of_double_3d

 !-----------------------------------------------------------------------
  subroutine write_double_4d(ncid, iodesc, varname, data4d)
    !
    use pio

    implicit none

    type(file_desc_t)     , intent(inout)        :: ncid
    type(io_desc_t)       , intent (in)          :: iodesc
    character(len=*)      , intent(in)           :: varname
    real(r8)              , pointer, intent(in)  :: data4d(:,:,:,:)

    type(var_desc_t)      :: varid
    integer               :: ier
    character(len=32)                            :: subname = 'write_double_2d'

    call check_ret(PIO_inq_varid(ncid, varname, varid), subname)
    call PIO_write_darray(ncid, varid, iodesc, data4d, ier)
    call PIO_syncfile(ncid)

  end subroutine write_double_4d

  !-----------------------------------------------------------------------
  subroutine get_dimlen(fname, dimname, dimlen)
    !
    use pio
    !
    implicit none
    !
    character(len=*) , intent(in)  :: fname
    character(len=*) , intent(in)  :: dimname
    integer          , intent(out) :: dimlen
    !
    type(file_desc_t)              :: ncid
    type(iosystem_desc_t)          :: pioIoSystem
    integer                        :: dimid
    integer                        :: ierr
    character(len=32)              :: subname = 'get_dimlen'

    ! open the netcdf file
    call OpenFilePIO(fname, pioIoSystem, ncid, PIO_NOWRITE)

    ! get dimid
    call check_ret(PIO_inq_dimid(ncid, dimname, dimid), subname)

    ! get dimlen
    call check_ret(PIO_inq_dimlen(ncid, dimid, dimlen), subname)

    ! clean up
    call PIO_closefile(ncid)
    call PIO_finalize(pioIoSystem, ierr)

  end subroutine get_dimlen

end module piofileutils
