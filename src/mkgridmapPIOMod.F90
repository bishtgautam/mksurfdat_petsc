module mkgridmapPIOMod

#include "petsc/finclude/petscsys.h"
#include <petsc/finclude/petscmat.h>
#include <petsc/finclude/petscvec.h>

  use petsc
  use pio
  use shr_kind_mod , only : r8 => shr_kind_r8
  use mkvarpar     , only : re
  use nanMod       , only : nan, bigint
  use spmdMod      , only : iam, npes, masterproc, mpicom
  use piofileutils , only : check_ret

  implicit none
  private

  type grid_type
     real(r8), pointer :: xc(:)
     real(r8), pointer :: yc(:)
     real(r8), pointer :: area(:)
     real(r8), pointer :: frac(:)
     integer , pointer :: mask(:)
  end type grid_type

  type gridmap_dim_type
     integer          :: nglb
     integer          :: nloc
     integer          :: begd, endd
     integer, pointer :: compdof(:)
  end type gridmap_dim_type

  type gridmap_pio_type
     character(len=32) :: set
     character(len=32) :: name
     
     type(grid_type) :: src
     type(grid_type) :: dst

     type(gridmap_dim_type) :: dim_na
     type(gridmap_dim_type) :: dim_nb
     type(gridmap_dim_type) :: dim_ns

     integer, pointer  :: row(:)  ! corresponding row index
     integer, pointer  :: col(:)  ! corresponding column index
     real(r8), pointer :: wovr(:) ! wt of overlapped input cell
     
     ! PETSc
     Mat :: map_mat, map_frac_mat
     Vec :: src_vec
     Vec :: dst_vec

  end type gridmap_pio_type

  public :: gridmap_pio_type

  public :: gridmap_mapread_pio
  public :: gridmap_clean_pio
  public :: gridmap_areaave_pio
  public :: gridmap_dominant_value_pio

  interface gridmap_areaave_pio
     module procedure gridmap_areaave_default_pio
     module procedure gridmap_areaave_srcmask_pio
  end interface gridmap_areaave_pio

  character(len=32), parameter :: isSet = "gridmap_IsSet"

contains

  !------------------------------------------------------------------------------
  subroutine gridmap_grid_setup(pioIoSystem, ncid, dim, grid_suffix, grid)

    implicit none

    type(iosystem_desc_t) , intent(in)    :: pioIoSystem
    type(file_desc_t)     , intent(in)    :: ncid
    type(gridmap_dim_type) , intent(in)    :: dim
    character(len=2)       , intent(in)    :: grid_suffix
    type(grid_type)        , intent(inout) :: grid

    character(len=24)                      :: var_name

    ! allocate memory
    allocate(grid%xc(dim%begd:dim%endd))
    allocate(grid%yc(dim%begd:dim%endd))
    allocate(grid%area(dim%begd:dim%endd))
    allocate(grid%frac(dim%begd:dim%endd))
    allocate(grid%mask(dim%begd:dim%endd))

    var_name = 'xc' // grid_suffix
    call gridmap_read_float_or_double_1d(pioIoSystem, ncid, var_name, dim%begd, dim%endd, dim%nglb, dim%compdof, grid%xc)

    var_name = 'yc' // grid_suffix
    call gridmap_read_float_or_double_1d(pioIoSystem, ncid, var_name, dim%begd, dim%endd, dim%nglb, dim%compdof, grid%yc)

    var_name = 'area' // grid_suffix
    call gridmap_read_float_or_double_1d(pioIoSystem, ncid, var_name, dim%begd, dim%endd, dim%nglb, dim%compdof, grid%area)

    var_name = 'frac' // grid_suffix
    call gridmap_read_float_or_double_1d(pioIoSystem, ncid, var_name, dim%begd, dim%endd, dim%nglb, dim%compdof, grid%frac)

    var_name = 'mask' // grid_suffix
    call gridmap_read_integer_1d(pioIoSystem, ncid, var_name, dim%begd, dim%endd, dim%nglb, dim%compdof, grid%mask)

  end subroutine gridmap_grid_setup

  !------------------------------------------------------------------------------
  subroutine gridmap_grid_clean(grid)

    implicit none

    type(grid_type), intent(inout) :: grid

    deallocate(grid%xc)
    deallocate(grid%yc)
    deallocate(grid%area)
    deallocate(grid%frac)
    deallocate(grid%mask)

  end subroutine gridmap_grid_clean

  !------------------------------------------------------------------------------
  subroutine gridmap_read_dimsize_pio(pioIoSystem, ncid, dim_name, ndim)

    implicit none

    type(iosystem_desc_t) , intent(in)    :: pioIoSystem
    type(file_desc_t)     , intent(in)    :: ncid
    character(len=*)      , intent(in)    :: dim_name
    integer               , intent(out)   :: ndim

    integer               :: dimid

    character(*), parameter :: subname = 'gridmap_read_dimsize_pio'

    call check_ret(PIO_inq_dimid(ncid, dim_name, dimid), subname)
    call check_ret(PIO_inq_dimid  (ncid, dim_name, dimid), subname)
    call check_ret(PIO_inq_dimlen (ncid, dimid, ndim), subname)

  end subroutine gridmap_read_dimsize_pio

  !------------------------------------------------------------------------------
  subroutine gridmap_compute_dimsize_local(iam, npes, dim_glb, dim_loc)

    implicit none

    integer, intent(in)  :: iam, npes, dim_glb
    integer, intent(out) :: dim_loc

    dim_loc = dim_glb/npes
    if (iam < dim_glb - dim_loc * npes) then
       dim_loc = dim_loc + 1
    end if

  end subroutine gridmap_compute_dimsize_local

  !------------------------------------------------------------------------------
  subroutine gridmap_dim_setup(pioIoSystem, ncid, dim_name, dim)

    implicit none

    type(iosystem_desc_t) , intent(in)    :: pioIoSystem
    type(file_desc_t)     , intent(in)    :: ncid
    character(len=*)      , intent(in)    :: dim_name
    type(gridmap_dim_type), intent(out)   :: dim

    integer                               :: i
    integer                               :: ier
    
    call gridmap_read_dimsize_pio(pioIoSystem, ncid, dim_name, dim%nglb)
    call gridmap_compute_dimsize_local(iam, npes, dim%nglb, dim%nloc)

    call MPI_Scan(dim%nloc, dim%endd, 1, MPI_INTEGER, MPI_SUM, mpicom, ier)
    dim%begd = dim%endd - dim%nloc + 1

    allocate(dim%compdof(dim%begd:dim%endd))

    dim%compdof(dim%begd:dim%endd) = (/(i, i=dim%begd, dim%endd)/)

  end subroutine gridmap_dim_setup

  !------------------------------------------------------------------------------
  subroutine gridmap_read_float_or_double_1d(pioIoSystem, ncid, var_name, begd, endd, nglb, compdof, dataBuffer1d)
    implicit none
    type(iosystem_desc_t) , intent(in)           :: pioIoSystem
    type(file_desc_t)     , intent(in)           :: ncid
    integer               , intent(in)           :: begd, endd
    character(len=*)      , intent(in)           :: var_name
    integer               , intent(in)           :: nglb
    integer               , intent(in), pointer  :: compdof(:)
    real(r8)              , intent(out), pointer :: dataBuffer1d(:)

    integer, pointer   :: dim_glb(:)
    real,     pointer  :: dataBuffer1dFloat(:)
    type(var_desc_t)   :: varid
    type(io_desc_t)    :: iodescNCells
    integer            :: vartype
    integer            :: ier
    character(len= 32) :: subname = 'read_float_or_double_1d'

    ! Inquiry the variable ID and variable type
    call check_ret(PIO_inq_varid(ncid, var_name, varid), subname)
    call check_ret(PIO_inq_vartype(ncid, varid, vartype), subname)

    allocate(dim_glb(1))
    dim_glb(1) = nglb
    call PIO_initdecomp(pioIoSystem, vartype, dim_glb, compdof, iodescNCells)
    deallocate(dim_glb)

    ! Read the variable
    if (vartype == PIO_DOUBLE) then
       call PIO_read_darray(ncid, varid, iodescNCells, dataBuffer1d, ier)
    else
       allocate(dataBuffer1dFloat(begd:endd))
       call PIO_read_darray(ncid, varid, iodescNCells, dataBuffer1dFloat, ier)
       dataBuffer1d = dble(dataBuffer1dFloat)
       deallocate(dataBuffer1dFloat)
    end if

    ! Delete decomposition
    call PIO_freedecomp(pioIoSystem, iodescNCells)

  end subroutine gridmap_read_float_or_double_1d

  !------------------------------------------------------------------------------
  subroutine gridmap_read_integer_1d(pioIoSystem, ncid, var_name, begd, endd, nglb, compdof, dataBuffer1d)
    implicit none
    type(iosystem_desc_t) , intent(in)           :: pioIoSystem
    type(file_desc_t)     , intent(in)           :: ncid
    integer               , intent(in)           :: begd, endd
    character(len=*)      , intent(in)           :: var_name
    integer               , intent(in)           :: nglb
    integer               , intent(in), pointer  :: compdof(:)
    integer               , intent(out), pointer :: dataBuffer1d(:)

    integer, pointer   :: dim_glb(:)
    type(var_desc_t)   :: varid
    type(io_desc_t)    :: iodescNCells
    integer            :: vartype
    integer            :: ier
    character(len= 32) :: subname = 'read_float_or_double_1d'

    ! Inquiry the variable ID and variable type
    call check_ret(PIO_inq_varid(ncid, var_name, varid), subname)
    call check_ret(PIO_inq_vartype(ncid, varid, vartype), subname)

    allocate(dim_glb(1))
    dim_glb(1) = nglb
    call PIO_initdecomp(pioIoSystem, vartype, dim_glb, compdof, iodescNCells)
    deallocate(dim_glb)

    ! Read the variable
    call PIO_read_darray(ncid, varid, iodescNCells, dataBuffer1d, ier)

    ! Delete decomposition
    call PIO_freedecomp(pioIoSystem, iodescNCells)

  end subroutine gridmap_read_integer_1d

  !------------------------------------------------------------------------------
  subroutine gridmap_mapread_pio(gridmap_pio, fileName)

    implicit none

    type(gridmap_pio_type), intent(out) :: gridmap_pio   ! mapping data
    character(len=*)      , intent(in)  :: filename  ! netCDF file to read

    integer                :: niotasks
    integer                :: numAggregator
    integer                :: stride
    integer                :: optBase
    integer                :: iotype
    integer                :: no, ii, n

    type(file_desc_t)      :: ncid
    type(iosystem_desc_t)  :: pioIoSystem

    Mat                    :: temp_mat
    PetscReal, pointer     :: src_p(:), dst_p(:)
    PetscErrorCode         :: ierr
    
    character(*),parameter :: subname = 'gridmap_mapread_pio'
    character(*),parameter :: F00 = '("(gridmap_mapread_pio) ",4a)'
    character(*),parameter :: F01 = '("(gridmap_mapread_pio) ",2(a,i7))'

    if (masterproc) write(6,F00) "reading mapping matrix data using PIO..."

    ! open & read the file
    if (masterproc) write(6,F00) "* file name                  : ",trim(fileName)

    stride        = 1
    numAggregator = 0
    iotype        = PIO_iotype_pnetcdf
    niotasks      = npes

    if (npes == 1) then
       optBase = 0
    else
       optBase = 1
    end if

    ! Initialize the PIO System
    call PIO_init(iam, mpicom, niotasks, numAggregator, stride, PIO_rearr_subset, pioIoSystem, optBase)

    ! Open the file
    call check_ret(PIO_openfile(pioIoSystem, ncid, iotype, trim(fileName), PIO_NOWRITE), subname)

    ! Read information about dimensions
    call gridmap_dim_setup(pioIoSystem, ncid, 'n_a', gridmap_pio%dim_na)
    call gridmap_dim_setup(pioIoSystem, ncid, 'n_b', gridmap_pio%dim_nb)
    call gridmap_dim_setup(pioIoSystem, ncid, 'n_s', gridmap_pio%dim_ns)

    call gridmap_grid_setup(pioIoSystem, ncid, gridmap_pio%dim_na, '_a', gridmap_pio%src)
    call gridmap_grid_setup(pioIoSystem, ncid, gridmap_pio%dim_nb, '_b', gridmap_pio%dst)

    allocate(gridmap_pio%wovr(gridmap_pio%dim_ns%begd:gridmap_pio%dim_ns%endd))
    allocate(gridmap_pio%row(gridmap_pio%dim_ns%begd:gridmap_pio%dim_ns%endd))
    allocate(gridmap_pio%col(gridmap_pio%dim_ns%begd:gridmap_pio%dim_ns%endd))

    call gridmap_read_float_or_double_1d(pioIoSystem, ncid, 'S', gridmap_pio%dim_ns%begd, &
         gridmap_pio%dim_ns%endd, gridmap_pio%dim_ns%nglb, gridmap_pio%dim_ns%compdof, gridmap_pio%wovr)

    call gridmap_read_integer_1d(pioIoSystem, ncid, 'row', gridmap_pio%dim_ns%begd, &
         gridmap_pio%dim_ns%endd, gridmap_pio%dim_ns%nglb, gridmap_pio%dim_ns%compdof, gridmap_pio%row)

    call gridmap_read_integer_1d(pioIoSystem, ncid, 'col', gridmap_pio%dim_ns%begd, &
         gridmap_pio%dim_ns%endd, gridmap_pio%dim_ns%nglb, gridmap_pio%dim_ns%compdof, gridmap_pio%col)

    ! Create PETSc-related data structure
    PetscCallA(VecCreate(PETSC_COMM_WORLD, gridmap_pio%src_vec, ierr))
    PetscCallA(VecSetSizes(gridmap_pio%src_vec, PETSC_DECIDE, gridmap_pio%dim_na%nglb, ierr))
    PetscCallA(VecSetFromOptions(gridmap_pio%src_vec, ierr))

    PetscCallA(VecCreate(PETSC_COMM_WORLD, gridmap_pio%dst_vec, ierr))
    PetscCallA(VecSetSizes(gridmap_pio%dst_vec, PETSC_DECIDE, gridmap_pio%dim_nb%nglb, ierr))
    PetscCallA(VecSetFromOptions(gridmap_pio%dst_vec, ierr))

    !
    ! map_mat: Create sparse matrix based on 'row', 'col', and 'S' values
    !
    PetscCallA(MatCreate(PETSC_COMM_WORLD, gridmap_pio%map_mat, ierr))
    PetscCall(MatSetsizes(gridmap_pio%map_mat, PETSC_DECIDE, PETSC_DECIDE, gridmap_pio%dim_nb%nglb, gridmap_pio%dim_na%nglb, ierr))
    PetscCall(MatSetFromOptions(gridmap_pio%map_mat, ierr))

    do n = gridmap_pio%dim_ns%begd, gridmap_pio%dim_ns%endd
       PetscCallA(MatSetValues(gridmap_pio%map_mat, 1, [gridmap_pio%row(n) - 1], 1, [gridmap_pio%col(n) - 1], [gridmap_pio%wovr(n)], INSERT_VALUES, ierr))
    end do
    PetscCallA(MatAssemblyBegin(gridmap_pio%map_mat, MAT_FINAL_ASSEMBLY, ierr))
    PetscCallA(MatAssemblyEnd(gridmap_pio%map_mat, MAT_FINAL_ASSEMBLY, ierr))

    !
    ! map_frac_mat: Matrix in which each row of map_mat is normalized with frac_b(row)
    !
    PetscCallA(MatCreate(PETSC_COMM_WORLD, temp_mat, ierr))
    PetscCall(MatSetsizes(temp_mat, PETSC_DECIDE, PETSC_DECIDE, gridmap_pio%dim_nb%nglb, gridmap_pio%dim_nb%nglb, ierr))
    PetscCall(MatSetFromOptions(temp_mat, ierr))

    PetscCallA(MatCreate(PETSC_COMM_WORLD, gridmap_pio%map_frac_mat, ierr))
    PetscCall(MatSetsizes(gridmap_pio%map_frac_mat, PETSC_DECIDE, PETSC_DECIDE, gridmap_pio%dim_nb%nglb, gridmap_pio%dim_na%nglb, ierr))
    PetscCall(MatSetFromOptions(gridmap_pio%map_frac_mat, ierr))

    ! Put 1.0 on src_vec
    PetscCallA(VecGetArray(gridmap_pio%src_vec, src_p, ierr))
    do ii = 1, gridmap_pio%dim_na%nloc
       src_p(ii) = 1._r8
    end do
    PetscCallA(VecRestoreArray(gridmap_pio%src_vec, src_p, ierr))

    ! dst_vec = map_mat * src_vec
    PetscCallA(MatMult(gridmap_pio%map_mat, gridmap_pio%src_vec, gridmap_pio%dst_vec, ierr))

    ! Now, create a temporary matrix that has 1/dst_vec(:) as the diagonal
    PetscCallA(VecGetArray(gridmap_pio%dst_vec, dst_p, ierr))
    do no = gridmap_pio%dim_nb%begd, gridmap_pio%dim_nb%endd
       ii = no - gridmap_pio%dim_nb%begd + 1
       if (dst_p(ii) > 0._r8) then
          PetscCallA(MatSetValues(temp_mat, 1, [no - 1], 1, [no - 1], [1._r8/dst_p(ii)], INSERT_VALUES, ierr))
       end if
    end do
    PetscCallA(VecRestoreArray(gridmap_pio%dst_vec, dst_p, ierr))

    PetscCallA(MatAssemblyBegin(temp_mat, MAT_FINAL_ASSEMBLY, ierr))
    PetscCallA(MatAssemblyEnd(temp_mat, MAT_FINAL_ASSEMBLY, ierr))

    ! map_frac_mat = temp_mat * map_mat
    PetscCallA(MatMatMult(temp_mat, gridmap_pio%map_mat, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, gridmap_pio%map_frac_mat, ierr))

    PetscCallA(MatDestroy(temp_mat, ierr))

    ! map is set
    gridmap_pio%set = isSet

  end subroutine gridmap_mapread_pio

  !------------------------------------------------------------------------------
  subroutine gridmap_areaave_default_pio(gridmap_pio, nrow, row_indices, src_array, dst_array, nodata)

    implicit none

    type(gridmap_pio_type)           , intent(in)  :: gridmap_pio
    integer                          , intent(in)  :: nrow
    integer                , pointer , intent (in) :: row_indices(:)
    real(r8)                         , intent(in)  :: src_array(:)
    real(r8)                         , intent(out) :: dst_array(:)
    real(r8)                         , intent(in)  :: nodata               ! value to apply where there are no input data

    integer                              :: no, idx
    PetscReal              , pointer     :: dst_p(:)
    PetscErrorCode                       :: ierr
    character(*)           , parameter   :: subName = '(gridmap_areaave_default_pio) '

    call gridmap_checkifset_pio(gridmap_pio, subname)

    ! fill the source Vec
    PetscCallA(VecSetValues(gridmap_pio%src_vec, nrow, row_indices, src_array, INSERT_VALUES, ierr))
    PetscCallA(VecAssemblyBegin(gridmap_pio%src_vec, ierr))
    PetscCallA(VecAssemblyEnd(gridmap_pio%src_vec, ierr))

    ! do matrix-vec multiplication
    PetscCallA(MatMult(gridmap_pio%map_frac_mat, gridmap_pio%src_vec, gridmap_pio%dst_vec, ierr))

    ! fill the destination array
    PetscCallA(VecGetArray(gridmap_pio%dst_vec, dst_p, ierr))
    do no = 1, gridmap_pio%dim_nb%nloc
       idx = gridmap_pio%dim_nb%begd + no - 1
       if (gridmap_pio%dst%frac(idx) <= 0._r8) then
          dst_array(no) = nodata
       else
          dst_array(no) = dst_p(no)
       end if
    end do
    PetscCallA(VecRestoreArray(gridmap_pio%dst_vec, dst_p, ierr))

  end subroutine gridmap_areaave_default_pio

  !------------------------------------------------------------------------------
  subroutine gridmap_areaave_srcmask_pio(gridmap_pio, nrow, row_indices, src_array, dst_array, nodata, mask_src)

    implicit none

    type(gridmap_pio_type)           , intent(in)  :: gridmap_pio
    integer                          , intent(in)  :: nrow
    integer                , pointer , intent(in)  :: row_indices(:)
    real(r8)                         , intent(in)  :: src_array(:)
    real(r8)                         , intent(out) :: dst_array(:)
    real(r8)                         , intent(in)  :: nodata               ! value to apply where there are no input data
    real(r8)                         , intent(in)  :: mask_src(:)

    integer                              :: no, idx
    PetscReal              , pointer     :: dst_p(:), tmp_p(:)
    Vec                                  :: tmp_src_vec, tmp_dst_vec
    PetscErrorCode                       :: ierr
    character(*)           , parameter   :: subName = '(gridmap_areaave_srcmask_pio) '

    call gridmap_checkifset_pio(gridmap_pio, subname)

    ! create a temporary vectors
    PetscCallA(VecDuplicate(gridmap_pio%src_vec, tmp_src_vec, ierr))
    PetscCallA(VecDuplicate(gridmap_pio%dst_vec, tmp_dst_vec, ierr))

    ! fill the source and temporary Vec
    PetscCallA(VecSetValues(gridmap_pio%src_vec, nrow, row_indices, src_array, INSERT_VALUES, ierr))
    PetscCallA(VecAssemblyBegin(gridmap_pio%src_vec, ierr))
    PetscCallA(VecAssemblyEnd(gridmap_pio%src_vec, ierr))

    PetscCallA(VecSetValues(tmp_src_vec, nrow, row_indices, mask_src, INSERT_VALUES, ierr))
    PetscCallA(VecAssemblyBegin(tmp_src_vec, ierr))
    PetscCallA(VecAssemblyEnd(tmp_src_vec, ierr))

    ! do matrix-vec multiplication
    PetscCallA(MatMult(gridmap_pio%map_frac_mat, gridmap_pio%src_vec, gridmap_pio%dst_vec, ierr))
    PetscCallA(MatMult(gridmap_pio%map_frac_mat, tmp_src_vec        , tmp_dst_vec        , ierr))

    ! fill the destination array
    PetscCallA(VecGetArray(gridmap_pio%dst_vec, dst_p, ierr))
    PetscCallA(VecGetArray(tmp_dst_vec        , tmp_p, ierr))
    do no = 1, gridmap_pio%dim_nb%nloc
       idx = gridmap_pio%dim_nb%begd + no - 1
       if (tmp_p(no) <= 0._r8) then
          dst_array(no) = nodata
       else
          dst_array(no) = dst_p(no)/tmp_p(no)
       end if
    end do
    PetscCallA(VecRestoreArray(tmp_dst_vec        , tmp_p, ierr))
    PetscCallA(VecRestoreArray(gridmap_pio%dst_vec, dst_p, ierr))

  end subroutine gridmap_areaave_srcmask_pio

  !------------------------------------------------------------------------------
  subroutine gridmap_dominant_value_pio(gridmap_pio, row_indices, src_array, minval_int, maxval_int, nodata, dst_array)

    implicit none

    type(gridmap_pio_type)           , intent(in)  :: gridmap_pio
    integer                , pointer , intent(in)  :: row_indices(:)
    integer                          , intent(in)  :: src_array(:)
    integer                          , intent(in)  :: minval_int
    integer                          , intent(in)  :: maxval_int
    integer                          , intent(in)  :: nodata               ! value to apply where there are no input data
    integer                          , intent(out) :: dst_array(:)

    integer                              :: ni, no, idx, ival
    PetscReal              , pointer     :: dst_p(:)
    real(r8)               , pointer     :: max_wts(:)
    real(r8)               , pointer     :: src_values(:)
    PetscErrorCode                       :: ierr
    character(*)           , parameter   :: subName = '(gridmap_domainant_value_pio) '

    ! check the map
    call gridmap_checkifset_pio(gridmap_pio, subname)

    ! allocate memory to save the maximum wt
    allocate(max_wts(gridmap_pio%dim_nb%nloc))
    allocate(src_values(gridmap_pio%dim_na%nloc))

    max_wts(:) = 0._r8
    dst_array(:) = nodata

    do ival = minval_int, maxval_int

       ! fill the value only for the cells that have values corresponding to 'ival'
       do ni = 1, gridmap_pio%dim_na%nloc
          if (src_array(ni) == ival) then
             src_values(ni) = 1._r8
          else
             src_values(ni) = 0._r8
          end if
       end do
       PetscCallA(VecSetValues(gridmap_pio%src_vec, gridmap_pio%dim_na%nloc, row_indices, src_values, INSERT_VALUES, ierr))
       PetscCallA(VecAssemblyBegin(gridmap_pio%src_vec, ierr))
       PetscCallA(VecAssemblyEnd(gridmap_pio%src_vec, ierr))

       ! do matrix-vec multiplication
       PetscCallA(MatMult(gridmap_pio%map_frac_mat, gridmap_pio%src_vec, gridmap_pio%dst_vec, ierr))

       ! fill the destination array if the weight corresponding to the 'ival' is max
       PetscCallA(VecGetArray(gridmap_pio%dst_vec, dst_p, ierr))
       do no = 1, gridmap_pio%dim_nb%nloc
          if (dst_p(no) > max_wts(no)) then
             dst_array(no) = ival
             max_wts(no)   = dst_p(no)
          end if
       end do
       PetscCallA(VecRestoreArray(gridmap_pio%dst_vec, dst_p, ierr))

    end do

    ! set value to nodata in destination array for which fraction is or less than zero
    do no = 1, gridmap_pio%dim_nb%nloc
       idx = gridmap_pio%dim_nb%begd + no - 1
       if (gridmap_pio%dst%frac(idx) <= 0._r8) then
          dst_array(no) = nodata
       end if
    end do

    ! free up memory
    deallocate(max_wts)
    deallocate(src_values)

  end subroutine gridmap_dominant_value_pio

  !------------------------------------------------------------------------------
  subroutine gridmap_dim_clean(dim)
    implicit none
    type(gridmap_dim_type), intent(inout) :: dim

    deallocate(dim%compdof)

  end subroutine gridmap_dim_clean

  !------------------------------------------------------------------------------
  subroutine gridmap_clean_pio(gridmap_pio)

    implicit none
    type(gridmap_pio_type), intent(inout) :: gridmap_pio

    character(len=*), parameter :: subName = "gridmap_pio_clean"
    PetscErrorCode              :: ierr

    if (gridmap_pio%set .eq. isSet) then

       call gridmap_dim_clean(gridmap_pio%dim_na)
       call gridmap_dim_clean(gridmap_pio%dim_nb)
       call gridmap_dim_clean(gridmap_pio%dim_ns)

       call gridmap_grid_clean(gridmap_pio%src)
       call gridmap_grid_clean(gridmap_pio%dst)

       deallocate(gridmap_pio%wovr)
       deallocate(gridmap_pio%row)
       deallocate(gridmap_pio%col)

       PetscCallA(VecDestroy(gridmap_pio%src_vec      , ierr))
       PetscCallA(VecDestroy(gridmap_pio%dst_vec      , ierr))
       PetscCallA(MatDestroy(gridmap_pio%map_mat      , ierr))
       PetscCallA(MatDestroy(gridmap_pio%map_frac_mat , ierr))
    else
       if (iam == 0) then
          write(6,*) subName//' Warning: claling '//trim(subName)//' on unallocated gridmap_pio'
       end if
    end if

    gridmap_pio%set = "NOT-set"

  end subroutine gridmap_clean_pio

  !------------------------------------------------------------------------------
  subroutine gridmap_checkifset_pio( gridmap_pio, subname )

    implicit none
    type(gridmap_pio_type), intent(in) :: gridmap_pio
    character(len=*)      , intent(in) :: subname

    if ( gridmap_pio%set .ne. IsSet )then
       write(6,*) SubName//' ERROR: gridmap_pio NOT set yet, run gridmap_mapread_pio first'
       call abort()
    end if

  end subroutine gridmap_checkifset_pio

  
end module mkgridmapPIOMod
