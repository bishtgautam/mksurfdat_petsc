module mkdomainPIOMod

  use petsc
  use shr_kind_mod , only : r8 => shr_kind_r8
  use mkvarpar     , only : re
  use nanMod       , only : nan, bigint
  use pio          , only : PIO_init, PIO_rearr_subset, iosystem_desc_t, file_desc_t
  use pio          , only : PIO_finalize, PIO_noerr, PIO_iotype_netcdf, PIO_createfile
  use pio          , only : PIO_iotype_pnetcdf
  use pio          , only : PIO_int,var_desc_t, PIO_redef, PIO_def_dim, PIO_def_var, PIO_enddef
  use pio          , only : PIO_closefile, io_desc_t, PIO_initdecomp, PIO_write_darray
  use pio          , only : PIO_freedecomp, PIO_clobber, PIO_read_darray, PIO_syncfile, PIO_OFFSET_KIND
  use pio          , only : PIO_nowrite, PIO_openfile
  use pio
  use spmdMod      , only : iam, npes, masterproc, mpicom

  implicit none
  private

  public :: domain_pio_type

  type domain_pio_type
     character*16     :: set           ! flag to check if domain is set
     integer          :: ns_glb        ! local size of domain
     integer          :: ns            ! local size of domain
     integer          :: begs, ends
     integer          :: ni_glb,nj_glb ! global size of domain for 2d domains only
     integer          :: ni,nj         ! local size of domain for 2d domains only
     integer, pointer :: dim_glb(:)    ! global dimension
     real(r8)         :: edgen         ! lsmedge north
     real(r8)         :: edgee         ! lsmedge east
     real(r8)         :: edges         ! lsmedge south
     real(r8)         :: edgew         ! lsmedge west
     integer ,pointer :: mask(:)       ! land mask: 1 = land, 0 = ocean
     real(r8),pointer :: frac(:)       ! fractional land
     real(r8),pointer :: latc(:)       ! latitude of grid cell (deg)
     real(r8),pointer :: lonc(:)       ! longitude of grid cell (deg)
     real(r8),pointer :: lats(:)       ! grid cell latitude, S edge (deg)
     real(r8),pointer :: latn(:)       ! grid cell latitude, N edge (deg)
     real(r8),pointer :: lonw(:)       ! grid cell longitude, W edge (deg)
     real(r8),pointer :: lone(:)       ! grid cell longitude, E edge (deg)
     real(r8),pointer :: area(:)       ! grid cell area (km**2) (only used for output grid)
     logical          :: is_2d         ! if this is a 2-d domain
     logical          :: fracset       ! if frac is set
     logical          :: maskset       ! if mask is set
     integer,pointer :: compdof(:)
  end type domain_pio_type

#include "petsc/finclude/petscsys.h"

  public domain_read_map_pio

contains

  !-----------------------------------------------------------------------
  !BOP
  !
  ! !IROUTINE: check_ret
  !
  ! !INTERFACE:
  subroutine check_ret(ret, calling)
    !
    ! !DESCRIPTION:
    ! Check return status from netcdf call
    !
    ! !ARGUMENTS:
    implicit none
    integer, intent(in) :: ret
    character(len=*)    :: calling
    !
    integer             :: ier
    character(len=256)  :: errmsg
    !
    ! !REVISION HISTORY:
    !
    !EOP
    !-----------------------------------------------------------------------

    if (ret /= PIO_NOERR) then
       ier = PIO_STRERROR(ret, errmsg)
       write(6,*)'netcdf error from ',trim(calling), ' rcode = ', ret, &
            ' error = ', trim(errmsg)
       call abort()
    end if

  end subroutine check_ret

  !----------------------------------------------------------------------------
  subroutine read_float_or_double(domain, pioIoSystem, ncid, varname, dataBuffer)

    implicit none

    type(domain_pio_type) , intent(inout)        :: domain
    type(file_desc_t)     , intent(in)           :: ncid
    type(iosystem_desc_t) , intent(in)           :: pioIoSystem
    character(len=*)      , intent(in)           :: varname
    real(r8)              , intent(out), pointer :: dataBuffer(:)

    type(var_desc_t)   :: varid
    type(io_desc_t)    :: iodescNCells
    integer            :: vartype
    integer            :: ier
    character(len= 32) :: subname = 'read_float_or_double'

    ! Inquiry the variable ID and variable type
    call check_ret(PIO_inq_varid(ncid, varname, varid), subname)
    call check_ret(PIO_inq_vartype(ncid, varid, vartype), subname)

    ! Define a decomposition
    call PIO_initdecomp(pioIoSystem, vartype, domain%dim_glb, domain%compdof(domain%begs:domain%ends), iodescNCells)

    ! Read the variable
    call PIO_read_darray(ncid, varid, iodescNCells, dataBuffer, ier)

    ! Delete decomposition
    call PIO_freedecomp(pioIoSystem, iodescNCells)

  end subroutine read_float_or_double

  !----------------------------------------------------------------------------

  subroutine read_integer(domain, pioIoSystem, ncid, varname, dataBuffer)

    implicit none

    type(domain_pio_type) , intent(inout)        :: domain
    type(file_desc_t)     , intent(in)           :: ncid
    type(iosystem_desc_t) , intent(in)           :: pioIoSystem
    character(len=*)      , intent(in)           :: varname
    integer               , intent(out), pointer :: dataBuffer(:)

    type(var_desc_t)   :: varid
    type(io_desc_t)    :: iodescNCells
    integer            :: vartype
    integer            :: ier
    character(len= 32) :: subname = 'read_float_or_double'

    ! Inquiry the variable ID and variable type
    call check_ret(PIO_inq_varid(ncid, varname, varid), subname)

    ! Define a decomposition
    call PIO_initdecomp(pioIoSystem, PIO_Int, domain%dim_glb, domain%compdof(domain%begs:domain%ends), iodescNCells)

    ! Read the variable
    call PIO_read_darray(ncid, varid, iodescNCells, dataBuffer, ier)

    ! Delete decomposition
    call PIO_freedecomp(pioIoSystem, iodescNCells)

  end subroutine read_integer

  !----------------------------------------------------------------------------
  logical function domain_read_map_pio(domain, fname)

    use mkutilsMod, only : convert_latlon

    implicit none

    type(domain_pio_type),intent(inout) :: domain
    character(len=*)      ,intent(in)   :: fname

    integer               :: niotasks
    integer               :: numAggregator
    integer               :: stride
    integer               :: optBase
    integer               :: iotype
    integer               :: retVal, ier
    integer               :: dimid
    integer               :: grid_rank                       ! rank of domain grid 
    integer               :: ns                              ! size of domain grid
    integer               :: vartype
    integer               :: i,j

    type(file_desc_t)     :: ncid
    type(var_desc_t)      :: varid
    type(iosystem_desc_t) :: pioIoSystem
    type(io_desc_t)       :: iodescNCells

    integer :: grid_dims(2)
    character(len= 32)    :: subname = 'domain_read_map_pio'

    domain_read_map_pio = .true.

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

    ! Assume unstructured grid

    domain%ni = -9999
    domain%nj = -9999
    domain%ni_glb = -9999
    domain%nj_glb = -9999
    domain%is_2d = .false.

    call check_ret(PIO_openfile(pioIoSystem, ncid, iotype, trim(fname), PIO_NOWRITE), subname)

    ier = PIO_inq_dimid(ncid, 'n_b', dimid)

    if ( ier /= PIO_NOERR) then
       domain_read_map_pio = .false.
    else
       call check_ret(PIO_inq_dimlen(ncid, dimid, domain%ns_glb), subname)

       call check_ret(PIO_inq_dimid  (ncid, 'dst_grid_rank', dimid), subname)
       call check_ret(PIO_inq_dimlen (ncid, dimid, grid_rank), subname)

       if (grid_rank == 2) then
          allocate(domain%dim_glb(2))
          write(6,*)'add code to support 2D domain grid'
          call abort()
       else
          allocate(domain%dim_glb(1))
          domain%dim_glb(1) = domain%ns_glb
       end if

       domain%ns = domain%ns_glb/npes
       if (iam < domain%ns_glb - domain%ns * npes) then
          domain%ns = domain%ns + 1
       end if

       call MPI_Scan(domain%ns, domain%ends, 1, MPI_INTEGER, MPI_SUM, mpicom, ier)
       domain%begs = domain%ends - domain%ns + 1

       allocate(domain%mask(domain%begs:domain%ends), &
            domain%frac(domain%begs:domain%ends), &
            domain%latc(domain%begs:domain%ends), &
            domain%lonc(domain%begs:domain%ends), &
            domain%lats(domain%begs:domain%ends), &
            domain%latn(domain%begs:domain%ends), &
            domain%lonw(domain%begs:domain%ends), &
            domain%lone(domain%begs:domain%ends), &
            domain%area(domain%begs:domain%ends), &
            domain%compdof(domain%begs:domain%ends), &
            stat=ier)

       domain%compdof(domain%begs:domain%ends) = (/(i, i=domain%begs,domain%ends, 1)/)

       call read_float_or_double(domain, pioIoSystem, ncid, 'xc_b', domain%lonc)
       call read_float_or_double(domain, pioIoSystem, ncid, 'yc_b', domain%latc)
       call read_float_or_double(domain, pioIoSystem, ncid, 'frac_b', domain%frac)
       call read_float_or_double(domain, pioIoSystem, ncid, 'area_b', domain%area)
       domain%area = domain%area * re**2._r8

       call read_integer(domain, pioIoSystem, ncid, 'mask_b', domain%mask)

       call PIO_closefile(ncid)
       call PIO_finalize(pioIoSystem, ier)
    end if

  end function domain_read_map_pio
  !----------------------------------------------------------------------------

end module mkdomainPIOMod
