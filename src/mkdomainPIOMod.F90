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
     integer          :: ns_loc        ! local size of domain
     integer          :: begs, ends

     integer          :: begi, endi
     integer          :: begj, endj
     integer          :: ni_glb,nj_glb ! global size of domain for 2d domains only
     integer          :: ni_loc,nj_loc ! local size of domain for 2d domains only

     integer, pointer :: dim_glb(:)    ! global dimension

     real(r8)         :: edgen         ! lsmedge north
     real(r8)         :: edgee         ! lsmedge east
     real(r8)         :: edges         ! lsmedge south
     real(r8)         :: edgew         ! lsmedge west
     integer ,pointer :: mask1d(:)       ! land mask: 1 = land, 0 = ocean
     real(r8),pointer :: frac1d(:)       ! fractional land
     real(r8),pointer :: latc1d(:)       ! latitude of grid cell (deg)
     real(r8),pointer :: lonc1d(:)       ! longitude of grid cell (deg)
     real(r8),pointer :: lats1d(:)       ! grid cell latitude, S edge (deg)
     real(r8),pointer :: latn1d(:)       ! grid cell latitude, N edge (deg)
     real(r8),pointer :: lonw1d(:)       ! grid cell longitude, W edge (deg)
     real(r8),pointer :: lone1d(:)       ! grid cell longitude, E edge (deg)
     real(r8),pointer :: area1d(:)       ! grid cell area (km**2) (only used for output grid)

     integer ,pointer :: mask2d(:,:)       ! land mask: 1 = land, 0 = ocean
     real(r8),pointer :: frac2d(:,:)       ! fractional land
     real(r8),pointer :: latc2d(:,:)       ! latitude of grid cell (deg)
     real(r8),pointer :: lonc2d(:,:)       ! longitude of grid cell (deg)
     real(r8),pointer :: lats2d(:,:)       ! grid cell latitude, S edge (deg)
     real(r8),pointer :: latn2d(:,:)       ! grid cell latitude, N edge (deg)
     real(r8),pointer :: lonw2d(:,:)       ! grid cell longitude, W edge (deg)
     real(r8),pointer :: lone2d(:,:)       ! grid cell longitude, E edge (deg)
     real(r8),pointer :: area2d(:,:)       ! grid cell area (km**2) (only used for output grid)

     logical          :: is_2d         ! if this is a 2-d domain
     logical          :: fracset       ! if frac is set
     logical          :: maskset       ! if mask is set
     integer,pointer :: compdof(:)
  end type domain_pio_type

#include "petsc/finclude/petscsys.h"

  public domain_read_pio
  public domain_read_map_pio

  real(r8) :: flandmin = 0.001            !minimum land frac for land cell

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
  subroutine read_float_or_double(domain, pioIoSystem, ncid, varname, dataBuffer1d, dataBuffer2d)

    use pio, only : PIO_DOUBLE
    !
    implicit none

    type(domain_pio_type) , intent(inout)        :: domain
    type(file_desc_t)     , intent(in)           :: ncid
    type(iosystem_desc_t) , intent(in)           :: pioIoSystem
    character(len=*)      , intent(in)           :: varname
    real(r8)              , intent(out), pointer :: dataBuffer1d(:)
    real(r8)              , intent(out), pointer :: dataBuffer2d(:,:)

    real, pointer      :: dataBuffer1DFloat(:), dataBuffer2DFloat(:,:)
    type(var_desc_t)   :: varid
    type(io_desc_t)    :: iodescNCells
    integer            :: vartype
    integer            :: ier
    character(len= 32) :: subname = 'read_float_or_double'

    ! Inquiry the variable ID and variable type
    call check_ret(PIO_inq_varid(ncid, varname, varid), subname)
    call check_ret(PIO_inq_vartype(ncid, varid, vartype), subname)

    call PIO_initdecomp(pioIoSystem, vartype, domain%dim_glb, domain%compdof, iodescNCells)

    ! Read the variable
    if (domain%is_2d) then

       if (vartype == PIO_DOUBLE) then
          call PIO_read_darray(ncid, varid, iodescNCells, dataBuffer2d, ier)
       else
          allocate(dataBuffer2DFloat(domain%begi:domain%endi, domain%begj:domain%endj))
          call PIO_read_darray(ncid, varid, iodescNCells, dataBuffer2DFloat, ier)
          dataBuffer2D = dble(dataBuffer2DFloat)
          deallocate(dataBuffer2DFloat)
       end if
    else

       if (vartype == PIO_DOUBLE) then
          call PIO_read_darray(ncid, varid, iodescNCells, dataBuffer1d, ier)
          allocate(dataBuffer2DFloat(domain%begi:domain%endi, domain%begj:domain%endj))
       else
          allocate(dataBuffer1DFloat(domain%begs:domain%ends))
          call PIO_read_darray(ncid, varid, iodescNCells, dataBuffer1DFloat, ier)
          dataBuffer1D = dble(dataBuffer1DFloat)
          deallocate(dataBuffer1DFloat)

       end if
    end if

    ! Delete decomposition
    call PIO_freedecomp(pioIoSystem, iodescNCells)

  end subroutine read_float_or_double

  !----------------------------------------------------------------------------
  subroutine check_and_read_float_or_double(domain, pioIoSystem, ncid, varname, dataRead, dataBuffer1d, dataBuffer2d)

    implicit none

    type(domain_pio_type) , intent(inout)        :: domain
    type(file_desc_t)     , intent(in)           :: ncid
    type(iosystem_desc_t) , intent(in)           :: pioIoSystem
    character(len=*)      , intent(in)           :: varname
    logical               , intent(out)          :: dataRead
    real(r8)              , intent(out), pointer :: dataBuffer1d(:)
    real(r8)              , intent(out), pointer :: dataBuffer2d(:,:)

    integer                                      :: varid
    integer                                      :: ier

    dataRead = .false.

    ier = PIO_inq_varid(ncid, varname, varid)
    if (ier == PIO_NOERR) then
       dataRead = .true.
       call read_float_or_double(domain, pioIoSystem, ncid, varname, dataBuffer1d, dataBuffer2d)
    end if

  end subroutine check_and_read_float_or_double

  !----------------------------------------------------------------------------

  subroutine read_integer(domain, pioIoSystem, ncid, varname, dataBuffer1d, dataBuffer2d)

    implicit none

    type(domain_pio_type) , intent(inout)        :: domain
    type(file_desc_t)     , intent(in)           :: ncid
    type(iosystem_desc_t) , intent(in)           :: pioIoSystem
    character(len=*)      , intent(in)           :: varname
    integer               , intent(out), pointer :: dataBuffer1d(:)
    integer               , intent(out), pointer :: dataBuffer2d(:,:)

    type(var_desc_t)   :: varid
    type(io_desc_t)    :: iodescNCells
    integer            :: vartype
    integer            :: ier
    character(len= 32) :: subname = 'read_float_or_double'

    ! Inquiry the variable ID and variable type
    call check_ret(PIO_inq_varid(ncid, varname, varid), subname)

    ! Define a decomposition
    call PIO_initdecomp(pioIoSystem, PIO_Int, domain%dim_glb, domain%compdof, iodescNCells)

    ! Read the variable
    if (domain%is_2d) then
       call PIO_read_darray(ncid, varid, iodescNCells, dataBuffer2d, ier)
    else
       call PIO_read_darray(ncid, varid, iodescNCells, dataBuffer1d, ier)
    end if

    ! Delete decomposition
    call PIO_freedecomp(pioIoSystem, iodescNCells)

  end subroutine read_integer

  !----------------------------------------------------------------------------
  subroutine check_and_read_integer(domain, pioIoSystem, ncid, varname, dataRead, dataBuffer1d, dataBuffer2d)

    implicit none

    type(domain_pio_type) , intent(inout)        :: domain
    type(file_desc_t)     , intent(in)           :: ncid
    type(iosystem_desc_t) , intent(in)           :: pioIoSystem
    character(len=*)      , intent(in)           :: varname
    logical               , intent(out)          :: dataRead
    integer               , intent(out), pointer :: dataBuffer1d(:)
    integer               , intent(out), pointer :: dataBuffer2d(:,:)

    integer                                      :: varid
    integer                                      :: ier

    dataRead = .false.

    ier = PIO_inq_varid(ncid, varname, varid)
    if (ier == PIO_NOERR) then
       dataRead = .true.
       call read_integer(domain, pioIoSystem, ncid, varname, dataBuffer1d, dataBuffer2d)
    end if

  end subroutine check_and_read_integer

  !----------------------------------------------------------------------------
  subroutine domain_read_dims_pio(domain, pioIoSystem, ncid)

    implicit none

    type(domain_pio_type), intent(inout) :: domain
    type(iosystem_desc_t) :: pioIoSystem
    type(file_desc_t)     :: ncid

    logical :: dimset
    character(len= 32) :: subname = 'domain_read_dims_pio'

    dimset = .false.

    call domain_read_dims_2d_pio(domain, pioIoSystem, ncid, 'lsmlon', 'lsmlat', dimset)
    call domain_read_dims_2d_pio(domain, pioIoSystem, ncid, 'ni', 'nj', dimset)
    call domain_read_dims_2d_pio(domain, pioIoSystem, ncid, 'lon', 'lat', dimset)

    ! ----- If we haven't found any info, abort -----

    if (.not. dimset) then
       write(6,*) trim(subname),' ERROR: dims not set'
       call abort()
    endif

  end subroutine domain_read_dims_pio

  !----------------------------------------------------------------------------
  subroutine domain_read_dims_2d_pio(domain, pioIoSystem, ncid, lon_name, lat_name, dimset)

    implicit none

    type(domain_pio_type) , intent(inout) :: domain
    type(iosystem_desc_t) , intent(in)    :: pioIoSystem
    type(file_desc_t)     , intent(in)    :: ncid
    character(len=*)      , intent(in)    :: lon_name
    character(len=*)      , intent(in)    :: lat_name
    logical               , intent(inout) :: dimset

    integer               :: dimid
    integer               :: nlon,nlat
    integer               :: ier
    type(var_desc_t)      :: varid
    type(io_desc_t)       :: iodescNCells

    character(len= 32) :: subname = 'domain_read_dims_2d_pio'

    ier = PIO_inq_dimid(ncid, lon_name, dimid)

    if (ier == PIO_NOERR) then
       if (dimset) then
          if (masterproc) write(6,*) trim(subname),' WARNING: dimension sizes already set; skipping ', &
               lon_name, '/', lat_name
       else
          if (masterproc) write(6,*) trim(subname),' read lon and lat dims from ', lon_name, '/', lat_name

          call check_ret(PIO_inq_dimid  (ncid, lon_name, dimid), subname)
          call check_ret(PIO_inq_dimlen (ncid, dimid, nlon), subname)

          call check_ret(PIO_inq_dimid  (ncid, lat_name, dimid), subname)
          call check_ret(PIO_inq_dimlen (ncid, dimid, nlat), subname)

          domain%ni_glb = nlon
          domain%nj_glb = nlat
          domain%ns_glb = nlon * nlat
          domain%is_2d = .true.
          dimset       = .true.

          allocate(domain%dim_glb(2))
          domain%dim_glb(1) = nlon
          domain%dim_glb(2) = nlat
          if (masterproc) write(6,*) trim(subname),' nlon ',nlon,'nlat ',nlat

       end if
    end if


  end subroutine domain_read_dims_2d_pio

  !----------------------------------------------------------------------------
  subroutine domain_read_pio(domain, fname)

    implicit none

    type(domain_pio_type), intent(inout) :: domain
    character(len=*)     , intent(in)    :: fname

    integer               :: niotasks
    integer               :: numAggregator
    integer               :: stride
    integer               :: optBase
    integer               :: iotype
    integer               :: ier

    integer               :: i,j,n
    type(file_desc_t)     :: ncid
    type(var_desc_t)      :: varid
    type(iosystem_desc_t) :: pioIoSystem
    type(io_desc_t)       :: iodescNCells
    logical :: edgeNESWset                     ! local EDGE[NESW]
    logical :: lonlatset                       ! local lon(:,:), lat(:,:)
    logical :: llneswset                       ! local lat[ns],lon[we]
    logical :: landfracset                     ! local landfrac
    logical :: maskset                         ! local mask
    logical :: lreadmask                       ! local readmask
    character(len= 32) :: lonvar               ! name of 2-d longitude variable
    character(len= 32) :: latvar               ! name of 2-d latitude variable

    character(len= 32)    :: subname = 'domain_read_pio'

    edgeNESWset = .false. 
    llneswset   = .false. 
    lreadmask   = .false.

    !if (present(readmask)) then
    !   lreadmask = readmask
    !end if

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
    call check_ret(PIO_openfile(pioIoSystem, ncid, iotype, trim(fname), PIO_NOWRITE), subname)

    !
    call domain_read_dims_pio(domain, pioIoSystem, ncid)
    call domain_init_pio(domain)
    if (masterproc) write(6,*) trim(subname),' initialized domain'


    ! ----- Set lat/lon 

    lonlatset   = .false. 
    call check_and_read_float_or_double(domain, pioIoSystem, ncid, 'LONGXY', lonlatset, domain%lonc1d, domain%lonc2d)
    if (lonlatset) then
       ! 'LONGXY' was present, so let's read 'LATIXY'
       call check_and_read_float_or_double(domain, pioIoSystem, ncid, 'LATIXY', lonlatset, domain%latc1d, domain%latc2d)
    else
       ! 'LONGXY' was not present, so check if 'lon' is present in the file
       call check_and_read_float_or_double(domain, pioIoSystem, ncid, 'lon', lonlatset, domain%lonc1d, domain%lonc2d)
       if (lonlatset) then
          ! 'lon' was present, so let's read 'lat'
          call check_and_read_float_or_double(domain, pioIoSystem, ncid, 'lat', lonlatset, domain%latc1d, domain%latc2d)
       else
          ! 'lon' was also not present, so check if 'LATITUDE' is present in the file
          call check_and_read_float_or_double(domain, pioIoSystem, ncid, 'LATITUDE', lonlatset, domain%lonc1d, domain%lonc2d)
          if (lonlatset) then
             ! 'LATITUDE' was present, so let's read 'LONGITUDE'
             call check_and_read_float_or_double(domain, pioIoSystem, ncid, 'LONGITUDE', lonlatset, domain%latc1d, domain%lonc2d)
          else
             write(6,*)'lon/lat values not set' 
             write(6,*)'currently assume either that lon/lat or LONGXY/LATIXY', &
                  ' or LONGITUDE/LATITUDE variables are on input dataset'
             call abort()
          end if
       end if
    end if

    ! Read 'frac' or 'LANDFRAC'
    landfracset = .false. 
    call check_and_read_float_or_double(domain, pioIoSystem, ncid, 'frac', landfracset, domain%frac1d, domain%frac2d)
    if (.not. landfracset) then
       call check_and_read_float_or_double(domain, pioIoSystem, ncid, 'LANDFRAC', landfracset, domain%frac1d, domain%frac2d)
    end if

    ! Read the 'mask' or 'LANDMASK'
    maskset = .false.
    if (lreadmask) then
       call check_and_read_integer(domain, pioIoSystem, ncid, 'mask', maskset, domain%mask1d, domain%mask2d)
    else
       call check_and_read_integer(domain, pioIoSystem, ncid, 'mask', maskset, domain%mask1d, domain%mask2d)
       if (.not. maskset) then
          call check_and_read_integer(domain, pioIoSystem, ncid, 'LANDMASK', maskset, domain%mask1d, domain%mask2d)
       end if
    end if
    
    ! clean up
    call PIO_closefile(ncid)
    call PIO_finalize(pioIoSystem, ier)

    ! ----- set derived variables ----
    
    if (.not.maskset.and.landfracset) then
       maskset = .true.
       if (domain%is_2d) then
          do i = domain%begi, domain%endi
             do j = domain%begj, domain%endj
                if (domain%frac2d(i,j) < flandmin) then
                   domain%mask2d(i,j) = 0 ! ocean
                else
                   domain%mask2d(i,j) = 1 ! land
                end if
             end do
          end do
       else
          do n = domain%begs,domain%ends
             if (domain%frac1d(n) < flandmin) then
                domain%mask1d(n) = 0 ! ocean
             else
                domain%mask1d(n) = 1 ! land
             end if
          end do
       end if
    endif

    if (.not.landfracset.and.maskset) then
       landfracset = .true.

       if (domain%is_2d) then
          do i = domain%begi, domain%endi
             do j = domain%begj, domain%endj
                if (domain%mask2d(i,j) == 0) then
                   domain%frac2d(i,j) = 0._r8 ! ocean
                else
                   domain%frac2d(i,j) = 1._r8 ! land
                end if
             end do
          end do
       else

          do n = domain%begs,domain%ends
             if ( domain%mask1d(n) == 0 )then
                domain%frac1d(n) = 0._r8     !ocean
             else
                domain%frac1d(n) = 1._r8     !land
             end if
          end do
       end if
    endif
    domain%maskset = maskset
    domain%fracset = landfracset

  end subroutine domain_read_pio

!----------------------------------------------------------------------------
  subroutine domain_init_pio (domain)

    implicit none

    type(domain_pio_type),intent(inout) :: domain

    if (.not. domain%is_2d) then
       call domain_init_pio_1d (domain)
    else
       call domain_init_pio_2d (domain)
    end if
    
  end subroutine domain_init_pio

  !----------------------------------------------------------------------------
  subroutine find_start_and_end_indices(len, myrank, ntasks, ista, isto )
    implicit none
    integer :: len, myrank, ntasks
    integer :: ista, isto

    integer :: arrIdxPerPe, local_arrIdxPerPe, ierr

    arrIdxPerPe = len/ ntasks
    local_arrIdxPerPe = arrIdxPerPe
    if (myRank < LEN - arrIdxPerPe * ntasks) then
       local_arrIdxPerPe = local_arrIdxPerPe + 1
    end if
    call MPI_Scan(local_arrIdxPerPe, isto, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)

    if (arrIdxPerPe < 1) then
       write(*,*) "Not enough work to distribute among pes"
       call exit(0)
    endif

    ista = isto - local_arrIdxPerPe + 1
    isto = isto 

  end subroutine find_start_and_end_indices
  
  !----------------------------------------------------------------------------
  subroutine domain_init_pio_1d (domain)

    implicit none

    type(domain_pio_type),intent(inout) :: domain

    integer :: i, ier

    domain%ns_glb = domain%dim_glb(1)

    domain%ns_loc = domain%ns_glb/npes

    if (iam < domain%ns_glb - domain%ns_loc * npes) then
       domain%ns_loc = domain%ns_loc + 1
    end if

    call MPI_Scan(domain%ns_loc, domain%ends, 1, MPI_INTEGER, MPI_SUM, mpicom, ier)
    domain%begs = domain%ends - domain%ns_loc + 1

    allocate(domain%mask1d(domain%begs:domain%ends), &
         domain%frac1d(domain%begs:domain%ends), &
         domain%latc1d(domain%begs:domain%ends), &
         domain%lonc1d(domain%begs:domain%ends), &
         domain%lats1d(domain%begs:domain%ends), &
         domain%latn1d(domain%begs:domain%ends), &
         domain%lonw1d(domain%begs:domain%ends), &
         domain%lone1d(domain%begs:domain%ends), &
         domain%area1d(domain%begs:domain%ends), &
         domain%compdof(domain%begs:domain%ends), &
         stat=ier)

    domain%compdof(domain%begs:domain%ends) = (/(i, i=domain%begs,domain%ends, 1)/)

  end subroutine domain_init_pio_1d
  
  !----------------------------------------------------------------------------
  subroutine domain_init_pio_2d (domain)

    implicit none

    type(domain_pio_type),intent(inout) :: domain

    integer :: i,j,count,ier

    domain%ni_glb = domain%dim_glb(1) ! nlon
    domain%nj_glb = domain%dim_glb(2) ! nlat

    call find_start_and_end_indices(domain%ni_glb, iam, npes, domain%begi, domain%endi)
    domain%begj = 1
    domain%endj = domain%nj_glb

    domain%ns_loc = (domain%endi - domain%begi + 1)*(domain%endj - domain%begj + 1)

    allocate(domain%mask2d(domain%begi:domain%endi, domain%begj:domain%endj), &
         domain%frac2d(domain%begi:domain%endi, domain%begj:domain%endj), &
         domain%latc2d(domain%begi:domain%endi, domain%begj:domain%endj), &
         domain%lonc2d(domain%begi:domain%endi, domain%begj:domain%endj), &
         domain%lats2d(domain%begi:domain%endi, domain%begj:domain%endj), &
         domain%latn2d(domain%begi:domain%endi, domain%begj:domain%endj), &
         domain%lonw2d(domain%begi:domain%endi, domain%begj:domain%endj), &
         domain%lone2d(domain%begi:domain%endi, domain%begj:domain%endj), &
         domain%area2d(domain%begi:domain%endi, domain%begj:domain%endj), &
         domain%compdof(1:domain%ns_loc), &
         stat=ier)

    count = 0
    do j = 1,(domain%endj - domain%begj + 1)
       do i = 1, (domain%endi - domain%begi + 1)
          count = count + 1
          domain%compdof(count) = i + (j-1)*domain%ni_glb + (domain%begi - 1)
       end do
    end do

  end subroutine domain_init_pio_2d

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

    domain%ni_loc = -9999
    domain%nj_loc = -9999
    domain%ns_loc = -9999
    domain%ni_glb = -9999
    domain%nj_glb = -9999
    domain%ns_glb = -9999
    domain%begi   = 0
    domain%endi   = 0
    domain%begj   = 0
    domain%endj   = 0
    domain%begs   = 0
    domain%ends   = 0
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
          write(6,*)'add code to support a map for a 2D domain grid'
          call abort()
       else
          allocate(domain%dim_glb(1))
          domain%dim_glb(1) = domain%ns_glb
       end if

       call domain_init_pio(domain)

       call read_float_or_double(domain, pioIoSystem, ncid, 'xc_b', domain%lonc1d, domain%lonc2d)
       call read_float_or_double(domain, pioIoSystem, ncid, 'yc_b', domain%latc1d, domain%latc2d)
       call read_float_or_double(domain, pioIoSystem, ncid, 'frac_b', domain%frac1d, domain%frac2d)
       call read_float_or_double(domain, pioIoSystem, ncid, 'area_b', domain%area1d, domain%area2d)
       if (domain%is_2d) then
          domain%area2d = domain%area2d * re**2._r8
       else
          domain%area1d = domain%area1d * re**2._r8
       end if

       call read_integer(domain, pioIoSystem, ncid, 'mask_b', domain%mask1d, domain%mask2d)

    end if

    call PIO_closefile(ncid)
    call PIO_finalize(pioIoSystem, ier)

  end function domain_read_map_pio

  !----------------------------------------------------------------------------

end module mkdomainPIOMod
