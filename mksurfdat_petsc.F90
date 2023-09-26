program mksurfdat_petsc

  use shr_kind_mod , only : r8 => shr_kind_r8, r4 => shr_kind_r4
  use mkdomainMod  , only : domain_type, domain_read_map, domain_read, &
       domain_write
  use mkpftMod     , only : mkpft
  use mkglcmecMod
  use mksoilMod
  use mkpftMod
  use mkvarctl
  use petsc
  use fileutils

  implicit none

#include "petsc/finclude/petscsys.h"

  integer               :: ier              ! error status
  character(len=256)    :: fsurlog          ! output surface log file name
  character(len=256)    :: fsurdat          ! output surface data file name
  character(len=256)    :: fdyndat          ! dynamic landuse data file name
  PetscErrorCode        :: ierr             ! PETSc error status
  integer               :: ndiag
  real(r8), allocatable :: pctlnd_pft(:)    ! PFT data: % of gridcell for PFTs
  real(r8), pointer     :: pctpft_full(:,:) ! PFT data: % cover of each pft and cft on the vegetated landunits
                                            ! ('full' denotes inclusion of CFTs as well as natural PFTs in this array)

  type(domain_type) :: ldomain

  namelist /elmexp/              &
       mksrf_fgrid,              &	
       mksrf_gridtype,           &
       mksrf_fvegtyp,            &
       mksrf_fsoitex,            &
       mksrf_forganic,           &
       mksrf_fsoicol,            &
       mksrf_fsoiord,            &
       mksrf_fvocef,             &
       mksrf_flakwat,            &
       mksrf_fwetlnd,            &
       mksrf_fglacier,           &
       mksrf_furbtopo,           &
       mksrf_flndtopo,           &
       mksrf_fmax,               &
       mksrf_furban,             &
       mksrf_flai,               &
       mksrf_fdynuse,            &
       mksrf_fgdp,               &
       mksrf_fpeat,              &
       mksrf_fabm,               &
       mksrf_ftopostats,         &
       mksrf_fvic,               &
       mksrf_fch4,               &
       mksrf_fphosphorus,        &
       mksrf_fgrvl,              &
       mksrf_fslp10,             &
       mksrf_fero,               &
       nglcec,                   &
       numpft,                   &
       soil_color,               &
       soil_order,               &
       soil_sand,                &
       soil_fmax,                &
       soil_clay,                &
       pft_idx,                  &
       pft_frc,                  &
       all_urban,                &
       no_inlandwet,             &
       map_fpft,                 &
       map_flakwat,              &
       map_fwetlnd,              &
       map_fglacier,             &
       map_fsoitex,              &
       map_fsoicol,              &
       map_fsoiord,              &
       map_furban,               &
       map_furbtopo,             &
       map_flndtopo,             &
       map_fmax,                 &
       map_forganic,             &
       map_fvocef,               &
       map_flai,                 &
       map_fharvest,             &
       map_fgdp,                 &
       map_fpeat,                &
       map_fabm,                 &
       map_ftopostats,           &
       map_fvic,                 &
       map_fch4,                 &
       map_fphosphorus,          &
       map_fgrvl,                &
       map_fslp10,               &
       map_fero,                 &
       outnc_large_files,        &
       outnc_double,             &
       outnc_dims,               &
       fsurdat,                  &
       fdyndat,                  &   
       fsurlog

  ! Initialize PETSc
  PetscCallA(PetscInitialize(ierr))

  write(6,*) 'Attempting to initialize control settings .....'

  call setup_namelist()

  write(6,*)'calling domain_read'
  if ( .not. domain_read_map(ldomain, mksrf_fgrid) )then
     call domain_read(ldomain, mksrf_fgrid)
  end if
  write(6,*)'finished domain_read'  

  ! Determine if will have 1d output

  if (ldomain%ni /= -9999 .and. ldomain%nj /= -9999) then
     write(6,*)'fsurdat is 2d lat/lon grid'
     write(6,*)'nlon= ',ldomain%ni,' nlat= ',ldomain%nj
     if (outnc_dims == 1) then
        write(6,*)' writing output file in 1d gridcell format'
     end if
  else
     write(6,*)'fsurdat is 1d gridcell grid'
     outnc_dims = 1
  end if

  outnc_1d = .false.
  if ((ldomain%ni == -9999 .and. ldomain%nj == -9999) .or. outnc_dims==1) then
     outnc_1d = .true.
     write(6,*)'output file will be 1d'
  end if

  ! allocate memory for all variables
  call allocate_memory(ldomain)

  ! ----------------------------------------------------------------------
  ! Make surface dataset fields
  ! ----------------------------------------------------------------------

  ! Make PFTs [pctpft_full] from dataset [fvegtyp]

  call mkpft(ldomain, mapfname=map_fpft, fpft=mksrf_fvegtyp, &
       ndiag=ndiag, pctlnd_o=pctlnd_pft, pctpft_o=pctpft_full )


  ! deallocate memory for all variables
  call deallocate_memory()

  ! Finalize PETSc
  PetscCallA(PetscFinalize(ierr))

contains

  !-----------------------------------------------------------------------
  subroutine setup_namelist()

    implicit none

    ! Default settings
    mksrf_gridtype = 'global'
    mksrf_fgrid    = ' '
    fsurlog        = ' '
    map_fpft       = ' '
    mksrf_fvegtyp  = ' '

    read(5, elmexp, iostat=ier)
    if (ier /= 0) then
       write(6,*)'error: namelist input resulted in error code ',ier
       call abort()
    endif

    call check_namelist_variable(mksrf_fgrid   , 'mksrf_fgrid')
    call check_namelist_variable(map_fpft      , 'map_fpft')
    call check_namelist_variable(mksrf_fvegtyp , 'mksrf_fvegtyp')

    call check_namelist_variable(fsurlog, 'fsurlog')
    ndiag = getavu()
    call opnfil(fsurlog, ndiag, 'f')
    
  end subroutine setup_namelist

  !-----------------------------------------------------------------------
  subroutine check_namelist_variable(var_value, var_name)

    implicit none
    character (len=*), intent(in) :: var_value
    character (len=*), intent(in) :: var_name

    if (var_value == ' ') then
       write(6,*)'must specify a value for ' // trim(var_name)
       call abort()
    end if

  end subroutine check_namelist_variable

  !-----------------------------------------------------------------------
  subroutine allocate_memory(ldomain)

    implicit none

    type(domain_type) :: ldomain
    
    integer :: ns_o
    integer :: ier

    ns_o = ldomain%ns

    allocate(pctlnd_pft(ns_o))
    allocate(pctpft_full(ns_o,0:numpft))
    
  end subroutine allocate_memory

  !-----------------------------------------------------------------------
  subroutine deallocate_memory()

    implicit none

    deallocate(pctlnd_pft)
    deallocate(pctpft_full)
    
  end subroutine deallocate_memory

end program mksurfdat_petsc


