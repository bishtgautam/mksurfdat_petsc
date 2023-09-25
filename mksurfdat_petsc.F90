program mksurfdat_petsc

  use shr_kind_mod, only : r8 => shr_kind_r8, r4 => shr_kind_r4
  use mkdomainMod , only : domain_type, domain_read_map, domain_read, &
                           domain_write
  use mkvarctl
  use petsc

  implicit none

#include "petsc/finclude/petscsys.h"

  integer  :: ier                         ! error status
  character(len=256) :: fsurlog           ! output surface log file name
  PetscErrorCode :: ierr                  ! PETSc error status

  type(domain_type) :: ldomain

  namelist /elmexp/              &
       mksrf_fgrid,              &	
       mksrf_gridtype,           &
       outnc_large_files,        &
       outnc_double,             &
       outnc_dims,               &
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

  ! Finalize PETSc
  PetscCallA(PetscFinalize(ierr))

contains

  !-----------------------------------------------------------------------
  subroutine setup_namelist()

    implicit none

    ! Default settings
    mksrf_gridtype = 'global'
    mksrf_fgrid    = ' '
    fsurlog        = ''

    read(5, elmexp, iostat=ier)
    if (ier /= 0) then
       write(6,*)'error: namelist input resulted in error code ',ier
       call abort()
    endif

    call check_namelist_variable(mksrf_fgrid, 'mksrf_fgrid')
    call check_namelist_variable(fsurlog, 'fsurlog')

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

end program mksurfdat_petsc


