program mksurfdat_petsc

  use shr_kind_mod , only : r8 => shr_kind_r8, r4 => shr_kind_r4
  use mkdomainMod  , only : domain_type, domain_read_map, domain_read, &
       domain_write
  use mkpftMod     , only : mkpft
  use mkglcmecMod
  use mksoilMod
  use mkpftMod
  use mkvarctl
  use mkurbanparMod
  use mkvarpar
  use mklanwatMod
  use petsc
  use fileutils

  implicit none

#include "petsc/finclude/petscsys.h"

  integer  :: nsoicol                     ! number of model color classes
  integer  :: nsoiord                     ! number of model order classes

  integer               :: ier              ! error status
  character(len=256)    :: fsurlog          ! output surface log file name
  character(len=256)    :: fsurdat          ! output surface data file name
  character(len=256)    :: fdyndat          ! dynamic landuse data file name

  logical :: all_veg ! if gridcell will be 100% vegetation land-cover

  PetscErrorCode        :: ierr             ! PETSc error status
  integer               :: ndiag

  real(r8), allocatable  :: landfrac_pft(:)    ! PFT data: % land per gridcell
  real(r8), allocatable  :: pctlnd_pft(:)      ! PFT data: % of gridcell for PFTs
  real(r8), allocatable  :: pctlnd_pft_dyn(:)  ! PFT data: % of gridcell for dyn landuse PFTs
  integer , allocatable  :: pftdata_mask(:)    ! mask indicating real or fake land type
  real(r8), pointer      :: pctpft_full(:,:)   ! PFT data: % cover of each pft and cft on the vegetated landunits
                                               ! ('full' denotes inclusion of CFTs as well as natural PFTs in this array)
  real(r8), allocatable  :: pctnatveg(:)       ! percent of grid cell that is natural veg landunit
  real(r8), allocatable  :: pctcrop(:)         ! percent of grid cell that is crop landunit
  real(r8), allocatable  :: pctnatpft(:,:)     ! % of each pft on the natural veg landunit (adds to 100%)
  real(r8), allocatable  :: pctcft(:,:)        ! % of each cft on the crop landunit (adds to 100%)
  real(r8), pointer      :: harvest(:,:)       ! harvest data: normalized harvesting
  real(r8), allocatable  :: pctgla(:)          ! percent of grid cell that is glacier
  real(r8), allocatable  :: pctglc_gic(:)      ! percent of grid cell that is gic (% of glc landunit)
  real(r8), allocatable  :: pctglc_icesheet(:) ! percent of grid cell that is ice sheet (% of glc landunit)
  real(r8), allocatable  :: pctglcmec(:,:)     ! glacier_mec pct coverage in each class (% of landunit)
  real(r8), allocatable  :: topoglcmec(:,:)    ! glacier_mec sfc elevation in each gridcell and class
  real(r8), allocatable  :: pctglcmec_gic(:,:) ! GIC pct coverage in each class (% of landunit)
  real(r8), allocatable  :: pctglcmec_icesheet(:,:) ! icesheet pct coverage in each class (% of landunit)
  real(r8), allocatable  :: elevclass(:)       ! glacier_mec elevation classes
  real(r8), allocatable  :: pctlak(:)          ! percent of grid cell that is lake
  real(r8), allocatable  :: pctwet(:)          ! percent of grid cell that is wetland
  real(r8), allocatable  :: pcturb(:)          ! percent of grid cell that is urbanized (total across all urban classes)
  real(r8), allocatable  :: urbn_classes(:,:)  ! percent cover of each urban class, as % of total urban area
  real(r8), allocatable  :: urbn_classes_g(:,:)! percent cover of each urban class, as % of grid cell
  real(r8), allocatable  :: elev(:)            ! glc elevation (m)
  real(r8), allocatable  :: topo(:)            ! land elevation (m)
  real(r8), allocatable  :: fmax(:)            ! fractional saturated area
  integer , allocatable  :: soicol(:)          ! soil color
  integer , allocatable  :: soiord(:)          ! soil order
  real(r8), allocatable  :: pctsand(:,:)       ! soil texture: percent sand
  real(r8), allocatable  :: pctclay(:,:)       ! soil texture: percent clay
  real(r8), allocatable  :: ef1_btr(:)         ! Isoprene emission factor for broadleaf
  real(r8), allocatable  :: ef1_fet(:)         ! Isoprene emission factor for fine/everg
  real(r8), allocatable  :: ef1_fdt(:)         ! Isoprene emission factor for fine/dec
  real(r8), allocatable  :: ef1_shr(:)         ! Isoprene emission factor for shrubs
  real(r8), allocatable  :: ef1_grs(:)         ! Isoprene emission factor for grasses
  real(r8), allocatable  :: ef1_crp(:)         ! Isoprene emission factor for crops
  real(r8), allocatable  :: organic(:,:)       ! organic matter density (kg/m3)
  real(r8), allocatable  :: gdp(:)             ! GDP (x1000 1995 US$/capita)
  real(r8), allocatable  :: fpeat(:)           ! peatland fraction of gridcell
  integer , allocatable  :: agfirepkmon(:)     ! agricultural fire peak month
  integer , allocatable  :: urban_region(:)    ! urban region ID
  real(r8), allocatable  :: topo_stddev(:)     ! standard deviation of elevation (m)
  real(r8), allocatable  :: slope(:)           ! topographic slope (degrees)
  real(r8), allocatable  :: vic_binfl(:)       ! VIC b parameter (unitless)
  real(r8), allocatable  :: vic_ws(:)          ! VIC Ws parameter (unitless)
  real(r8), allocatable  :: vic_dsmax(:)       ! VIC Dsmax parameter (mm/day)
  real(r8), allocatable  :: vic_ds(:)          ! VIC Ds parameter (unitless)
  real(r8), allocatable  :: lakedepth(:)       ! lake depth (m)
  real(r8), allocatable  :: f0(:)              ! max fractional inundated area (unitless)
  real(r8), allocatable  :: p3(:)              ! coefficient for qflx_surf_lag for finundated (s/mm)
  real(r8), allocatable  :: zwt0(:)            ! decay factor for finundated (m)
  real(r8), allocatable  :: apatiteP(:)        ! apptite phosphorus
  real(r8), allocatable  :: labileP(:)         ! labile phosphorus
  real(r8), allocatable  :: occludedP(:)       ! occluded phosphorus
  real(r8), allocatable  :: secondaryP(:)      ! secondaryP phosphorus
  real(r8), allocatable  :: grvl(:,:)          ! soil gravel content (percent)
  real(r8), allocatable  :: slp10(:,:)         ! slope percentile (km/km)
  real(r8), allocatable  :: ero_c1(:)          ! ELM-Erosion c1 parameter (unitless)
  real(r8), allocatable  :: ero_c2(:)          ! ELM-Erosion c2 parameter (unitless)
  real(r8), allocatable  :: ero_c3(:)          ! ELM-Erosion c3 parameter (unitless)
  real(r8), allocatable  :: tillage(:)         ! conserved tillage fraction (fraction)
  real(r8), allocatable  :: litho(:)           ! lithology erodiblity index (unitless)

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

  ! Reads namelist
  call setup_namelist()

  ! Initializes modules
  call mksoilInit()
  call mkpftInit(all_urban, all_veg)
  allocate ( elevclass(nglcec+1) )
  call mkglcmecInit (elevclass)
  call mkurbanInit (mksrf_furban)
  
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

  ! Make inland water [pctlak, pctwet] [flakwat] [fwetlnd]

  call mklakwat (ldomain, mapfname=map_flakwat, datfname=mksrf_flakwat, &
       ndiag=ndiag, zero_out=all_urban.or.all_veg, lake_o=pctlak)

  call mkwetlnd (ldomain, mapfname=map_fwetlnd, datfname=mksrf_fwetlnd, &
       ndiag=ndiag, zero_out=all_urban.or.all_veg.or.no_inlandwet, swmp_o=pctwet)

  ! Make glacier fraction [pctgla] from [fglacier] dataset

  call mkglacier (ldomain, mapfname=map_fglacier, datfname=mksrf_fglacier, &
       ndiag=ndiag, zero_out=all_urban.or.all_veg, glac_o=pctgla)

  ! Make soil texture [pctsand, pctclay]  [fsoitex]

  call mksoiltex (ldomain, mapfname=map_fsoitex, datfname=mksrf_fsoitex, &
       ndiag=ndiag, sand_o=pctsand, clay_o=pctclay)

  ! Make soil color classes [soicol] [fsoicol]

  call mksoilcol (ldomain, mapfname=map_fsoicol, datfname=mksrf_fsoicol, &
       ndiag=ndiag, soil_color_o=soicol, nsoicol=nsoicol)

  ! Make soil order classes [soiord] [fsoiord]

  call mksoilord (ldomain, mapfname=map_fsoiord, datfname=mksrf_fsoiord, &
       ndiag=ndiag, pctglac_o=pctgla, soil_order_o=soiord, nsoiord=nsoiord)

  ! Make fmax [fmax] from [fmax] dataset

  call mkfmax (ldomain, mapfname=map_fmax, datfname=mksrf_fmax, &
       ndiag=ndiag, fmax_o=fmax)

  ! ----------------------------------------------------------------------
  ! deallocate memory for all variables
  ! ----------------------------------------------------------------------
  call deallocate_memory()

  ! Finalize PETSc
  PetscCallA(PetscFinalize(ierr))

contains

  !-----------------------------------------------------------------------
  subroutine setup_namelist()

    implicit none

    ! Default settings
    mksrf_gridtype    = 'global'

    mksrf_fgrid       = ' '
    mksrf_fvegtyp     = ' '
    mksrf_fsoitex     = ' '
    mksrf_forganic    = ' '
    mksrf_fsoicol     = ' '
    mksrf_fsoiord     = ' '
    mksrf_fvocef      = ' '
    mksrf_flakwat     = ' '
    mksrf_fwetlnd     = ' '
    mksrf_fglacier    = ' '
    mksrf_furbtopo    = ' '
    mksrf_flndtopo    = ' '
    mksrf_fmax        = ' '
    mksrf_furban      = ' '
    mksrf_flai        = ' '
    mksrf_fdynuse     = ' '
    mksrf_fgdp        = ' '
    mksrf_fpeat       = ' '
    mksrf_fabm        = ' '
    mksrf_ftopostats  = ' '
    mksrf_fvic        = ' '
    mksrf_fch4        = ' '
    mksrf_fphosphorus = ' '
    mksrf_fgrvl       = ' '
    mksrf_fslp10      = ' '
    mksrf_fero        = ' '

    map_flakwat     = ' '
    map_fwetlnd     = ' '
    map_fglacier    = ' '
    map_fsoitex     = ' '
    map_fsoicol     = ' '
    map_fsoiord     = ' '
    map_furban      = ' '
    map_furbtopo    = ' '
    map_flndtopo    = ' '
    map_fmax        = ' '
    map_forganic    = ' '
    map_fvocef      = ' '
    map_flai        = ' '
    map_fharvest    = ' '
    map_fgdp        = ' '
    map_fpeat       = ' '
    map_fabm        = ' '
    map_ftopostats  = ' '
    map_fvic        = ' '
    map_fch4        = ' '
    map_fphosphorus = ' '
    map_fgrvl       = ' '
    map_fslp10      = ' '
    map_fero        = ' '

    fsurlog        = ' '
    map_fpft       = ' '
    mksrf_fvegtyp  = ' '

    outnc_large_files = .false.
    outnc_double      = .true.
    all_urban         = .false.
    no_inlandwet      = .true.

    read(5, elmexp, iostat=ier)
    if (ier /= 0) then
       write(6,*)'error: namelist input resulted in error code ',ier
       call abort()
    endif

    call check_namelist_variable(mksrf_fgrid       ,'mksrf_fgrid'       )
    call check_namelist_variable(mksrf_fvegtyp     ,'mksrf_fvegtyp'     )
    call check_namelist_variable(mksrf_fsoitex     ,'mksrf_fsoitex'     )
    call check_namelist_variable(mksrf_forganic    ,'mksrf_forganic'    )
    call check_namelist_variable(mksrf_fsoicol     ,'mksrf_fsoicol'     )
    call check_namelist_variable(mksrf_fsoiord     ,'mksrf_fsoiord'     )
    call check_namelist_variable(mksrf_fvocef      ,'mksrf_fvocef'      )
    call check_namelist_variable(mksrf_flakwat     ,'mksrf_flakwat'     )
    call check_namelist_variable(mksrf_fwetlnd     ,'mksrf_fwetlnd'     )
    call check_namelist_variable(mksrf_fglacier    ,'mksrf_fglacier'    )
    call check_namelist_variable(mksrf_furbtopo    ,'mksrf_furbtopo'    )
    call check_namelist_variable(mksrf_flndtopo    ,'mksrf_flndtopo'    )
    call check_namelist_variable(mksrf_fmax        ,'mksrf_fmax'        )
    call check_namelist_variable(mksrf_furban      ,'mksrf_furban'      )
    call check_namelist_variable(mksrf_flai        ,'mksrf_flai'        )
    call check_namelist_variable(mksrf_fgdp        ,'mksrf_fgdp'        )
    call check_namelist_variable(mksrf_fpeat       ,'mksrf_fpeat'       )
    call check_namelist_variable(mksrf_fabm        ,'mksrf_fabm'        )
    call check_namelist_variable(mksrf_ftopostats  ,'mksrf_ftopostats'  )
    call check_namelist_variable(mksrf_fvic        ,'mksrf_fvic'        )
    call check_namelist_variable(mksrf_fch4        ,'mksrf_fch4'        )
    call check_namelist_variable(mksrf_fphosphorus ,'mksrf_fphosphorus' )
    call check_namelist_variable(mksrf_fgrvl       ,'mksrf_fgrvl'       )
    call check_namelist_variable(mksrf_fslp10      ,'mksrf_fslp10'      )
    call check_namelist_variable(mksrf_fero        ,'mksrf_fero'        )

    call check_namelist_variable(map_flakwat       ,'map_flakwat'       )
    call check_namelist_variable(map_fwetlnd       ,'map_fwetlnd'       )
    call check_namelist_variable(map_fglacier      ,'map_fglacier'      )
    call check_namelist_variable(map_fsoitex       ,'map_fsoitex'       )
    call check_namelist_variable(map_fsoicol       ,'map_fsoicol'       )
    call check_namelist_variable(map_fsoiord       ,'map_fsoiord'       )
    call check_namelist_variable(map_furban        ,'map_furban'        )
    call check_namelist_variable(map_furbtopo      ,'map_furbtopo'      )
    call check_namelist_variable(map_flndtopo      ,'map_flndtopo'      )
    call check_namelist_variable(map_fmax          ,'map_fmax'          )
    call check_namelist_variable(map_forganic      ,'map_forganic'      )
    call check_namelist_variable(map_fvocef        ,'map_fvocef'        )
    call check_namelist_variable(map_flai          ,'map_flai'          )
    call check_namelist_variable(map_fharvest      ,'map_fharvest'      )
    call check_namelist_variable(map_fgdp          ,'map_fgdp'          )
    call check_namelist_variable(map_fpeat         ,'map_fpeat'         )
    call check_namelist_variable(map_fabm          ,'map_fabm'          )
    call check_namelist_variable(map_ftopostats    ,'map_ftopostats'    )
    call check_namelist_variable(map_fvic          ,'map_fvic'          )
    call check_namelist_variable(map_fch4          ,'map_fch4'          )
    call check_namelist_variable(map_fphosphorus   ,'map_fphosphorus'   )
    call check_namelist_variable(map_fgrvl         ,'map_fgrvl'         )
    call check_namelist_variable(map_fslp10        ,'map_fslp10'        )
    call check_namelist_variable(map_fero          ,'map_fero'          )

    call check_namelist_variable(fsurlog, 'fsurlog')
    ndiag = getavu()
    call opnfil(fsurlog, ndiag, 'f')


    if (trim(mksrf_gridtype) == 'global' .or. &
         trim(mksrf_gridtype) == 'regional') then
       write(6,*)'mksrf_gridtype = ',trim(mksrf_gridtype)
    else
       write(6,*)'mksrf_gridtype = ',trim(mksrf_gridtype)
       write (6,*)'illegal mksrf_gridtype, must be global or regional '
       call abort()
    endif

    if ( outnc_large_files )then
       write(6,*)'Output file in NetCDF 64-bit large_files format'
    end if
    if ( outnc_double )then
       write(6,*)'Output ALL data in file as 64-bit'
    end if
    if ( all_urban )then
       write(6,*) 'Output ALL data in file as 100% urban'
    end if
    if ( no_inlandwet )then
       write(6,*) 'Set wetland to 0% over land'
    end if

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

    allocate ( landfrac_pft(ns_o)                 , &
               pctlnd_pft(ns_o)                   , &
               pftdata_mask(ns_o)                 , &
               pctpft_full(ns_o,0:numpft)         , &
               pctnatveg(ns_o)                    , &
               pctcrop(ns_o)                      , &
               pctnatpft(ns_o,natpft_lb:natpft_ub), &
               pctcft(ns_o,cft_lb:cft_ub)         , &
               pctgla(ns_o)                       , &
               pctlak(ns_o)                       , &
               pctwet(ns_o)                       , &
               pcturb(ns_o)                       , &
               urban_region(ns_o)                 , &
               urbn_classes(ns_o,numurbl)         , &
               urbn_classes_g(ns_o,numurbl)       , &
               pctsand(ns_o,nlevsoi)              , &
               pctclay(ns_o,nlevsoi)              , &
               soicol(ns_o)                       , &
               soiord(ns_o)                       , &
               gdp(ns_o)                          , &
               fpeat(ns_o)                        , &
               agfirepkmon(ns_o)                  , &
               topo_stddev(ns_o)                  , &
               slope(ns_o)                        , &
               vic_binfl(ns_o)                    , &
               vic_ws(ns_o)                       , &
               vic_dsmax(ns_o)                    , &
               vic_ds(ns_o)                       , &
               lakedepth(ns_o)                    , &
               f0(ns_o)                           , &
               p3(ns_o)                           , &
               zwt0(ns_o)                         , &
               apatiteP(ns_o)                     , &
               labileP(ns_o)                      , &
               occludedP(ns_o)                    , &
               secondaryP(ns_o)                   , &
               grvl(ns_o,nlevsoi)                 , &
               slp10(ns_o,nlevslp)                , &
               ero_c1(ns_o)                       , &
               ero_c2(ns_o)                       , &
               ero_c3(ns_o)                       , &
               tillage(ns_o)                      , &
               litho(ns_o)                        , &
               fmax(ns_o)                           &
               )

    landfrac_pft(:)       = spval
    pctlnd_pft(:)         = spval
    pftdata_mask(:)       = -999
    pctpft_full(:,:)      = spval
    pctnatveg(:)          = spval
    pctcrop(:)            = spval
    pctnatpft(:,:)        = spval
    pctcft(:,:)           = spval
    pctgla(:)             = spval
    pctlak(:)             = spval
    pctwet(:)             = spval
    pcturb(:)             = spval
    urban_region(:)       = -999
    urbn_classes(:,:)     = spval
    urbn_classes_g(:,:)   = spval
    pctsand(:,:)          = spval
    pctclay(:,:)          = spval
    soicol(:)             = -999
    soiord(:)             = -999
    gdp(:)                = spval
    fpeat(:)              = spval
    agfirepkmon(:)        = -999
    topo_stddev(:)        = spval
    slope(:)              = spval
    vic_binfl(:)          = spval
    vic_ws(:)             = spval
    vic_dsmax(:)          = spval
    vic_ds(:)             = spval
    lakedepth(:)          = spval
    f0(:)                 = spval
    p3(:)                 = spval
    zwt0(:)               = spval
    apatiteP(:)           = spval
    labileP(:)            = spval
    occludedP(:)          = spval
    secondaryP(:)         = spval
    grvl(:,:)             = spval
    slp10(:,:)            = spval
    ero_c1(:)             = spval
    ero_c2(:)             = spval
    ero_c3(:)             = spval
    tillage(:)            = spval
    litho(:)              = spval
    fmax(:)               = spval

  end subroutine allocate_memory

  !-----------------------------------------------------------------------
  subroutine deallocate_memory()

    implicit none

    deallocate(pctlnd_pft)
    deallocate(pctpft_full)

  end subroutine deallocate_memory

end program mksurfdat_petsc


