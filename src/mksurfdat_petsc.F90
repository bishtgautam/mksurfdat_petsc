program mksurfdat_petsc

  use shr_kind_mod , only : r8 => shr_kind_r8, r4 => shr_kind_r4
  use mkdomainMod  , only : domain_type, domain_read_map, domain_read, &
       domain_write
  use mkpftMod     , only : mkpft
  use mkdomainPIOMod, only : domain_pio_type, domain_read_map_pio, domain_read_pio
  use mkglcmecMod
  use mksoilMod
  use mkpftMod
  use mkpftPIOMod, only : mkpft_pio
  use mkvarctl
  use mkurbanparMod
  use mkurbanparCommonMod
  use mkvarpar
  use mklanwatMod
  use mkgdpMod
  use mkpeatMod
  use mkagfirepkmonthMod
  use mkVICparamsMod
  use mkCH4inversionMod
  use mkvocefMod
  use mkSedMod
  use mksoilphosphorusMod
  use mkutilsMod
  use petsc
  use fileutils
  use spmdMod
  use mkdataPIOMod
  use mkfileMod

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

  real(r8), pointer  :: landfrac_pft(:)    ! PFT data: % land per gridcell
  real(r8), pointer  :: pctlnd_pft(:)      ! PFT data: % of gridcell for PFTs
  !real(r8), pointer  :: pctlnd_pft_dyn(:)  ! PFT data: % of gridcell for dyn landuse PFTs
  integer , pointer  :: pftdata_mask(:)    ! mask indicating real or fake land type
  real(r8), pointer      :: pctpft_full(:,:)   ! PFT data: % cover of each pft and cft on the vegetated landunits
  ! ('full' denotes inclusion of CFTs as well as natural PFTs in this array)
  real(r8), pointer  :: pctnatveg(:)       ! percent of grid cell that is natural veg landunit
  real(r8), pointer  :: pctcrop(:)         ! percent of grid cell that is crop landunit
  real(r8), pointer  :: pctnatpft(:,:)     ! % of each pft on the natural veg landunit (adds to 100%)
  real(r8), pointer  :: pctcft(:,:)        ! % of each cft on the crop landunit (adds to 100%)
  !real(r8), pointer      :: harvest(:,:)       ! harvest data: normalized harvesting
  real(r8), pointer  :: pctgla(:)          ! percent of grid cell that is glacier
  real(r8), pointer  :: pctglc_gic(:)      ! percent of grid cell that is gic (% of glc landunit)
  real(r8), pointer  :: pctglc_icesheet(:) ! percent of grid cell that is ice sheet (% of glc landunit)
  real(r8), pointer  :: pctglcmec(:,:)     ! glacier_mec pct coverage in each class (% of landunit)
  real(r8), pointer  :: topoglcmec(:,:)    ! glacier_mec sfc elevation in each gridcell and class
  !real(r8), pointer  :: pctglcmec_gic(:,:) ! GIC pct coverage in each class (% of landunit)
  !real(r8), pointer  :: pctglcmec_icesheet(:,:) ! icesheet pct coverage in each class (% of landunit)
  real(r8), pointer  :: elevclass(:)       ! glacier_mec elevation classes
  real(r8), pointer  :: pctlak(:)          ! percent of grid cell that is lake
  real(r8), pointer  :: pctwet(:)          ! percent of grid cell that is wetland
  real(r8), pointer  :: pcturb(:)          ! percent of grid cell that is urbanized (total across all urban classes)
  real(r8), pointer  :: urbn_classes(:,:)  ! percent cover of each urban class, as % of total urban area
  real(r8), pointer  :: urbn_classes_g(:,:)! percent cover of each urban class, as % of grid cell
  real(r8), pointer  :: elev(:)            ! glc elevation (m)
  real(r8), pointer  :: topo(:)            ! land elevation (m)
  real(r8), pointer  :: fmax(:)            ! fractional saturated area
  integer , pointer  :: soicol(:)          ! soil color
  integer , pointer  :: soiord(:)          ! soil order
  real(r8), pointer  :: pctsand(:,:)       ! soil texture: percent sand
  real(r8), pointer  :: pctclay(:,:)       ! soil texture: percent clay
  real(r8), pointer  :: ef1_btr(:)         ! Isoprene emission factor for broadleaf
  real(r8), pointer  :: ef1_fet(:)         ! Isoprene emission factor for fine/everg
  real(r8), pointer  :: ef1_fdt(:)         ! Isoprene emission factor for fine/dec
  real(r8), pointer  :: ef1_shr(:)         ! Isoprene emission factor for shrubs
  real(r8), pointer  :: ef1_grs(:)         ! Isoprene emission factor for grasses
  real(r8), pointer  :: ef1_crp(:)         ! Isoprene emission factor for crops
  real(r8), pointer  :: organic(:,:)       ! organic matter density (kg/m3)
  real(r8), pointer  :: gdp(:)             ! GDP (x1000 1995 US$/capita)
  real(r8), pointer  :: fpeat(:)           ! peatland fraction of gridcell
  integer , pointer  :: agfirepkmon(:)     ! agricultural fire peak month
  integer , pointer  :: urban_region(:)    ! urban region ID
  real(r8), pointer  :: topo_stddev(:)     ! standard deviation of elevation (m)
  real(r8), pointer  :: slope(:)           ! topographic slope (degrees)
  real(r8), pointer  :: vic_binfl(:)       ! VIC b parameter (unitless)
  real(r8), pointer  :: vic_ws(:)          ! VIC Ws parameter (unitless)
  real(r8), pointer  :: vic_dsmax(:)       ! VIC Dsmax parameter (mm/day)
  real(r8), pointer  :: vic_ds(:)          ! VIC Ds parameter (unitless)
  real(r8), pointer  :: lakedepth(:)       ! lake depth (m)
  real(r8), pointer  :: f0(:)              ! max fractional inundated area (unitless)
  real(r8), pointer  :: p3(:)              ! coefficient for qflx_surf_lag for finundated (s/mm)
  real(r8), pointer  :: zwt0(:)            ! decay factor for finundated (m)
  real(r8), pointer  :: apatiteP(:)        ! apptite phosphorus
  real(r8), pointer  :: labileP(:)         ! labile phosphorus
  real(r8), pointer  :: occludedP(:)       ! occluded phosphorus
  real(r8), pointer  :: secondaryP(:)      ! secondaryP phosphorus
  real(r8), pointer  :: grvl(:,:)          ! soil gravel content (percent)
  real(r8), pointer  :: slp10(:,:)         ! slope percentile (km/km)
  real(r8), pointer  :: ero_c1(:)          ! ELM-Erosion c1 parameter (unitless)
  real(r8), pointer  :: ero_c2(:)          ! ELM-Erosion c2 parameter (unitless)
  real(r8), pointer  :: ero_c3(:)          ! ELM-Erosion c3 parameter (unitless)
  real(r8), pointer  :: tillage(:)         ! conserved tillage fraction (fraction)
  real(r8), pointer  :: litho(:)           ! lithology erodiblity index (unitless)

  type(domain_type) :: ldomain
  type(domain_pio_type) :: ldomain_pio

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

  mpicom = PETSC_COMM_WORLD

  PetscCallMPIA(MPI_Comm_size(PETSC_COMM_WORLD,npes,ierr))
  PetscCallMPIA(MPI_Comm_rank(PETSC_COMM_WORLD,iam,ierr))
  masterproc = (iam == 0)

  if (masterproc) write(6,*) 'Attempting to initialize control settings .....'

  ! Reads namelist
  call setup_namelist()

  ! Initializes modules
  call mksoilInit()
  call mkpftInit(all_urban, all_veg)
  allocate ( elevclass(nglcec+1) )
  call mkglcmecInit (elevclass)
  call mkurbanInit (mksrf_furban)

  if (.not. domain_read_map_pio(ldomain_pio, mksrf_fgrid)) then
     if (masterproc) then
        write(6,*)'domain_read_map_pio returned FALSE. Add code to support this case'
        call abort()
     end if
  endif

  if (masterproc) write(6,*)'calling domain_read'
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
  call allocate_memory(ldomain_pio%ns_loc)

  ! ----------------------------------------------------------------------
  ! Make surface dataset fields
  ! ----------------------------------------------------------------------

  ! Make PFTs [pctpft_full] from dataset [fvegtyp]

  call mkpft_pio(ldomain_pio, mapfname=map_fpft, fpft=mksrf_fvegtyp, &
       ndiag=ndiag, pctlnd_o=pctlnd_pft, pctpft_o=pctpft_full )

  ! Make inland water [pctlak, pctwet] [flakwat] [fwetlnd]

  call mklakwat_pio (ldomain_pio, mapfname=map_flakwat, datfname=mksrf_flakwat, &
       ndiag=ndiag, zero_out=all_urban.or.all_veg, lake_o=pctlak)

  call mkwetlnd_pio (ldomain_pio, mapfname=map_fwetlnd, datfname=mksrf_fwetlnd, &
       ndiag=ndiag, zero_out=all_urban.or.all_veg.or.no_inlandwet, swmp_o=pctwet)

  ! Make glacier fraction [pctgla] from [fglacier] dataset

  call mkglacier_pio (ldomain_pio, mapfname=map_fglacier, datfname=mksrf_fglacier, &
       ndiag=ndiag, zero_out=all_urban.or.all_veg, glac_o=pctgla)

  ! Make soil texture [pctsand, pctclay]  [fsoitex]

  if (npes == 1) then
     call mksoiltex (ldomain, mapfname=map_fsoitex, datfname=mksrf_fsoitex, &
          ndiag=ndiag, sand_o=pctsand, clay_o=pctclay)
  endif
  !TODO: call mksoiltex_pio (ldomain_pio, mapfname=map_fsoitex, datfname=mksrf_fsoitex, &
  !     ndiag=ndiag, sand_o=pctsand, clay_o=pctclay)

  ! Make soil color classes [soicol] [fsoicol]

  call mksoilcol_pio (ldomain_pio, mapfname=map_fsoicol, datfname=mksrf_fsoicol, &
       ndiag=ndiag, soil_color_o=soicol, nsoicol=nsoicol)

  ! Make soil order classes [soiord] [fsoiord]

  call mksoilord_pio (ldomain_pio, mapfname=map_fsoiord, datfname=mksrf_fsoiord, &
       ndiag=ndiag, pctglac_o=pctgla, soil_order_o=soiord, nsoiord=nsoiord)

  ! Make fmax [fmax] from [fmax] dataset

  call mkfmax_pio (ldomain_pio, mapfname=map_fmax, datfname=mksrf_fmax, &
       ndiag=ndiag, fmax_o=fmax)

  ! Make GDP data [gdp] from [gdp]

  call mkgdp_pio (ldomain_pio, mapfname=map_fgdp, datfname=mksrf_fgdp, &
       ndiag=ndiag, gdp_o=gdp)

  ! Make peat data [fpeat] from [peatf]

  call mkpeat_pio (ldomain_pio, mapfname=map_fpeat, datfname=mksrf_fpeat, &
       ndiag=ndiag, peat_o=fpeat)

  ! Make agricultural fire peak month data [abm] from [abm]

  call mkagfirepkmon_pio (ldomain_pio, mapfname=map_fabm, datfname=mksrf_fabm, &
       ndiag=ndiag, agfirepkmon_o=agfirepkmon)

  ! Make urban fraction [pcturb] from [furban] dataset

  call mkurban_pio (ldomain_pio, mapfname=map_furban, datfname=mksrf_furban, &
       ndiag=ndiag, zero_out=all_veg, urbn_o=pcturb, urbn_classes_o=urbn_classes, &
       region_o=urban_region)

  if ( .not. all_urban .and. .not. all_veg )then
     call mkelev_pio (ldomain_pio, mapfname=map_furbtopo, datfname=mksrf_furbtopo, &
          varname='TOPO_ICE', ndiag=ndiag, elev_o=elev);

     where (elev .gt. elev_thresh)
        pcturb = 0._r8
     end where
     deallocate(elev)
  end if

  ! Determine topography

  call mkelev_pio (ldomain_pio, mapfname=map_flndtopo, datfname=mksrf_flndtopo, &
       varname='TOPO', ndiag=ndiag, elev_o=topo);

  ! Make VIC parameters [binfl, ws, dsmax, ds] from [fvic]

  call mkVICparams_pio (ldomain_pio, mapfname=map_fvic, datfname=mksrf_fvic, ndiag=ndiag, &
       binfl_o=vic_binfl, ws_o=vic_ws, dsmax_o=vic_dsmax, ds_o=vic_ds)

  ! Make lake depth [lakedepth] from [flakwat]

  call mklakparams_pio (ldomain_pio, mapfname=map_flakwat, datfname=mksrf_flakwat, ndiag=ndiag, &
       lakedepth_o=lakedepth)

  ! Make inversion-derived CH4 parameters [f0, p3, zwt0] from [fch4]

  call mkCH4inversion_pio (ldomain_pio, mapfname=map_fch4, datfname=mksrf_fch4, ndiag=ndiag, &
       f0_o=f0, p3_o=p3, zwt0_o=zwt0)

  ! Make organic matter density [organic] [forganic]

  call mkorganic_pio (ldomain_pio, mapfname=map_forganic, datfname=mksrf_forganic, &
       ndiag=ndiag, organic_o=organic)

  ! Make VOC emission factors for isoprene &
  ! [ef1_btr,ef1_fet,ef1_fdt,ef1_shr,ef1_grs,ef1_crp]

  call mkvocef_pio (ldomain_pio, mapfname=map_fvocef, datfname=mksrf_fvocef, ndiag=ndiag, &
       ef_btr_o=ef1_btr, ef_fet_o=ef1_fet, ef_fdt_o=ef1_fdt,  &
       ef_shr_o=ef1_shr, ef_grs_o=ef1_grs, ef_crp_o=ef1_crp)

  call mksoilphosphorus_pio (ldomain_pio, mapfname=map_fphosphorus, datfname=mksrf_fphosphorus, &
       ndiag=ndiag, apatiteP_o=apatiteP, labileP_o=labileP, occludedP_o=occludedP, &
       secondaryP_o=secondaryP)

  call mkgrvl_pio(ldomain_pio, mapfname=map_fgrvl, datfname=mksrf_fgrvl, ndiag=ndiag, grvl_o=grvl)

  call mkslp10_pio(ldomain_pio, mapfname=map_fslp10, datfname=mksrf_fslp10, ndiag=ndiag, slp10_o=slp10)

  call mkEROparams_pio(ldomain_pio, mapfname=map_fero, datfname=mksrf_fero, ndiag=ndiag, &
       ero_c1_o=ero_c1, ero_c2_o=ero_c2, ero_c3_o=ero_c3, tillage_o=tillage, &
       litho_o=litho)

  ! Do landuse changes such as for the poles, etc.

  call change_landuse( ldomain_pio%ns_loc, dynpft=.false. )

  ! Perform few sanity checks and update the values, if needed
  call check_and_update_values()

  ! ----------------------------------------------------------------------
  ! write out the netcdf file
  ! ----------------------------------------------------------------------

  if (npes == 1) then
     call write_surface_dataset()
  endif
  call write_surface_dataset_pio()

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

    if (masterproc) then
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

    end if

    call mpi_bcast(mksrf_fgrid       , len(mksrf_fgrid)       , MPI_CHARACTER, 0, PETSC_COMM_WORLD, ier) 
    call mpi_bcast(mksrf_fvegtyp     , len(mksrf_fvegtyp)     , MPI_CHARACTER, 0, PETSC_COMM_WORLD, ier) 
    call mpi_bcast(mksrf_fsoitex     , len(mksrf_fsoitex)     , MPI_CHARACTER, 0, PETSC_COMM_WORLD, ier) 
    call mpi_bcast(mksrf_forganic    , len(mksrf_forganic)    , MPI_CHARACTER, 0, PETSC_COMM_WORLD, ier) 
    call mpi_bcast(mksrf_fsoicol     , len(mksrf_fsoicol)     , MPI_CHARACTER, 0, PETSC_COMM_WORLD, ier) 
    call mpi_bcast(mksrf_fsoiord     , len(mksrf_fsoiord)     , MPI_CHARACTER, 0, PETSC_COMM_WORLD, ier) 
    call mpi_bcast(mksrf_fvocef      , len(mksrf_fvocef)      , MPI_CHARACTER, 0, PETSC_COMM_WORLD, ier) 
    call mpi_bcast(mksrf_flakwat     , len(mksrf_flakwat)     , MPI_CHARACTER, 0, PETSC_COMM_WORLD, ier) 
    call mpi_bcast(mksrf_fwetlnd     , len(mksrf_fwetlnd)     , MPI_CHARACTER, 0, PETSC_COMM_WORLD, ier) 
    call mpi_bcast(mksrf_fglacier    , len(mksrf_fglacier)    , MPI_CHARACTER, 0, PETSC_COMM_WORLD, ier) 
    call mpi_bcast(mksrf_furbtopo    , len(mksrf_furbtopo)    , MPI_CHARACTER, 0, PETSC_COMM_WORLD, ier) 
    call mpi_bcast(mksrf_flndtopo    , len(mksrf_flndtopo)    , MPI_CHARACTER, 0, PETSC_COMM_WORLD, ier) 
    call mpi_bcast(mksrf_fmax        , len(mksrf_fmax)        , MPI_CHARACTER, 0, PETSC_COMM_WORLD, ier) 
    call mpi_bcast(mksrf_furban      , len(mksrf_furban)      , MPI_CHARACTER, 0, PETSC_COMM_WORLD, ier) 
    call mpi_bcast(mksrf_flai        , len(mksrf_flai)        , MPI_CHARACTER, 0, PETSC_COMM_WORLD, ier) 
    call mpi_bcast(mksrf_fgdp        , len(mksrf_fgdp)        , MPI_CHARACTER, 0, PETSC_COMM_WORLD, ier) 
    call mpi_bcast(mksrf_fpeat       , len(mksrf_fpeat)       , MPI_CHARACTER, 0, PETSC_COMM_WORLD, ier) 
    call mpi_bcast(mksrf_fabm        , len(mksrf_fabm)        , MPI_CHARACTER, 0, PETSC_COMM_WORLD, ier) 
    call mpi_bcast(mksrf_ftopostats  , len(mksrf_ftopostats)  , MPI_CHARACTER, 0, PETSC_COMM_WORLD, ier) 
    call mpi_bcast(mksrf_fvic        , len(mksrf_fvic)        , MPI_CHARACTER, 0, PETSC_COMM_WORLD, ier) 
    call mpi_bcast(mksrf_fch4        , len(mksrf_fch4)        , MPI_CHARACTER, 0, PETSC_COMM_WORLD, ier) 
    call mpi_bcast(mksrf_fphosphorus , len(mksrf_fphosphorus) , MPI_CHARACTER, 0, PETSC_COMM_WORLD, ier) 
    call mpi_bcast(mksrf_fgrvl       , len(mksrf_fgrvl)       , MPI_CHARACTER, 0, PETSC_COMM_WORLD, ier) 
    call mpi_bcast(mksrf_fslp10      , len(mksrf_fslp10)      , MPI_CHARACTER, 0, PETSC_COMM_WORLD, ier) 
    call mpi_bcast(mksrf_fero        , len(mksrf_fero)        , MPI_CHARACTER, 0, PETSC_COMM_WORLD, ier) 
    call mpi_bcast(map_fpft          , len(map_fpft)          , MPI_CHARACTER, 0, PETSC_COMM_WORLD, ier)
    call mpi_bcast(map_flakwat       , len(map_flakwat)       , MPI_CHARACTER, 0, PETSC_COMM_WORLD, ier) 
    call mpi_bcast(map_fwetlnd       , len(map_fwetlnd)       , MPI_CHARACTER, 0, PETSC_COMM_WORLD, ier) 
    call mpi_bcast(map_fglacier      , len(map_fglacier)      , MPI_CHARACTER, 0, PETSC_COMM_WORLD, ier) 
    call mpi_bcast(map_fsoitex       , len(map_fsoitex)       , MPI_CHARACTER, 0, PETSC_COMM_WORLD, ier) 
    call mpi_bcast(map_fsoicol       , len(map_fsoicol)       , MPI_CHARACTER, 0, PETSC_COMM_WORLD, ier) 
    call mpi_bcast(map_fsoiord       , len(map_fsoiord)       , MPI_CHARACTER, 0, PETSC_COMM_WORLD, ier) 
    call mpi_bcast(map_furban        , len(map_furban)        , MPI_CHARACTER, 0, PETSC_COMM_WORLD, ier) 
    call mpi_bcast(map_furbtopo      , len(map_furbtopo)      , MPI_CHARACTER, 0, PETSC_COMM_WORLD, ier) 
    call mpi_bcast(map_flndtopo      , len(map_flndtopo)      , MPI_CHARACTER, 0, PETSC_COMM_WORLD, ier) 
    call mpi_bcast(map_fmax          , len(map_fmax)          , MPI_CHARACTER, 0, PETSC_COMM_WORLD, ier) 
    call mpi_bcast(map_forganic      , len(map_forganic)      , MPI_CHARACTER, 0, PETSC_COMM_WORLD, ier) 
    call mpi_bcast(map_fvocef        , len(map_fvocef)        , MPI_CHARACTER, 0, PETSC_COMM_WORLD, ier) 
    call mpi_bcast(map_flai          , len(map_flai)          , MPI_CHARACTER, 0, PETSC_COMM_WORLD, ier) 
    call mpi_bcast(map_fharvest      , len(map_fharvest)      , MPI_CHARACTER, 0, PETSC_COMM_WORLD, ier) 
    call mpi_bcast(map_fgdp          , len(map_fgdp)          , MPI_CHARACTER, 0, PETSC_COMM_WORLD, ier) 
    call mpi_bcast(map_fpeat         , len(map_fpeat)         , MPI_CHARACTER, 0, PETSC_COMM_WORLD, ier) 
    call mpi_bcast(map_fabm          , len(map_fabm)          , MPI_CHARACTER, 0, PETSC_COMM_WORLD, ier) 
    call mpi_bcast(map_ftopostats    , len(map_ftopostats)    , MPI_CHARACTER, 0, PETSC_COMM_WORLD, ier) 
    call mpi_bcast(map_fvic          , len(map_fvic)          , MPI_CHARACTER, 0, PETSC_COMM_WORLD, ier) 
    call mpi_bcast(map_fch4          , len(map_fch4)          , MPI_CHARACTER, 0, PETSC_COMM_WORLD, ier) 
    call mpi_bcast(map_fphosphorus   , len(map_fphosphorus)   , MPI_CHARACTER, 0, PETSC_COMM_WORLD, ier) 
    call mpi_bcast(map_fgrvl         , len(map_fgrvl)         , MPI_CHARACTER, 0, PETSC_COMM_WORLD, ier) 
    call mpi_bcast(map_fslp10        , len(map_fslp10)        , MPI_CHARACTER, 0, PETSC_COMM_WORLD, ier) 
    call mpi_bcast(map_fero          , len(map_fero)          , MPI_CHARACTER, 0, PETSC_COMM_WORLD, ier) 

    call mpi_bcast(fsurdat           , len(fsurdat)           , MPI_CHARACTER, 0, PETSC_COMM_WORLD, ier)

    call mpi_bcast(nglcec            , 1                      , MPI_INT      , 0, PETSC_COMM_WORLD, ier)
    call mpi_bcast(numpft            , 1                      , MPI_INT      , 0, PETSC_COMM_WORLD, ier)
    call mpi_bcast(soil_color        , 1                      , MPI_INT      , 0, PETSC_COMM_WORLD, ier)
    call mpi_bcast(soil_order        , 1                      , MPI_INT      , 0, PETSC_COMM_WORLD, ier)

    call mpi_bcast(soil_sand         , 1                      , MPI_DOUBLE   , 0, PETSC_COMM_WORLD, ier)

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
  subroutine allocate_memory(ns_o)

    implicit none

    integer, intent (in) :: ns_o

    allocate ( landfrac_pft(ns_o)                  , &
         pctlnd_pft(ns_o)                    , &
         pftdata_mask(ns_o)                  , &
         pctpft_full(ns_o,0:numpft)          , &
         pctnatveg(ns_o)                     , &
         pctcrop(ns_o)                       , &
         pctnatpft(ns_o,natpft_lb:natpft_ub) , &
         pctcft(ns_o,cft_lb:cft_ub)          , &
         pctgla(ns_o)                        , &
         pctlak(ns_o)                        , &
         pctwet(ns_o)                        , &
         pcturb(ns_o)                        , &
         urban_region(ns_o)                  , &
         urbn_classes(ns_o,numurbl)          , &
         urbn_classes_g(ns_o,numurbl)        , &
         pctsand(ns_o,nlevsoi)               , &
         pctclay(ns_o,nlevsoi)               , &
         soicol(ns_o)                        , &
         soiord(ns_o)                        , &
         gdp(ns_o)                           , &
         fpeat(ns_o)                         , &
         agfirepkmon(ns_o)                   , &
         topo_stddev(ns_o)                   , &
         slope(ns_o)                         , &
         vic_binfl(ns_o)                     , &
         vic_ws(ns_o)                        , &
         vic_dsmax(ns_o)                     , &
         vic_ds(ns_o)                        , &
         lakedepth(ns_o)                     , &
         f0(ns_o)                            , &
         p3(ns_o)                            , &
         zwt0(ns_o)                          , &
         apatiteP(ns_o)                      , &
         labileP(ns_o)                       , &
         occludedP(ns_o)                     , &
         secondaryP(ns_o)                    , &
         grvl(ns_o,nlevsoi)                  , &
         slp10(ns_o,nlevslp)                 , &
         ero_c1(ns_o)                        , &
         ero_c2(ns_o)                        , &
         ero_c3(ns_o)                        , &
         tillage(ns_o)                       , &
         litho(ns_o)                         , &
         fmax(ns_o)                          , &
         topo(ns_o)                          , &
         organic(ns_o,nlevsoi)               , &
         ef1_btr(ns_o)                       , & 
         ef1_fet(ns_o)                       , & 
         ef1_fdt(ns_o)                       , & 
         ef1_shr(ns_o)                       , & 
         ef1_grs(ns_o)                       , & 
         ef1_crp(ns_o)                         &
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
    topo(:)               = spval
    organic(:,:)          = spval
    ef1_btr(:) = 0._r8
    ef1_fet(:) = 0._r8
    ef1_fdt(:) = 0._r8
    ef1_shr(:) = 0._r8
    ef1_grs(:) = 0._r8
    ef1_crp(:) = 0._r8

    if ( .not. all_urban .and. .not. all_veg )then
       allocate(elev(ns_o))
       elev(:) = spval
    end if

  end subroutine allocate_memory

  !-----------------------------------------------------------------------
  subroutine deallocate_memory()

    implicit none

    deallocate(pctlnd_pft)
    deallocate ( organic )
    deallocate ( ef1_btr, ef1_fet, ef1_fdt, ef1_shr, ef1_grs, ef1_crp )
    if ( nglcec > 0 ) deallocate ( pctglcmec, topoglcmec)
    if ( nglcec > 0 ) deallocate ( pctglc_gic, pctglc_icesheet)
    deallocate ( elevclass )
    deallocate ( fmax )
    deallocate ( pctsand, pctclay )
    deallocate ( soicol )
    deallocate ( soiord )
    deallocate ( gdp, fpeat, agfirepkmon )
    deallocate ( topo_stddev, slope )
    deallocate ( vic_binfl, vic_ws, vic_dsmax, vic_ds )
    deallocate ( lakedepth )
    deallocate ( f0, p3, zwt0 )
    deallocate ( apatiteP, labileP, occludedP, secondaryP )
    deallocate ( grvl, slp10 )
    deallocate ( ero_c1, ero_c2, ero_c3, tillage, litho )

  end subroutine deallocate_memory

  !-----------------------------------------------------------------------

  subroutine check_and_update_values()

    implicit none

    integer :: k, n, ns_o
    real(r8) :: suma
    character(len=32) :: subname = 'perform_sanity_check'

    ns_o = ldomain_pio%ns_loc

    do n = 1,ns_o

       ! Assume wetland and/or lake when dataset landmask implies ocean 
       ! (assume medium soil color (15), soil order(15) and loamy texture).

       ! Also set pftdata_mask here
       ! Note that pctpft_full is NOT adjusted here, so that we still have information
       ! about the landunit breakdown into PFTs (pctnatveg and pctcrop will later become 0
       ! due to wetland and lake adding to 100%).

       if (pctlnd_pft(n) < 1.e-6_r8) then
          pftdata_mask(n)  = 0
          soicol(n)        = 15
          soiord(n)        = 15
          pctwet(n)        = 100._r8 - pctlak(n)
          pcturb(n)        = 0._r8
          pctgla(n)        = 0._r8
          pctsand(n,:)     = 43._r8
          pctclay(n,:)     = 18._r8
          organic(n,:)   = 0._r8
          grvl(n,:)        = 0._r8
          slp10(n,:)       = 0._r8
          ero_c1(n)        = 0._r8
          ero_c2(n)        = 0._r8
          ero_c3(n)        = 0._r8
          tillage(n)       = 0._r8
          litho(n)         = 0._r8
       else
          pftdata_mask(n) = 1
       end if

       ! Truncate all percentage fields on output grid. This is needed to
       ! insure that wt is zero (not a very small number such as
       ! 1e-16) where it really should be zero

       do k = 1,nlevsoi
          pctsand(n,k) = float(nint(pctsand(n,k)))
          pctclay(n,k) = float(nint(pctclay(n,k)))
          grvl(n,k)    = float(nint(grvl(n,k)))
       end do
       pctlak(n) = float(nint(pctlak(n)))
       pctwet(n) = float(nint(pctwet(n)))
       pctgla(n) = float(nint(pctgla(n)))

       ! Make sure sum of land cover types does not exceed 100. If it does,
       ! subtract excess from most dominant land cover.

       suma = pctlak(n) + pctwet(n) + pcturb(n) + pctgla(n)
       if (suma > 250._r4) then
          write (6,*) subname, ' error: sum of pctlak, pctwet,', &
               'pcturb and pctgla is greater than 250%'
          write (6,*)'n,pctlak,pctwet,pcturb,pctgla= ', &
               n,pctlak(n),pctwet(n),pcturb(n),pctgla(n)
          call abort()
       else if (suma > 100._r4) then
          pctlak(n) = pctlak(n) * 100._r8/suma
          pctwet(n) = pctwet(n) * 100._r8/suma
          pcturb(n) = pcturb(n) * 100._r8/suma
          pctgla(n) = pctgla(n) * 100._r8/suma
       end if

    end do

    call normalizencheck_landuse(ldomain_pio%ns_loc)

    ! Write out sum of PFT's

    do k = natpft_lb,natpft_ub
       suma = 0._r8
       do n = 1,ns_o
          suma = suma + pctnatpft(n,k)*pctnatveg(n)/100._r8
       enddo
       write(6,*) 'sum over domain of pft ',k,suma
    enddo
    write(6,*)

    do k = cft_lb,cft_ub
       suma = 0._r8
       do n = 1,ns_o
          suma = suma + pctcft(n,k)*pctcrop(n)/100._r8
       enddo
       write(6,*) 'sum over domain of cft ',k,suma
    enddo
    write(6,*)

    ! Make final values of percent urban by class
    ! This call needs to occur after all corrections are made to pcturb

    call normalize_classes_by_gcell(urbn_classes, pcturb, urbn_classes_g)

    ! Make glacier multiple elevation classes [pctglcmec,topoglcmec] from [fglacier] dataset
    ! This call needs to occur after pctgla has been adjusted for the final time

    if ( nglcec > 0 )then
       write(6,*) 'error: add new code to extend for the case nglec > 0'
       call abort()
    end if

    ! Determine fractional land from pft dataset

    do n = 1,ns_o
       landfrac_pft(n) = pctlnd_pft(n)/100._r8
    end do

  end subroutine check_and_update_values

  !-----------------------------------------------------------------------
  !BOP
  !
  ! !IROUTINE: change_landuse
  !
  ! !INTERFACE:
  subroutine change_landuse( ns_o, dynpft )
    !
    ! !DESCRIPTION:
    !
    ! Do landuse changes such as for the poles, etc.
    !
    ! !USES:
    implicit none
    !
    ! !ARGUMENTS:
    integer, intent(in)  :: ns_o
    logical, intent(in)  :: dynpft   ! if part of the dynpft section of code

    !
    ! !REVISION HISTORY:
    ! 9/10/09: Erik Kluzek spin off subroutine from original embedded code
    !
    !EOP
    !
    ! !LOCAL VARIABLES:
    logical  :: first_time = .true.         ! flag if this is the first pass through or not
    integer ,parameter :: bdtroptree = 6    ! Index for broadleaf decidious tropical tree
    integer ,parameter :: bdtemptree = 7    ! Index for broadleaf decidious temperate tree
    integer ,parameter :: bdtempshrub = 10  ! Index for broadleaf decidious temperate shrub
    real(r8),parameter :: troplat = 23.5_r8 ! Latitude to define as tropical
    integer  :: n                           ! indices
    character(len=32) :: subname = 'change_landuse'  ! subroutine name
    !-----------------------------------------------------------------------

    do n = 1,ns_o

       ! Set pfts 7 and 10 to 6 in the tropics to avoid lais > 1000
       ! Using P. Thornton's method found in surfrdMod.F90 in clm3.5

       if (abs(ldomain%latc(n))<troplat .and. pctpft_full(n,bdtemptree)>0._r8) then
          pctpft_full(n,bdtroptree) = pctpft_full(n,bdtroptree) + pctpft_full(n,bdtemptree)
          pctpft_full(n,bdtemptree) = 0._r8
          if ( first_time ) write (6,*) subname, ' Warning: all wgt of pft ', &
               bdtemptree, ' now added to pft ', bdtroptree
       end if
       if (abs(ldomain%latc(n))<troplat .and. pctpft_full(n,bdtempshrub)>0._r8) then
          pctpft_full(n,bdtroptree) = pctpft_full(n,bdtroptree) + pctpft_full(n,bdtempshrub)
          pctpft_full(n,bdtempshrub) = 0._r8
          if ( first_time ) write (6,*) subname, ' Warning: all wgt of pft ', &
               bdtempshrub, ' now added to pft ', bdtroptree
       end if
       first_time = .false.

       ! If have pole points on grid - set south pole to glacier
       ! north pole is assumed as non-land
       ! Note that pctpft_full is NOT adjusted here, so that we still have information
       ! about the landunit breakdown into PFTs (pctnatveg and pctcrop will later become 0
       ! due to glacier adding to 100%).

       if (abs((ldomain%latc(n) - 90._r8)) < 1.e-6_r8) then
          pctlak(n)   = 0._r8
          pctwet(n)   = 0._r8
          pcturb(n)   = 0._r8
          pctgla(n)   = 100._r8
          if ( .not. dynpft )then
             organic(n,:)   = 0._r8
             ef1_btr(n)     = 0._r8
             ef1_fet(n)     = 0._r8
             ef1_fdt(n)     = 0._r8
             ef1_shr(n)     = 0._r8
             ef1_grs(n)     = 0._r8
             ef1_crp(n)     = 0._r8
          end if
       end if

    end do

  end subroutine change_landuse

  !-----------------------------------------------------------------------
  !BOP
  !
  ! !IROUTINE: normalizencheck_landuse
  !
  ! !INTERFACE:
  subroutine normalizencheck_landuse(ns_o)
    !
    ! !DESCRIPTION:
    !
    ! Normalize land use and make sure things add up to 100% as well as
    ! checking that things are as they should be.
    !
    ! Note: this deallocates pctpft_full
    !
    ! !USES:
    use mkutilsMod, only : remove_small_cover
    use mkpftMod  , only : mkpft_normalize
    implicit none
    ! !ARGUMENTS:
    integer :: ns_o
    !
    ! !REVISION HISTORY:
    ! 9/10/09: Erik Kluzek spin off subroutine from original embedded code
    !
    !EOP
    !
    ! !LOCAL VARIABLES:
    integer  :: m,k,n                       ! indices
    integer  :: nsmall                      ! number of small PFT values for a single check
    integer  :: nsmall_tot                  ! total number of small PFT values in all grid cells
    real(r8) :: suma                        ! sum for error check
    real(r8) :: suma2                       ! another sum for error check
    real(r8) :: pcturb_rescaled             ! % urban normalized as if it were on the veg landunit
    real(r8) :: bare_urb_diff               ! difference between bare soil and urban %
    real(r8) :: pcturb_excess               ! excess urban % not accounted for by bare soil
    real(r8) :: sumpft                      ! sum of non-baresoil pfts
    real(r8) :: sum8, sum8a                 ! sum for error check
    real(r4) :: sum4a                       ! sum for error check
    real(r8), parameter :: toosmallPFT = 1.e-10_r8            ! tolerance for PFT's to ignore
    character(len=32) :: subname = 'normalizencheck_landuse'  ! subroutine name
    !-----------------------------------------------------------------------

    do n = 1,ns_o

       ! The following corrects pctpft_full to account for urban area. Urban needs to be
       ! handled specially because we replace bare soil preferentially with urban, rather
       ! than rescaling all PFTs equally. To accomplish this, we first rescale urban to
       ! the value it would have if it were in the pctpft_full array (pcturb_rescaled), which
       ! gives the fraction of each pft on the (vegetated + crop) landunits. Before
       ! (conceptually) adding urban to this array, pctpft_full added up to 100% of the
       ! (vegetated + crop) landunits, so we preserve this concept by adjusting pctpft_full so
       ! that sum(pctpft_full) + pcturb_rescaled = 100%. After doing this, we can rescale
       ! pctpft_full so it accounts for the other special landunits - but since we have already
       ! accounted for urban, note that that final rescaling only considers the sum of the
       ! non-urban special landunits.

       if (pcturb(n) .gt. 0._r8) then

          suma = pctlak(n)+pctwet(n)+pctgla(n)

          if (suma >= 100._r8) then
             write(6,*) 'ERROR: pcturb > 0 and other special landunits add to >= 100'
             write(6,*) 'n, pcturb, suma = ', n, pcturb(n), suma
             call abort()
          end if

          ! Rescale pcturb as if it were going on the (vegetated + crop) landunits
          pcturb_rescaled = pcturb(n)/(0.01_r8*(100._r8 - suma))

          ! Replace bare soil preferentially with urban
          bare_urb_diff = pctpft_full(n,0) - pcturb_rescaled
          pctpft_full(n,0) = max(0._r8,bare_urb_diff)
          pcturb_excess = abs(min(0._r8,bare_urb_diff))

          ! For any urban not accounted for by bare soil, replace other PFTs
          ! proportionally
          sumpft = sum(pctpft_full(n,1:numpft))
          if (pcturb_excess > 0._r8 .and. sumpft > 0._r8) then
             do m = 1, numpft
                pctpft_full(n,m) = pctpft_full(n,m) - pcturb_excess*pctpft_full(n,m)/sumpft
                if ( pctpft_full(n,m) < 0.0_r8 )then
                   write (6,*)'pctpft_full < 0.0 = ', pctpft_full(n,m), &
                        ' n, m, suma, pcturb_excess, sumpft = ',  n, m, suma, pcturb_excess, sumpft
                   ! Note that the tolerance here (0.00001_r8) matches the
                   ! tolerance for the error check on the sum in mkpftMod: mkpft
                   if ( abs(pctpft_full(n,m)) > 1.e-5_r8 )then
                      call abort()
                   end if
                   pctpft_full(n,m) = 0.0_r8
                end if
             end do
          end if

          ! Confirm that we have done the rescaling correctly: 
          ! make sure sum(pctpft_full) + pcturb_rescaled = 100%
          ! Note that the tolerance here (0.00001_r8) matches the tolerance for the error
          ! check on the sum in mkpftMod: mkpft
          if (abs( (sum(pctpft_full(n,0:numpft)) + pcturb_rescaled) - 100._r8 ) > 0.00001_r8) then
             write(6,*) 'ERROR: pctpft_full + pcturb_rescaled not equal to 100'
             write(6,*) 'n, sum(pctpft_full), pcturb_rescaled, diff = '
             write(6,*) n, sum(pctpft_full(n,0:numpft)), pcturb_rescaled, abs( (sum(pctpft_full(n,0:numpft)) + pcturb_rescaled) - 100._r8 )
             call abort()
          end if
       end if

       ! Normalize pctpft_full to be the remainder of [100 - (special landunits)]
       ! Note that we have already accounted for pcturb, above, so here we just account
       ! for the other special landunits
       suma = pctlak(n) + pctwet(n) + pctgla(n)
       call mkpft_normalize (pctpft_full(n,:), suma, pctnatveg(n), pctcrop(n), pctnatpft(n,:), pctcft(n,:))
    end do

    ! Deallocate no-longer-needed pctpft_full so it won't accidentally be updated further
    deallocate(pctpft_full)

    nsmall_tot = 0

    do n = 1,ns_o

       suma = pctlak(n) + pctwet(n) + pcturb(n) + pctgla(n) + pctnatveg(n) + pctcrop(n)

       if (suma < 90._r8) then
          write (6,*) subname, ' error: sum of pctlak, pctwet,', &
               'pcturb, pctgla, pctnatveg and pctcrop is less than 90'
          write (6,*)'n,pctlak,pctwet,pcturb,pctgla,pctnatveg,pctcrop= ', &
               n,pctlak(n),pctwet(n),pcturb(n),pctgla(n),pctnatveg(n),pctcrop(n)
          call abort()
       else if (suma > 100._r8 + 1.e-4_r8) then
          write (6,*) subname, ' error: sum of pctlak, pctwet,', &
               'pcturb, pctgla, pctnatveg and pctcrop is greater than 100'
          write (6,*)'n,pctlak,pctwet,pcturb,pctgla,pctnatveg,pctcrop,sum= ', &
               n,pctlak(n),pctwet(n),pcturb(n),pctgla(n),pctnatveg(n),pctcrop(n),suma
          call abort()
       else
          ! If the coverage of any PFT or CFT is too small at the gridcell level, set its
          ! % cover to 0, then renormalize everything else as needed
          call remove_small_cover(pctnatveg(n), pctnatpft(n,:), nsmall, toosmallPFT, suma)
          nsmall_tot = nsmall_tot + nsmall
          call remove_small_cover(pctcrop(n), pctcft(n,:), nsmall, toosmallPFT, suma)
          nsmall_tot = nsmall_tot + nsmall

          suma = pctlak(n) + pctwet(n) + pcturb(n) + pctgla(n) + pctnatveg(n) + pctcrop(n)
          if ( abs(suma - 100.0_r8) > 2.0*epsilon(suma) )then
             pctlak(n)    = pctlak(n)    * 100._r8/suma
             pctwet(n)    = pctwet(n)    * 100._r8/suma
             pcturb(n)    = pcturb(n)    * 100._r8/suma
             pctgla(n)    = pctgla(n)    * 100._r8/suma
             pctnatveg(n) = pctnatveg(n) * 100._r8/suma
             pctcrop(n)   = pctcrop(n)   * 100._r8/suma
          end if
       end if

       ! Roundoff error fix
       suma = pctlak(n) + pctwet(n) + pcturb(n) + pctgla(n)
       suma2 = pctnatveg(n) + pctcrop(n)
       if ( (suma < 100._r8 .and. suma > (100._r8 - 1.e-6_r8)) .or. &
            (suma2 > 0.0_r8 .and. suma2 <  1.e-6_r8) ) then
          write (6,*) 'Special land units near 100%, but not quite for n,suma =',n,suma
          write (6,*) 'Adjusting special land units to 100%'
          if (pctlak(n) >= 25._r8) then
             pctlak(n) = 100._r8 - (pctwet(n) + pcturb(n) + pctgla(n))
          else if (pctwet(n) >= 25._r8) then
             pctwet(n) = 100._r8 - (pctlak(n) + pcturb(n) + pctgla(n))
          else if (pcturb(n) >= 25._r8) then
             pcturb(n) = 100._r8 - (pctlak(n) + pctwet(n) + pctgla(n))
          else if (pctgla(n) >= 25._r8) then
             pctgla(n) = 100._r8 - (pctlak(n) + pctwet(n) + pcturb(n))
          else
             write (6,*) subname, 'Error: sum of special land units nearly 100% but none is >= 25% at ', &
                  'n,pctlak(n),pctwet(n),pcturb(n),pctgla(n),pctnatveg(n),pctcrop(n),suma = ', &
                  n,pctlak(n),pctwet(n),pcturb(n),pctgla(n),pctnatveg(n),pctcrop(n),suma
             call abort()
          end if
          pctnatveg(n) = 0._r8
          pctcrop(n) = 0._r8
       end if
       if ( any(pctnatpft(n,:)*pctnatveg(n) > 0.0_r8 .and. pctnatpft(n,:)*pctnatveg(n) < toosmallPFT ) .or. &
            any(pctcft(n,:)*pctcrop(n)   > 0.0_r8 .and. pctcft(n,:)*pctcrop(n)   < toosmallPFT )) then
          write (6,*) 'pctnatpft or pctcft is small at n=', n
          write (6,*) 'pctnatpft = ', pctnatpft(n,:)
          write (6,*) 'pctcft = ', pctcft(n,:)
          write (6,*) 'pctnatveg = ', pctnatveg
          write (6,*) 'pctcrop = ', pctcrop
          call abort()
       end if

       suma = pctlak(n) + pctwet(n) + pcturb(n) + pctgla(n)
       if (suma < 100._r8-epsilon(suma) .and. suma > (100._r8 - 4._r8*epsilon(suma))) then
          write (6,*) subname, 'n,pctlak,pctwet,pcturb,pctgla,pctnatveg,pctcrop= ', &
               n,pctlak(n),pctwet(n),pcturb(n),pctgla(n),pctnatveg(n),pctcrop(n)
          call abort()
       end if
       suma = suma + pctnatveg(n) + pctcrop(n)
       if ( abs(suma-100._r8) > 1.e-10_r8) then
          write (6,*) subname, ' error: sum of pctlak, pctwet,', &
               'pcturb, pctgla, pctnatveg and pctcrop is NOT equal to 100'
          write (6,*)'n,pctlak,pctwet,pcturb,pctgla,pctnatveg,pctcrop,sum= ', &
               n,pctlak(n),pctwet(n),pcturb(n),pctgla(n),&
               pctnatveg(n),pctcrop(n), suma
          call abort()
       end if

    end do

    ! Check that when pctnatveg+pctcrop identically zero, sum of special landunits is identically 100%

    if ( .not. outnc_double )then
       do n = 1,ns_o
          sum8  =         real(pctlak(n),r4)
          sum8  = sum8  + real(pctwet(n),r4)
          sum8  = sum8  + real(pcturb(n),r4)
          sum8  = sum8  + real(pctgla(n),r4)
          sum4a =         real(pctnatveg(n),r4)
          sum4a = sum4a + real(pctcrop(n),r4)
          if ( sum4a==0.0_r4 .and. sum8 < 100._r4-2._r4*epsilon(sum4a) )then
             write (6,*) subname, ' error: sum of pctlak, pctwet,', &
                  'pcturb, pctgla is < 100% when pctnatveg+pctcrop==0 sum = ', sum8
             write (6,*)'n,pctlak,pctwet,pcturb,pctgla,pctnatveg,pctcrop= ', &
                  n,pctlak(n),pctwet(n),pcturb(n),pctgla(n), pctnatveg(n),pctcrop(n)
             call abort()
          end if
       end do
    else
       do n = 1,ns_o
          sum8  =         pctlak(n)
          sum8  = sum8  + pctwet(n)
          sum8  = sum8  + pcturb(n)
          sum8  = sum8  + pctgla(n)
          sum8a =         pctnatveg(n)
          sum8a = sum8a + pctcrop(n)
          if ( sum8a==0._r8 .and. sum8 < (100._r8-4._r8*epsilon(sum8)) )then
             write (6,*) subname, ' error: sum of pctlak, pctwet,', &
                  'pcturb, pctgla is < 100% when pctnatveg+pctcrop==0 sum = ', sum8
             write (6,*) 'Total error, error/epsilon = ',100._r8-sum8, ((100._r8-sum8)/epsilon(sum8))
             write (6,*)'n,pctlak,pctwet,pcturb,pctgla,pctnatveg,pctcrop,epsilon= ', &
                  n,pctlak(n),pctwet(n),pcturb(n),pctgla(n),pctnatveg(n),pctcrop(n), epsilon(sum8)
             call abort()
          end if
       end do
    end if
    do n = 1,ns_o
       do k = natpft_lb,natpft_ub
          if ( pctnatpft(n,k)*pctnatveg(n) < 0.0_r8 )then
             write (6,*)'pctnatpft*pctnatveg < 0.0 = ', pctnatpft(n,k)*pctnatveg(n)
             call abort()
          end if
       end do
       do k = cft_lb,cft_ub
          if ( pctcft(n,k)*pctcrop(n) < 0.0_r8 )then
             write (6,*)'pctcft*pctcrop < 0.0 = ', pctcft(n,k)*pctcrop(n)
             call abort()
          end if
       end do
    end do

    ! Make sure that there is no vegetation outside the pft mask
    do n = 1,ns_o
       if (pftdata_mask(n) == 0 .and. (pctnatveg(n) > 0 .or. pctcrop(n) > 0)) then
          write (6,*)'vegetation found outside the pft mask at n=',n
          write (6,*)'pctnatveg,pctcrop=', pctnatveg, pctcrop
          call abort()
       end if
    end do

    ! Make sure that sums at the landunit level all add to 100%
    ! (Note that we don't check pctglcmec here, because it isn't computed at the point
    ! that this subroutine is called -- but the check of sum(pctglcmec) is done in
    ! mkglcmecMod)
    do n = 1,ns_o
       if (abs(sum(pctnatpft(n,:)) - 100._r8) > 1.e-12_r8) then
          write(6,*) 'sum(pctnatpft(n,:)) != 100: ', n, sum(pctnatpft(n,:))
          call abort()
       end if
       if (size(pctcft,2) > 0) then
          if (abs(sum(pctcft(n,:)) - 100._r8) > 1.e-12_r8) then
             write(6,*) 'sum(pctcft(n,:)) != 100: ', n, sum(pctcft(n,:))
             call abort()
          end if
       end if
       if (abs(sum(urbn_classes(n,:)) - 100._r8) > 1.e-12_r8) then
          write(6,*) 'sum(urbn_classes(n,:)) != 100: ', n, sum(urbn_classes(n,:))
          call abort()
       end if
    end do

    if ( nsmall_tot > 0 )then
       write (6,*)'number of small pft = ', nsmall_tot
    end if

  end subroutine normalizencheck_landuse


  ! ----------------------------------------------------------------------
  subroutine write_surface_dataset()

    use mkfileMod
    use mklaiMod
    use mkncdio            , only : check_ret

    implicit none

    include 'netcdf.inc'

    integer :: n
    integer :: ns_o
    integer  :: ncid                        ! netCDF id
    integer  :: omode                       ! netCDF output mode
    integer  :: varid                       ! netCDF variable id
    character(len=32) :: subname = 'write_surface_dataset'

    ns_o = ldomain%ns

    ! ----------------------------------------------------------------------
    ! Create surface dataset
    ! ----------------------------------------------------------------------

    ! Create netCDF surface dataset.  

    if (fsurdat == ' ') then
       write(6,*)' must specify fsurdat in namelist'
       stop
    end if

    call mkfile(ldomain, trim(fsurdat), dynlanduse = .false.)

    call domain_write(ldomain, fsurdat)

    call check_ret(nf_open(trim(fsurdat), nf_write, ncid), subname)
    call check_ret(nf_set_fill (ncid, nf_nofill, omode), subname)

    ! Write fields OTHER THAN lai, sai, heights, and urban parameters to netcdf surface dataset

    call check_ret(nf_inq_varid(ncid, 'natpft', varid), subname)
    call check_ret(nf_put_var_int(ncid, varid, (/(n,n=natpft_lb,natpft_ub)/)), subname)

    if (num_cft > 0) then
       call check_ret(nf_inq_varid(ncid, 'cft', varid), subname)
       call check_ret(nf_put_var_int(ncid, varid, (/(n,n=cft_lb,cft_ub)/)), subname)
    end if

    call check_ret(nf_inq_varid(ncid, 'PFTDATA_MASK', varid), subname)
    call check_ret(nf_put_var_int(ncid, varid, pftdata_mask), subname)

    call check_ret(nf_inq_varid(ncid, 'LANDFRAC_PFT', varid), subname)
    call check_ret(nf_put_var_double(ncid, varid, landfrac_pft), subname)

    call check_ret(nf_inq_varid(ncid, 'mxsoil_color', varid), subname)
    call check_ret(nf_put_var_int(ncid, varid, nsoicol), subname)

    call check_ret(nf_inq_varid(ncid, 'SOIL_COLOR', varid), subname)
    call check_ret(nf_put_var_int(ncid, varid, soicol), subname)

    call check_ret(nf_inq_varid(ncid, 'mxsoil_order', varid), subname)
    call check_ret(nf_put_var_int(ncid, varid, nsoiord), subname)

    call check_ret(nf_inq_varid(ncid, 'SOIL_ORDER', varid), subname)
    call check_ret(nf_put_var_int(ncid, varid, soiord), subname) 

    call check_ret(nf_inq_varid(ncid, 'PCT_SAND', varid), subname)
    call check_ret(nf_put_var_double(ncid, varid, pctsand), subname)

    call check_ret(nf_inq_varid(ncid, 'PCT_CLAY', varid), subname)
    call check_ret(nf_put_var_double(ncid, varid, pctclay), subname)

    call check_ret(nf_inq_varid(ncid, 'PCT_WETLAND', varid), subname)
    call check_ret(nf_put_var_double(ncid, varid, pctwet), subname)

    call check_ret(nf_inq_varid(ncid, 'PCT_LAKE', varid), subname)
    call check_ret(nf_put_var_double(ncid, varid, pctlak), subname)

    call check_ret(nf_inq_varid(ncid, 'PCT_GLACIER', varid), subname)
    call check_ret(nf_put_var_double(ncid, varid, pctgla), subname)

    if ( nglcec > 0 )then
       write(6,*)'error: add new code to support nglec > 0'
       call abort()
    end if

    call check_ret(nf_inq_varid(ncid, 'TOPO', varid), subname)
    call check_ret(nf_put_var_double(ncid, varid, topo), subname)

    call check_ret(nf_inq_varid(ncid, 'PCT_URBAN', varid), subname)
    call check_ret(nf_put_var_double(ncid, varid, urbn_classes_g), subname)

    call check_ret(nf_inq_varid(ncid, 'PCT_NATVEG', varid), subname)
    call check_ret(nf_put_var_double(ncid, varid, pctnatveg), subname)

    call check_ret(nf_inq_varid(ncid, 'PCT_CROP', varid), subname)
    call check_ret(nf_put_var_double(ncid, varid, pctcrop), subname)

    call check_ret(nf_inq_varid(ncid, 'PCT_NAT_PFT', varid), subname)
    call check_ret(nf_put_var_double(ncid, varid, pctnatpft), subname)

    if (num_cft > 0) then
       call check_ret(nf_inq_varid(ncid, 'PCT_CFT', varid), subname)
       call check_ret(nf_put_var_double(ncid, varid, pctcft), subname)
    end if

    call check_ret(nf_inq_varid(ncid, 'FMAX', varid), subname)
    call check_ret(nf_put_var_double(ncid, varid, fmax), subname)

    call check_ret(nf_inq_varid(ncid, 'gdp', varid), subname)
    call check_ret(nf_put_var_double(ncid, varid, gdp), subname)

    call check_ret(nf_inq_varid(ncid, 'peatf', varid), subname)
    call check_ret(nf_put_var_double(ncid, varid, fpeat), subname)

    call check_ret(nf_inq_varid(ncid, 'abm', varid), subname)
    call check_ret(nf_put_var_int(ncid, varid, agfirepkmon), subname)

    call check_ret(nf_inq_varid(ncid, 'SLOPE', varid), subname)
    call check_ret(nf_put_var_double(ncid, varid, slope), subname)

    call check_ret(nf_inq_varid(ncid, 'STD_ELEV', varid), subname)
    call check_ret(nf_put_var_double(ncid, varid, topo_stddev), subname)

    call check_ret(nf_inq_varid(ncid, 'binfl', varid), subname)
    call check_ret(nf_put_var_double(ncid, varid, vic_binfl), subname)

    call check_ret(nf_inq_varid(ncid, 'Ws', varid), subname)
    call check_ret(nf_put_var_double(ncid, varid, vic_ws), subname)

    call check_ret(nf_inq_varid(ncid, 'Dsmax', varid), subname)
    call check_ret(nf_put_var_double(ncid, varid, vic_dsmax), subname)

    call check_ret(nf_inq_varid(ncid, 'Ds', varid), subname)
    call check_ret(nf_put_var_double(ncid, varid, vic_ds), subname)

    call check_ret(nf_inq_varid(ncid, 'LAKEDEPTH', varid), subname)
    call check_ret(nf_put_var_double(ncid, varid, lakedepth), subname)

    call check_ret(nf_inq_varid(ncid, 'F0', varid), subname)
    call check_ret(nf_put_var_double(ncid, varid, f0), subname)

    call check_ret(nf_inq_varid(ncid, 'P3', varid), subname)
    call check_ret(nf_put_var_double(ncid, varid, p3), subname)

    call check_ret(nf_inq_varid(ncid, 'ZWT0', varid), subname)
    call check_ret(nf_put_var_double(ncid, varid, zwt0), subname)

    call check_ret(nf_inq_varid(ncid, 'EF1_BTR', varid), subname)
    call check_ret(nf_put_var_double(ncid, varid, ef1_btr), subname)

    call check_ret(nf_inq_varid(ncid, 'EF1_FET', varid), subname)
    call check_ret(nf_put_var_double(ncid, varid, ef1_fet), subname)

    call check_ret(nf_inq_varid(ncid, 'EF1_FDT', varid), subname)
    call check_ret(nf_put_var_double(ncid, varid, ef1_fdt), subname)

    call check_ret(nf_inq_varid(ncid, 'EF1_SHR', varid), subname)
    call check_ret(nf_put_var_double(ncid, varid, ef1_shr), subname)

    call check_ret(nf_inq_varid(ncid, 'EF1_GRS', varid), subname)
    call check_ret(nf_put_var_double(ncid, varid, ef1_grs), subname)

    call check_ret(nf_inq_varid(ncid, 'EF1_CRP', varid), subname)
    call check_ret(nf_put_var_double(ncid, varid, ef1_crp), subname)

    call check_ret(nf_inq_varid(ncid, 'ORGANIC', varid), subname)
    call check_ret(nf_put_var_double(ncid, varid, organic), subname)

    call check_ret(nf_inq_varid(ncid, 'URBAN_REGION_ID', varid), subname)
    call check_ret(nf_put_var_int(ncid, varid, urban_region), subname)

    call check_ret(nf_inq_varid(ncid, 'APATITE_P', varid), subname)
    call check_ret(nf_put_var_double(ncid, varid, apatiteP), subname)

    call check_ret(nf_inq_varid(ncid, 'LABILE_P', varid), subname)
    call check_ret(nf_put_var_double(ncid, varid, labileP), subname)

    call check_ret(nf_inq_varid(ncid, 'OCCLUDED_P', varid), subname)
    call check_ret(nf_put_var_double(ncid, varid, occludedP), subname)

    call check_ret(nf_inq_varid(ncid, 'SECONDARY_P', varid), subname)
    call check_ret(nf_put_var_double(ncid, varid, secondaryP), subname)

    call check_ret(nf_inq_varid(ncid, 'PCT_GRVL', varid), subname)
    call check_ret(nf_put_var_double(ncid, varid, grvl), subname) 

    call check_ret(nf_inq_varid(ncid, 'SLP_P10', varid), subname)
    call check_ret(nf_put_var_double(ncid, varid, slp10), subname)

    call check_ret(nf_inq_varid(ncid, 'parEro_c1', varid), subname)
    call check_ret(nf_put_var_double(ncid, varid, ero_c1), subname)

    call check_ret(nf_inq_varid(ncid, 'parEro_c2', varid), subname)
    call check_ret(nf_put_var_double(ncid, varid, ero_c2), subname)

    call check_ret(nf_inq_varid(ncid, 'parEro_c3', varid), subname)
    call check_ret(nf_put_var_double(ncid, varid, ero_c3), subname)

    call check_ret(nf_inq_varid(ncid, 'Tillage', varid), subname)
    call check_ret(nf_put_var_double(ncid, varid, tillage), subname)

    call check_ret(nf_inq_varid(ncid, 'Litho', varid), subname)
    call check_ret(nf_put_var_double(ncid, varid, litho), subname)

    ! Synchronize the disk copy of a netCDF dataset with in-memory buffers

    call check_ret(nf_sync(ncid), subname)

    ! ----------------------------------------------------------------------
    ! Make Urban Parameters from raw input data and write to surface dataset 
    ! Write to netcdf file is done inside mkurbanpar routine
    ! ----------------------------------------------------------------------

    write(6,*)'calling mkurbanpar'
    call mkurbanpar(datfname=mksrf_furban, ncido=ncid, region_o=urban_region, &
         urbn_classes_gcell_o=urbn_classes_g)

    ! ----------------------------------------------------------------------
    ! Make LAI and SAI from 1/2 degree data and write to surface dataset 
    ! Write to netcdf file is done inside mklai routine
    ! ----------------------------------------------------------------------

    write(6,*)'calling mklai'
    call mklai(ldomain, mapfname=map_flai, datfname=mksrf_flai, &
         ndiag=ndiag, ncido=ncid )

    ! Close surface dataset

    call check_ret(nf_close(ncid), subname)

    write (6,'(72a1)') ("-",n=1,60)
    write (6,*)' land model surface data set successfully created for ', &
         'grid of size ',ns_o

  end subroutine write_surface_dataset

  ! ----------------------------------------------------------------------
  subroutine write_surface_dataset_pio()

    use mkfileMod
    use mklaiMod
    use pio
    use piofileutils

    implicit none

    type(iosystem_desc_t) :: pioIoSystem
    type(file_desc_t)     :: ncid
    type(io_desc_t)       :: iodesc
    type(io_desc_t)       :: iodesc_scalar, iodesc_rad, iodesc_urb
    integer               :: i, j, k, count
    integer               :: ier
    integer, pointer      :: compdof(:)
    integer               :: dim2d(2), dim3d(3)
    integer               :: ncid_dummy

    if (fsurdat == ' ') then
       write(6,*)' must specify fsurdat in namelist'
       stop
    end if

    !call mkfile_pio(ldomain, trim(fsurdat), dynlanduse = .false.)
    call mkfile_pio(ldomain_pio, trim(fsurdat) // '.pio', .false.)
    call domain_write_pio(ldomain_pio, trim(fsurdat) // '.pio' )

    ! Open the file for writing

    call OpenFilePIO(trim(fsurdat) // '.pio' , pioIoSystem, ncid, PIO_WRITE)

    ! TOWRITE:
    !   natpft
    !   cft
    !   mxsoil_color
    !   mxsoil_order

    ! Write PIO_INT data with dim of (gridcell) or (lsmlat, lsmlon)

    call PIO_initdecomp(pioIoSystem, PIO_INT, ldomain_pio%dim_glb, ldomain_pio%compdof, iodesc)

    call write_integer_1d(ncid, iodesc, 'PFTDATA_MASK', pftdata_mask )
    call write_integer_1d(ncid, iodesc, 'SOIL_COLOR'  , soicol       )
    call write_integer_1d(ncid, iodesc, 'SOIL_ORDER'  , soiord       )
    call write_integer_1d(ncid, iodesc, 'abm'         , agfirepkmon  )
    call write_integer_1d(ncid, iodesc, 'URBAN_REGION_ID', urban_region  )

    call PIO_freedecomp(pioIoSystem, iodesc)

    ! Write PIO_DOUBLE data with dim of (gridcell) or (lsmlat, lsmlon)

    call PIO_initdecomp(pioIoSystem, PIO_DOUBLE, ldomain_pio%dim_glb, ldomain_pio%compdof, iodesc)

    call write_double_1d(ncid, iodesc, 'LANDFRAC_PFT', landfrac_pft)
    call write_double_1d(ncid, iodesc, 'PCT_WETLAND' , pctwet)
    call write_double_1d(ncid, iodesc, 'PCT_LAKE'    , pctlak)
    call write_double_1d(ncid, iodesc, 'PCT_GLACIER' , pctgla)
    call write_double_1d(ncid, iodesc, 'PCT_NATVEG'  , pctnatveg)
    call write_double_1d(ncid, iodesc, 'PCT_CROP'    , pctcrop)
    call write_double_1d(ncid, iodesc, 'parEro_c1'   , ero_c1)
    call write_double_1d(ncid, iodesc, 'parEro_c2'   , ero_c2)
    call write_double_1d(ncid, iodesc, 'parEro_c3'   , ero_c3)
    call write_double_1d(ncid, iodesc, 'Tillage'     , tillage)
    call write_double_1d(ncid, iodesc, 'Litho'       , litho)

    if (num_cft > 0) then
       write(6,*)'error: add code to support output of CFT'
       call abort()
    end if

    call write_double_1d(ncid, iodesc, 'FMAX', fmax)
    call write_double_1d(ncid, iodesc, 'gdp', gdp)
    call write_double_1d(ncid, iodesc, 'peatf', fpeat)
    call write_double_1d(ncid, iodesc, 'SLOPE', slope)
    call write_double_1d(ncid, iodesc, 'STD_ELEV', topo_stddev)
    call write_double_1d(ncid, iodesc, 'binfl', vic_binfl)
    call write_double_1d(ncid, iodesc, 'Ws'    , vic_ws)
    call write_double_1d(ncid, iodesc, 'Dsmax'    , vic_dsmax)
    call write_double_1d(ncid, iodesc, 'Ds'    , vic_ds)
    call write_double_1d(ncid, iodesc, 'LAKEDEPTH'    , lakedepth)
    call write_double_1d(ncid, iodesc, 'F0'    , f0)
    call write_double_1d(ncid, iodesc, 'P3'    , p3)
    call write_double_1d(ncid, iodesc, 'ZWT0'    , zwt0)
    call write_double_1d(ncid, iodesc, 'EF1_BTR'    , ef1_btr)
    call write_double_1d(ncid, iodesc, 'EF1_FET'    , ef1_fet)
    call write_double_1d(ncid, iodesc, 'EF1_FDT'    , ef1_fdt)
    call write_double_1d(ncid, iodesc, 'EF1_SHR'    , ef1_shr)
    call write_double_1d(ncid, iodesc, 'EF1_GRS'    , ef1_grs)
    call write_double_1d(ncid, iodesc, 'EF1_CRP'    , ef1_crp)
    call write_double_1d(ncid, iodesc, 'TOPO'       , topo)
    call write_double_1d(ncid, iodesc, 'APATITE_P'  , apatiteP)
    call write_double_1d(ncid, iodesc, 'LABILE_P'   , labileP)
    call write_double_1d(ncid, iodesc, 'OCCLUDED_P' , occludedP)
    call write_double_1d(ncid, iodesc, 'SECONDARY_P', secondaryP)

    call PIO_freedecomp(pioIoSystem, iodesc)

    ! Write PIO_DOUBLE data with dim of (nlevsoi, gridcell) or (nlevsoi, lsmlat, lsmlon)

    allocate(compdof(nlevsoi * ldomain_pio%ns_loc))
    count = 0
    do j = 1, nlevsoi
       do i = 1, ldomain_pio%ns_loc
          count = count + 1
          compdof(count) = i + (j-1)*ldomain_pio%ns_glb + (ldomain_pio%begs - 1)
       end do
    end do
    dim2d(1) = ldomain_pio%ns_glb; dim2d(2) = nlevsoi;
    call PIO_initdecomp(pioIoSystem, PIO_DOUBLE, dim2d, compdof, iodesc)

    call write_double_2d(ncid, iodesc, 'PCT_SAND', pctsand)
    call write_double_2d(ncid, iodesc, 'PCT_CLAY', pctclay)
    call write_double_2d(ncid, iodesc, 'ORGANIC' , organic)
    call write_double_2d(ncid, iodesc, 'PCT_GRVL', grvl)

    call PIO_freedecomp(pioIoSystem, iodesc)
    deallocate(compdof)

    ! Write PIO_DOUBLE data with dim of (nlevslp, gridcell) or (nlevslp, lsmlat, lsmlon)

    allocate(compdof((natpft_ub - natpft_lb + 1) * ldomain_pio%ns_loc))

    count = 0
    do j = natpft_lb, natpft_ub
      do i = 1, ldomain_pio%ns_loc
         count = count + 1
         compdof(count) = i + (j-1)*ldomain_pio%ns_glb + (ldomain_pio%begs - 1)
      end do
   end do

   dim2d(1) = ldomain_pio%ns_glb; dim2d(2) = natpft_ub - natpft_lb + 1;
   call PIO_initdecomp(pioIoSystem, PIO_DOUBLE, dim2d, compdof, iodesc)

   call write_double_2d(ncid, iodesc, 'PCT_NAT_PFT', pctnatpft)

   call PIO_freedecomp(pioIoSystem, iodesc)
   deallocate(compdof)

    ! Write PIO_DOUBLE data with dim of (nlevslp, gridcell) or (nlevslp, lsmlat, lsmlon)

    allocate(compdof(nlevslp * ldomain_pio%ns_loc))
    count = 0
    do j = 1, nlevslp
       do i = 1, ldomain_pio%ns_loc
          count = count + 1
          compdof(count) = i + (j-1)*ldomain_pio%ns_glb + (ldomain_pio%begs - 1)
       end do
    end do
    dim2d(1) = ldomain_pio%ns_glb; dim2d(2) = nlevslp;
    call PIO_initdecomp(pioIoSystem, PIO_DOUBLE, dim2d, compdof, iodesc)

    call write_double_2d(ncid, iodesc, 'SLP_P10', slp10)

    call PIO_freedecomp(pioIoSystem, iodesc)
    deallocate(compdof)

    !
    ! Write urban parameter
    !

    ! create IO-descripto of urban parameters with no extra dimensions
    allocate(compdof(numurbl * ldomain_pio%ns_loc))
    count = 0
    do j = 1, numurbl
       do i = 1, ldomain_pio%ns_loc
          count = count + 1
          compdof(count) = i + (j-1)*ldomain_pio%ns_glb + (ldomain_pio%begs - 1)
       end do
    end do
    dim2d(1) = ldomain_pio%ns_glb; dim2d(2) = numurbl
    call PIO_initdecomp(pioIoSystem, PIO_DOUBLE, dim2d, compdof, iodesc_scalar)

    ! create IO-descripto of urban parameters dimensioned by numrad and numurbl
    allocate(compdof(numrad * numurbl * ldomain_pio%ns_loc))
    count = 0
    do k = 1, numrad
       do j = 1, numurbl
          do i = 1, ldomain_pio%ns_loc
             count = count + 1
             compdof(count) = i + (j-1)*ldomain_pio%ns_glb + (k-1)*ldomain_pio%ns_glb * numurbl + (ldomain_pio%begs - 1)
          end do
       end do
    end do
    dim3d(1) = ldomain_pio%ns_glb; dim3d(2) = numurbl; dim3d(3) = numrad
    call PIO_initdecomp(pioIoSystem, PIO_DOUBLE, dim3d, compdof, iodesc_rad)

    ! create IO-descripto of urban parameters dimensioned by nlevurb
    allocate(compdof(nlevurb * numurbl * ldomain_pio%ns_loc))
    count = 0
    do k = 1, nlevurb
       do j = 1, numurbl
          do i = 1, ldomain_pio%ns_loc
             count = count + 1
             compdof(count) = i + (j-1)*ldomain_pio%ns_glb + (k-1)*ldomain_pio%ns_glb * numurbl + (ldomain_pio%begs - 1)
          end do
       end do
    end do
    dim3d(1) = ldomain_pio%ns_glb; dim3d(2) = numurbl; dim3d(3) = nlevurb
    call PIO_initdecomp(pioIoSystem, PIO_DOUBLE, dim3d, compdof, iodesc_urb)

    call mkurbanpar(datfname=mksrf_furban, ncid_pio=ncid, region_o=urban_region, &
         urbn_classes_gcell_o=urbn_classes_g, iodesc_scalar=iodesc_scalar, &
         iodesc_rad=iodesc_rad, iodesc_urb=iodesc_urb)

    !
    ! Write LAI
    !

    call mklai_pio(ldomain_pio, map_flai, mksrf_flai, ndiag, pioIoSystem, ncid)

    ! Clean up

    call PIO_closefile(ncid)
    call PIO_finalize(pioIoSystem, ier)

  end subroutine write_surface_dataset_pio

  !------------------------------------------------------------------------------
  subroutine domain_write_pio(domain_pio, fname)
    !
    use pio
    use piofileutils
    !
    implicit none
    !
    type(domain_pio_type)        :: domain_pio
    character(len=*) ,intent(in) :: fname
    !
    type(iosystem_desc_t) :: pioIoSystem
    type(file_desc_t)     :: ncid
    type(var_desc_t)      :: varid
    type(io_desc_t)       :: iodescNCells
    integer               :: vartype
    integer               :: ier
    character(len=32)     :: subname = 'domain_write_pio'

    ! Open the file for writing

    call OpenFilePIO(trim(fname), pioIoSystem, ncid, PIO_WRITE)

    ! Write variables

    call check_ret(PIO_inq_varid(ncid, 'AREA', varid), subname)
    call check_ret(PIO_inq_vartype(ncid, varid, vartype), subname)

    call PIO_initdecomp(pioIoSystem, vartype, domain_pio%dim_glb, domain_pio%compdof, iodescNCells)
    call PIO_write_darray(ncid, varid, iodescNCells, domain_pio%area1d, ier)
    call PIO_syncfile(ncid)

    call check_ret(PIO_inq_varid(ncid, 'LONGXY', varid), subname)
    call PIO_write_darray(ncid, varid, iodescNCells, domain_pio%lonc1d, ier)
    call PIO_syncfile(ncid)

    call check_ret(PIO_inq_varid(ncid, 'LATIXY', varid), subname)
    call PIO_write_darray(ncid, varid, iodescNCells, domain_pio%latc1d, ier)
    call PIO_syncfile(ncid)

    ! Clean up

    call PIO_freedecomp(pioIoSystem, iodescNCells)
    call PIO_closefile(ncid)
    call PIO_finalize(pioIoSystem, ier)

  end subroutine domain_write_pio

  end program mksurfdat_petsc


