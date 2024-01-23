module mkdataPIOMod

  use shr_kind_mod, only : r8 => shr_kind_r8, r4=>shr_kind_r4
  use shr_sys_mod , only : shr_sys_flush
  use mkchecksMod, only : min_bad, max_bad

  implicit none

  public mkdata_double_2d_pio
  public mkdata_double_3d_pio

contains

  !-----------------------------------------------------------------------
  subroutine mkdata_double_2d_pio(ldomain_pio, mapfname, datfname, varname, data_descrip, &
       ndiag, zero_out, nodata_value, data_o, threshold_o, min_valid_value, max_valid_value)
    !
    ! !DESCRIPTION:
    !
    ! !USES:
    use mkdomainPIOMod, only : domain_pio_type, domain_clean_pio, domain_read_pio
    use mkgridmapMod
    use mkvarpar
    use mkvarctl
    use mkncdio
    use pio
    use piofileutils
    use utils
    !
    ! !ARGUMENTS:

    implicit none
    type(domain_pio_type) , intent(in)           :: ldomain_pio
    character(len=*)      , intent(in)           :: mapfname          ! input mapping file name
    character(len=*)      , intent(in)           :: datfname          ! input data file name
    character(len=*)      , intent(in)           :: varname           ! input data file name
    character(len=*)      , intent(in)           :: data_descrip      ! data description
    integer               , intent(in)           :: ndiag             ! variable name
    real(r8)              , intent(in)           :: nodata_value      !
    logical               , intent(in)           :: zero_out          ! if should zero glacier out
    real(r8)              , intent(out)          :: data_o(:)         ! output grid: %lake
    real(r8)              , intent(in), optional :: threshold_o       !
    real(r8)              , intent(in), optional :: min_valid_value   !
    real(r8)              , intent(in), optional :: max_valid_value
                                                                      !
    type(gridmap_type)                           :: tgridmap
    type(domain_pio_type)                        :: tdomain_pio       ! local domain
    integer                                      :: no
    integer                                      :: ns_loc_i,ns_loc_o !  indices
    integer                                      :: ierr              ! error status

    type(file_desc_t)                            :: ncid
    type(iosystem_desc_t)                        :: pioIoSystem
    real(r8) , pointer                           :: data2d_i(:,:)
    real(r8) , pointer                           :: data1d_i(:)
    integer                                      :: dim_idx(2,2)
    logical                                      :: threshold_specified
    logical                                      :: min_valid_specified
    logical                                      :: max_valid_specified

    !-----------------------------------------------------------------------

    write (6,*) 'Attempting to make ' // trim(data_descrip) // ' .....'
    call shr_sys_flush(6)

    threshold_specified = .false.
    min_valid_specified = .false.
    max_valid_specified = .false.

    if (present(threshold_o    )) threshold_specified = .true.
    if (present(min_valid_value)) min_valid_specified = .true.
    if (present(max_valid_value)) max_valid_specified = .true.

    ! -----------------------------------------------------------------
    ! Read input file
    ! -----------------------------------------------------------------

    if ( .not. zero_out )then

       write(6,*)'Open file: ', trim(datfname)

       ! Obtain input grid info, read local fields
       
       call domain_read_pio(tdomain_pio, datfname)

       ! Open the netcdf file
       call OpenFilePIO(datfname, pioIoSystem, ncid)

       ! Read the variable
       call read_float_or_double_2d(tdomain_pio, pioIoSystem, ncid, varname, dim_idx, data2d_i)

       call PIO_closefile(ncid)
       call PIO_finalize(pioIoSystem, ierr)

       ! Read the map
       call gridmap_mapread(tgridmap, mapfname )

       ! Convert 2D vector to 1D vector
       ns_loc_i = (dim_idx(1,2) - dim_idx(1,1) + 1) * (dim_idx(2,2) - dim_idx(2,1) + 1)
       allocate(data1d_i(ns_loc_i))
       call convert_2d_to_1d_array(dim_idx(1,1), dim_idx(1,2), dim_idx(2,1), dim_idx(2,2), data2d_i, data1d_i)

       ! Determine data_o on output grid
       call gridmap_areaave(tgridmap, data1d_i(:), data_o, nodata=nodata_value)

       if (min_valid_specified) then
          if (min_bad(data_o, min_valid_value, data_descrip)) then
             write (6,*) 'Value on output grid is below the valid value of ', min_valid_value
             call abort()
          end if
       end if

       if (max_valid_specified) then
          if (max_bad(data_o, max_valid_value, data_descrip)) then
             write (6,*) 'Value on output grid is above the valid value of ', max_valid_value
             call abort()
          end if
       end if

       if (threshold_specified) then
          ns_loc_o = ldomain_pio%ns_loc
          do no = 1, ns_loc_o
             if (data_o(no) < threshold_o) data_o(no) = 0._r8
          end do
       end if

    else
       data_o(:) = 0._r8
    end if

    ! Deallocate dynamic memory

    if ( .not. zero_out )then
       call domain_clean_pio(tdomain_pio)

       call gridmap_clean(tgridmap)
       deallocate (data2d_i)
       deallocate (data1d_i)
    end if

    write (6,*) 'Successfully made ' // trim(data_descrip)
    write (6,*)
    call shr_sys_flush(6)

  end subroutine mkdata_double_2d_pio

  !-----------------------------------------------------------------------
  subroutine mkdata_double_3d_pio(ldomain_pio, mapfname, datfname, varname, data_descrip, &
       ndiag, zero_out, nodata_value, start_id_for_dim3, data_o, threshold_o, &
       min_valid_value, max_valid_value)
    !
    ! !DESCRIPTION:
    !
    ! !USES:
    use mkdomainPIOMod, only : domain_pio_type, domain_clean_pio, domain_read_pio
    use mkgridmapMod
    use mkvarpar
    use mkvarctl
    use mkncdio
    use pio
    use piofileutils
    use utils
    !
    ! !ARGUMENTS:

    implicit none
    type(domain_pio_type) , intent(in)           :: ldomain_pio
    character(len=*)      , intent(in)           :: mapfname          ! input mapping file name
    character(len=*)      , intent(in)           :: datfname          ! input data file name
    character(len=*)      , intent(in)           :: varname           ! input data file name
    character(len=*)      , intent(in)           :: data_descrip      ! data description
    integer               , intent(in)           :: ndiag             ! variable name
    real(r8)              , intent(in)           :: nodata_value      !
    integer               , intent(in)           :: start_id_for_dim3
    logical               , intent(in)           :: zero_out          ! if should zero glacier out
    real(r8)              , intent(out)          :: data_o(:,:)       ! output grid: %lake
    real(r8)              , intent(in), optional :: threshold_o       !
    real(r8)              , intent(in), optional :: min_valid_value   !
    real(r8)              , intent(in), optional :: max_valid_value
                                                                      !
    type(gridmap_type)                           :: tgridmap
    type(domain_pio_type)                        :: tdomain_pio       ! local domain
    integer                                      :: no,m
    integer                                      :: ns_loc_i,ns_loc_o !  indices
    integer                                      :: ierr              ! error status

    type(file_desc_t)                            :: ncid
    type(iosystem_desc_t)                        :: pioIoSystem
    real(r8) , pointer                           :: data3d_i(:,:,:)
    real(r8) , pointer                           :: data1d_i(:)
    integer                                      :: dim_idx(3,2)
    logical                                      :: threshold_specified
    logical                                      :: min_valid_specified
    logical                                      :: max_valid_specified

    !-----------------------------------------------------------------------

    write (6,*) 'Attempting to make ' // trim(data_descrip) // ' .....'
    call shr_sys_flush(6)

    threshold_specified = .false.
    min_valid_specified = .false.
    max_valid_specified = .false.

    if (present(threshold_o    )) threshold_specified = .true.
    if (present(min_valid_value)) min_valid_specified = .true.
    if (present(max_valid_value)) max_valid_specified = .true.

    ! -----------------------------------------------------------------
    ! Read input file
    ! -----------------------------------------------------------------

    if ( .not. zero_out )then

       write(6,*)'Open file: ', trim(datfname)

       ! Obtain input grid info, read local fields
       
       write(*,*)'  datfname: ',trim(datfname)
       call domain_read_pio(tdomain_pio, datfname)

       ! Open the netcdf file
       write(*,*)'  call OpenFilePIO: '
       call OpenFilePIO(datfname, pioIoSystem, ncid)

       ! Read the variable
       write(*,*)'  call read_float_or_double_3d(): ',trim(varname)
       call read_float_or_double_3d(tdomain_pio, pioIoSystem, ncid, varname, start_id_for_dim3, dim_idx, data3d_i)

       call PIO_closefile(ncid)
       call PIO_finalize(pioIoSystem, ierr)

       ! Read the map
       call gridmap_mapread(tgridmap, mapfname )

       ! Convert 2D vector to 1D vector
       ns_loc_i = (dim_idx(1,2) - dim_idx(1,1) + 1) * (dim_idx(2,2) - dim_idx(2,1) + 1)
       allocate(data1d_i(ns_loc_i))

       do m = start_id_for_dim3, start_id_for_dim3 + (dim_idx(3,2) - dim_idx(3,1))

          ! Convert data from 2D to 1D
          call convert_2d_to_1d_array(dim_idx(1,1), dim_idx(1,2), dim_idx(2,1), dim_idx(2,2), data3d_i(:,:,m), data1d_i)

          ! Determine data_o on output grid
          call gridmap_areaave(tgridmap, data1d_i(:), data_o(:,m), nodata=nodata_value)

          if (min_valid_specified) then
             if (min_bad(data_o(:,m), min_valid_value, data_descrip)) then
                write (6,*) 'Value on output grid is below the valid value of ', min_valid_value
                call abort()
             end if
          end if

          if (max_valid_specified) then
             if (max_bad(data_o(:,m), max_valid_value, data_descrip)) then
                write (6,*) 'Value on output grid is above the valid value of ', max_valid_value
                call abort()
             end if
          end if

          if (threshold_specified) then
             ns_loc_o = ldomain_pio%ns_loc
             do no = 1, ns_loc_o
                if (data_o(no,m) < threshold_o) data_o(no,m) = 0._r8
             end do
          end if
       end do

    else
       data_o(:,:) = 0._r8
    end if

    ! Deallocate dynamic memory

    if ( .not. zero_out )then
       call domain_clean_pio(tdomain_pio)

       call gridmap_clean(tgridmap)
       deallocate (data3d_i)
       deallocate (data1d_i)
    end if

    write (6,*) 'Successfully made ' // trim(data_descrip)
    write (6,*)
    call shr_sys_flush(6)

  end subroutine mkdata_double_3d_pio

end module mkdataPIOMod
