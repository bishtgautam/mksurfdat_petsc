module utils

  use shr_kind_mod, only : r8 => shr_kind_r8

  implicit none

  public convert_2d_to_1d_array

contains

  subroutine convert_2d_to_1d_array(dim_idx, data2d, data1d)
    !
    implicit none
    !
    integer          , intent(in)  :: dim_idx(2,2)
    real(r8), pointer, intent(in)  :: data2d(:,:)
    real(r8), pointer, intent(out) :: data1d(:)
    !
    integer                        :: i, j, count
    
    count = 0
    do j = dim_idx(2,1), dim_idx(2,2)
       do i = dim_idx(1,1), dim_idx(1,2)
          count = count + 1
          data1d(count) = data2d(i,j)
       end do
    end do
    
  end subroutine convert_2d_to_1d_array
  
end module utils
