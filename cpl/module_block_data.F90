!> @file
!> @brief Copies block data containing real*4, real*8, or integer into
!> ESMF_KIND_R8 arrays, with an optional scaling factor. Can also fill
!> ESMF_KIND_R8 arrays with a constant value.
!>
!> @author Raffaele Montuoro @date 7/1/21
module module_block_data


  use ESMF,              only: ESMF_KIND_R8, ESMF_SUCCESS, &
                               ESMF_RC_PTR_NOTALLOC, ESMF_RC_VAL_OUTOFRANGE
  use block_control_mod, only: block_control_type

  implicit none

  interface block_data_copy
    module procedure block_copy_1d_i4_to_2d_r8
    module procedure block_copy_1d_r8_to_2d_r8
    module procedure block_copy_spval_1d_r8_to_2d_r8
    module procedure block_copy_2d_r8_to_2d_r8
    module procedure block_copy_2d_r8_to_3d_r8
    module procedure block_copy_3d_r8_to_3d_r8
    module procedure block_copy_1dslice_r8_to_2d_r8
    module procedure block_copy_1dslice2_r8_to_2d_r8
    module procedure block_copy_3dslice_r8_to_3d_r8
    module procedure block_copy_1d_r4_to_2d_r8
    module procedure block_copy_spval_1d_r4_to_2d_r8
    module procedure block_copy_2d_r4_to_2d_r8
    module procedure block_copy_2d_r4_to_3d_r8
    module procedure block_copy_3d_r4_to_3d_r8
    module procedure block_copy_1dslice_r4_to_2d_r8
    module procedure block_copy_1dslice2_r4_to_2d_r8
    module procedure block_copy_3dslice_r4_to_3d_r8
  end interface block_data_copy

  interface block_data_fill
    module procedure block_fill_2d_r8
    module procedure block_fill_3d_r8
  end interface block_data_fill

  interface block_data_copy_or_fill
    module procedure block_copy_or_fill_1d_r8_to_2d_r8
    module procedure block_copy_or_fill_2d_r8_to_3d_r8
    module procedure block_copy_or_fill_1dslice_r8_to_2d_r8
    module procedure block_copy_or_fill_1dslice2_r8_to_2d_r8
    module procedure block_copy_or_fill_1d_r4_to_2d_r8
    module procedure block_copy_or_fill_2d_r4_to_3d_r8
    module procedure block_copy_or_fill_1dslice_r4_to_2d_r8
    module procedure block_copy_or_fill_1dslice2_r4_to_2d_r8
  end interface block_data_copy_or_fill

  interface block_data_combine_fractions
    module procedure block_combine_frac_1d_r8_to_2d_r8
    module procedure block_combine_frac_1d_r4_to_2d_r8
  end interface block_data_combine_fractions

  interface block_atmos_copy
    module procedure block_array_copy_2d_r8_to_2d_r8
    module procedure block_array_copy_3d_r8_to_3d_r8
    module procedure block_array_copy_3dslice_r8_to_3d_r8
    module procedure block_array_copy_2d_r4_to_2d_r8
    module procedure block_array_copy_3d_r4_to_3d_r8
    module procedure block_array_copy_3dslice_r4_to_3d_r8
  end interface block_atmos_copy

  private

  public :: block_atmos_copy

  public :: block_data_copy
  public :: block_data_fill
  public :: block_data_copy_or_fill
  public :: block_data_combine_fractions

contains

  !> Copies a 1D array of real*8 values to a 2D array of real*8 values.
  !>
  !> This subroutine copies values from a 1D source array to a 2D destination
  !> array, applying a scale factor and handling special values. The copy
  !> operation is performed for a specific block and block index.
  !>
  !> @param[out] destin_ptr Pointer to the destination 2D array.
  !> @param[in]  source_ptr Pointer to the source 1D array.
  !> @param[in]  block      The block structure containing information about
  !>                        the block dimensions.
  !> @param[in]  block_index The index of the block to copy.
  !> @param[in]  scale_factor The factor by which to scale the source values.
  !> @param[in]  special_value The special value to handle in the source array.
  !> @param[in]  offset     The offset to apply to the destination indices.
  !> @param[out] rc         Return code indicating success or failure of the
  !>                        operation.
  !>
  !> @author Raffaele Montuoro @date 7/1/21
  subroutine block_copy_spval_1d_r8_to_2d_r8(destin_ptr, source_ptr, block, block_index, scale_factor, special_value, offset, rc)

    ! -- arguments
    real(ESMF_KIND_R8),        pointer     :: destin_ptr(:,:)
    real(kind=8),              pointer     :: source_ptr(:)
    type(block_control_type),  intent(in)  :: block
    integer,                   intent(in)  :: block_index
    real(kind=8),              intent(in)  :: scale_factor
    real(kind=8),              intent(in)  :: special_value
    integer,                   intent(in)  :: offset
    integer,         optional, intent(out) :: rc

    ! -- local variables
    integer      :: localrc
    integer      :: i, ib, ix, im, j, jb

    ! -- begin
    localrc = ESMF_RC_PTR_NOTALLOC
    if (associated(destin_ptr) .and. associated(source_ptr)) then
!$omp parallel do private(ix,im,ib,jb,i,j)
       do ix = 1, block%blksz(block_index)
          im = offset + ix - 1
          ib = block%index(block_index)%ii(ix)
          jb = block%index(block_index)%jj(ix)
          i = ib - block%isc + 1
          j = jb - block%jsc + 1
          if (source_ptr(im) .ne. special_value) then
             destin_ptr(i,j) = scale_factor * source_ptr(im)
          else
             destin_ptr(i,j) = special_value
          end if
       enddo
       localrc = ESMF_SUCCESS
    end if

    if (present(rc)) rc = localrc

  end subroutine block_copy_spval_1d_r8_to_2d_r8

  !> Copies a 1D slice of real*8 data to a 2D real*8 array.
  !>
  !> This subroutine copies a 1D slice of real*8 data from the source
  !> array to a 2D real*8 destination array. The slice is specified by
  !> the block and block_index parameters. The data can be scaled and
  !> offset during the copy process.
  !>
  !> @param[out] destin_ptr Pointer to the destination 2D real*8 array.
  !> @param[in]  source_ptr Pointer to the source 1D real*8 array.
  !> @param[in]  slice      Integer specifying the slice to be copied.
  !> @param[in]  block      Integer specifying the block to be copied.
  !> @param[in]  block_index Integer specifying the index within the block.
  !> @param[in]  scale_factor Real*8 value to scale the data during the copy.
  !> @param[in]  offset      Real*8 value to offset the data during the copy.
  !> @param[out] rc          Integer return code (0 for success, non-zero for error).
  !>
  !> @author Raffaele Montuoro @date 7/1/21
  subroutine block_copy_1dslice_r8_to_2d_r8(destin_ptr, source_ptr, slice, block, block_index, scale_factor, offset, rc)

    ! -- arguments
    real(ESMF_KIND_R8),        pointer     :: destin_ptr(:,:)
    real(kind=8),              pointer     :: source_ptr(:,:)
    integer,                   intent(in)  :: slice
    type (block_control_type), intent(in)  :: block
    integer,                   intent(in)  :: block_index
    real(kind=8),    optional, intent(in)  :: scale_factor
    integer,                   intent(in)  :: offset
    integer,         optional, intent(out) :: rc

    ! -- local variables
    integer      :: localrc
    integer      :: i, ib, ix, im, j, jb
    real(kind=8) :: factor

    ! -- begin
    localrc = ESMF_RC_PTR_NOTALLOC
    if (associated(destin_ptr) .and. associated(source_ptr)) then
      localrc = ESMF_RC_VAL_OUTOFRANGE
      if (slice > 0 .and. slice <= size(source_ptr, dim=2)) then
        factor = 1._8
        if (present(scale_factor)) factor = scale_factor
!$omp parallel do private(ix,im,ib,jb,i,j)
        do ix = 1, block%blksz(block_index)
          im = offset + ix - 1
          ib = block%index(block_index)%ii(ix)
          jb = block%index(block_index)%jj(ix)
          i = ib - block%isc + 1
          j = jb - block%jsc + 1
          destin_ptr(i,j) = factor * source_ptr(im,slice)
        enddo
        localrc = ESMF_SUCCESS
      end if
    end if

    if (present(rc)) rc = localrc

  end subroutine block_copy_1dslice_r8_to_2d_r8

  !> Copies a 1D slice of real*8 data to a 2D real*8 array.
  !>
  !> This subroutine takes a 1D slice of real*8 data from the source array
  !> and copies it into a 2D real*8 destination array. The data can be
  !> scaled and offset during the copy process.
  !>
  !> @param[out] destin_ptr Pointer to the destination 2D real*8 array.
  !> @param[in]  source_ptr Pointer to the source 1D real*8 array.
  !> @param[in]  slice1     The first dimension index of the slice.
  !> @param[in]  slice2     The second dimension index of the slice.
  !> @param[in]  block      The block size for the copy operation.
  !> @param[in]  block_index The index of the block to be copied.
  !> @param[in]  scale_factor The factor by which to scale the source data.
  !> @param[in]  offset     The offset to be added to the scaled data.
  !> @param[out] rc         Return code indicating success or failure of the operation.
  !>
  !> @author Raffaele Montuoro @date 7/1/21
  subroutine block_copy_1dslice2_r8_to_2d_r8(destin_ptr, source_ptr, slice1, slice2, block, block_index, scale_factor, offset, rc)

    ! -- arguments
    real(ESMF_KIND_R8),        pointer     :: destin_ptr(:,:)
    real(kind=8),              pointer     :: source_ptr(:,:,:)
    integer,                   intent(in)  :: slice1
    integer,                   intent(in)  :: slice2
    type (block_control_type), intent(in)  :: block
    integer,                   intent(in)  :: block_index
    real(kind=8),    optional, intent(in)  :: scale_factor
    integer,                   intent(in)  :: offset
    integer,         optional, intent(out) :: rc

    ! -- local variables
    integer      :: localrc
    integer      :: i, ib, ix, im, j, jb
    real(kind=8) :: factor

    ! -- begin
    localrc = ESMF_RC_PTR_NOTALLOC
    if (associated(destin_ptr) .and. associated(source_ptr)) then
      localrc = ESMF_RC_VAL_OUTOFRANGE
      if (slice1 > 0 .and. slice1 <= size(source_ptr, dim=2) .and. slice2 > 0 .and. slice2 <= size(source_ptr, dim=3)) then
        factor = 1._8
        if (present(scale_factor)) factor = scale_factor
!$omp parallel do private(ix,im,ib,jb,i,j)
        do ix = 1, block%blksz(block_index)
          im = offset + ix - 1
          ib = block%index(block_index)%ii(ix)
          jb = block%index(block_index)%jj(ix)
          i = ib - block%isc + 1
          j = jb - block%jsc + 1
          destin_ptr(i,j) = factor * source_ptr(im,slice1,slice2)
        enddo
        localrc = ESMF_SUCCESS
      end if
    end if

    if (present(rc)) rc = localrc

  end subroutine block_copy_1dslice2_r8_to_2d_r8

  !> Copies a 2D real*8 block of data to a 3D real*8 block of data.
  !>
  !> This subroutine copies data from a 2D source array to a 3D destination
  !> array, applying a scale factor and an offset to each element during
  !> the copy process.
  !>
  !> @param[inout] destin_ptr Pointer to the destination 3D array.
  !> @param[in] source_ptr Pointer to the source 2D array.
  !> @param[in] block The block of data to be copied.
  !> @param[in] block_index The index of the block in the destination array.
  !> @param[in] scale_factor The factor by which to scale the source data.
  !> @param[in] offset The offset to be added to the scaled source data.
  !> @param[out] rc Return code indicating success or failure of the operation.
  !>
  !> @author Raffaele Montuoro @date 7/1/21
  subroutine block_copy_2d_r8_to_3d_r8(destin_ptr, source_ptr, block, block_index, scale_factor, offset, rc)

    ! -- arguments
    real(ESMF_KIND_R8),        pointer     :: destin_ptr(:,:,:)
    real(kind=8),              pointer     :: source_ptr(:,:)
    type (block_control_type), intent(in)  :: block
    integer,                   intent(in)  :: block_index
    real(kind=8),    optional, intent(in)  :: scale_factor
    integer,                   intent(in)  :: offset
    integer,         optional, intent(out) :: rc

    ! -- local variables
    integer      :: localrc
    integer      :: i, ib, ix, im, j, jb, k
    real(kind=8) :: factor

    ! -- begin
    localrc = ESMF_RC_PTR_NOTALLOC
    if (associated(destin_ptr) .and. associated(source_ptr)) then
      factor = 1._8
      if (present(scale_factor)) factor = scale_factor
      do k = 1, size(source_ptr, dim=2)
!$omp parallel do private(ix,im,ib,jb,i,j)
        do ix = 1, block%blksz(block_index)
          im = offset + ix - 1
          ib = block%index(block_index)%ii(ix)
          jb = block%index(block_index)%jj(ix)
          i = ib - block%isc + 1
          j = jb - block%jsc + 1
          destin_ptr(i,j,k) = factor * source_ptr(im,k)
        enddo
      enddo
      localrc = ESMF_SUCCESS
    end if

    if (present(rc)) rc = localrc

  end subroutine block_copy_2d_r8_to_3d_r8

  !> Copies a 3D block of real*8 data from source to destination.
  !>
  !> This subroutine copies a 3D block of real*8 (double precision) data
  !> from the source array to the destination array. The block to be copied
  !> is specified by the block and block_index parameters. The data can be
  !> scaled and offset during the copy operation.
  !>
  !> @param[out] destin_ptr Pointer to the destination array.
  !> @param[in] source_ptr Pointer to the source array.
  !> @param[in] block Specifies the block dimensions to be copied.
  !> @param[in] block_index Specifies the starting index of the block in the source array.
  !> @param[in] scale_factor Scaling factor to be applied to the source data.
  !> @param[in] offset Offset to be added to the scaled source data.
  !> @param[out] rc Return code indicating success or failure of the operation.
  !>
  !> @author Raffaele Montuoro @date 7/1/21
  subroutine block_copy_3d_r8_to_3d_r8(destin_ptr, source_ptr, block, block_index, scale_factor, offset, rc)

    ! -- arguments
    real(ESMF_KIND_R8),        pointer     :: destin_ptr(:,:,:)
    real(kind=8),              pointer     :: source_ptr(:,:,:)
    type (block_control_type), intent(in)  :: block
    integer,                   intent(in)  :: block_index
    real(kind=8),    optional, intent(in)  :: scale_factor
    integer,                   intent(in)  :: offset
    integer,         optional, intent(out) :: rc

    ! -- local variables
    integer      :: localrc
    integer      :: i, ib, ix, im, j, jb, k
    real(kind=8) :: factor

    ! -- begin
    localrc = ESMF_RC_PTR_NOTALLOC
    if (associated(destin_ptr) .and. associated(source_ptr)) then
      factor = 1._8
      if (present(scale_factor)) factor = scale_factor
      do k = 1, size(source_ptr, dim=3)
!$omp parallel do private(ix,im,ib,jb,i,j)
        do ix = 1, block%blksz(block_index)
          im = offset + ix - 1
          ib = block%index(block_index)%ii(ix)
          jb = block%index(block_index)%jj(ix)
          i = ib - block%isc + 1
          j = jb - block%jsc + 1
          destin_ptr(i,j,k) = factor * source_ptr(ib,jb,k)
        enddo
      enddo
      localrc = ESMF_SUCCESS
    end if

    if (present(rc)) rc = localrc

  end subroutine block_copy_3d_r8_to_3d_r8

  !> Copies a 3D real*8 array from source to destination with scaling and offset.
  !>
  !> This subroutine copies a 3D array of real*8 (double precision) values from the source array
  !> to the destination pointer, applying a scaling factor and an offset to each element.
  !>
  !> @param[out] destin_ptr Pointer to the destination 3D array of real*8 values.
  !> @param[in]  source_arr 3D array of real*8 values to be copied.
  !> @param[in]  block      Integer specifying the block size or index.
  !> @param[in]  block_index Integer specifying the index within the block.
  !> @param[in]  scale_factor Real*8 value to scale each element of the source array.
  !> @param[in]  offset      Real*8 value to add to each scaled element of the source array.
  !> @param[out] rc          Integer return code indicating success or failure of the operation.
  !>
  !> @author Raffaele Montuoro @date 7/1/21
  subroutine block_array_copy_3d_r8_to_3d_r8(destin_ptr, source_arr, block, block_index, scale_factor, offset, rc)

    ! -- arguments
    real(ESMF_KIND_R8),        pointer     :: destin_ptr(:,:,:)
    real(kind=8),              intent(in)  :: source_arr(:,:,:)
    type (block_control_type), intent(in)  :: block
    integer,                   intent(in)  :: block_index
    real(kind=8),    optional, intent(in)  :: scale_factor
    integer,                   intent(in)  :: offset
    integer,         optional, intent(out) :: rc

    ! -- local variables
    integer      :: localrc
    integer      :: i, ib, ix, im, j, jb, k
    real(kind=8) :: factor

    ! -- begin
    localrc = ESMF_RC_PTR_NOTALLOC
    if (associated(destin_ptr)) then
      factor = 1._8
      if (present(scale_factor)) factor = scale_factor
      do k = 1, size(source_arr, dim=3)
!$omp parallel do private(ix,im,ib,jb,i,j)
        do ix = 1, block%blksz(block_index)
          im = offset + ix - 1
          ib = block%index(block_index)%ii(ix)
          jb = block%index(block_index)%jj(ix)
          i = ib - block%isc + 1
          j = jb - block%jsc + 1
          destin_ptr(i,j,k) = factor * source_arr(i,j,k)
        enddo
      enddo
      localrc = ESMF_SUCCESS
    end if

    if (present(rc)) rc = localrc

  end subroutine block_array_copy_3d_r8_to_3d_r8

  !> copy: 3D slice to 3D ???
  !> 
  !> @param[in] destin_ptr ???
  !> @param[in] source_ptr ???
  !> @param[in] slice ???
  !> @param[in] block ???
  !> @param[in] block_index ???
  !> @param[in] scale_factor ???
  !> @param[in] offset ???
  !> @param[out] rc ???
  !>
  !> @author
  subroutine block_copy_3dslice_r8_to_3d_r8(destin_ptr, source_ptr, slice, block, block_index, scale_factor, offset, rc)

    ! -- arguments
    real(ESMF_KIND_R8),        pointer     :: destin_ptr(:,:,:)
    real(kind=8),              pointer     :: source_ptr(:,:,:,:)
    integer,                   intent(in)  :: slice
    type (block_control_type), intent(in)  :: block
    integer,                   intent(in)  :: block_index
    real(kind=8),    optional, intent(in)  :: scale_factor
    integer,                   intent(in)  :: offset
    integer,         optional, intent(out) :: rc

    ! -- local variables
    integer      :: localrc
    integer      :: i, ib, ix, im, j, jb, k
    real(kind=8) :: factor

    ! -- begin
    localrc = ESMF_RC_PTR_NOTALLOC
    if (associated(destin_ptr) .and. associated(source_ptr)) then
      localrc = ESMF_RC_VAL_OUTOFRANGE
      if (slice > 0 .and. slice <= size(source_ptr, dim=4)) then
        factor = 1._8
        if (present(scale_factor)) factor = scale_factor
        do k = 1, size(source_ptr, dim=3)
!$omp parallel do private(ix,im,ib,jb,i,j)
          do ix = 1, block%blksz(block_index)
            im = offset + ix - 1
            ib = block%index(block_index)%ii(ix)
            jb = block%index(block_index)%jj(ix)
            i = ib - block%isc + 1
            j = jb - block%jsc + 1
            destin_ptr(i,j,k) = factor * source_ptr(ib,jb,k,slice)
          enddo
        enddo
        localrc = ESMF_SUCCESS
      end if
    end if

    if (present(rc)) rc = localrc

  end subroutine block_copy_3dslice_r8_to_3d_r8

  !> ???
  !> 
  !> @param[in] destin_ptr ???
  !> @param[in] source_ptr ???
  !> @param[in] slice ???
  !> @param[in] block ???
  !> @param[in] block_index ???
  !> @param[in] scale_factor ???
  !> @param[in] offset ???
  !> @param[out] rc ???
  !>
  !> @author
  subroutine block_array_copy_3dslice_r8_to_3d_r8(destin_ptr, source_arr, slice, block, block_index, scale_factor, offset, rc)

    ! -- arguments
    real(ESMF_KIND_R8),        pointer     :: destin_ptr(:,:,:)
    real(kind=8),              intent(in)  :: source_arr(:,:,:,:)
    integer,                   intent(in)  :: slice
    type (block_control_type), intent(in)  :: block
    integer,                   intent(in)  :: block_index
    real(kind=8),    optional, intent(in)  :: scale_factor
    integer,                   intent(in)  :: offset
    integer,         optional, intent(out) :: rc

    ! -- local variables
    integer      :: localrc
    integer      :: i, ib, ix, im, j, jb, k
    real(kind=8) :: factor

    ! -- begin
    localrc = ESMF_RC_PTR_NOTALLOC
    if (associated(destin_ptr)) then
      localrc = ESMF_RC_VAL_OUTOFRANGE
      if (slice > 0 .and. slice <= size(source_arr, dim=4)) then
        factor = 1._8
        if (present(scale_factor)) factor = scale_factor
        do k = 1, size(source_arr, dim=3)
!$omp parallel do private(ix,im,ib,jb,i,j)
          do ix = 1, block%blksz(block_index)
            im = offset + ix - 1
            ib = block%index(block_index)%ii(ix)
            jb = block%index(block_index)%jj(ix)
            i = ib - block%isc + 1
            j = jb - block%jsc + 1
            destin_ptr(i,j,k) = factor * source_arr(i,j,k,slice)
          enddo
        enddo
        localrc = ESMF_SUCCESS
      end if
    end if

    if (present(rc)) rc = localrc

  end subroutine block_array_copy_3dslice_r8_to_3d_r8

  !> fill: 2D ???
  !> 
  !> @param[in] destin_ptr ???
  !> @param[in] fill_value ???
  !> @param[in] block ???
  !> @param[in] block_index ???
  !> @param[in] offset ???
  !> @param[out] rc ???
  !>
  !> @author
  subroutine block_fill_2d_r8(destin_ptr, fill_value, block, block_index, offset, rc)

    ! -- arguments
    real(ESMF_KIND_R8),        pointer     :: destin_ptr(:,:)
    real(ESMF_KIND_R8),        intent(in)  :: fill_value
    type (block_control_type), intent(in)  :: block
    integer,                   intent(in)  :: block_index
    integer,                   intent(in)  :: offset
    integer, optional,         intent(out) :: rc

    ! -- local variables
    integer :: localrc
    integer :: i, ib, ix, j, jb, im

    ! -- begin
    localrc = ESMF_RC_PTR_NOTALLOC
    if (associated(destin_ptr)) then
!$omp parallel do private(ix,im,ib,jb,i,j)
      do ix = 1, block%blksz(block_index)
        im = offset + ix - 1
        ib = block%index(block_index)%ii(ix)
        jb = block%index(block_index)%jj(ix)
        i = ib - block%isc + 1
        j = jb - block%jsc + 1
        destin_ptr(i,j) = fill_value
      enddo
      localrc = ESMF_SUCCESS
    end if

    if (present(rc)) rc = localrc

  end subroutine block_fill_2d_r8

  !> fill: 3D ???
  !> 
  !> @param[in] destin_ptr ???
  !> @param[in] fill_value ???
  !> @param[in] block ???
  !> @param[in] block_index ???
  !> @param[in] offset ???
  !> @param[out] rc ???
  !>
  !> @author
  subroutine block_fill_3d_r8(destin_ptr, fill_value, block, block_index, offset, rc)

    ! -- arguments
    real(ESMF_KIND_R8),        pointer     :: destin_ptr(:,:,:)
    real(ESMF_KIND_R8),        intent(in)  :: fill_value
    type (block_control_type), intent(in)  :: block
    integer,                   intent(in)  :: block_index
    integer,                   intent(in)  :: offset
    integer, optional,         intent(out) :: rc

    ! -- local variables
    integer :: localrc
    integer :: i, ib, ix, im, j, jb, k

    ! -- begin
    localrc = ESMF_RC_PTR_NOTALLOC
    if (associated(destin_ptr)) then
      do k = 1, size(destin_ptr, dim=3)
!$omp parallel do private(ix,im,ib,jb,i,j)
        do ix = 1, block%blksz(block_index)
          im = offset + ix - 1
          ib = block%index(block_index)%ii(ix)
          jb = block%index(block_index)%jj(ix)
          i = ib - block%isc + 1
          j = jb - block%jsc + 1
          destin_ptr(i,j,k) = fill_value
        enddo
      enddo
      localrc = ESMF_SUCCESS
    end if

    if (present(rc)) rc = localrc

  end subroutine block_fill_3d_r8

  !> copy/fill: 1D to 2D ???
  !> 
  !> @param[in] destin_ptr ???
  !> @param[in] source_ptr ???
  !> @param[in] fill_value ???
  !> @param[in] block ???
  !> @param[in] block_index ???
  !> @param[in] offset ???
  !> @param[out] rc ???
  !>
  !> @author
  subroutine block_copy_or_fill_1d_r8_to_2d_r8(destin_ptr, source_ptr, fill_value, block, block_index, offset, rc)

    ! -- arguments
    real(ESMF_KIND_R8),        pointer     :: destin_ptr(:,:)
    real(kind=8),              pointer     :: source_ptr(:)
    real(ESMF_KIND_R8),        intent(in)  :: fill_value
    type (block_control_type), intent(in)  :: block
    integer,                   intent(in)  :: block_index
    integer,                   intent(in)  :: offset
    integer, optional,         intent(out) :: rc

    ! -- begin
    if (present(rc)) rc = ESMF_RC_PTR_NOTALLOC

    if (associated(destin_ptr)) then
      if (associated(source_ptr)) then
        call block_copy_1d_r8_to_2d_r8(destin_ptr, source_ptr, block, block_index, offset=offset, rc=rc)
      else
        call block_fill_2d_r8(destin_ptr, fill_value, block, block_index, offset=offset, rc=rc)
      end if
    end if

  end subroutine block_copy_or_fill_1d_r8_to_2d_r8

  !> copy/fill: 1D slice to 2D ???
  !> 
  !> @param[in] destin_ptr ???
  !> @param[in] source_ptr ???
  !> @param[in] slice ???
  !> @param[in] fill_value ???
  !> @param[in] block ???
  !> @param[in] block_index ???
  !> @param[in] offset ???
  !> @param[out] rc ???
  !>
  !> @author
  subroutine block_copy_or_fill_1dslice_r8_to_2d_r8(destin_ptr, source_ptr, slice, fill_value, block, block_index, offset, rc)

    ! -- arguments
    real(ESMF_KIND_R8),        pointer     :: destin_ptr(:,:)
    real(kind=8),              pointer     :: source_ptr(:,:)
    integer,                   intent(in)  :: slice
    real(ESMF_KIND_R8),        intent(in)  :: fill_value
    type (block_control_type), intent(in)  :: block
    integer,                   intent(in)  :: block_index
    integer,                   intent(in)  :: offset
    integer, optional,         intent(out) :: rc

    ! -- begin
    if (present(rc)) rc = ESMF_RC_PTR_NOTALLOC

    if (associated(destin_ptr)) then
      if (associated(source_ptr)) then
        call block_copy_1dslice_r8_to_2d_r8(destin_ptr, source_ptr, slice, block, block_index, offset=offset, rc=rc)
      else
        call block_fill_2d_r8(destin_ptr, fill_value, block, block_index, offset=offset, rc=rc)
      end if
    end if

  end subroutine block_copy_or_fill_1dslice_r8_to_2d_r8

  !> copy/fill: 1D slice to 2D ???
  !> 
  !> @param[in] destin_ptr ???
  !> @param[in] source_ptr ???
  !> @param[in] slice1 ???
  !> @param[in] slice2 ???
  !> @param[in] fill_value ???
  !> @param[in] block ???
  !> @param[in] block_index ???
  !> @param[in] offset ???
  !> @param[out] rc ???
  !>
  !> @author
  subroutine block_copy_or_fill_1dslice2_r8_to_2d_r8(destin_ptr, source_ptr, slice1, slice2, fill_value, block, block_index, offset, rc)

    ! -- arguments
    real(ESMF_KIND_R8),        pointer     :: destin_ptr(:,:)
    real(kind=8),              pointer     :: source_ptr(:,:,:)
    integer,                   intent(in)  :: slice1
    integer,                   intent(in)  :: slice2
    real(ESMF_KIND_R8),        intent(in)  :: fill_value
    type (block_control_type), intent(in)  :: block
    integer,                   intent(in)  :: block_index
    integer,                   intent(in)  :: offset
    integer, optional,         intent(out) :: rc

    ! -- begin
    if (present(rc)) rc = ESMF_RC_PTR_NOTALLOC

    if (associated(destin_ptr)) then
      if (associated(source_ptr)) then
        call block_copy_1dslice2_r8_to_2d_r8(destin_ptr, source_ptr, slice1, slice2, block, block_index, offset=offset, rc=rc)
      else
        call block_fill_2d_r8(destin_ptr, fill_value, block, block_index, offset=offset, rc=rc)
      end if
    end if

  end subroutine block_copy_or_fill_1dslice2_r8_to_2d_r8

  !> copy/fill: 2D to 3D ???
  !> 
  !> @param[in] destin_ptr ???
  !> @param[in] source_ptr ???
  !> @param[in] fill_value ???
  !> @param[in] block ???
  !> @param[in] block_index ???
  !> @param[in] offset ???
  !> @param[out] rc ???
  !>
  !> @author
  subroutine block_copy_or_fill_2d_r8_to_3d_r8(destin_ptr, source_ptr, fill_value, block, block_index, offset, rc)

    ! -- arguments
    real(ESMF_KIND_R8),        pointer     :: destin_ptr(:,:,:)
    real(kind=8),              pointer     :: source_ptr(:,:)
    real(ESMF_KIND_R8),        intent(in)  :: fill_value
    type (block_control_type), intent(in)  :: block
    integer,                   intent(in)  :: block_index
    integer,                   intent(in)  :: offset
    integer, optional,         intent(out) :: rc

    ! -- begin
    if (present(rc)) rc = ESMF_RC_PTR_NOTALLOC

    if (associated(destin_ptr)) then
      if (associated(source_ptr)) then
        call block_copy_2d_r8_to_3d_r8(destin_ptr, source_ptr, block, block_index, offset=offset, rc=rc)
      else
        call block_fill_3d_r8(destin_ptr, fill_value, block, block_index, offset=offset, rc=rc)
      end if
    end if

  end subroutine block_copy_or_fill_2d_r8_to_3d_r8

  !> combine: 1D to 2D ???
  !> 
  !> @param[in] destin_ptr ???
  !> @param[in] fract1_ptr ???
  !> @param[in] fract2_ptr ???
  !> @param[in] block ???
  !> @param[in] block_index ???
  !> @param[in] offset ???
  !> @param[out] rc ???
  !>
  !> @author
  subroutine block_combine_frac_1d_r8_to_2d_r8(destin_ptr, fract1_ptr, fract2_ptr, block, block_index, offset, rc)

    ! -- arguments
    real(ESMF_KIND_R8),        pointer     :: destin_ptr(:,:)
    real(kind=8),              pointer     :: fract1_ptr(:)
    real(kind=8),              pointer     :: fract2_ptr(:)
    type(block_control_type),  intent(in)  :: block
    integer,                   intent(in)  :: block_index
    integer,                   intent(in)  :: offset
    integer, optional,         intent(out) :: rc

    ! -- local variables
    integer      :: localrc
    integer      :: i, ib, ix, im, j, jb

    ! -- begin
    localrc = ESMF_RC_PTR_NOTALLOC
    if (associated(destin_ptr) .and. &
        associated(fract1_ptr) .and. associated(fract2_ptr)) then
!$omp parallel do private(ix,im,ib,jb,i,j)
      do ix = 1, block%blksz(block_index)
        im = offset + ix - 1
        ib = block%index(block_index)%ii(ix)
        jb = block%index(block_index)%jj(ix)
        i = ib - block%isc + 1
        j = jb - block%jsc + 1
        destin_ptr(i,j) = fract1_ptr(im) * (1._8 - fract2_ptr(im))
      enddo
      localrc = ESMF_SUCCESS
    end if

    if (present(rc)) rc = localrc

  end subroutine block_combine_frac_1d_r8_to_2d_r8

  ! Real*4 Routines

  !> ???
  !> 
  !> @param[in] destin_ptr ???
  !> @param[in] source_ptr ???
  !> @param[in] block ???
  !> @param[in] block_index ???
  !> @param[in] scale_factor ???
  !> @param[in] offset ???
  !> @param[out] rc ???
  !>
  !> @author
  subroutine block_copy_1d_r4_to_2d_r8(destin_ptr, source_ptr, block, block_index, scale_factor, offset, rc)

    ! -- arguments
    real(ESMF_KIND_R8),        pointer     :: destin_ptr(:,:)
    real(kind=4),              pointer     :: source_ptr(:)
    type(block_control_type),  intent(in)  :: block
    integer,                   intent(in)  :: block_index
    real(kind=4),    optional, intent(in)  :: scale_factor
    integer,                   intent(in)  :: offset
    integer,         optional, intent(out) :: rc

    ! -- local variables
    integer      :: localrc
    integer      :: i, ib, ix, im, j, jb
    real(kind=4) :: factor

    ! -- begin
    localrc = ESMF_RC_PTR_NOTALLOC
    if (associated(destin_ptr) .and. associated(source_ptr)) then
      factor = 1._4
      if (present(scale_factor)) factor = scale_factor
!$omp parallel do private(ix,im,ib,jb,i,j)
      do ix = 1, block%blksz(block_index)
        im = offset + ix - 1
        ib = block%index(block_index)%ii(ix)
        jb = block%index(block_index)%jj(ix)
        i = ib - block%isc + 1
        j = jb - block%jsc + 1
        destin_ptr(i,j) = factor * source_ptr(im)
      enddo
      localrc = ESMF_SUCCESS
    end if

    if (present(rc)) rc = localrc

  end subroutine block_copy_1d_r4_to_2d_r8

  !> ???
  !> 
  !> @param[in] destin_ptr ???
  !> @param[in] source_ptr ???
  !> @param[in] block ???
  !> @param[in] block_index ???
  !> @param[in] scale_factor ???
  !> @param[in] special_value ???
  !> @param[in] offset ???
  !> @param[out] rc ???
  !>
  !> @author
  subroutine block_copy_spval_1d_r4_to_2d_r8(destin_ptr, source_ptr, block, block_index, scale_factor, special_value, offset, rc)

    ! -- arguments
    real(ESMF_KIND_R8),        pointer     :: destin_ptr(:,:)
    real(kind=4),              pointer     :: source_ptr(:)
    type(block_control_type),  intent(in)  :: block
    integer,                   intent(in)  :: block_index
    real(kind=4),              intent(in)  :: scale_factor
    real(kind=4),              intent(in)  :: special_value
    integer,                   intent(in)  :: offset
    integer,         optional, intent(out) :: rc

    ! -- local variables
    integer      :: localrc
    integer      :: i, ib, ix, im, j, jb

    ! -- begin
    localrc = ESMF_RC_PTR_NOTALLOC
    if (associated(destin_ptr) .and. associated(source_ptr)) then
!$omp parallel do private(ix,im,ib,jb,i,j)
       do ix = 1, block%blksz(block_index)
          im = offset + ix - 1
          ib = block%index(block_index)%ii(ix)
          jb = block%index(block_index)%jj(ix)
          i = ib - block%isc + 1
          j = jb - block%jsc + 1
          if (source_ptr(im) .ne. special_value) then
             destin_ptr(i,j) = scale_factor * source_ptr(im)
          else
             destin_ptr(i,j) = special_value
          end if
       enddo
       localrc = ESMF_SUCCESS
    end if

    if (present(rc)) rc = localrc

  end subroutine block_copy_spval_1d_r4_to_2d_r8

  !> copy: 1D slice to 2D ???
  !> 
  !> @param[in] destin_ptr ???
  !> @param[in] source_ptr ???
  !> @param[in] slice ???
  !> @param[in] block ???
  !> @param[in] block_index ???
  !> @param[in] scale_factor ???
  !> @param[in] offset ???
  !> @param[out] rc ???
  !>
  !> @author
  subroutine block_copy_1dslice_r4_to_2d_r8(destin_ptr, source_ptr, slice, block, block_index, scale_factor, offset, rc)

    ! -- arguments
    real(ESMF_KIND_R8),        pointer     :: destin_ptr(:,:)
    real(kind=4),              pointer     :: source_ptr(:,:)
    integer,                   intent(in)  :: slice
    type (block_control_type), intent(in)  :: block
    integer,                   intent(in)  :: block_index
    real(kind=4),    optional, intent(in)  :: scale_factor
    integer,                   intent(in)  :: offset
    integer,         optional, intent(out) :: rc

    ! -- local variables
    integer      :: localrc
    integer      :: i, ib, ix, im, j, jb
    real(kind=4) :: factor

    ! -- begin
    localrc = ESMF_RC_PTR_NOTALLOC
    if (associated(destin_ptr) .and. associated(source_ptr)) then
      localrc = ESMF_RC_VAL_OUTOFRANGE
      if (slice > 0 .and. slice <= size(source_ptr, dim=2)) then
        factor = 1._4
        if (present(scale_factor)) factor = scale_factor
!$omp parallel do private(ix,im,ib,jb,i,j)
        do ix = 1, block%blksz(block_index)
          im = offset + ix - 1
          ib = block%index(block_index)%ii(ix)
          jb = block%index(block_index)%jj(ix)
          i = ib - block%isc + 1
          j = jb - block%jsc + 1
          destin_ptr(i,j) = factor * source_ptr(im,slice)
        enddo
        localrc = ESMF_SUCCESS
      end if
    end if

    if (present(rc)) rc = localrc

  end subroutine block_copy_1dslice_r4_to_2d_r8

  !> copy: 1D slice to 2D ???
  !> 
  !> @param[in] destin_ptr ???
  !> @param[in] source_ptr ???
  !> @param[in] slice1 ???
  !> @param[in] slice2 ???
  !> @param[in] block ???
  !> @param[in] block_index ???
  !> @param[in] scale_factor ???
  !> @param[in] offset ???
  !> @param[out] rc ???
  !>
  !> @author
  subroutine block_copy_1dslice2_r4_to_2d_r8(destin_ptr, source_ptr, slice1, slice2, block, block_index, scale_factor, offset, rc)

    ! -- arguments
    real(ESMF_KIND_R8),        pointer     :: destin_ptr(:,:)
    real(kind=4),              pointer     :: source_ptr(:,:,:)
    integer,                   intent(in)  :: slice1
    integer,                   intent(in)  :: slice2
    type (block_control_type), intent(in)  :: block
    integer,                   intent(in)  :: block_index
    real(kind=4),    optional, intent(in)  :: scale_factor
    integer,                   intent(in)  :: offset
    integer,         optional, intent(out) :: rc

    ! -- local variables
    integer      :: localrc
    integer      :: i, ib, ix, im, j, jb
    real(kind=4) :: factor

    ! -- begin
    localrc = ESMF_RC_PTR_NOTALLOC
    if (associated(destin_ptr) .and. associated(source_ptr)) then
      localrc = ESMF_RC_VAL_OUTOFRANGE
      if (slice1 > 0 .and. slice1 <= size(source_ptr, dim=2) .and. slice2 > 0 .and. slice2 <= size(source_ptr, dim=3)) then
        factor = 1._4
        if (present(scale_factor)) factor = scale_factor
!$omp parallel do private(ix,im,ib,jb,i,j)
        do ix = 1, block%blksz(block_index)
          im = offset + ix - 1
          ib = block%index(block_index)%ii(ix)
          jb = block%index(block_index)%jj(ix)
          i = ib - block%isc + 1
          j = jb - block%jsc + 1
          destin_ptr(i,j) = factor * source_ptr(im,slice1,slice2)
        enddo
        localrc = ESMF_SUCCESS
      end if
    end if

    if (present(rc)) rc = localrc

  end subroutine block_copy_1dslice2_r4_to_2d_r8

  !> copy: 2D to 3D ???
  !> 
  !> @param[in] destin_ptr ???
  !> @param[in] source_ptr ???
  !> @param[in] block ???
  !> @param[in] block_index ???
  !> @param[in] scale_factor ???
  !> @param[in] offset ???
  !> @param[out] rc ???
  !>
  !> @author
  subroutine block_copy_2d_r4_to_3d_r8(destin_ptr, source_ptr, block, block_index, scale_factor, offset, rc)

    ! -- arguments
    real(ESMF_KIND_R8),        pointer     :: destin_ptr(:,:,:)
    real(kind=4),              pointer     :: source_ptr(:,:)
    type (block_control_type), intent(in)  :: block
    integer,                   intent(in)  :: block_index
    real(kind=4),    optional, intent(in)  :: scale_factor
    integer,                   intent(in)  :: offset
    integer,         optional, intent(out) :: rc

    ! -- local variables
    integer      :: localrc
    integer      :: i, ib, ix, im, j, jb, k
    real(kind=4) :: factor

    ! -- begin
    localrc = ESMF_RC_PTR_NOTALLOC
    if (associated(destin_ptr) .and. associated(source_ptr)) then
      factor = 1._4
      if (present(scale_factor)) factor = scale_factor
      do k = 1, size(source_ptr, dim=2)
!$omp parallel do private(ix,im,ib,jb,i,j)
        do ix = 1, block%blksz(block_index)
          im = offset + ix - 1
          ib = block%index(block_index)%ii(ix)
          jb = block%index(block_index)%jj(ix)
          i = ib - block%isc + 1
          j = jb - block%jsc + 1
          destin_ptr(i,j,k) = factor * source_ptr(im,k)
        enddo
      enddo
      localrc = ESMF_SUCCESS
    end if

    if (present(rc)) rc = localrc

  end subroutine block_copy_2d_r4_to_3d_r8

  !> copy: 2D to 2D ???
  !> 
  !> @param[in] destin_ptr ???
  !> @param[in] source_ptr ???
  !> @param[in] block ???
  !> @param[in] block_index ???
  !> @param[in] scale_factor ???
  !> @param[in] offset ???
  !> @param[out] rc ???
  !>
  !> @author
  subroutine block_copy_2d_r4_to_2d_r8(destin_ptr, source_ptr, block, block_index, scale_factor, offset, rc)

    ! -- arguments
    real(ESMF_KIND_R8),        pointer     :: destin_ptr(:,:)
    real(kind=4),              pointer     :: source_ptr(:,:)
    type (block_control_type), intent(in)  :: block
    integer,                   intent(in)  :: block_index
    real(kind=4),    optional, intent(in)  :: scale_factor
    integer,                   intent(in)  :: offset
    integer,         optional, intent(out) :: rc

    ! -- local variables
    integer      :: localrc
    integer      :: i, ib, ix, im, j, jb
    real(kind=4) :: factor

    ! -- begin
    localrc = ESMF_RC_PTR_NOTALLOC
    if (associated(destin_ptr) .and. associated(source_ptr)) then
      factor = 1._4
      if (present(scale_factor)) factor = scale_factor
!$omp parallel do private(ix,im,ib,jb,i,j)
      do ix = 1, block%blksz(block_index)
        im = offset + ix - 1
        ib = block%index(block_index)%ii(ix)
        jb = block%index(block_index)%jj(ix)
        i = ib - block%isc + 1
        j = jb - block%jsc + 1
        destin_ptr(i,j) = factor * source_ptr(ib,jb)
      enddo
      localrc = ESMF_SUCCESS
    end if

    if (present(rc)) rc = localrc

  end subroutine block_copy_2d_r4_to_2d_r8

  !> ???
  !> 
  !> @param[in] destin_ptr ???
  !> @param[in] source_arr ???
  !> @param[in] block ???
  !> @param[in] block_index ???
  !> @param[in] scale_factor ???
  !> @param[in] offset ???
  !> @param[out] rc ???
  !>
  !> @author
  subroutine block_array_copy_2d_r4_to_2d_r8(destin_ptr, source_arr, block, block_index, scale_factor, offset, rc)

    ! -- arguments
    real(ESMF_KIND_R8),        pointer     :: destin_ptr(:,:)
    real(kind=4),              intent(in)  :: source_arr(:,:)
    type (block_control_type), intent(in)  :: block
    integer,                   intent(in)  :: block_index
    real(kind=4),    optional, intent(in)  :: scale_factor
    integer,                   intent(in)  :: offset
    integer,         optional, intent(out) :: rc

    ! -- local variables
    integer      :: localrc
    integer      :: i, ib, ix, im, j, jb
    real(kind=4) :: factor

    ! -- begin
    localrc = ESMF_RC_PTR_NOTALLOC
    if (associated(destin_ptr)) then
      factor = 1._4
      if (present(scale_factor)) factor = scale_factor
!$omp parallel do private(ix,im,ib,jb,i,j)
      do ix = 1, block%blksz(block_index)
        im = offset + ix - 1
        ib = block%index(block_index)%ii(ix)
        jb = block%index(block_index)%jj(ix)
        i = ib - block%isc + 1
        j = jb - block%jsc + 1
        destin_ptr(i,j) = factor * source_arr(i,j)
      enddo
      localrc = ESMF_SUCCESS
    end if

    if (present(rc)) rc = localrc

  end subroutine block_array_copy_2d_r4_to_2d_r8

  !> copy: 3D to 3D ???
  !> 
  !> @param[in] destin_ptr ???
  !> @param[in] source_ptr ???
  !> @param[in] block ???
  !> @param[in] block_index ???
  !> @param[in] scale_factor ???
  !> @param[in] offset ???
  !> @param[out] rc ???
  !>
  !> @author
  subroutine block_copy_3d_r4_to_3d_r8(destin_ptr, source_ptr, block, block_index, scale_factor, offset, rc)

    ! -- arguments
    real(ESMF_KIND_R8),        pointer     :: destin_ptr(:,:,:)
    real(kind=4),              pointer     :: source_ptr(:,:,:)
    type (block_control_type), intent(in)  :: block
    integer,                   intent(in)  :: block_index
    real(kind=4),    optional, intent(in)  :: scale_factor
    integer,                   intent(in)  :: offset
    integer,         optional, intent(out) :: rc

    ! -- local variables
    integer      :: localrc
    integer      :: i, ib, ix, im, j, jb, k
    real(kind=4) :: factor

    ! -- begin
    localrc = ESMF_RC_PTR_NOTALLOC
    if (associated(destin_ptr) .and. associated(source_ptr)) then
      factor = 1._4
      if (present(scale_factor)) factor = scale_factor
      do k = 1, size(source_ptr, dim=3)
!$omp parallel do private(ix,im,ib,jb,i,j)
        do ix = 1, block%blksz(block_index)
          im = offset + ix - 1
          ib = block%index(block_index)%ii(ix)
          jb = block%index(block_index)%jj(ix)
          i = ib - block%isc + 1
          j = jb - block%jsc + 1
          destin_ptr(i,j,k) = factor * source_ptr(ib,jb,k)
        enddo
      enddo
      localrc = ESMF_SUCCESS
    end if

    if (present(rc)) rc = localrc

  end subroutine block_copy_3d_r4_to_3d_r8

  !> ???
  !> 
  !> @param[in] destin_ptr ???
  !> @param[in] source_arr ???
  !> @param[in] block ???
  !> @param[in] block_index ???
  !> @param[in] scale_factor ???
  !> @param[in] offset ???
  !> @param[out] rc ???
  !>
  !> @author
  subroutine block_array_copy_3d_r4_to_3d_r8(destin_ptr, source_arr, block, block_index, scale_factor, offset, rc)

    ! -- arguments
    real(ESMF_KIND_R8),        pointer     :: destin_ptr(:,:,:)
    real(kind=4),              intent(in)  :: source_arr(:,:,:)
    type (block_control_type), intent(in)  :: block
    integer,                   intent(in)  :: block_index
    real(kind=4), optional,    intent(in)  :: scale_factor
    integer,                   intent(in)  :: offset
    integer, optional,         intent(out) :: rc

    ! -- local variables
    integer      :: localrc
    integer      :: i, ib, ix, im, j, jb, k
    real(kind=4) :: factor

    ! -- begin
    localrc = ESMF_RC_PTR_NOTALLOC
    if (associated(destin_ptr)) then
      factor = 1._4
      if (present(scale_factor)) factor = scale_factor
      do k = 1, size(source_arr, dim=3)
!$omp parallel do private(ix,im,ib,jb,i,j)
        do ix = 1, block%blksz(block_index)
          im = offset + ix - 1
          ib = block%index(block_index)%ii(ix)
          jb = block%index(block_index)%jj(ix)
          i = ib - block%isc + 1
          j = jb - block%jsc + 1
          destin_ptr(i,j,k) = factor * source_arr(i,j,k)
        enddo
      enddo
      localrc = ESMF_SUCCESS
    end if

    if (present(rc)) rc = localrc

  end subroutine block_array_copy_3d_r4_to_3d_r8

  !> copy: 3D slice to 3D ???
  !> 
  !> @param[in] destin_ptr ???
  !> @param[in] source_ptr ???
  !> @param[in] slice ???
  !> @param[in] block ???
  !> @param[in] block_index ???
  !> @param[in] scale_factor ???
  !> @param[in] offset ???
  !> @param[out] rc ???
  !>
  !> @author
  subroutine block_copy_3dslice_r4_to_3d_r8(destin_ptr, source_ptr, slice, block, block_index, scale_factor, offset, rc)

    ! -- arguments
    real(ESMF_KIND_R8),        pointer     :: destin_ptr(:,:,:)
    real(kind=4),              pointer     :: source_ptr(:,:,:,:)
    integer,                   intent(in)  :: slice
    type (block_control_type), intent(in)  :: block
    integer,                   intent(in)  :: block_index
    real(kind=4),    optional, intent(in)  :: scale_factor
    integer,                   intent(in)  :: offset
    integer,         optional, intent(out) :: rc

    ! -- local variables
    integer      :: localrc
    integer      :: i, ib, ix, im, j, jb, k
    real(kind=4) :: factor

    ! -- begin
    localrc = ESMF_RC_PTR_NOTALLOC
    if (associated(destin_ptr) .and. associated(source_ptr)) then
      localrc = ESMF_RC_VAL_OUTOFRANGE
      if (slice > 0 .and. slice <= size(source_ptr, dim=4)) then
        factor = 1._4
        if (present(scale_factor)) factor = scale_factor
        do k = 1, size(source_ptr, dim=3)
!$omp parallel do private(ix,im,ib,jb,i,j)
          do ix = 1, block%blksz(block_index)
            im = offset + ix - 1
            ib = block%index(block_index)%ii(ix)
            jb = block%index(block_index)%jj(ix)
            i = ib - block%isc + 1
            j = jb - block%jsc + 1
            destin_ptr(i,j,k) = factor * source_ptr(ib,jb,k,slice)
          enddo
        enddo
        localrc = ESMF_SUCCESS
      end if
    end if

    if (present(rc)) rc = localrc

  end subroutine block_copy_3dslice_r4_to_3d_r8

  !> ???
  !> 
  !> @param[in] destin_ptr ???
  !> @param[in] source_arr ???
  !> @param[in] slice ???
  !> @param[in] block ???
  !> @param[in] block_index ???
  !> @param[in] scale_factor ???
  !> @param[in] offset ???
  !> @param[out] rc ???
  !>
  !> @author
  subroutine block_array_copy_3dslice_r4_to_3d_r8(destin_ptr, source_arr, slice, block, block_index, scale_factor, offset, rc)

    ! -- arguments
    real(ESMF_KIND_R8),        pointer     :: destin_ptr(:,:,:)
    real(kind=4),              intent(in)  :: source_arr(:,:,:,:)
    integer,                   intent(in)  :: slice
    type (block_control_type), intent(in)  :: block
    integer,                   intent(in)  :: block_index
    real(kind=4),    optional, intent(in)  :: scale_factor
    integer,                   intent(in)  :: offset
    integer,         optional, intent(out) :: rc

    ! -- local variables
    integer      :: localrc
    integer      :: i, ib, ix, im, j, jb, k
    real(kind=4) :: factor

    ! -- begin
    localrc = ESMF_RC_PTR_NOTALLOC
    if (associated(destin_ptr)) then
      localrc = ESMF_RC_VAL_OUTOFRANGE
      if (slice > 0 .and. slice <= size(source_arr, dim=4)) then
        factor = 1._4
        if (present(scale_factor)) factor = scale_factor
        do k = 1, size(source_arr, dim=3)
!$omp parallel do private(ix,im,ib,jb,i,j)
          do ix = 1, block%blksz(block_index)
            im = offset + ix - 1
            ib = block%index(block_index)%ii(ix)
            jb = block%index(block_index)%jj(ix)
            i = ib - block%isc + 1
            j = jb - block%jsc + 1
            destin_ptr(i,j,k) = factor * source_arr(i,j,k,slice)
          enddo
        enddo
        localrc = ESMF_SUCCESS
      end if
    end if

    if (present(rc)) rc = localrc

  end subroutine block_array_copy_3dslice_r4_to_3d_r8

  !> copy/fill: 1D to 2D ???
  !> 
  !> @param[in] destin_ptr ???
  !> @param[in] source_ptr ???
  !> @param[in] fill_value ???
  !> @param[in] block ???
  !> @param[in] block_index ???
  !> @param[in] offset ???
  !> @param[out] rc ???
  !>
  !> @author
  subroutine block_copy_or_fill_1d_r4_to_2d_r8(destin_ptr, source_ptr, fill_value, block, block_index, offset, rc)

    ! -- arguments
    real(ESMF_KIND_R8),        pointer     :: destin_ptr(:,:)
    real(kind=4),              pointer     :: source_ptr(:)
    real(ESMF_KIND_R8),        intent(in)  :: fill_value
    type (block_control_type), intent(in)  :: block
    integer,                   intent(in)  :: block_index
    integer,                   intent(in)  :: offset
    integer, optional,         intent(out) :: rc

    ! -- begin
    if (present(rc)) rc = ESMF_RC_PTR_NOTALLOC

    if (associated(destin_ptr)) then
      if (associated(source_ptr)) then
        call block_copy_1d_r4_to_2d_r8(destin_ptr, source_ptr, block, block_index, offset=offset, rc=rc)
      else
        call block_fill_2d_r8(destin_ptr, fill_value, block, block_index, offset=offset, rc=rc)
      end if
    end if

  end subroutine block_copy_or_fill_1d_r4_to_2d_r8

  !> copy/fill: 1D slice to 2D ???
  !> 
  !> @param[in] destin_ptr ???
  !> @param[in] source_ptr ???
  !> @param[in] slice ???
  !> @param[in] fill_value ???
  !> @param[in] block ???
  !> @param[in] block_index ???
  !> @param[in] offset ???
  !> @param[out] rc ???
  !>
  !> @author
  subroutine block_copy_or_fill_1dslice_r4_to_2d_r8(destin_ptr, source_ptr, slice, fill_value, block, block_index, offset, rc)

    ! -- arguments
    real(ESMF_KIND_R8),        pointer     :: destin_ptr(:,:)
    real(kind=4),              pointer     :: source_ptr(:,:)
    integer,                   intent(in)  :: slice
    real(ESMF_KIND_R8),        intent(in)  :: fill_value
    type (block_control_type), intent(in)  :: block
    integer,                   intent(in)  :: block_index
    integer,                   intent(in)  :: offset
    integer, optional,         intent(out) :: rc

    ! -- begin
    if (present(rc)) rc = ESMF_RC_PTR_NOTALLOC

    if (associated(destin_ptr)) then
      if (associated(source_ptr)) then
        call block_copy_1dslice_r4_to_2d_r8(destin_ptr, source_ptr, slice, block, block_index, offset=offset, rc=rc)
      else
        call block_fill_2d_r8(destin_ptr, fill_value, block, block_index, offset=offset, rc=rc)
      end if
    end if

  end subroutine block_copy_or_fill_1dslice_r4_to_2d_r8

  !> copy/fill: 1D slice to 2D ???
  !> 
  !> @param[in] destin_ptr ???
  !> @param[in] source_ptr ???
  !> @param[in] slice1 ???
  !> @param[in] slice2 ???
  !> @param[in] fill_value ???
  !> @param[in] block ???
  !> @param[in] block_index ???
  !> @param[in] offset ???
  !> @param[out] rc ???
  !>
  !> @author
  subroutine block_copy_or_fill_1dslice2_r4_to_2d_r8(destin_ptr, source_ptr, slice1, slice2, fill_value, block, block_index, offset, rc)

    ! -- arguments
    real(ESMF_KIND_R8),        pointer     :: destin_ptr(:,:)
    real(kind=4),              pointer     :: source_ptr(:,:,:)
    integer,                   intent(in)  :: slice1
    integer,                   intent(in)  :: slice2
    real(ESMF_KIND_R8),        intent(in)  :: fill_value
    type (block_control_type), intent(in)  :: block
    integer,                   intent(in)  :: block_index
    integer,                   intent(in)  :: offset
    integer, optional,         intent(out) :: rc

    ! -- begin
    if (present(rc)) rc = ESMF_RC_PTR_NOTALLOC

    if (associated(destin_ptr)) then
      if (associated(source_ptr)) then
        call block_copy_1dslice2_r4_to_2d_r8(destin_ptr, source_ptr, slice1, slice2, block, block_index, offset=offset, rc=rc)
      else
        call block_fill_2d_r8(destin_ptr, fill_value, block, block_index, offset=offset, rc=rc)
      end if
    end if

  end subroutine block_copy_or_fill_1dslice2_r4_to_2d_r8

  !> copy/fill: 2D to 3D ???
  !> 
  !> @param[in] destin_ptr ???
  !> @param[in] source_ptr ???
  !> @param[in] fill_value ???
  !> @param[in] block ???
  !> @param[in] block_index ???
  !> @param[in] offset ???
  !> @param[out] rc ???
  !>
  !> @author
  subroutine block_copy_or_fill_2d_r4_to_3d_r8(destin_ptr, source_ptr, fill_value, block, block_index, offset, rc)

    ! -- arguments
    real(ESMF_KIND_R8),        pointer     :: destin_ptr(:,:,:)
    real(kind=4),              pointer     :: source_ptr(:,:)
    real(ESMF_KIND_R8),        intent(in)  :: fill_value
    type (block_control_type), intent(in)  :: block
    integer,                   intent(in)  :: block_index
    integer,                   intent(in)  :: offset
    integer, optional,         intent(out) :: rc

    ! -- begin
    if (present(rc)) rc = ESMF_RC_PTR_NOTALLOC

    if (associated(destin_ptr)) then
      if (associated(source_ptr)) then
        call block_copy_2d_r4_to_3d_r8(destin_ptr, source_ptr, block, block_index, offset=offset, rc=rc)
      else
        call block_fill_3d_r8(destin_ptr, fill_value, block, block_index, offset=offset, rc=rc)
      end if
    end if

  end subroutine block_copy_or_fill_2d_r4_to_3d_r8

  !> combine: 1D to 2D ???
  !> 
  !> @param[in] destin_ptr ???
  !> @param[in] fract1_ptr ???
  !> @param[in] fract2_ptr ???
  !> @param[in] block ???
  !> @param[in] block_index ???
  !> @param[in] offset ???
  !> @param[out] rc ???
  !>
  !> @author
  subroutine block_combine_frac_1d_r4_to_2d_r8(destin_ptr, fract1_ptr, fract2_ptr, block, block_index, offset, rc)

    ! -- arguments
    real(ESMF_KIND_R8),        pointer     :: destin_ptr(:,:)
    real(kind=4),              pointer     :: fract1_ptr(:)
    real(kind=4),              pointer     :: fract2_ptr(:)
    type(block_control_type),  intent(in)  :: block
    integer,                   intent(in)  :: block_index
    integer,                   intent(in)  :: offset
    integer, optional,         intent(out) :: rc

    ! -- local variables
    integer :: localrc
    integer :: i, ib, ix, im, j, jb

    ! -- begin
    localrc = ESMF_RC_PTR_NOTALLOC
    if (associated(destin_ptr) .and. &
        associated(fract1_ptr) .and. associated(fract2_ptr)) then
!$omp parallel do private(ix,im,ib,jb,i,j)
      do ix = 1, block%blksz(block_index)
        im = offset + ix - 1
        ib = block%index(block_index)%ii(ix)
        jb = block%index(block_index)%jj(ix)
        i = ib - block%isc + 1
        j = jb - block%jsc + 1
        destin_ptr(i,j) = fract1_ptr(im) * (1._4 - fract2_ptr(im))
      enddo
      localrc = ESMF_SUCCESS
    end if

    if (present(rc)) rc = localrc

  end subroutine block_combine_frac_1d_r4_to_2d_r8

end module module_block_data
