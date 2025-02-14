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

  !> Copies a 3D slice from a source array to a destination array.
  !>
  !> This subroutine copies a 3D slice from the source array to the 
  !> destination array, applying a scale factor and an offset to the 
  !> values during the copy process.
  !>
  !> @param[inout] destin_ptr Pointer to the destination array.
  !> @param[in] source_ptr Pointer to the source array.
  !> @param[in] slice Integer specifying the slice to be copied.
  !> @param[in] block Integer specifying the block to be copied.
  !> @param[in] block_index Integer specifying the index of the block.
  !> @param[in] scale_factor Real(8) scale factor to be applied to the values.
  !> @param[in] offset Real(8) offset to be applied to the values.
  !> @param[out] rc Integer return code indicating success or failure.
  !>
  !> @author Raffaele Montuoro @date 7/1/21
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

  !> Copies a 3D slice from a source array to a destination pointer array.
  !>
  !> This subroutine copies a specified 3D slice from the source array to the 
  !> destination pointer array. The slice is determined by the provided slice 
  !> index and block information. The copied values can be scaled and offset 
  !> by the provided scale factor and offset values.
  !>
  !> @param[out] destin_ptr   Pointer to the destination array where the slice 
  !>                          will be copied.
  !> @param[in]  source_arr   Source array from which the slice will be copied.
  !> @param[in]  slice        Integer specifying the slice index to be copied.
  !> @param[in]  block        Integer specifying the block index.
  !> @param[in]  block_index  Integer specifying the block index within the slice.
  !> @param[in]  scale_factor Real(8) value to scale the copied values.
  !> @param[in]  offset       Real(8) value to offset the copied values.
  !> @param[out] rc           Integer return code indicating the success or 
  !>                          failure of the operation.
  !>
  !> @author Raffaele Montuoro @date 7/1/21
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

  !> Fills a 2D block of real*8 data with a specified fill value.
  !>
  !> @param[inout] destin_ptr  Pointer to the destination array.
  !> @param[in]    fill_value  The value to fill the block with.
  !> @param[in]    block       The block dimensions.
  !> @param[in]    block_index The index of the block to fill.
  !> @param[in]    offset      The offset to apply to the block index.
  !> @param[out]   rc          Return code indicating success or failure.
  !>
  !> @author Raffaele Montuoro @date 7/1/21
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

  !> Fills a 3D block of real*8 data with a specified fill value.
  !>
  !> This subroutine fills a 3D block of real*8 (double precision) data
  !> pointed to by `destin_ptr` with the value `fill_value`. The block
  !> is specified by `block` and `block_index`, and the filling starts
  !> at the given `offset`. The result code is returned in `rc`.
  !>
  !> @param[inout] destin_ptr Pointer to the destination 3D array to be filled.
  !> @param[in] fill_value The value to fill the block with.
  !> @param[in] block The block dimensions to be filled.
  !> @param[in] block_index The index of the block to be filled.
  !> @param[in] offset The starting offset for filling the block.
  !> @param[out] rc The result code indicating success or failure of the operation.
  !>
  !> @author Raffaele Montuoro @date 7/1/21
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

  !> Copies or fills a 1D real*8 array to a 2D real*8 array.
  !>
  !> This subroutine takes a 1D array and either copies its contents to a
  !> 2D array or fills the 2D array with a specified fill value. The
  !> operation is determined by the provided block and block_index.
  !>
  !> @param[inout] destin_ptr Pointer to the destination 2D real*8 array.
  !> @param[in] source_ptr Pointer to the source 1D real*8 array.
  !> @param[in] fill_value The value used to fill the destination array if
  !>                       the source array is not used.
  !> @param[in] block The block size or dimension for the operation.
  !> @param[in] block_index The index of the block to be copied or filled.
  !> @param[in] offset The offset to be applied during the copy or fill
  !>                   operation.
  !> @param[out] rc Return code indicating the success or failure of the
  !>                operation.
  !>
  !> @author Raffaele Montuoro @date 7/1/21
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

  !> Copies or fills a 1D slice of real*8 data to a 2D real*8 array.
  !>
  !> This subroutine either copies data from a 1D source array to a 2D
  !> destination array or fills the destination array with a specified
  !> fill value if the source array is not provided.
  !>
  !> @param[inout] destin_ptr Pointer to the destination 2D real*8 array.
  !> @param[in] source_ptr Pointer to the source 1D real*8 array.
  !> @param[in] slice Integer specifying the slice of the destination array to fill.
  !> @param[in] fill_value Real*8 value used to fill the destination array if source_ptr is not provided.
  !> @param[in] block Integer specifying the block size for the operation.
  !> @param[in] block_index Integer specifying the index of the block to operate on.
  !> @param[in] offset Integer specifying the offset within the block.
  !> @param[out] rc Integer return code indicating success or failure of the operation.
  !>
  !> @author Raffaele Montuoro @date 7/1/21
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

  !> Copies or fills a 1D slice of real*8 data to a 2D real*8 array.
  !>
  !> This subroutine either copies data from a 1D source array to a 2D 
  !> destination array or fills the destination array with a specified 
  !> fill value if the source data is not available.
  !>
  !> @param[inout] destin_ptr Pointer to the destination 2D real*8 array.
  !> @param[in] source_ptr Pointer to the source 1D real*8 array.
  !> @param[in] slice1 The first dimension of the slice to be copied or filled.
  !> @param[in] slice2 The second dimension of the slice to be copied or filled.
  !> @param[in] fill_value The value to fill the destination array if source data is not available.
  !> @param[in] block The block of data to be processed.
  !> @param[in] block_index The index of the block within the data.
  !> @param[in] offset The offset within the block where the data starts.
  !> @param[out] rc Return code indicating success or failure of the operation.
  !>
  !> @author Raffaele Montuoro @date 7/1/21
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

  !> Copies or fills a 2D real*8 array into a 3D real*8 array.
  !>
  !> This subroutine takes a 2D array and either copies its contents or 
  !> fills a 3D array with a specified fill value. The operation is 
  !> controlled by the provided block and block_index parameters.
  !>
  !> @param[out] destin_ptr Pointer to the destination 3D array.
  !> @param[in]  source_ptr Pointer to the source 2D array.
  !> @param[in]  fill_value Value used to fill the destination array if 
  !>                        the source array is not used.
  !> @param[in]  block      Specifies the block of the destination array 
  !>                        to be filled or copied into.
  !> @param[in]  block_index Index within the block where the operation 
  !>                         starts.
  !> @param[in]  offset     Offset to be applied during the copy or fill 
  !>                        operation.
  !> @param[out] rc         Return code indicating the success or failure 
  !>                        of the operation.
  !>
  !> @author Raffaele Montuoro @date 7/1/21
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

  !> Copies a 1D array of real(4) values to a 2D array of real(8) values.
  !>
  !> This subroutine performs a block copy from a 1D source array of 
  !> single-precision real numbers to a 2D destination array of 
  !> double-precision real numbers. The copy operation includes scaling 
  !> and offset adjustments.
  !>
  !> @param[out] destin_ptr Pointer to the destination 2D array (real(8)).
  !> @param[in]  source_ptr Pointer to the source 1D array (real(4)).
  !> @param[in]  block      The size of the block to be copied.
  !> @param[in]  block_index The starting index in the destination array.
  !> @param[in]  scale_factor The factor by which to scale the source values.
  !> @param[in]  offset      The offset to add to the scaled source values.
  !> @param[out] rc          Return code indicating success or failure of the operation.
  !>
  !> @author Raffaele Montuoro @date 7/1/21
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

  !> Copies a 1D array of single precision real values to a 2D array of 
  !> double precision real values, applying a scale factor and handling special 
  !> values.
  !>
  !> @param[out] destin_ptr Pointer to the destination 2D array of double precision real values.
  !> @param[in] source_ptr Pointer to the source 1D array of single precision real values.
  !> @param[in] block The block size for the copy operation.
  !> @param[in] block_index The index of the block to be copied.
  !> @param[in] scale_factor The factor by which to scale the source values.
  !> @param[in] special_value The special value in the source array to be handled.
  !> @param[in] offset The offset to be applied to the destination array.
  !> @param[out] rc Return code indicating the success or failure of the operation.
  !>
  !> @author Raffaele Montuoro @date 7/1/21
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

  !> Copies a 1D slice of real(4) data to a 2D real(8) array.
  !>
  !> This subroutine performs the operation of copying a 1-dimensional 
  !> slice of single precision real data to a 2-dimensional double 
  !> precision real array. The operation includes scaling and offsetting 
  !> the data as specified.
  !>
  !> @param[out] destin_ptr Pointer to the destination 2D real(8) array.
  !> @param[in]  source_ptr Pointer to the source 1D real(4) array.
  !> @param[in]  slice      Integer specifying the slice to be copied.
  !> @param[in]  block      Integer specifying the block of data.
  !> @param[in]  block_index Integer specifying the index within the block.
  !> @param[in]  scale_factor Real(8) value used to scale the source data.
  !> @param[in]  offset      Real(8) value used to offset the source data.
  !> @param[out] rc          Integer return code indicating success or failure.
  !>
  !> @author Raffaele Montuoro @date 7/1/21
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

  !> Copies a 1D slice of real(4) data to a 2D array of real(8) data.
  !>
  !> This subroutine performs a block copy from a 1D slice of single precision
  !> floating-point data to a 2D array of double precision floating-point data.
  !> The data is scaled and offset during the copy process.
  !>
  !> @param[out] destin_ptr Pointer to the destination 2D array (real(8)).
  !> @param[in]  source_ptr Pointer to the source 1D slice (real(4)).
  !> @param[in]  slice1     The first dimension of the slice to be copied.
  !> @param[in]  slice2     The second dimension of the slice to be copied.
  !> @param[in]  block      The block size for the copy operation.
  !> @param[in]  block_index The index of the block to be copied.
  !> @param[in]  scale_factor The factor by which to scale the source data.
  !> @param[in]  offset     The offset to be added to the scaled source data.
  !> @param[out] rc         Return code indicating success or failure of the operation.
  !>
  !> @author Raffaele Montuoro @date 7/1/21
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

  !> Copies a 2D real (kind=4) block to a 3D real (kind=8) block.
  !>
  !> This subroutine performs a copy operation from a 2D source array of 
  !> single precision real numbers to a 3D destination array of double 
  !> precision real numbers. The copy operation includes scaling and 
  !> offset adjustments.
  !>
  !> @param[out] destin_ptr Pointer to the destination 3D array (real, kind=8).
  !> @param[in]  source_ptr Pointer to the source 2D array (real, kind=4).
  !> @param[in]  block      Specifies the block dimensions and indices.
  !> @param[in]  block_index Index of the block to be copied.
  !> @param[in]  scale_factor Scaling factor to be applied during the copy.
  !> @param[in]  offset      Offset to be applied during the copy.
  !> @param[out] rc          Return code indicating success or failure of the operation.
  !>
  !> @author Raffaele Montuoro @date 7/1/21
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

  !> Copies a 2D block of real*4 data to a 2D block of real*8 data.
  !>
  !> This subroutine performs a copy operation from a source 2D array of 
  !> single precision real numbers (real*4) to a destination 2D array of 
  !> double precision real numbers (real*8). The copy operation includes 
  !> scaling and offsetting the source data.
  !>
  !> @param[out] destin_ptr Pointer to the destination 2D array (real*8).
  !> @param[in]  source_ptr Pointer to the source 2D array (real*4).
  !> @param[in]  block      The block size to be copied.
  !> @param[in]  block_index The index of the block to be copied.
  !> @param[in]  scale_factor The factor by which to scale the source data.
  !> @param[in]  offset     The offset to be added to the scaled source data.
  !> @param[out] rc         Return code indicating success or failure of the operation.
  !>
  !> @author Raffaele Montuoro @date 7/1/21
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

  !> Copies a 2D array of real(4) to a 2D array of real(8) with scaling and offset.
  !>
  !> This subroutine performs a copy operation from a 2D array of single precision
  !> floating-point numbers (real(4)) to a 2D array of double precision floating-point
  !> numbers (real(8)). The copy operation includes scaling and offset adjustments.
  !>
  !> @param[out] destin_ptr Pointer to the destination 2D array of real(8).
  !> @param[in]  source_arr 2D array of real(4) to be copied.
  !> @param[in]  block      Block size for the copy operation.
  !> @param[in]  block_index Index of the block to be copied.
  !> @param[in]  scale_factor Scaling factor to be applied to the source array elements.
  !> @param[in]  offset     Offset to be added to the scaled source array elements.
  !> @param[out] rc         Return code indicating the success or failure of the operation.
  !>
  !> @author Raffaele Montuoro @date 7/1/21
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

  !> Copies a 3D block of real*4 data to a 3D block of real*8 data.
  !>
  !> This subroutine performs a copy operation from a source 3D array of 
  !> single-precision real numbers (real*4) to a destination 3D array of 
  !> double-precision real numbers (real*8). The copy operation can include 
  !> scaling and offset adjustments.
  !>
  !> @param[out] destin_ptr Pointer to the destination 3D array (real*8).
  !> @param[in]  source_ptr Pointer to the source 3D array (real*4).
  !> @param[in]  block      The block size to be copied.
  !> @param[in]  block_index The index of the block to be copied.
  !> @param[in]  scale_factor The factor by which to scale the source data.
  !> @param[in]  offset     The offset to be added to the scaled source data.
  !> @param[out] rc         Return code indicating success or failure of the operation.
  !>
  !> @author Raffaele Montuoro @date 7/1/21
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

  !> Copies a 3D real(4) array to a 3D real(8) array with scaling and offset.
  !>
  !> This subroutine performs a copy operation from a source 3D array of 
  !> single precision real numbers (real(4)) to a destination 3D array of 
  !> double precision real numbers (real(8)). The copy operation includes 
  !> scaling and offset adjustments.
  !>
  !> @param[out] destin_ptr Pointer to the destination 3D array (real(8)).
  !> @param[in]  source_arr Source 3D array (real(4)).
  !> @param[in]  block      Block size or dimensions.
  !> @param[in]  block_index Index of the block to be copied.
  !> @param[in]  scale_factor Scaling factor to be applied during the copy.
  !> @param[in]  offset     Offset to be added during the copy.
  !> @param[out] rc         Return code indicating success or failure of the operation.
  !>
  !> @author Raffaele Montuoro @date 7/1/21
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

  !> Copies a 3D slice of real(4) data to a 3D real(8) array.
  !>
  !> This subroutine copies a 3D slice from a source array of real(4) 
  !> precision to a destination array of real(8) precision. The slice 
  !> is specified by the slice index and block information. The data 
  !> can be scaled and offset during the copy process.
  !>
  !> @param[out] destin_ptr Pointer to the destination array (real(8)).
  !> @param[in]  source_ptr Pointer to the source array (real(4)).
  !> @param[in]  slice      Index of the slice to be copied.
  !> @param[in]  block      Block information for the data.
  !> @param[in]  block_index Index of the block within the data.
  !> @param[in]  scale_factor Scaling factor to be applied to the data.
  !> @param[in]  offset     Offset to be applied to the data.
  !> @param[out] rc         Return code indicating success or failure.
  !>
  !> @author Raffaele Montuoro @date 7/1/21
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

  !> Copies a 3D slice of a real(4) array to a real(8) 3D array.
  !>
  !> This subroutine copies a specified 3D slice from a source array of 
  !> real(4) type to a destination array of real(8) type. The copy operation 
  !> can include scaling and offset adjustments.
  !>
  !> @param[out] destin_ptr Pointer to the destination array (real(8)).
  !> @param[in]  source_arr Source array (real(4)).
  !> @param[in]  slice      Index of the slice to be copied.
  !> @param[in]  block      Block size for the copy operation.
  !> @param[in]  block_index Index of the block to be copied.
  !> @param[in]  scale_factor Scaling factor to be applied during the copy.
  !> @param[in]  offset     Offset to be applied during the copy.
  !> @param[out] rc         Return code indicating success or failure of the operation.
  !>
  !> @author Raffaele Montuoro @date 7/1/21
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

  !> Copies or fills a 1D array of real(4) values to a 2D array of real(8) values.
  !>
  !> This subroutine takes a 1D array of real(4) values and either copies it or fills it into a 2D array of real(8) values.
  !> If the source array is smaller than the destination block, the remaining elements are filled with a specified fill value.
  !>
  !> @param[out] destin_ptr Pointer to the destination 2D array of real(8) values.
  !> @param[in]  source_ptr Pointer to the source 1D array of real(4) values.
  !> @param[in]  fill_value The value used to fill the remaining elements if the source array is smaller than the destination block.
  !> @param[in]  block      The size of the block to be copied or filled.
  !> @param[in]  block_index The index of the block in the destination array.
  !> @param[in]  offset     The offset in the destination array where the block starts.
  !> @param[out] rc         Return code indicating the success or failure of the operation.
  !>
  !> @author Raffaele Montuoro @date 7/1/21
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

  !> Copies or fills a 1D slice of real(4) data to a 2D real(8) array.
  !>
  !> This subroutine either copies data from a 1D slice of real(4) values
  !> to a 2D real(8) array or fills the 2D array with a specified fill value.
  !>
  !> @param[out] destin_ptr Pointer to the destination 2D real(8) array.
  !> @param[in]  source_ptr Pointer to the source 1D real(4) array.
  !> @param[in]  slice      Integer specifying the slice of the source array to copy.
  !> @param[in]  fill_value Real(8) value used to fill the destination array if copying is not performed.
  !> @param[in]  block      Integer specifying the block size.
  !> @param[in]  block_index Integer specifying the index of the block.
  !> @param[in]  offset     Integer specifying the offset for the destination array.
  !> @param[out] rc         Integer return code indicating success or failure of the operation.
  !>
  !> @author Raffaele Montuoro @date 7/1/21
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

  !> Copies or fills a 1D slice of real(4) data to a 2D array of real(8) data.
  !>
  !> This subroutine takes a 1D slice of real(4) data and either copies it
  !> or fills it into a 2D array of real(8) data. The operation is controlled
  !> by the provided parameters.
  !>
  !> @param[out] destin_ptr Pointer to the destination 2D array of real(8) data.
  !> @param[in]  source_ptr Pointer to the source 1D slice of real(4) data.
  !> @param[in]  slice1     The first dimension of the slice to be copied or filled.
  !> @param[in]  slice2     The second dimension of the slice to be copied or filled.
  !> @param[in]  fill_value The value to fill in the destination array if the source is not provided.
  !> @param[in]  block      The block of data to be processed.
  !> @param[in]  block_index The index of the block within the data.
  !> @param[in]  offset     The offset to be applied during the copy or fill operation.
  !> @param[out] rc         Return code indicating the success or failure of the operation.
  !>
  !> @author Raffaele Montuoro @date 7/1/21
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

  !> Copies or fills a 2D real(4) array into a 3D real(8) array.
  !>
  !> This subroutine takes a 2D array of real(4) values and either copies it
  !> into or fills a 3D array of real(8) values. The operation is controlled
  !> by the provided block and block_index parameters.
  !>
  !> @param[out] destin_ptr Pointer to the destination 3D real(8) array.
  !> @param[in]  source_ptr Pointer to the source 2D real(4) array.
  !> @param[in]  fill_value Value used to fill the destination array if needed.
  !> @param[in]  block      Specifies the block of data to be copied or filled.
  !> @param[in]  block_index Index of the block in the destination array.
  !> @param[in]  offset     Offset to be applied during the copy or fill operation.
  !> @param[out] rc         Return code indicating success or failure of the operation.
  !>
  !> @author Raffaele Montuoro @date 7/1/21
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

  !> Combines two 1D real(4) fractional arrays into a 2D real(8) array.
  !>
  !> @param destin_ptr Pointer to the destination 2D real(8) array.
  !> @param fract1_ptr Pointer to the first 1D real(4) fractional array.
  !> @param fract2_ptr Pointer to the second 1D real(4) fractional array.
  !> @param block Integer specifying the block size.
  !> @param block_index Integer specifying the block index.
  !> @param offset Integer specifying the offset.
  !> @param rc Integer return code.
  !>
  !> @author Raffaele Montuoro @date 7/1/21
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
