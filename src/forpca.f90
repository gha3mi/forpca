module forpca

   use kinds
   use foreig

   implicit none

   private
   public tpca

   !===============================================================================
   !> author: Seyed Ali Ghasemi
   type tpca
      integer                               :: ncol, nrow, npc
      real(rk), dimension(:,:), allocatable :: matrix
      real(rk), dimension(:,:), allocatable :: coeff
      real(rk), dimension(:,:), allocatable :: score
      real(rk), dimension(:,:), allocatable :: matrix_app
      real(rk), dimension(:,:), allocatable :: mean_data
      real(rk), dimension(:),   allocatable :: latent
      real(rk), dimension(:),   allocatable :: explained_variance
   contains
      procedure :: initialize
      procedure :: compute_coeff
      procedure :: compute_score
      procedure :: reconstruct_data
      procedure :: cmp_explained_variance
      procedure :: pca
      procedure :: dlloc => deallocate_tpca
   end type tpca
   !===============================================================================

contains

   !===============================================================================
   !> author: Seyed Ali Ghasemi
   pure subroutine initialize(this, matrix, npc)
      class(tpca),              intent(inout)        :: this
      real(rk), dimension(:,:), intent(in)           :: matrix
      integer,                  intent(in), optional :: npc

      this%matrix = matrix
      this%nrow = size(matrix,1)
      this%ncol = size(matrix,2)

      if (.not.present(npc)) then
         this%npc = this%ncol
      else
         this%npc = npc
      end if
   end subroutine initialize
   !===============================================================================


   !===============================================================================
   !> author: Seyed Ali Ghasemi
   pure subroutine compute_coeff(this)
      class(tpca), intent(inout)               :: this
      real(rk), dimension(this%ncol)           :: mean
      real(rk), dimension(this%ncol,this%ncol) :: cov
      logical,  dimension(this%ncol)           :: mask
      integer,  dimension(this%ncol)           :: order
      integer                                  :: i

      mean = sum(this%matrix, dim=1) / real(this%nrow, kind=rk)

      allocate(this%mean_data(this%nrow,this%ncol))
      do i = 1, this%nrow
         this%mean_data(i, :) = this%matrix(i, :) - mean
      end do

      cov = matmul(transpose(this%mean_data), this%mean_data)/real(this%nrow - 1, kind=rk)

      call eig(cov, this%coeff, this%latent)

      ! Sort
      mask = .true.
      do i = lbound(this%latent,1), ubound(this%latent,1)
         order(i) = maxloc(this%latent,1,mask)
         mask(order(i)) = .false.
      end do

      this%latent = this%latent(order)
      this%coeff  = this%coeff(:,order)

   end subroutine compute_coeff
   !===============================================================================


   !===============================================================================
   !> author: Seyed Ali Ghasemi
   pure subroutine compute_score(this)
      class(tpca), intent(inout) :: this
      integer                    :: i, j

      allocate(this%score(this%nrow,this%ncol))
      do i = 1, this%nrow
         do j = 1, this%npc
            this%score(i, j) = dot_product(this%mean_data(i, :), this%coeff(:, j))
         end do
      end do
   end subroutine compute_score
   !===============================================================================


   !===============================================================================
   !> author: Seyed Ali Ghasemi
   pure subroutine reconstruct_data(this)
      class(tpca), intent(inout)                :: this
      real(rk), dimension(this%nrow, this%ncol) :: X_centered
      real(rk), dimension(this%nrow, this%npc)  :: pca_X
      integer                                   :: i, j

      X_centered = this%matrix - this%mean_data
      pca_X = matmul(X_centered, this%coeff)
      this%matrix_app = matmul(pca_X, transpose(this%coeff)) + this%mean_data
   end subroutine reconstruct_data
   !===============================================================================


   !===============================================================================
   !> author: Seyed Ali Ghasemi
   pure subroutine cmp_explained_variance(this)
      class(tpca), intent(inout) :: this
      real(rk)                   :: sum_latent

      sum_latent = sum(this%latent(1:this%npc))
      this%explained_variance = this%latent / sum_latent
   end subroutine cmp_explained_variance
   !===============================================================================


   !===============================================================================
   !> author: Seyed Ali Ghasemi
   pure subroutine pca(this, matrix, npc, coeff, score, latent, explained, matrix_app)
      class(tpca),                           intent(inout)         :: this
      real(rk), dimension(:,:),              intent(in)            :: matrix
      integer,                               intent(in),  optional :: npc
      real(rk), dimension(:,:), allocatable, intent(out)           :: coeff
      real(rk), dimension(:,:), allocatable, intent(out), optional :: score
      real(rk), dimension(:),   allocatable, intent(out), optional :: latent
      real(rk), dimension(:),   allocatable, intent(out), optional :: explained
      real(rk), dimension(:,:), allocatable, intent(out), optional :: matrix_app

      call this%initialize(matrix, npc)

      call this%compute_coeff()
      coeff = this%coeff

      if(present(score)) then
         call this%compute_score()
         score = this%score
      end if

      if(present(latent)) then
         latent = this%latent
      end if

      if(present(explained)) then
         call this%cmp_explained_variance()
         explained = this%explained_variance*100.0_rk
      end if

      if(present(matrix_app)) then
         call this%reconstruct_data()
         matrix_app = this%matrix_app
      end if

   end subroutine pca
   !===============================================================================


   !===============================================================================
   !> author: Seyed Ali Ghasemi
   elemental pure subroutine deallocate_tpca(this)
      class(tpca), intent(inout) :: this

      if (allocated(this%matrix))     deallocate(this%matrix)
      if (allocated(this%coeff))      deallocate(this%coeff)
      if (allocated(this%mean_data))  deallocate(this%mean_data)
      if (allocated(this%score))      deallocate(this%score)
      if (allocated(this%matrix_app)) deallocate(this%matrix_app)
   end subroutine deallocate_tpca
   !===============================================================================

end module forpca
