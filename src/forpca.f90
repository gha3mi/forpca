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
      character(:),             allocatable :: method
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
      procedure :: finalize => deallocate_tpca
   end type tpca
   !===============================================================================

contains

   !===============================================================================
   !> author: Seyed Ali Ghasemi
#if defined(USE_COARRAY)
   impure subroutine initialize(this, matrix, npc, method)
#else
   pure subroutine initialize(this, matrix, npc, method)
#endif
      class(tpca),              intent(inout)        :: this
      real(rk), dimension(:,:), intent(in)           :: matrix
      integer,                  intent(in), optional :: npc
      character(*),             intent(in), optional :: method


      allocate(this%matrix(size(matrix,1),size(matrix,2)))

#if defined(USE_COARRAY)
      if (this_image() == 1) then
#endif
         this%matrix = matrix
         this%nrow = size(matrix,1)
         this%ncol = size(matrix,2)

         if (.not.present(npc)) then
            this%npc = this%ncol
         else
            this%npc = npc
         end if
#if defined(USE_COARRAY)
      end if
#endif

      if (.not.present(method)) then
         this%method = 'svd'
      else
         this%method = method
      end if

#if defined(USE_COARRAY)
      call co_broadcast(this%matrix, source_image=1)
      call co_broadcast(this%nrow, source_image=1)
      call co_broadcast(this%ncol, source_image=1)
      call co_broadcast(this%npc, source_image=1)
#endif
   end subroutine initialize
   !===============================================================================


   !===============================================================================
   !> author: Seyed Ali Ghasemi
#if defined(USE_COARRAY)
   impure subroutine compute_coeff(this)
#else
   pure subroutine compute_coeff(this)
#endif
      use forsvd
      use formatmul, only: matmul

      class(tpca), intent(inout)               :: this
      real(rk), dimension(this%ncol)           :: mean
      real(rk), dimension(this%ncol,this%ncol) :: cov
      integer                                  :: i
      real(rk), dimension(:,:), allocatable    :: mean_data_T

      allocate(this%mean_data(this%nrow,this%ncol))
      allocate(mean_data_T(this%ncol,this%nrow))

#if defined(USE_COARRAY)
      if (this_image() == 1) then
#endif
         mean = sum(this%matrix, dim=1) / real(this%nrow, kind=rk)

         do i = 1, this%nrow
            this%mean_data(i, :) = this%matrix(i, :) - mean
         end do
         mean_data_T = transpose(this%mean_data)
#if defined(USE_COARRAY)
      end if
#endif

#if defined(USE_COARRAY)
      call co_broadcast(mean_data_T, source_image=1)
      call co_broadcast(this%mean_data, source_image=1)
      cov = matmul(mean_data_T, this%mean_data, method='coarray', option='m2')/real(this%nrow - 1, kind=rk)
#else
      cov = matmul(mean_data_T, this%mean_data, method='default', option='m2')/real(this%nrow - 1, kind=rk)
#endif

      allocate(this%latent(this%ncol))
      allocate(this%coeff(this%ncol,this%ncol))
#if defined(USE_COARRAY)
      if (this_image() == 1) then
#endif
         select case(this%method)
          case('svd')
            block
               real(rk), dimension(this%ncol,this%ncol) :: U
               real(rk), dimension(this%ncol,this%ncol) :: VT
               real(rk), dimension(this%ncol)           :: S

               call svd(cov, U,S,VT)
               this%latent = S**2/(this%ncol-1)
               this%coeff  = transpose(VT)
            end block
          case('eig')
            block
               logical,  dimension(this%ncol) :: mask
               integer,  dimension(this%ncol) :: order

               call eig(cov, this%coeff, this%latent)

               ! Sort
               mask = .true.
               do i = lbound(this%latent,1), ubound(this%latent,1)
                  order(i) = maxloc(this%latent,1,mask)
                  mask(order(i)) = .false.
               end do

               this%latent = this%latent(order)
               this%coeff  = this%coeff(:,order)
            end block
         end select

#if defined(USE_COARRAY)
      end if
#endif

   end subroutine compute_coeff
   !===============================================================================


   !===============================================================================
   !> author: Seyed Ali Ghasemi
#if defined(USE_COARRAY)
   impure subroutine compute_score(this)
#else
   pure subroutine compute_score(this)
#endif
      class(tpca), intent(inout) :: this
      integer                    :: i, j

#if defined(USE_COARRAY)
      if (this_image() == 1) then
#endif
         allocate(this%score(this%nrow,this%ncol))
         do i = 1, this%nrow
            do j = 1, this%npc
               this%score(i, j) = dot_product(this%mean_data(i, :), this%coeff(:, j))
            end do
         end do
#if defined(USE_COARRAY)
      end if
#endif

   end subroutine compute_score
   !===============================================================================


   !===============================================================================
   !> author: Seyed Ali Ghasemi
#if defined(USE_COARRAY)
   impure subroutine reconstruct_data(this)
#else
   pure subroutine reconstruct_data(this)
#endif
      use formatmul, only: matmul

      class(tpca), intent(inout)             :: this
      real(rk), dimension(:, :), allocatable :: X_centered
      real(rk), dimension(:, :), allocatable :: pca_X
      real(rk), dimension(:, :), allocatable :: coeff_T
      integer                                :: i, j

      X_centered = this%matrix - this%mean_data

#if defined(USE_COARRAY)
      call co_broadcast(this%coeff, source_image=1)
      pca_X = matmul(X_centered, this%coeff,method='coarray', option='m2')
      call co_broadcast(pca_X, source_image=1)
      coeff_T = transpose(this%coeff)
      this%matrix_app = matmul(pca_X, coeff_T, method='coarray', option='m2') + this%mean_data
#else
      coeff_T = transpose(this%coeff)
      pca_X = matmul(X_centered, this%coeff, method='default', option='m2')
      this%matrix_app = matmul(pca_X, coeff_T, method='default', option='m2') + this%mean_data
#endif
   end subroutine reconstruct_data
   !===============================================================================


   !===============================================================================
   !> author: Seyed Ali Ghasemi
#if defined(USE_COARRAY)
   impure subroutine cmp_explained_variance(this)
#else
   pure subroutine cmp_explained_variance(this)
#endif
      class(tpca), intent(inout) :: this
      real(rk)                   :: sum_latent

#if defined(USE_COARRAY)
      if (this_image() == 1) then
#endif
         sum_latent = sum(this%latent(1:this%npc))
         this%explained_variance = this%latent / sum_latent
#if defined(USE_COARRAY)
      end if
#endif

   end subroutine cmp_explained_variance
   !===============================================================================


   !===============================================================================
   !> author: Seyed Ali Ghasemi
#if defined(USE_COARRAY)
   impure subroutine pca(this, matrix, npc, method, coeff, score, latent, explained, matrix_app)
#else
   pure subroutine pca(this, matrix, npc, method, coeff, score, latent, explained, matrix_app)
#endif
      class(tpca),                           intent(inout)         :: this
      real(rk), dimension(:,:),              intent(in)            :: matrix
      integer,                               intent(in),  optional :: npc
      character(*),                          intent(in),  optional :: method
      real(rk), dimension(:,:), allocatable, intent(out)           :: coeff
      real(rk), dimension(:,:), allocatable, intent(out), optional :: score
      real(rk), dimension(:),   allocatable, intent(out), optional :: latent
      real(rk), dimension(:),   allocatable, intent(out), optional :: explained
      real(rk), dimension(:,:), allocatable, intent(out), optional :: matrix_app

      call this%initialize(matrix, npc, method)

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
