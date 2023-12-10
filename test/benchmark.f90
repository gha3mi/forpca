program benchmark

   use kinds
   use forpca, only: tpca
   use fortime, only: timer

   implicit none

   real(rk), dimension(:,:), allocatable :: matrix
   real(rk), dimension(:,:), allocatable :: matrix_app
   real(rk), dimension(:,:), allocatable :: coeff
   real(rk), dimension(:,:), allocatable :: score
   real(rk), dimension(:),   allocatable :: latent
   real(rk), dimension(:),   allocatable :: explained
   type(tpca)                            :: p
   type(timer)                           :: t

   allocate(matrix(100,100))

#if defined(USE_COARRAY)
   sync all
   call t%timer_start()
   call p%pca(matrix, 50, 'svd', coeff, score, latent, explained, matrix_app)
   call t%timer_stop(message=' Elapsed time (benchmark: pca, coarray):')
   sync all
   call p%finalize()
#else
   call t%timer_start()
   call p%pca(matrix, 50, 'svd', coeff, score, latent, explained, matrix_app)
   call t%timer_stop(message=' Elapsed time (benchmark: pca):')   
#endif

end program benchmark
