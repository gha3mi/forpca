program test1

   use kinds
   use forpca

   implicit none

   real(rk), dimension(:,:), allocatable :: matrix
   real(rk), dimension(:,:), allocatable :: matrix_app
   real(rk), dimension(:,:), allocatable :: coeff
   real(rk), dimension(:,:), allocatable :: score
   real(rk), dimension(:),   allocatable :: latent
   real(rk), dimension(:),   allocatable :: explained
   type(tpca)                            :: p

   allocate(matrix(5,2))
   matrix(1, 1) = 0.123_rk
   matrix(1, 2) = 0.456_rk
   matrix(2, 1) = 0.789_rk
   matrix(2, 2) = 0.234_rk
   matrix(3, 1) = 0.567_rk
   matrix(3, 2) = 0.890_rk
   matrix(4, 1) = 0.123_rk
   matrix(4, 2) = 0.678_rk
   matrix(5, 1) = 0.901_rk
   matrix(5, 2) = 0.345_rk

   call p%pca(matrix, 2, coeff, score, latent, explained, matrix_app)
   print*, p%matrix_app(1, 1)
   print*, p%matrix_app(1, 2)
   print*, p%matrix_app(2, 1)
   print*, p%matrix_app(2, 2)
   print*, p%matrix_app(3, 1)
   print*, p%matrix_app(3, 2)
   print*, p%matrix_app(4, 1)
   print*, p%matrix_app(4, 2)
   print*, p%matrix_app(5, 1)
   print*, p%matrix_app(5, 2)
   call p%dlloc()

   deallocate(matrix, matrix_app, coeff, score, latent, explained)

end program test1
