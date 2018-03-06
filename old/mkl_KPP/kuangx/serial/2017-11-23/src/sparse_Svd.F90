      subroutine sparse_Svd(ncup, N, temp_hat1, UP_hat, S, B)
      integer::ncup, N, i, j, k, lwork, info
      real*8 :: S(1:N*N)
      complex*16:: UP_hat(1:(N/2 + 1)*N, 1:ncup), temp_hat1(1:N*N,1:ncup),B(1:N*N,1:ncup)
      complex *16:: work_temp, CZERO, CONE
      complex*16, allocatable ::CT(:,:), work(:)
      real*8, allocatable::rwork(:)


     do k = 1, ncup
       do j = 1, N
         do i = 1, N/2 + 1
           temp_hat1(i+(j-1)*N,k) = UP_hat(i+(j-1)*(N/2+1),k)
         end do
       end do
       ! Symmetrize
       do j = 2, N
         do i = N/2 + 2, N
           temp_hat1(i+(j-1)*N,k)=conjg(UP_hat(N+2-i+(N+1-j)*(N/2+1),k))
         end do
       end do
       do i = N/2 + 2, N
         temp_hat1(i,k)=conjg(UP_hat(N+2-i,k))
       end do
     end do
     allocate(CT(1:ncup, 1:ncup))
     allocate(rwork(1:5*min(N*N, ncup)))
     S(1:N*N) = 0.0D+0
     lwork = -1
     call zgesvd('S', 'N', N*N, ncup, temp_hat1, N*N, S, B, &
                 N*N, CT, ncup, work_temp, lwork, rwork, info)

     ! Compute the SVD
     lwork = int(work_temp)
     allocate(work(1:lwork))
     call zgesvd('S', 'N', N*N, ncup, temp_hat1, N*N, S, B, &
                 N*N, CT, ncup, work, lwork, rwork, info)
     deallocate(work)
     deallocate(rwork)
     deallocate(CT)
      end
