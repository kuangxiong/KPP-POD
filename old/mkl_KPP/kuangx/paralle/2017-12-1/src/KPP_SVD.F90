      subroutine sparse_Svd(ncup, N, temp_hat1, UP_hat, S, B)
      integer::ncup, N, i, j, k, lwork, info
      real*8 :: S(1:N*N)
      complex*16:: UP_hat(1:(N/2 + 1)*N, 1:ncup), &
                   temp_hat1(1:N*N,1:ncup),B(1:N*N,1:ncup)
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






      subroutine sparse_mpi_Svd(ncup, N, mpitmp_hat1, S, mpiB, localN,&
          myid, numprocs)
      include'mpif.h'
      integer::ncup, N, lwork, info, myid, numprocs,&
               status(MPI_STATUS_SIZE), myrow, mycol, npcol,&
               nprow, ONE 
      real*8 :: S(1: ncup)
      complex*16:: mpitmp_hat1(1:N*localN, 1:ncup),&
                    work_temp, mpiB(localN*N, localN*N)
      
      complex*16, allocatable :: work(:)
      real*8, allocatable:: rwork(:)
      integer :: context, descA(1:9), descU(1:9),&
              rwkoptint, wkopt, rwkopt
     
      nprow = numprocs
      npcol = 1

      call blacs_pinfo(myid, numprocs)
      call blacs_get(-1, 0, context)
      call blacs_gridinit(context, "Row", nprow, npcol)
      call blacs_gridinfo(context, nprow, npcol, myrow, mycol)

      descA=(/1, 0, N*N, ncup, 1, 1, 0, 0, localN*N/)
      descU=(/1, 0, N*N, N*N,  1, 1, 0, 0, localN*N/)
      
!     descA(1) = 1, descA(2) = 0, descA(3) = N*N, descA(4) = ncup 
!     descA(5) = 1, descA(6) = 1, descA(7) = 0, descA(8) = 0
!     descA(9) = localN *N
!     descU(1) = 1, descU(2) = 0, descU(3) = N*N, descU(4) = ncup,
!     descU(5) = 1, descU(6) = 1, descU(7) = 0, descU(8) = 0
!     descU(9) = localN *N
!      lwork = -1
!      info = 0
!     call pzgesvd('V', 'N', N*N, ncup, mpitmp_hat1, ONE, ONE, descA, S,&
!      mpiB, ONE, ONE, descU, NULL, NULL, NULL, NULL, wkopt, lwork,&
!      rwkopt, info)
!     rwkoptint = int(rwkopt)
!     lwork = int(wkopt)
      ONE = 1
      lwork = N*N*ncup*ncup  
      rwkoptint = N*N*ncup*ncup 
!      write(*,*), rwkoptint, lwork
      allocate(work(1:lwork))
      allocate(rwork(1:rwkoptint))
      call pzgesvd('V', 'N', N*N, ncup, mpitmp_hat1, ONE, ONE, descA, S,&
      mpiB, ONE, ONE, descU, NULL, NULL, NULL, NULL, work, lwork,&
      rwork, info)
!#endif 
!      

!      do k = 1, ncup
!       do j = 1, N
!         do i = 1, N/2 + 1
!           temp_hat1(i+(j-1)*N,k) = UP_hat(i+(j-1)*(N/2+1),k)
!         end do
!       end do
!       ! Symmetrize
!       do j = 2, N
!         do i = N/2 + 2, N
!           temp_hat1(i+(j-1)*N,k)=conjg(UP_hat(N+2-i+(N+1-j)*(N/2+1),k))
!         end do
!       end do
!       do i = N/2 + 2, N
!         temp_hat1(i,k)=conjg(UP_hat(N+2-i,k))
!       end do
!      end do
!      allocate(CT(1:ncup, 1:ncup))
!      allocate(rwork(1:5*min(N*N, ncup)))
!      S(1:N*N) = 0.0D+0
!      lwork = -1
!      call zgesvd('S', 'N', N*N, ncup, temp_hat1, N*N, S, B, &
!                 N*N, CT, ncup, work_temp, lwork, rwork, info)

     ! Compute the SVD
!      lwork = int(work_temp)
!      allocate(work(1:lwork))
!      call zgesvd('S', 'N', N*N, ncup, temp_hat1, N*N, S, B, &
!                 N*N, CT, ncup, work, lwork, rwork, info)
!      deallocate(work)
!      deallocate(rwork)
!      deallocate(CT)
!      deallocate(tmp_v)
      end
