       subroutine kpp_build_fftmatrix(ncup, localN, N, mpiUP_hat,&
                                      tmp_hat1, myrank, nprocs)
       use mpi
       integer:: ncup, localN, N, i, j, k, myrank, nprocs, dest1, dest2,&
                 status(MPI_STATUS_SIZE), ierr
       complex*16 :: mpiUP_hat(1:(N/2+1)*localN,1:ncup),&
             tmp_hat1(1:N*localN,1:ncup), tmp_v(1:(N/2-1)*localN,1:ncup)

       !if(myrank .eq. 0) then
      !do i=1, ncup
      !     write(*,*), mpiUP_hat(2,i)
      !end do 
      !end if
       !build tmp_hat1 for SVD
       do k=1, ncup
          do j=1, localN
            do i=1, N/2+1
                tmp_hat1(i + j * N, k) = mpiUP_hat(i + j * (N/2+1), k)
            end do
          end do 
       end do

       !get conjugate matrix of mpiUP_hat
       do k=1, ncup
         do j=1, localN
            do i=1, N/2-1
        tmp_v(i+ j*(N/2-1), k) = conjg(mpiUP_hat(N/2+1-i+ j*(N/2+1), k))
            end do 
         end do 
       end do 
        
       dest1 = nprocs - myrank -1
       dest2 = nprocs - myrank

       if(myrank .eq. 0 .or. (myrank .eq. nprocs/2))then
            dest2 = myrank
       end if 

       do k=1, ncup
         do j=2, localN
      call MPI_Sendrecv(tmp_v(j*(N/2-1), k), N/2-1, MPI_COMPLEX, dest1,&
            990, tmp_hat1(N/2+2 +(localN+1-j)*N, k), N/2-1, MPI_COMPLEX,&
            dest1,990, MPI_COMM_WORLD, status, ierr)   
         end do 
       end do 

       do k=1, ncup
         call MPI_Sendrecv(tmp_v(1, k), N/2-1, MPI_COMPLEX, dest1, 991,&
             tmp_hat1(N/2+2, k), N/2-1, MPI_COMPLEX, dest1,991,&
             MPI_COMM_WORLD, status, ierr)  
       end do 
      end 
        
