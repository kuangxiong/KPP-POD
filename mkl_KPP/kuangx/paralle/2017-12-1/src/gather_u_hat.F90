subroutine gather_u_hat(N, ierr, u_hat, uu_hat,  local_M)
   use mpi
   integer :: ierr, num
   integer*8 :: N, i, j, local_M
   complex*16:: u_hat(1:N,1:local_M), uu_hat(1:N/2+1, 1:N)
   complex*16:: tmpu_hat(1:N/2+1, 1:local_M)
   do j=1, local_M
       do i = 1, N/2+1
           tmpu_hat(i, j) = u_hat(i, j)
       end do
   end do 
   num = (N/2+1) * local_M
   call MPI_Gather(tmpu_hat, num , MPI_COMPLEX, uu_hat, num , MPI_COMPLEX, 0, MPI_COMM_WORLD, ierr)
end


