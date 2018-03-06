   program main
  !  include"mpif.h" 
   use mpi 
   implicit none
   integer :: ierr, mycomm, i, j, k, kmax, iter, conv_flg, inniter
   real*8 :: eps, tau, theta, dt, laml, lamr, laml1, lamr1, &
             vaeps, y1, y2, Amp, t, cM

   real*8 :: to_etime
   integer :: myrank, nprocs, tostime, toetime, clock_rate
   
   integer*8 ::N, localN, local_0_start, alloc_local
   complex*16, allocatable::u_hat(:,:), uu_hat(:,:), Nu_hat(:,:)
     
   call MPI_Init(ierr)
!     call MPI_Comm_dup(MPI_COMM_WORLD, mycomm, ierr)
   call MPI_COMM_RANK(MPI_COMM_WORLD, myrank, ierr)
   call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)


   call system_clock(count_rate=clock_rate)
   call system_clock(count=tostime)

   N = 160 
   laml = 0.0D+0
   lamr = 10.0D+0

     eps = 1.0D+0
     tau = 1.0D+0
     theta = 1.0D+0
     dt = 1.0e-6
     Amp = 1.0D+0
     vaeps = 5e-3

     if(mod(N, nprocs) .ne. 0 ) then
         write(*,*)'err0r:nprocs shoule be divided by N'
         stop
     end if
     if(myrank .eq. 0) then
     print*, 'step size', dt, 'theta', theta, 'Amp', Amp
     print*, 'Parameter tau', tau
     print*, 'The convergence criteria for the outer loop', vaeps
     end if
     laml1 = lamr - (dsqrt(5.0D+0) - 1)/2*(lamr - laml)
     cM = eps*laml1*laml1 + laml1*Amp*sqrt(theta*theta + 1.0D+0) &
          + 1.0D+0/tau
     if(myrank .eq. 0) then
     print*, 'For lambda =', laml1, 'eps', eps, 'cM', cM
     end if
     call getLocalInf(N, ierr, alloc_local, localN, local_0_start)
   
     allocate(uu_hat(1:N, 1:localN))
     allocate(u_hat(1:N/2+1, 1:N))
     allocate(Nu_hat(1:N/2+1, 1:localN))
     do j=1, localN
         do i=1, N
            uu_hat(i, j) = 1
         end do
     end do
     call fftw_mpi_2d(N, ierr, uu_hat, Nu_hat, alloc_local, localN, local_0_start) 
     
!     call gather_u_hat(N, ierr, uu_hat, u_hat, localN)
 !    if (myrank==0) then
      call KPP_POD(N, Nu_hat, laml1, eps, tau, theta, dt, t, Amp, y1, localN, local_0_start, myrank, nprocs,  ierr)
!     end if
 !    if (myrank==0) then
 !        do j=1, localN
 !            do i=1, N
 !               write(*,*)'hahaha',N, localN, i, j, u_hat(i,j) 
 !            end do
 !        end do
 !   end if
    ! call sparse_run(N, laml1, eps, tau, theta, dt, t, Amp, y1)

    ! print*, 'Stopping time:', t
    ! print*, 'left point C*', y1

    ! lamr1 = laml + (dsqrt(5.0D+0) - 1)/2*(lamr - laml)
    ! cM = eps*lamr1*lamr1 + lamr1*Amp*sqrt(theta*theta + 1.0D+0) &
    !      + 1.0D+0/tau
    ! print*, 'For lambda =', lamr1, 'eps', eps, 'cM', cM
!
    ! call sparse_run(N, lamr1, eps, tau, theta, dt, t, Amp, y2)

    ! print*, 'Stopping time:', t
    ! print*, 'right point C*', y2

     ! The outer loop
    ! iter = 0
     !do while( abs(y1 - y2) .gt. vaeps .or. abs(lamr1-laml1) .ge. 0.6)
 
     !  do while( abs(y1 - y2) .gt. vaeps)
     !    iter = iter + 1
     !    print*, 'For step', iter
     !    if( y1 .gt. y2) then
     !      laml = laml1
     !      laml1 = lamr1
     !      y1 = y2

     !      lamr1 = laml + (dsqrt(5.0D+0) - 1)/2*(lamr - laml)
     !      cM = eps*lamr1*lamr1 + lamr1*Amp*sqrt(theta*theta + 1.0D+0) &
     !           + 1.0D+0/tau

     !      print*, 'For lambda =', lamr1, 'eps', eps, 'cM', cM
     !      call sparse_run(N, lamr1, eps, tau, theta, dt, t, Amp, y2)

     !      print*, 'Stopping time:', t
     !      print*, 'right point C*', y2
     !    else
     !      lamr = lamr1 
     !      lamr1 = laml1
     !      y2 = y1

     !      laml1 = lamr - (dsqrt(5.0D+0) - 1)/2*(lamr - laml)
     !      cM = eps*laml1*laml1 + laml1*Amp*sqrt(theta*theta + 1.0D+0) &
     !           + 1.0D+0/tau

     !      print*, 'For lambda =', laml1, 'eps', eps, 'cM', cM
     !      call sparse_run(N, laml1, eps, tau, theta, dt, t, Amp, y1)

     !      print*, 'Stopping time:', t
     !      print*, 'left point C*', y1
     !    end if
     !  end do

     if (myrank .eq. 1) then
     print*, 'Final result: ', 'left', laml1, y1, 'right', &
                   lamr1, y2
     print*, 'Final C value', (y1 + y2)/2

     call system_clock(count=toetime)
     to_etime = dble(toetime - tostime)/dble(clock_rate)
     print*, 'Total elasped time is', to_etime, ' seconds.'
     end if
     call MPI_FINALIZE(ierr)
   end program
