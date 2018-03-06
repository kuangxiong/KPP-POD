
   program main
     implicit none
     integer :: i, j, k, kmax, N, iter, conv_flg, inniter
     real*8 :: eps, tau, theta, dt, laml, lamr, laml1, lamr1, &
               vaeps, y1, y2, Amp, t, cM

     real*8 :: to_etime
     integer :: tostime, toetime, clock_rate

     call system_clock(count_rate=clock_rate)
     call system_clock(count=tostime)

     N = 8 
     laml = 0.0D+0
     lamr = 10.0D+0

     eps = 1.0D+0
     tau = 1.0D+0
     theta = 1.0D+0
     dt = 1.0e-6
     Amp = 1.0D+0
     vaeps = 5e-3

     print*, 'step size', dt, 'theta', theta, 'Amp', Amp
     print*, 'Parameter tau', tau
     print*, 'The convergence criteria for the outer loop', vaeps

     laml1 = lamr - (dsqrt(5.0D+0) - 1)/2*(lamr - laml)
     cM = eps*laml1*laml1 + laml1*Amp*sqrt(theta*theta + 1.0D+0) &
          + 1.0D+0/tau
     print*, 'For lambda =', laml1, 'eps', eps, 'cM', cM
     call sparse_run(N, laml1, eps, tau, theta, dt, t, Amp, y1)

     print*, 'Stopping time:', t
     print*, 'left point C*', y1

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

     print*, 'Final result: ', 'left', laml1, y1, 'right', &
                   lamr1, y2
     print*, 'Final C value', (y1 + y2)/2

     call system_clock(count=toetime)
     to_etime = dble(toetime - tostime)/dble(clock_rate)
     print*, 'Total elasped time is', to_etime, ' seconds.'

   end program
