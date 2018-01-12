! Full K space sample with err indicator
   subroutine KPP_run(N, lam, eps, tau, theta, dt, t, Amp, y)
#include "/opt/fftw/fftw-3/include/fftw3.f"
     ! Arguments
     integer :: N
     real*8 :: lam, eps, tau, theta, dt, t, Amp, y
     ! Variables
     integer :: i, j, k, kmax, niter, nc, ncup, m, exiter, ts, initer,&
                exflg, m2
     real*8 :: pi, cM, h, C, tmax, cost, ineps, res, gamma,&
               lim, limold, sums, sacc, pc, t2, lims, limsold, pc2, pc3
     real*8, allocatable :: u(:, :), cosx(:), sinx(:), rhs(:, :), &
                            S(:), Bx(:,:), By(:,:), rwork(:), St(:),&
                            RBNx(:,:,:), RBNy(:,:,:), RBBNx(:,:,:), &
                            RBBNy(:,:,:), RBN(:,:,:), RBBN(:,:,:),& 
                            conv(:,:), convs(:,:), S2(:), rmtr(:,:),&
                            ovlap(:,:), rorth(:,:), rct(:,:)

     complex*16, allocatable :: u_hat(:, :),ux_hat(:, :), u_hat1(:),&
                                uy_hat(:,:),rhs_hat(:,:),ut_hat(:, :),&
                                UP_hat(:,:), BN(:,:), B(:,:), CT(:,:),&
                                temp_hat1(:,:),BN1(:,:), su_hat1(:),&
                                a_hat(:), su_hat(:,:), BNlap(:,:), &
                                BNx(:,:), BNy(:,:), BBNx(:,:),&
                                BBNy(:,:),BBN(:,:),BNt(:,:),Fixm(:,:),&
                                Varm(:,:), Varm1(:,:), at_hat(:), &
                                conv_hat(:,:),convs_hat(:,:), &
                                BBNyt(:,:),anew_hat(:),sunew_hat(:,:),&
                                BN2(:,:), C2(:,:), B2(:,:), tauq(:), &
                                taup(:), mtr(:,:), suma_hat(:)

     character( len = 512 ) :: Fileout, Sigout, Errout
     integer :: info, lwork
     complex*16, allocatable :: work(:)
     complex*16 :: work_temp, CZERO, CONE
     !real*8, allocatable :: work(:)
     !real*8 :: work_temp

     integer*8 :: p1, p2
     ! Declearization for external functions
     real*8, external :: dasum, dnrm2
     complex*16, external :: zdotc
     real :: to_etime
     integer :: tostime, toetime, clock_rate, svdstime, svdetime

     pi=dacos(-1.0D+0)
     t = 0.0D+0
     tmax = 0.01D+0
     t2 = 20.0D+0 - tmax  
     h = 1.0D+0/N
     cM = eps*lam*lam + lam*Amp*sqrt(theta*theta + 1.0D+0) + 1.0D+0/tau
     C = -lam*Amp*sqrt(theta*theta + 1.0D+0)
     ineps = 6e-4
     CZERO = dcmplx(0.0D+0, 0.0D+0)
     CONE = dcmplx(1.0D+0, 0.0D+0)

     !kmax = int(tmax)
     kmax = 1000
     nc = 100 ! Take the snapshots of u every nc steps
     ncup = int(tmax/(nc*dt)) ! Number of columns of UP
     pc = 0.9999999999d0
     pc2 = 0.9999999999d0
     pc3 = 0.999d0
     gamma = 0.5D+0
     initer = 10000

     call system_clock(count_rate=clock_rate)
     call system_clock(count=tostime)

     allocate(u(1:N, 1:N))
     allocate(cosx(1:N))
     allocate(sinx(1:N))
     allocate(u_hat(1:N/2 + 1, 1:N))
     allocate(u_hat1(1:N*N))
     allocate(ux_hat(1:N/2 + 1, 1:N))
     allocate(uy_hat(1:N/2 + 1, 1:N))
     allocate(rhs(1:N, 1:N))
     allocate(rhs_hat(1:N/2 + 1, 1:N))
     allocate(ut_hat(1:N/2 + 1, 1:N))
     allocate(UP_hat(1:(N/2 + 1)*N, 1:ncup))
     allocate(temp_hat1(1:N*N, 1:ncup))
     allocate(S(1:N*N))
     allocate(B(1:N*N, 1:ncup))
     allocate(su_hat(1:N/2 + 1, 1:N))
     allocate(su_hat1(1:N*N))
     allocate(Bx(1:N, 1:N))
     allocate(By(1:N, 1:N))
     allocate(conv_hat(1:N/2 + 1, 1:N))
     allocate(conv(1:N, 1:N))
     allocate(convs_hat(1:N/2 + 1, 1:N))
     allocate(convs(1:N, 1:N))
     allocate(sunew_hat(1:N/2 + 1, 1:N))

     forall(i = 1:N)
       cosx(i) = cos(2*pi*(i-1)*h)
     end forall

     forall(i = 1:N)
       sinx(i) = sin(2*pi*(i-1)*h)
     end forall 

     u(1:N, 1:N) = 1.0D+0

     call dfftw_plan_dft_r2c_2d(p1, N, N, u, u_hat, FFTW_ESTIMATE)
     call dfftw_plan_dft_c2r_2d(p2, N, N, u_hat, u, FFTW_ESTIMATE)

     write(Fileout, * ) lam
     Fileout = 'Lam' // Trim(AdjustL(Fileout)) // '.out'
     open(unit = 100, file = Fileout, form='formatted')
     write(100, *) 'stepsize', dt, 'lam', lam, 'eps', eps, 'cM', cM,&
                   'theta', theta, 'tau', tau, 'Amp', Amp
     print*, 'Initial time to get the sparse basis', tmax
     print*, 'Take the snapshots of u every', nc, ' steps'
     print*, 'The parameter for the error indicator', gamma 
     print*, 'pc', pc, 'pc2', pc2, 'pc3', pc3

     !write(100, *) 'Initial u:'
     !do i = 1, N
     !  do j = 1, N
     !    write(100, '(f16.8)', advance='no') u(i, j)
     !  end do
     !  write(100, *)
     !end do

     call dfftw_execute_dft_r2c(p1, u, u_hat)

     lim = 1.0D+0
     limold = 0.0D+0

     niter = 0
     exiter = 0
!     do while( t .le. tmax)
!     do while( niter .le. 400)
       niter = niter + 1
       cost = cos(t)
       t = t + dt
 !      print*,'u_hat', u_hat(10,10) 
       call KPP_ComputePlane(pi, cost, dt, N, C, Amp, lam,  theta, eps, u_hat, ux_hat, uy_hat, rhs_hat)

!     print*,'u_hat', u_hat(10, 10)
!      if(mod(niter, nc) .eq. 0) then
!        call zcopy((N/2 + 1)*N, u_hat, 1, UP_hat(:, niter/nc), 1)
!      end if
!
!      if(mod(exiter, 100) .eq. 0) then
!        call zcopy((N/2 + 1)*N, u_hat, 1, conv_hat, 1)
!      else
!        call zaxpy((N/2 + 1)*N, dcmplx(1.0D+0,0.0D+0), u_hat, 1, &
!                   conv_hat, 1)
!      end if
!
!   if(mod(exiter, 100) .eq. 99) then
!       limold = lim
!       lim = dlog(real(conv_hat(1, 1))/(100*N*N))/(t - 50*dt)
!
!        write(100, *) 'time', t, real(conv_hat(1, 1))/(100*N*N), &
!                      'lim', lim, &
!                      'diff', abs(lim - limold)/(100*dt), 'C*',&
!                      (cM + dlog(real(u_hat(1, 1))/(N*N))/t)/lam
!        call dfftw_execute_dft_c2r(p2, conv_hat, conv)
!        
!      !  write(100, *) dasum(N*N, conv, 1)*h*h/(100*N*N)
!       if(abs(lim - limold)/(100*dt) .lt. ineps) exit
!      end if

       exiter = exiter + 1
!      end do
      do j=1, N
          do i=1, N/2+1
            write(*,*) 'test:', u_hat(i, 1)
          end do 
       end do 

!!     print*, 'ncup', ncup
!     call system_clock(count=svdstime)
!!!     ! Copy 
!     call sparse_Svd(ncup, N, temp_hat1, UP_hat, S, B)
!     do i=1,ncup 
!        write(*,*)'S',S(i)
!     end do
!
!!     deallocate(UP_hat)
!!
!!     call system_clock(count=svdetime)
!!     to_etime = real(svdetime - svdstime)/real(clock_rate)
!!     print*, 'Elasped time for SVD ', to_etime, ' seconds.'
!!     ! The singular value
!!     write(Sigout, * ) lam
!!     Sigout = 'Lam' // Trim(AdjustL(Sigout)) // 'sig'
!!     open(unit=11, file= Sigout, form='formatted')
!!     write(11, *) 'for lambda =', lam
!!     do i = 1, ncup
!!       write(11, '(f20.12)') S(i)
!!     end do
!!     close(11)
!!
!!     sums = sum(S(1:ncup))
!!     sacc = 0.d0
!!     do i = 1, N
!!       sacc = sacc + S(i)
!!       if( sacc .gt. pc*sums) exit
!!     end do
!!
!!     m = min(i, ncup)
!!
!!     allocate(ovlap(1:m, 1:m))
!!     allocate(Fixm(1:m, 1:m))
!!     call dgemm('C', 'N', m, m, 2*N*N, 1.0D+0, B,&
!!                2*N*N, B, 2*N*N, 0.0D+0, ovlap, m)
!!     deallocate(ovlap)
!!     deallocate(Fixm)
!!     
!!     print*, 'Number of basis', m
!!
!!     allocate(BN(1:(N/2+1)*N, 1:m))
!!     allocate(BN1(1:N*N, 1:m))
!!     allocate(BNlap(1:(N/2+1)*N, 1:m))
!!     allocate(BNx(1:(N/2+1)*N, 1:m))
!!     allocate(BNy(1:(N/2+1)*N, 1:m))
!!     allocate(RBNx(1:N, 1:N, 1:m))
!!     allocate(RBNy(1:N, 1:N, 1:m))
!!     allocate(RBBNx(1:N, 1:N, 1:m))
!!     allocate(RBBNy(1:N, 1:N, 1:m))
!!     allocate(BBNx(1:(N/2+1)*N, 1:m))
!!     allocate(BBNy(1:(N/2+1)*N, 1:m))
!!     allocate(a_hat(1:m))
!!     allocate(RBN(1:N, 1:N, 1:m))
!!     allocate(RBBN(1:N, 1:N, 1:m))
!!     allocate(BBN(1:(N/2+1)*N, 1:m))
!!     allocate(BNt(1:(N/2+1)*N, 1:m))
!!     allocate(Fixm(1:m, 1:m))
!!     allocate(Varm(1:m, 1:m))
!!     allocate(Varm1(1:m, 1:m))
!!     allocate(at_hat(1:m))
!!     allocate(BBNyt(1:(N/2+1)*N, 1:m))
!!     allocate(anew_hat(1:m))
!!     allocate(suma_hat(1:m))
!!
!!     !call dcopy(2*N*N*m, temp_hat1, 1, BN1, 1)
!!     call dcopy(2*N*N*m, B, 1, BN1, 1)
!!
!!     ! Copy 
!!     deallocate(B)
!!     do k = 1, m
!!       do j = 1, N
!!         do i = 1, N/2 + 1
!!           BN(i+(j-1)*(N/2+1), k) = BN1(i+(j-1)*N,k)
!!         end do
!!       end do
!!     end do
!!
!!     call zgemm('C', 'N', m, m, N*N, dcmplx(1.0D+0,0.0D+0), BN1,&
!!                N*N, BN1, N*N, dcmplx(0.0D+0,0.0D+0), &
!!                Fixm, m)
!!   do j = 1, N
!!      do i = 1, N/2 + 1
!!        u_hat1(i+(j-1)*N) = u_hat(i,j)
!!      end do
!!    end do
!!    ! Symmetrize
!!    do j = 2, N
!!      do i = N/2 + 2, N
!!        u_hat1(i+(j-1)*N)=conjg(u_hat(N+2-i, N+2-j))
!!      end do
!!    end do
!!    do i = N/2 + 2, N
!!      u_hat1(i)=conjg(u_hat(N+2-i,1))
!!    end do
!!
!!     call zgemv('C', N*N, m, dcmplx(1.0D+0,0.0D+0), BN1, N*N,&
!!                u_hat1, 1, dcmplx(0.0D+0,0.0D+0), a_hat, 1)
!!
!!     call zgemv('N', N*N, m, dcmplx(1.0D+0,0.0D+0), BN1, N*N,&
!!                a_hat, 1, dcmplx(0.0D+0,0.0D+0), su_hat1, 1)
!!
!!     print*, 'difference', dnrm2(N*N*2, su_hat1-u_hat1, 1)
!!     ! Form the matrices
!!     ! Fixm = eps*BN'BNlap + 2*lam*eps*BN'BNx + C*I
!!
!!     ! BNlap(1:(N/2+1)*N, 1:m) is the laplacian of BN
!!     call sparse_BuildBNlap(m, N, pi, BNlap, BN)
!!     call sparse_CopyMatrix_Symmetrize(m, N, temp_hat1, BNlap)
!!
!!     call zgemm('C', 'N', m, m, N*N, dcmplx(1.0D+0,0.0D+0), BN1,&
!!                N*N, temp_hat1, N*N, dcmplx(0.0D+0,0.0D+0), &
!!                Fixm, m)
!!
!!     call sparse_BuildBNx(m, N, pi, BNx, BN)
!!
!!     call sparse_CopyMatrix_Symmetrize(m, N, temp_hat1, BNx)
!!!
!!     call zgemm('C', 'N', m, m, N*N, dcmplx(2*eps*lam,0.0D+0), &
!!                BN1, N*N, temp_hat1, N*N, dcmplx(eps,0.0D+0), &
!!                Fixm, m)
!!
!!     do i = 1, m
!!       Fixm(i, i) = Fixm(i, i) + C
!!     end do
!!
!!     ! BNy(1:(N/2+1)*N, 1:m) is the derivatives of BN to y
!!     call sparse_BuildBNy(m, N, pi, BNy, BN)
!!
!!     ! RBNx(1:N, 1:N, 1:m) is the derivatives of BN to x in real space
!!     call zcopy((N/2 + 1)*N*m, BNx, 1, BNt, 1)
!!     do k = 1, m
!!       call dfftw_execute_dft_c2r(p2, BNt(:,k), RBNx(:,:,k))
!!       call dscal(N*N, 1.0D+0/(N*N), RBNx(:,:,k), 1)
!!     end do
!!
!!     ! RBNy(1:N, 1:N, 1:m) is the derivatives of BN to y in real space
!!     do k = 1, m
!!       call dfftw_execute_dft_c2r(p2, BNy(:,k), RBNy(:,:,k))
!!       call dscal(N*N, 1.0D+0/(N*N), RBNy(:,:,k), 1)
!!     end do
!!      
!!     call sparse_BuildRBBNxRBBNy1(m, N, cosx, Amp, Bx, By, RBNx, RBNy, RBBNx, RBBNy)
!!     ! The part of Bx without time t
!!     do k = 1, m
!!       call dfftw_execute_dft_r2c(p1, RBBNx(:,:,k), BBNx(:,k))
!!     end do
!!
!!     ! BBNy(1:(N/2+1)*N, 1:m) is RBBNy without time t in k space
!!     do k = 1, m
!!       call dfftw_execute_dft_r2c(p1, RBBNy(:,:,k), BBNy(:,k))
!!     end do
!!
!!     ! BBNy = BBNy + BBNx
!!     call zaxpy((N/2+1)*N*m, dcmplx(1.0D+0,0.0D+0), BBNx, 1, BBNy, 1)
!!
!!     ! RBN(1:N, 1:N, 1:m) is BN in real space
!!     call zcopy((N/2 + 1)*N*m, BN, 1, BNt, 1)
!!     do k = 1, m
!!       call dfftw_execute_dft_c2r(p2, BNt(:,k), RBN(:,:,k))
!!       call dscal(N*N, 1.0D+0/(N*N), RBN(:,:,k), 1)
!!     end do
!!
!!     ! RBBN(1:N, 1:N, 1:m) is Bx*RBN in real space without time t
!!     do k = 1, m
!!       forall(j = 1:N)
!!         forall(i = 1:N)
!!           RBBN(i, j, k) = Bx(i, j)*RBN(i, j, k)
!!         end forall
!!       end forall
!!     end do
!!
!!     ! BBN(1:(N/2+1)*N, 1:m) is RBBN in k space
!!     do k = 1, m
!!       call dfftw_execute_dft_r2c(p1, RBBN(:,:,k), BBN(:,k))
!!     end do
!!
!!     ! BBNy = BBNy + lam*BBN
!!     call zaxpy((N/2+1)*N*m, dcmplx(lam,0.0D+0), BBN, 1, BBNy, 1)
!!
!!     call sparse_CopyMatrix_Symmetrize(m, N, temp_hat1, BBNy)
!!
!!     call zgemm('C', 'N', m, m, N*N, dcmplx(1.0D+0,0.0D+0), &
!!                BN1,N*N, temp_hat1, N*N, dcmplx(1.0D+0,0.0D+0), &
!!                Fixm, m)
!!
!!     call sparse_BuildRBBNxRBBNy2(m, N, sinx, Amp, theta, Bx, By, RBNx, RBNy,&
!!     RBBNx, RBBNy)
!!
!!     ! BBNx(1:(N/2+1)*N, 1:m) is RBBNx with time t in k space
!!     do k = 1, m
!!       call dfftw_execute_dft_r2c(p1, RBBNx(:,:,k), BBNx(:,k))
!!     end do
!!
!!     ! BBNy(1:(N/2+1)*N, 1:m) is RBBNy with time t in k space
!!     do k = 1, m
!!       !call dfftw_execute_dft_r2c(p1, RBBNy(:,:,k), BBNy(:,k))
!!       call dfftw_execute_dft_r2c(p1, RBBNy(:,:,k), BBNyt(:,k))
!!     end do
!!
!!     ! BBNy = BBNy + BBNx
!!     !call zaxpy((N/2+1)*N*m, dcmplx(1.0D+0,0.0D+0), BBNx, 1, BBNy, 1)
!!     call zaxpy((N/2+1)*N*m, dcmplx(1.0D+0,0.0D+0), BBNx, 1, BBNyt, 1)
!!
!!     ! RBBN(1:N, 1:N, 1:m) is Bx*RBN in real space with time t
!!     do k = 1, m
!!       forall(j = 1:N)
!!         forall(i = 1:N)
!!           RBBN(i, j, k) = Bx(i, j)*RBN(i, j, k)
!!         end forall
!!       end forall
!!     end do
!!
!!     ! BBN(1:(N/2+1)*N, 1:m) is RBBN in k space
!!     do k = 1, m
!!       call dfftw_execute_dft_r2c(p1, RBBN(:,:,k), BBN(:,k))
!!     end do
!!
!!     ! BBNy = BBNy + lam*BBN
!!     !call zaxpy((N/2+1)*N*m, dcmplx(lam,0.0D+0), BBN, 1, BBNy, 1)
!!     call zaxpy((N/2+1)*N*m, dcmplx(lam,0.0D+0), BBN, 1, BBNyt, 1)
!!     call sparse_CopyMatrix_Symmetrize(m, N, temp_hat1, BBNyt)
!!
!!     call zgemm('C', 'N', m, m, N*N, dcmplx(1.0D+0,0.0D+0), &
!!                BN1, N*N, temp_hat1, N*N, dcmplx(0.0D+0,0.0D+0), &
!!                Varm1, m)
!!
!!     deallocate(temp_hat1)
!!
!!     lims = lim
!!     limsold = limold
!!
!!     write(Errout, * ) lam
!!     Errout = 'Erridx' // Trim(AdjustL(Errout)) // '.out'
!!     open(unit = 22, file = Errout, form='formatted')
!!
!!     exiter = 0
!!     ts = 0
!!     exflg = 0
!!     do while( t .le. tmax + t2)
!!       cost = cos(t)
!!       t = t + dt
!!
!!       !call sparse_ComputePlane(pi, cost, dt, N, C, Amp, lam,  theta, eps, u_hat, ux_hat, uy_hat, rhs_hat)
!!       forall(j = 1:m)
!!         forall(i = 1:m)
!!           Varm(i, j) = Fixm(i, j) + cost*Varm1(i, j)
!!         end forall
!!       end forall
!!
!!       ! update a_hat
!!       call zgemv('N', m, m, dcmplx(dt,0.0D+0), Varm, m, a_hat, 1, &
!!                  dcmplx(0.0D+0,0.0D+0), at_hat, 1)
!!       do i = 1, m
!!         !a_hat(i) = a_hat(i) + at_hat(i)
!!         anew_hat(i) = a_hat(i) + at_hat(i)
!!       end do
!!
!!       ! The error indicator
!!       if(mod(exiter, 100) .eq. 0) then
!!         call zgemv('N', (N/2+1)*N, m, dcmplx(1.0D+0,0.0D+0), BN, &
!!                (N/2+1)*N, a_hat, 1, dcmplx(0.0D+0,0.0D+0), su_hat, 1)
!!         call zgemv('N', (N/2+1)*N, m, dcmplx(1.0D+0,0.0D+0), BN, &
!!           (N/2+1)*N, anew_hat, 1, dcmplx(0.0D+0,0.0D+0), sunew_hat, 1)
!!
!!         call zgemv('N', (N/2+1)*N, m, dcmplx(eps,0.0D+0), BNlap, &
!!           (N/2+1)*N, a_hat, 1, dcmplx(0.0D+0,0.0D+0), rhs_hat, 1)
!!         call zgemv('N', (N/2+1)*N, m, dcmplx(2*eps*lam,0.0D+0), BNx, &
!!           (N/2+1)*N, a_hat, 1, dcmplx(1.0D+0,0.0D+0), rhs_hat, 1)
!!         call zgemv('N', (N/2+1)*N, m, dcmplx(1.0D+0,0.0D+0), BBNy, &
!!           (N/2+1)*N, a_hat, 1, dcmplx(1.0D+0,0.0D+0), rhs_hat, 1)
!!         call zgemv('N', (N/2+1)*N, m, dcmplx(cost,0.0D+0), BBNyt, &
!!           (N/2+1)*N, a_hat, 1, dcmplx(1.0D+0,0.0D+0), rhs_hat, 1)
!!         call zgemv('N', (N/2+1)*N, m, dcmplx(C,0.0D+0), BN, &
!!           (N/2+1)*N, a_hat, 1, dcmplx(1.0D+0,0.0D+0), rhs_hat, 1)
!!
!!         do j = 1, N
!!           do i = 1, N/2 + 1
!!             rhs_hat(i, j)=(sunew_hat(i,j)-su_hat(i,j))/dt-rhs_hat(i,j)
!!           end do
!!         end do
!!         res = dnrm2((N/2+1)*N*2,rhs_hat,1)/dnrm2((N/2+1)*N*2,su_hat,1)
!!         write(22, *), 'time ', t, 'err_ind', res, &
!!           'err',  dnrm2((N/2+1)*N*2, rhs_hat,1), &
!!           'su', dnrm2((N/2+1)*N*2, su_hat,1)
!!
!!         !if ( res .ge. gamma .and. ts .le. 5) then
!!#if 0
!!         if ( res .ge. gamma ) then
!!           ts = ts + 1
!!           print*, 'ts =', ts, 'res', res, 'time', t
!!           su_hat = sunew_hat
!!           allocate(UP_hat(1:(N/2 + 1)*N, 1:initer/nc))
!!           allocate(temp_hat1(1:N*N, 1:initer/nc))
!!           write(22, *), 'time ', t, 'u', dnrm2((N/2+1)*N*2, u_hat,1),&
!!             'rel_diff', dnrm2((N/2+1)*N*2, su_hat-u_hat,1)/&
!!             dnrm2((N/2+1)*N*2, u_hat,1)
!!           call zcopy((N/2 + 1)*N, su_hat, 1, convs_hat, 1)
!!           exiter = exiter + 1
!!
!!           do niter = 1, initer
!!             cost = cos(t)
!!             t = t + dt
!!             !! u
!!             !call sparse_ComputePlane(pi, cost, dt, N, C, Amp, lam,  theta, eps, u_hat, ux_hat, uy_hat, rhs_hat)
!!
!!             ! su
!!             call sparse_ComputePlane(pi, cost, dt, N, C, Amp, lam,  theta, eps, su_hat, ux_hat, uy_hat, rhs_hat)
!!
!!             ! Copy
!!             if(mod(niter, nc) .eq. 0) then
!!               call zcopy((N/2 + 1)*N, su_hat, 1, UP_hat(:, niter/nc), 1)
!!             end if
!!
!!             if(mod(exiter, 100) .eq. 0) then
!!               call zcopy((N/2 + 1)*N, su_hat, 1, convs_hat, 1)
!!             else
!!               call zaxpy((N/2 + 1)*N, dcmplx(1.0D+0,0.0D+0), su_hat, 1, &
!!                          convs_hat, 1)
!!             end if
!!
!!             !if(mod(exiter, 100) .eq. 0) then
!!             !  write(22, *), 'time ', t, 'u', dnrm2((N/2+1)*N*2, u_hat,1),&
!!             !    'rel_diff', dnrm2((N/2+1)*N*2, su_hat-u_hat,1)/&
!!             !    dnrm2((N/2+1)*N*2, u_hat,1)
!!             !end if
!!             if(mod(exiter, 100) .eq. 99) then
!!               limsold = lims
!!               lims = dlog(real(convs_hat(1, 1))/(100*N*N))/(t - 50*dt)
!!               y = (cM + dlog(real(su_hat(1, 1))/(N*N))/t)/lam
!!
!!               write(100, *) 'time', t, real(convs_hat(1, 1))/(100*N*N), &
!!                             'lims', lims, &
!!                             'diffs', abs(lims - limsold)/(100*dt), 'Cs*', y
!!               call dfftw_execute_dft_c2r(p2, convs_hat, convs)
!!               write(100, *) dasum(N*N, convs, 1)*h*h/(100*N*N)
!!
!!               if(abs(lims - limsold)/(100*dt) .lt. ineps) then
!!                 exflg = 1
!!                 exit
!!               end if
!!             end if
!!
!!             exiter = exiter + 1
!!           end do
!!
!!           if (exflg .eq. 0) then
!!
!!             deallocate(BNlap)
!!             deallocate(BNx)
!!             deallocate(BNy)
!!             deallocate(RBNx)
!!             deallocate(RBNy)
!!             deallocate(RBBNx)
!!             deallocate(RBBNy)
!!             deallocate(BBNx)
!!             deallocate(BBNy)
!!             deallocate(a_hat)
!!             deallocate(RBN)
!!             deallocate(RBBN)
!!             deallocate(BBN)
!!             deallocate(BNt)
!!             deallocate(Fixm)
!!             deallocate(Varm)
!!             deallocate(Varm1)
!!             deallocate(at_hat)
!!             deallocate(BBNyt)
!!             deallocate(anew_hat)
!!             deallocate(suma_hat)
!!             call sparse_CopyMatrix_Symmetrize(m, N, temp_hat1, UP_hat)
!!
!!               temp_hat1(1,k) = dcmplx(dble(temp_hat1(1,k)), 0.0D+0)
!!               temp_hat1(1+N/2*N,k)=dcmplx(dble(temp_hat1(1+N/2*N,k)),&
!!                 0.0D+0)
!!               temp_hat1(1+N/2,k)=dcmplx(dble(temp_hat1(1+N/2,k)),&
!!                 0.0D+0)
!!               temp_hat1(1+N/2+N/2*N,k)= & 
!!                 dcmplx(dble(temp_hat1(1+N/2+N/2*N,k)), 0.0D+0)
!!             end do
!!
!!
!!             deallocate(UP_hat)
!!             allocate(B(1:N*N, 1:initer/nc))
!!             allocate(CT(1:initer/nc, 1:initer/nc))
!!             allocate(rwork(1:5*min(N*N, initer/nc)))
!!             S(1:N*N) = 0.0D+0
!!             call system_clock(count=svdstime)
!!             lwork = -1
!!             call zgesvd('S', 'N', N*N, initer/nc, temp_hat1, N*N,&
!!               S, B, N*N, CT, initer/nc, work_temp, lwork, rwork, info)
!!
!!             ! Compute the SVD
!!             lwork = int(work_temp)
!!             allocate(work(1:lwork))
!!             call zgesvd('S', 'N', N*N, initer/nc, temp_hat1, N*N,&
!!               S, B, N*N, CT, initer/nc, work, lwork, rwork, info)
!!             deallocate(work)
!!             deallocate(rwork)
!!             deallocate(CT)
!!             deallocate(temp_hat1)
!!
!!             call system_clock(count=svdetime)
!!             to_etime = real(svdetime - svdstime)/real(clock_rate)
!!             print*, 'Elasped time for SVD ', to_etime, ' seconds.'
!!
!!             sums = sum(S(1:initer/nc))
!!             sacc = 0.d0
!!             do i = 1, initer/nc
!!               sacc = sacc + S(i)
!!               if( sacc .gt. pc2*sums) exit
!!             end do
!!             m2 = min(i, initer/nc)
!!             !m2 = max(m2, 20)
!!            ! Copy 
!!            do j = 1, N
!!              do i = 1, N/2 + 1
!!                u_hat1(i+(j-1)*N) = su_hat(i,j)
!!              end do
!!            end do
!!            ! Symmetrize
!!            do j = 2, N
!!              do i = N/2 + 2, N
!!                u_hat1(i+(j-1)*N)=conjg(su_hat(N+2-i, N+2-j))
!!              end do
!!            end do
!!            do i = N/2 + 2, N
!!              u_hat1(i)=conjg(su_hat(N+2-i,1))
!!            end do
!!
!!             allocate(a_hat(1:m2))
!!             call zgemv('C', N*N, m2, dcmplx(1.0D+0,0.0D+0), B, N*N,&
!!                        u_hat1, 1, dcmplx(0.0D+0,0.0D+0), a_hat, 1)
!!
!!             call zgemv('N', N*N, m2, dcmplx(1.0D+0,0.0D+0), B, N*N,&
!!                        a_hat, 1, dcmplx(0.0D+0,0.0D+0), su_hat1, 1)
!!
!!             print*, 'difference', dnrm2(N*N*2, su_hat1-u_hat1, 1)
!!             deallocate(a_hat)
!!
!!             ! Copy
!!             allocate(BN2(1:N*N, 1:(m + m2)))
!!             call  sparse_CopyMatrix_Symmetrize(m, N, BN2, BN)
!!             do k = m + 1, m + m2
!!               do j = 1, N
!!                 do i = 1, N/2 + 1
!!                   BN2(i+(j-1)*N,k) = B(i+(j-1)*N,k-m)
!!                 end do
!!               end do
!!               ! Symmetrize
!!               do j = 2, N
!!                 do i = N/2 + 2, N
!!                   BN2(i+(j-1)*N,k)=conjg(B(N+2-i+(N+1-j)*N,k-m))
!!                 end do
!!               end do
!!               do i = N/2 + 2, N
!!                 BN2(i,k)=conjg(B(N+2-i,k-m))
!!               end do
!!             end do
!!             deallocate(B)
!!             m = m + m2
!!
!!             print*, 'm =', m, ' m2 =', m2
!!             allocate(tauq(1:m))
!!             lwork = -1
!!             call zgeqrf(N*N, m, BN2, N*N, tauq, work_temp, &
!!                         lwork, info)
!!             lwork = int(work_temp)
!!             allocate(work(1:lwork))
!!             call zgeqrf(N*N, m, BN2, N*N, tauq, work, &
!!                         lwork, info)
!!             deallocate(work)
!!
!!             ! mtr is the R matrix of the QR factorization of BN2
!!             allocate(mtr(1:m, 1:m))
!!             call zlacpy('U', m, m, BN2, N*N, mtr, m)
!!             call zlaset('L', m-1, m-1, CZERO, CZERO, mtr(2, 1), m)
!!
!!             ! Form the Q factor in BN2
!!             lwork = -1
!!             call zungqr(N*N, m, m, BN2, N*N, tauq, work_temp, &
!!                         lwork, info)
!!             lwork = int(work_temp)
!!             allocate(work(1:lwork))
!!             call zungqr(N*N, m, m, BN2, N*N, tauq, work, &
!!                         lwork, info)
!!             deallocate(work)
!!             deallocate(tauq)
!!
!!             ! Perform SVD factorization for the real mtr
!!             allocate(rmtr(1:m, 1:m))
!!             do j = 1, m
!!               do i = 1, m
!!                 rmtr(i, j) = dble(mtr(i, j))
!!               end do
!!             end do
!!             allocate(rorth(1:m, 1:m))
!!             allocate(rct(1:m, 1:m))
!!             allocate(S2(1:m))
!!             S2(1:m) = 0.0D+0
!!             lwork = -1
!!             call dgesvd('S', 'N', m, m, rmtr, m,&
!!               S2, rorth, m, rct, m, work_temp, lwork, info)
!!
!!             ! Compute the SVD
!!             lwork = int(work_temp)
!!             allocate(work(1:lwork))
!!             call dgesvd('S', 'N', m, m, rmtr, m,&
!!               S2, rorth, m, rct, m, work, lwork, info)
!!             deallocate(work)
!!             deallocate(rct)
!!             deallocate(rmtr)
!!
!!             do j = 1, m
!!               do i = 1, m
!!                 mtr(i, j) = dcmplx(rorth(i, j), 0.0D+0)
!!               end do
!!             end do
!!             deallocate(rorth)
!!
!!             allocate(B2(1:N*N, 1:m))
!!             call zgemm('N', 'N', N*N, m, m, CONE, BN2, N*N, mtr, &
!!                        m, CZERO, B2, N*N)
!!             deallocate(mtr)
!!             sums = sum(S2(1:m))
!!             sacc = 0.d0
!!             do i = 1, m
!!               sacc = sacc + S2(i)
!!               if( sacc .gt. pc3*sums) exit
!!             end do
!!             m = min(i, m)
!!             print*, "m =", m
!!             deallocate(S2)
!!
!!             allocate(Fixm(1:m, 1:m))
!!             call zgemm('C', 'N', m, m, N*N, dcmplx(1.0D+0,0.0D+0),B2,&
!!                        N*N, B2, N*N, dcmplx(0.0D+0,0.0D+0), &
!!                        Fixm, m)
!!             deallocate(Fixm)
!!
!!             deallocate(BN1)
!!             deallocate(BN)
!!             allocate(BN1(1:N*N, 1:m))
!!             allocate(BN(1:(N/2+1)*N, 1:m))
!!
!!             call dcopy(2*N*N*m, B2, 1, BN1, 1)
!!             do k = 1, m
!!               do j = 1, N
!!                 do i = 1, N/2 + 1
!!                   BN(i+(j-1)*(N/2+1),k)=BN1(i+(j-1)*N,k)
!!                 end do
!!               end do
!!             end do
!!             deallocate(B2)
!!
!!             allocate(BNlap(1:(N/2+1)*N, 1:m))
!!             allocate(BNx(1:(N/2+1)*N, 1:m))
!!             allocate(BNy(1:(N/2+1)*N, 1:m))
!!             allocate(RBNx(1:N, 1:N, 1:m))
!!             allocate(RBNy(1:N, 1:N, 1:m))
!!             allocate(RBBNx(1:N, 1:N, 1:m))
!!             allocate(RBBNy(1:N, 1:N, 1:m))
!!             allocate(BBNx(1:(N/2+1)*N, 1:m))
!!             allocate(BBNy(1:(N/2+1)*N, 1:m))
!!             allocate(RBN(1:N, 1:N, 1:m))
!!             allocate(RBBN(1:N, 1:N, 1:m))
!!             allocate(BBN(1:(N/2+1)*N, 1:m))
!!             allocate(BNt(1:(N/2+1)*N, 1:m))
!!             allocate(a_hat(1:m))
!!             allocate(Fixm(1:m, 1:m))
!!             allocate(Varm(1:m, 1:m))
!!             allocate(Varm1(1:m, 1:m))
!!             allocate(at_hat(1:m))
!!             allocate(BBNyt(1:(N/2+1)*N, 1:m))
!!             allocate(anew_hat(1:m))
!!             allocate(suma_hat(1:m))
!!             
!!          !   call sparse_CopyMatrix_Symmetrize(m, N, u_hat1, su_hat)
!!            ! Copy 
!!            do j = 1, N
!!              do i = 1, N/2 + 1
!!                u_hat1(i+(j-1)*N) = su_hat(i,j)
!!              end do
!!            end do
!!            ! Symmetrize
!!            do j = 2, N
!!              do i = N/2 + 2, N
!!                u_hat1(i+(j-1)*N)=conjg(su_hat(N+2-i, N+2-j))
!!              end do
!!            end do
!!            do i = N/2 + 2, N
!!              u_hat1(i)=conjg(su_hat(N+2-i,1))
!!            end do
!!
!!             call zgemv('C', N*N, m, dcmplx(1.0D+0,0.0D+0), BN1, N*N,&
!!                        u_hat1, 1, dcmplx(0.0D+0,0.0D+0), a_hat, 1)
!!
!!             call zgemv('N', N*N, m, dcmplx(1.0D+0,0.0D+0), BN1, N*N,&
!!                        a_hat, 1, dcmplx(0.0D+0,0.0D+0), su_hat1, 1)
!!
!!             print*, 'difference', dnrm2(N*N*2, su_hat1-u_hat1, 1)
!!
!!             ! Form the matrices
!!             ! Fixm = eps*BN'BNlap + 2*lam*eps*BN'BNx + C*I
!!
!!             ! BNlap(1:(N/2+1)*N, 1:m) is the laplacian of BN
!!             call sparse_BuildBNlap(m, N, pi, BNlap, BN)
!!
!!            call sparse_CopyMatrix_Symmetrize(m , N, BN2, BNlap)
!!
!!             call zgemm('C', 'N', m, m, N*N, dcmplx(1.0D+0,0.0D+0), BN1,&
!!                        N*N, BN2, N*N, dcmplx(0.0D+0,0.0D+0), &
!!                        Fixm, m)
!!
!!             ! BNx(1:(N/2+1)*N, 1:m) is the derivatives of BN to x
!!             call sparse_BuildBNx(m, N, pi, BNx, BN)
!!
!!             call sparse_CopyMatrix_Symmetrize(m, N, BN2, BNx)
!!
!!             call zgemm('C', 'N', m, m, N*N, dcmplx(2*eps*lam,0.0D+0), &
!!                        BN1, N*N, BN2, N*N, dcmplx(eps,0.0D+0), &
!!                        Fixm, m)
!!
!!             do i = 1, m
!!               Fixm(i, i) = Fixm(i, i) + C
!!             end do
!!
!!             ! BNy(1:(N/2+1)*N, 1:m) is the derivatives of BN to y
!!            call sparse_BuildBNy(m, N, pi, BNy, BN)
!!
!!             ! RBNx(1:N, 1:N, 1:m) is the derivatives of BN to x in real space
!!             call zcopy((N/2 + 1)*N*m, BNx, 1, BNt, 1)
!!             do k = 1, m
!!               call dfftw_execute_dft_c2r(p2, BNt(:,k), RBNx(:,:,k))
!!               call dscal(N*N, 1.0D+0/(N*N), RBNx(:,:,k), 1)
!!             end do
!!
!!             ! RBNy(1:N, 1:N, 1:m) is the derivatives of BN to y in real space
!!             do k = 1, m
!!               call dfftw_execute_dft_c2r(p2, BNy(:,k), RBNy(:,:,k))
!!               call dscal(N*N, 1.0D+0/(N*N), RBNy(:,:,k), 1)
!!             end do
!!
!!            call sparse_BuildRBBNxRBBNy2(m, N, sinx, Amp, theta, Bx, By, RBNx, RBNy,&
!!            RBBNx, RBBNy)
!!
!!             ! BBNx(1:(N/2+1)*N, 1:m) is RBBNx without time t in k space
!!             do k = 1, m
!!               call dfftw_execute_dft_r2c(p1, RBBNx(:,:,k), BBNx(:,k))
!!             end do
!!
!!             ! BBNy(1:(N/2+1)*N, 1:m) is RBBNy without time t in k space
!!             do k = 1, m
!!               call dfftw_execute_dft_r2c(p1, RBBNy(:,:,k), BBNy(:,k))
!!             end do
!!
!!             ! BBNy = BBNy + BBNx
!!             call zaxpy((N/2+1)*N*m, dcmplx(1.0D+0,0.0D+0), BBNx, 1, BBNy, 1)
!!
!!             ! RBN(1:N, 1:N, 1:m) is BN in real space
!!             call zcopy((N/2 + 1)*N*m, BN, 1, BNt, 1)
!!             do k = 1, m
!!               call dfftw_execute_dft_c2r(p2, BNt(:,k), RBN(:,:,k))
!!               call dscal(N*N, 1.0D+0/(N*N), RBN(:,:,k), 1)
!!             end do
!!
!!             ! RBBN(1:N, 1:N, 1:m) is Bx*RBN in real space without time t
!!             do k = 1, m
!!               forall(j = 1:N)
!!                 forall(i = 1:N)
!!                   RBBN(i, j, k) = Bx(i, j)*RBN(i, j, k)
!!                 end forall
!!               end forall
!!             end do
!!
!!             ! BBN(1:(N/2+1)*N, 1:m) is RBBN in k space
!!             do k = 1, m
!!               call dfftw_execute_dft_r2c(p1, RBBN(:,:,k), BBN(:,k))
!!             end do
!!
!!             ! BBNy = BBNy + lam*BBN
!!             call zaxpy((N/2+1)*N*m, dcmplx(lam,0.0D+0), BBN, 1, BBNy, 1)
!! 
!!             call sparse_CopyMatrix_Symmetrize(m, N, BN2, BBNy)
!!
!!             call zgemm('C', 'N', m, m, N*N, dcmplx(1.0D+0,0.0D+0), &
!!                        BN1,N*N, BN2, N*N, dcmplx(1.0D+0,0.0D+0), &
!!                        Fixm, m)
!!
!!
!!            call sparse_BuildRBBNxRBBNy2(m, N, sinx, Amp, theta, Bx, By, RBNx, RBNy,&
!!            RBBNx, RBBNy)
!!             do k = 1, m
!!               call dfftw_execute_dft_r2c(p1, RBBNx(:,:,k), BBNx(:,k))
!!             end do 
!!             
!!             ! BBNy(1:(N/2+1)*N, 1:m) is RBBNy with time t in k space
!!             do k = 1, m
!!               !call dfftw_execute_dft_r2c(p1, RBBNy(:,:,k), BBNy(:,k))
!!               call dfftw_execute_dft_r2c(p1, RBBNy(:,:,k), BBNyt(:,k))
!!             end do 
!!             
!!             ! BBNy = BBNy + BBNx
!!             !call zaxpy((N/2+1)*N*m, dcmplx(1.0D+0,0.0D+0), BBNx, 1, BBNy, 1)
!!             call zaxpy((N/2+1)*N*m, dcmplx(1.0D+0,0.0D+0), BBNx, 1, BBNyt, 1)
!!             
!!             ! RBBN(1:N, 1:N, 1:m) is Bx*RBN in real space with time t
!!             do k = 1, m
!!               forall(j = 1:N)
!!                 forall(i = 1:N)
!!                   RBBN(i, j, k) = Bx(i, j)*RBN(i, j, k)
!!                 end forall
!!               end forall
!!             end do
!!
!!             ! BBN(1:(N/2+1)*N, 1:m) is RBBN in k space
!!             do k = 1, m
!!               call dfftw_execute_dft_r2c(p1, RBBN(:,:,k), BBN(:,k))
!!             end do
!!
!!             ! BBNy = BBNy + lam*BBN
!!             !call zaxpy((N/2+1)*N*m, dcmplx(lam,0.0D+0), BBN, 1, BBNy, 1)
!!             call zaxpy((N/2+1)*N*m, dcmplx(lam,0.0D+0), BBN, 1, BBNyt, 1)
!!                
!!             call sparse_CopyMatrix_Symmetrize(m, N, BN2, BBNyt)
!!
!!             call zgemm('C', 'N', m, m, N*N, dcmplx(1.0D+0,0.0D+0), &
!!                        BN1, N*N, BN2, N*N, dcmplx(0.0D+0,0.0D+0), &
!!                        Varm1, m)
!!
!!             do i = 1, m
!!               anew_hat(i) = a_hat(i)
!!             end do
!!             deallocate(BN2)
!!           end if
!!         end if
!!#endif
!!       end if
!!
!!       do i = 1, m
!!         a_hat(i) = anew_hat(i)
!!       end do
!!
!!       if (exflg .eq. 1) exit
!!       !! The error
!!       !if(mod(exiter, 100) .eq. 0) then
!!       !  write(22, *), 'time ', t, 'u', dnrm2((N/2+1)*N*2, u_hat,1),&
!!       !    'rel_diff', &
!!       ! dnrm2((N/2+1)*N*2, su_hat-u_hat,1)/dnrm2((N/2+1)*N*2, u_hat,1)
!!       !end if
!!
!!       if(mod(exiter, 100) .eq. 0) then
!!         call zcopy(m, a_hat, 1, suma_hat, 1)
!!       else
!!         call zaxpy(m, dcmplx(1.0D+0,0.0D+0), a_hat, 1, suma_hat, 1)
!!       end if
!!
!!       if(mod(exiter, 100) .eq. 99) then
!!         call zgemv('N', (N/2+1)*N, m, dcmplx(1.0D+0,0.0D+0), BN, &
!!           (N/2+1)*N, suma_hat, 1, dcmplx(0.0D+0,0.0D+0), convs_hat, 1)
!!         call zgemv('N', (N/2+1)*N, m, dcmplx(1.0D+0,0.0D+0), BN, &
!!           (N/2+1)*N, a_hat, 1, dcmplx(0.0D+0,0.0D+0), su_hat, 1)
!!
!!         limsold = lims
!!         lims = dlog(real(convs_hat(1, 1))/(100*N*N))/(t - 50*dt)
!!         y = (cM + dlog(real(su_hat(1, 1))/(N*N))/t)/lam
!!
!!         write(100, *) 'time', t, real(convs_hat(1, 1))/(100*N*N), &
!!                       'lims', lims, &
!!                       'diffs', abs(lims - limsold)/(100*dt), 'Cs*', y
!!         call dfftw_execute_dft_c2r(p2, convs_hat, convs)
!!         write(100, *) dasum(N*N, convs, 1)*h*h/(100*N*N)
!!
!!         if(abs(lims - limsold)/(100*dt) .lt. ineps) exit
!!       end if
!!
!!       exiter = exiter + 1
!!
!!       !limold = lim
!!       !lim = dlog(real(u_hat(1, 1))/(N*N))/t
!!       !limsold = lims
!!       !lims = dlog(real(su_hat(1, 1))/(N*N))/t
!!       !do k = 1, kmax
!!       !  if((t .ge. k*1.0e-4*t2) .and. (t .lt. k*1.0e-4*t2 + dt)) then
!!       !    write(100, *) 'time', t, real(u_hat(1, 1))/(N*N), &
!!       !                  lim, 'rel_diff', (lim - limold)/abs(limold),&
!!       !                  'C*', (cM + lim)/lam
!!       !    !if(real(u_hat(1, 1)) .le. 0.D+0) then
!!       !      ut_hat(1:N/2 + 1, 1:N) = u_hat(1:N/2 + 1, 1:N)
!!       !      call dfftw_execute_dft_c2r(p2, ut_hat, u)
!!       !      write(100, *) dasum(N*N, u, 1)*h*h/(N*N)
!!       !    write(100, *) 'time', t, real(su_hat(1, 1))/(N*N), &
!!       !                  lims,'rels_diff',(lims-limsold)/abs(limsold),&
!!       !                  'Cs*', (cM + lims)/lam
!!       !    !end if
!!       !  end if
!!       !end do
!!
!!       !write(100, *) 'norm u_hat', dznrm2((N/2+1)*N, u_hat, 1), &
!!       !            'norm su_hat', dznrm2((N/2+1)*N, su_hat, 1), &
!!       !            'error', dznrm2((N/2+1)*N, su_hat-u_hat, 1), &
!!       ! 'relerror',dznrm2((N/2+1)*N, su_hat-u_hat, 1)/dznrm2((N/2+1)*N,&
!!       !            u_hat, 1)
!!     end do
!!
!!     !write(100, *) 'norm u_hat', dznrm2((N/2+1)*N, u_hat, 1), &
!!     !              'norm su_hat', dznrm2((N/2+1)*N, su_hat, 1), &
!!     !              dznrm2((N/2+1)*N, su_hat-u_hat, 1), &
!!     !              dznrm2((N/2+1)*N, su_hat-u_hat, 1)/dznrm2((N/2+1)*N,&
!!     !              u_hat, 1)
!!
!!     write(100, *) 'Number of basis', m
!!     write(100, *) 'time', t, real(u_hat(1, 1))/(N*N), &
!!                          dlog(real(u_hat(1, 1))/(N*N))/t
!!
!!     call dfftw_execute_dft_c2r(p2, u_hat, u)
!!     call dscal(N*N, 1.0D+0/(N*N), u, 1)
!!
!!     write(100, *) 'Final time t =', t
!!
!!     call dfftw_destroy_plan(p1)
!!     call dfftw_destroy_plan(p2)
!!
!!     call system_clock(count=toetime)
!!     to_etime = real(toetime - tostime)/real(clock_rate)
!!     print*, 'Total elasped time is ', to_etime, &
!!                               ' seconds.'
!!
!!     close(100)
!!     close(22)
!!     deallocate(u)
!!     deallocate(u_hat)
!!     deallocate(cosx)
!!     deallocate(sinx)
!!     deallocate(ux_hat)
!!     deallocate(uy_hat)
!!     deallocate(rhs)
!!     deallocate(rhs_hat)
!!     deallocate(ut_hat)
!!     deallocate(S)
!!     deallocate(BN)
!!     deallocate(BNlap)
!!     deallocate(BNx)
!!     deallocate(BNy)
!!     deallocate(a_hat)
!!     deallocate(su_hat)
!!     deallocate(Bx)
!!     deallocate(By)
!!     deallocate(RBNx)
!!     deallocate(RBNy)
!!     deallocate(RBBNx)
!!     deallocate(RBBNy)
!!     deallocate(BBNx)
!!     deallocate(BBNy)
!!     deallocate(RBN)
!!     deallocate(RBBN)
!!     deallocate(BBN)
!!     deallocate(BNt)
!!     deallocate(Fixm)
!!     deallocate(Varm)
!!     deallocate(Varm1)
!!     deallocate(at_hat)
!!     deallocate(conv)
!!     deallocate(conv_hat)
!!     deallocate(convs)
!!     deallocate(convs_hat)
!!     deallocate(u_hat1)
!!     deallocate(su_hat1)
!!     deallocate(BBNyt)
!!     deallocate(anew_hat)
!!     deallocate(sunew_hat)
!!     deallocate(BN1)
!!     deallocate(suma_hat)
   end
