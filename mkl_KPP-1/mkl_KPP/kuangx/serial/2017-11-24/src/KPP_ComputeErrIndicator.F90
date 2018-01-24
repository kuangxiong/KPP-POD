  subroutine KPP_ComputeErrIndicator(N, m, BN, BNx, BNlap, BBNy, BBNyt,&
       rhs_hat , a_hat, su_hat, anew_hat, sunew_hat, res)
    integer:: N, m
    real*8::res
    complex*16::BN(1:(N/2+1)*N, 1:m), BNlap(1:(N/2+1)*N, 1:m),&
               BNx(1:(N/2+1)*N, 1:m), BBNy(1:(N/2+1)*N, 1:m),&
               BBNyt(1:(N/2+1)*N, 1:m), rhs_hat(1:N/2+1, 1:N),&
               a_hat(1:m), su_hat(1:N/2+1, 1:N), anew_hat(1:m),&
               sunew_hat(1:N/2+1, 1:N)

    call zgemv('N', (N/2+1)*N, m, dcmplx(1.0D+0,0.0D+0), BN, &
           (N/2+1)*N, a_hat, 1, dcmplx(0.0D+0,0.0D+0), su_hat, 1)
    call zgemv('N', (N/2+1)*N, m, dcmplx(1.0D+0,0.0D+0), BN, &
           (N/2+1)*N, anew_hat, 1, dcmplx(0.0D+0,0.0D+0), sunew_hat, 1)

    call zgemv('N', (N/2+1)*N, m, dcmplx(eps,0.0D+0), BNlap, &
           (N/2+1)*N, a_hat, 1, dcmplx(0.0D+0,0.0D+0), rhs_hat, 1)
    call zgemv('N', (N/2+1)*N, m, dcmplx(2*eps*lam,0.0D+0), BNx, &
           (N/2+1)*N, a_hat, 1, dcmplx(1.0D+0,0.0D+0), rhs_hat, 1)
    call zgemv('N', (N/2+1)*N, m, dcmplx(1.0D+0,0.0D+0), BBNy, &
           (N/2+1)*N, a_hat, 1, dcmplx(1.0D+0,0.0D+0), rhs_hat, 1)
    call zgemv('N', (N/2+1)*N, m, dcmplx(cost,0.0D+0), BBNyt, &
           (N/2+1)*N, a_hat, 1, dcmplx(1.0D+0,0.0D+0), rhs_hat, 1)
    call zgemv('N', (N/2+1)*N, m, dcmplx(C,0.0D+0), BN, &
           (N/2+1)*N, a_hat, 1, dcmplx(1.0D+0,0.0D+0), rhs_hat, 1)

    do j = 1, N
        do i = 1, N/2 + 1
           rhs_hat(i, j)=(sunew_hat(i,j)-su_hat(i,j))/dt-rhs_hat(i,j)
        end do
    end do
         res = dnrm2((N/2+1)*N*2,rhs_hat,1)/dnrm2((N/2+1)*N*2,su_hat,1)
      end 
