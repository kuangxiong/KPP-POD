subroutine KPP_ComputePlane(pi, cost, dt, N, C, Amp, lam, theta, eps, u_hat,&
    ux_hat, uy_hat,  rhs_hat)
    integer :: N
    real*8 :: pi,  cost, dt, C,  Amp, lam, theta, eps
    complex*16:: u_hat(1:N/2+1, 1:N), ux_hat(1:N/2+1, 1:N), uy_hat(1:N/2+1, 1:N), rhs_hat(1:N/2+1, 1:N)
    integer :: i, j

    do j = 2, N-1
      do i = 1, N/2
        ux_hat(i,j)= Amp*((pi*(i-1)*dcmplx(theta*cost, 1.0D+0) + lam&
                     /2*dcmplx(1.0D+0,-theta*cost))*u_hat(i, j-1)&
                     + (pi*(i - 1)*dcmplx(-theta*cost, 1.0D+0)&
                     +lam/2*dcmplx(1.0D+0, theta * cost)) * u_hat(i, j+1))
      end do
    end do

    do i = 1, N/2
      ux_hat(i,1)=Amp*((pi*(i - 1)*dcmplx(theta*cost, 1.0D+0) + lam&
                  /2*dcmplx(1.0D+0,-theta*cost))*u_hat(i, N) &
                  +(pi*(i - 1)*dcmplx(-theta*cost, 1.0D+0)&
                  +lam/2*dcmplx(1.0D+0,theta*cost))*u_hat(i, 2))
    end do

    do i = 1, N/2
      ux_hat(i,N)=Amp*((pi*(i - 1)*dcmplx(theta*cost, 1.0D+0) + lam&
                  /2*dcmplx(1.0D+0,-theta*cost))*u_hat(i, N-1) &
                  +(pi*(i - 1)*dcmplx(-theta*cost, 1.0D+0)&
                  +lam/2*dcmplx(1.0D+0,theta*cost))*u_hat(i, 1))
    end do
    i = N/2 + 1
    do j = 2, N-1
      ux_hat(i, j) = Amp*lam/2*(dcmplx(1.0D+0,-theta*cost)*&
                           u_hat(i, j-1)+dcmplx(1.0D+0,theta*cost)&
                           *u_hat(i,j+1))
    end do
    ux_hat(i, 1) = Amp*lam/2*(dcmplx(1.0D+0,-theta*cost)*&
                         u_hat(i, N)+dcmplx(1.0D+0,theta*cost)&
                         *u_hat(i,2))
    ux_hat(i, N) = Amp*lam/2*(dcmplx(1.0D+0,-theta*cost)*&
                         u_hat(i, N-1)+dcmplx(1.0D+0,theta*cost)&
                         *u_hat(i,1))

    do j = 1, N/2
      do i = 2, N/2
        uy_hat(i,j)=Amp*pi*(j - 1)*(dcmplx(theta*cost, 1.0D+0)*&
                    u_hat(i - 1, j) + dcmplx(-theta*cost, 1.0D+0)&
                    *u_hat(i + 1, j))
      end do
    end do

    do j = N/2 + 2, N
      do i = 2, N/2
        uy_hat(i, j)=Amp*pi*(j - 1 - N)*(dcmplx(theta*cost, 1.0D+0)*&
                       u_hat(i - 1, j) + dcmplx(-theta*cost, 1.0D+0)&
                       *u_hat(i + 1, j))
      end do
    end do

    uy_hat(1, 1) = dcmplx(0.0D+0, 0.0D+0)
    do j = 2, N/2
      uy_hat(1, j) = Amp*pi*(j - 1)*(dcmplx(theta*cost, 1.0D+0)*&
               conjg(u_hat(2,N-j+2)) + dcmplx(-theta*cost, 1.0D+0)&
                     *u_hat(2, j))
    end do

    do j = N/2 + 2, N
      uy_hat(1,j)=Amp*pi*(j - 1 - N)*(dcmplx(theta*cost, 1.0D+0)*&
               conjg(u_hat(2,N-j+2)) + dcmplx(-theta*cost, 1.0D+0)&
                     *u_hat(2, j))
    end do

    uy_hat(N/2+1, 1) = dcmplx(0.0D+0, 0.0D+0)
    do j = 2, N/2
      uy_hat(N/2+1,j)=Amp*pi*(j - 1)*(dcmplx(theta*cost, 1.0D+0)*&
                        u_hat(N/2, j) + dcmplx(-theta*cost, 1.0D+0)&
                        *conjg(u_hat(N/2, N - j + 2)))
    end do

    do j = N/2 + 2, N
      uy_hat(N/2+1,j)=Amp*pi*(j - 1 - N)*(dcmplx(theta*cost, 1.0D+0)*&
                        u_hat(N/2, j) + dcmplx(-theta*cost, 1.0D+0)&
                        *conjg(u_hat(N/2, N - j + 2)))
    end do

    do i = 1, N/2 + 1
      uy_hat(i, N/2 + 1) = dcmplx(0.0D+0, 0.0D+0)
    end do

    do j = 1, N
      do i = 1, N/2
        rhs_hat(i, j) = ux_hat(i, j) + uy_hat(i, j) + &
                      dcmplx(C, 4*pi*lam*eps*(i-1))*u_hat(i, j)
 !       rhs_hat(i, j) =  uy_hat(i, j) 
      end do
    end do
    i = N/2 + 1
    do j = 1, N
      rhs_hat(i, j) = ux_hat(i, j) + uy_hat(i, j) + &
                    dcmplx(C, 0.0D+0)*u_hat(i, j)
!      rhs_hat(i, j) =  uy_hat(i, j) 
    end do

    forall(j = 1:N/2 + 1)
      forall(i = 1:N/2 + 1)
        u_hat(i, j) = 1.0D+0/(4*pi*pi*((i-1)*(i-1)+(j-1)*(j-1))&
                      *eps*dt+1)*(u_hat(i,j) + dt*rhs_hat(i,j))
      end forall
    end forall

    forall(j = N/2 + 2:N)
      forall(i = 1:N/2 + 1)
        u_hat(i,j)= 1.0D+0/(4*pi*pi*((i-1)*(i-1)+(j-1-N)*(j-1-N))&
                    *eps*dt+1)*(u_hat(i,j) + dt*rhs_hat(i,j))
      end forall
    end forall
    end 
      
