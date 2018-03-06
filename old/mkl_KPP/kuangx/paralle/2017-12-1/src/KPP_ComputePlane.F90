subroutine KPP_ComputePlane(pi, cost, dt, N, C, Amp, lam, theta, eps, u_hat,&
    rhs_hat)
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
      end do
    end do
    i = N/2 + 1
    do j = 1, N
      rhs_hat(i, j) = ux_hat(i, j) + uy_hat(i, j) + &
                    dcmplx(C, 0.0D+0)*u_hat(i, j)
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
    
subroutine KPP_mpi_ComputePlane(pi, cost, dt, N, C, Amp, lam, theta, eps, u_hat,&
    rhs_hat,localN, local_0_start, ierr, myrank, nprocs)
    include'mpif.h'
    integer :: N, nprocs
    real*8 :: pi,  cost, dt, C,  Amp, lam, theta, eps
    complex*16 :: u_hat(1:N/2+1, 1:localN), ux_hat(1:N/2+1, 1:localN), uy_hat(1:N/2+1, 1:localN), rhs_hat(1:N/2+1, 1:localN)
    integer :: i, j, indexj
    complex*16 :: leftu_hat(1:localN), rightu_hat(1:localN), gatherLu_hat(1:N), gatherRu_hat(1:N) 
    complex*16, allocatable ::  tmpv_up(:),tmpv_dn(:),tmpv(:)
    integer :: ierr, status(MPI_STATUS_SIZE), myrank, preRank, nextRank
    do i=1, localN
        leftu_hat(i) = conjg(u_hat(2, i))
        rightu_hat(i) = conjg(u_hat(N/2, i))
    end do
    call MPI_Allgather(leftu_hat, localN, MPI_COMPLEX16, gatherLu_hat, localN, MPI_COMPLEX16, MPI_COMM_WORLD, ierr)
    call MPI_Allgather(rightu_hat, localN, MPI_COMPLEX16, gatherRu_hat, localN, MPI_COMPLEX16,  MPI_COMM_WORLD, ierr)
     
    preRank=mod(myrank-1+nprocs, nprocs)
    nextRank=mod(myrank+1+nprocs, nprocs)
    
    allocate(tmpv_up(1:N/2+1))
    allocate(tmpv_dn(1:N/2+1))
    allocate(tmpv(1:N/2+1))

!    call MPI_Barrier(MPI_COMM_WORLD)
    do i=1,N/2+1
       tmpv_up(i) = u_hat(i, 1)
       tmpv_dn(i) = u_hat(i, localN)
       tmpv(i) = u_hat(i, localN)
    end do

    call MPI_SendRecv(tmpv_up, N/2+1, MPI_COMPLEX16, preRank, 99, tmpv_dn, N/2+1, MPI_COMPLEX16, nextRank,&
                99, MPI_COMM_WORLD, status, ierr)
    call MPI_SendRecv(tmpv, N/2+1, MPI_COMPLEX16, nextRank, 991, tmpv_up, N/2+1, MPI_COMPLEX16, preRank, &
                991, MPI_COMM_WORLD, status, ierr)
!    
!    write(*,*)'i am here111'
!
    do j = 2, localN-1
      do i = 1, N/2
        ux_hat(i,j)= Amp*((pi*(i-1)*dcmplx(theta*cost, 1.0D+0) + lam&
                     /2*dcmplx(1.0D+0,-theta*cost))*u_hat(i, j-1)&
                     + (pi*(i - 1)*dcmplx(-theta*cost, 1.0D+0)&
                     +lam/2*dcmplx(1.0D+0, theta * cost)) * u_hat(i, j+1))
      end do
    end do

    do i = 1, N/2
      ux_hat(i,1)=Amp*((pi*(i - 1)*dcmplx(theta*cost, 1.0D+0) + lam&
                  /2*dcmplx(1.0D+0,-theta*cost))*tmpv_up(i) &
                  +(pi*(i - 1)*dcmplx(-theta*cost, 1.0D+0)&
                  +lam/2*dcmplx(1.0D+0,theta*cost))*u_hat(i, 2))
    end do

    do i = 1, N/2
      ux_hat(i,localN)=Amp*((pi*(i - 1)*dcmplx(theta*cost, 1.0D+0) + lam&
                  /2*dcmplx(1.0D+0,-theta*cost))*u_hat(i, localN-1) &
                  +(pi*(i - 1)*dcmplx(-theta*cost, 1.0D+0)&
                  +lam/2*dcmplx(1.0D+0,theta*cost))*tmpv_dn(i))
    end do

    i = N/2 + 1
    do j = 2, localN-1
      ux_hat(i, j) = Amp*lam/2*(dcmplx(1.0D+0, -theta*cost)*&
                           u_hat(i, j-1)+dcmplx(1.0D+0, theta*cost)&
                           *u_hat(i,j+1))
    end do
    ux_hat(i, 1) = Amp*lam/2*(dcmplx(1.0D+0,-theta*cost)*&
                         tmpv_up(i)+dcmplx(1.0D+0, theta*cost)&
                         *u_hat(i,2))
    ux_hat(i, localN) = Amp*lam/2*(dcmplx(1.0D+0,-theta*cost)*&
                         u_hat(i, localN-1)+dcmplx(1.0D+0, theta*cost)&
                         *tmpv_dn(i))


    do j = 1, localN
      indexj = j+ local_0_start
      do i = 2, N/2
        if(indexj .lt. N/2+1) then
            uy_hat(i,j)=Amp*pi*(j - 1)*(dcmplx(theta*cost, 1.0D+0)*&
                        u_hat(i - 1, j) + dcmplx(-theta*cost, 1.0D+0)&
                        *u_hat(i + 1, j))

            uy_hat(1, j) = Amp*pi*(j - 1)*(dcmplx(theta*cost, 1.0D+0)*&
                         gatherLu_hat(N-indexj+2) + dcmplx(-theta*cost, 1.0D+0)&
                         *u_hat(2, j))

            uy_hat(N/2+1,j)=Amp*pi*(j - 1)*(dcmplx(theta*cost, 1.0D+0)*&
                          u_hat(N/2, j) + dcmplx(-theta*cost, 1.0D+0)&
                          *gatherRu_hat(N - indexj + 2))

        else if(indexj .eq.  N/2+1) then
            uy_hat(i, j) = dcmplx(0.0D+0, 0.0D+0)
            uy_hat(1, j) = dcmplx(0.0D+0, 0.0D+0)
            uy_hat(N/2+1, j) = dcmplx(0.0D+0, 0.0D+0)

        else 
            uy_hat(i, j)=Amp*pi*(indexj - 1 - N)*(dcmplx(theta*cost, 1.0D+0)*&
                       u_hat(i - 1, j) + dcmplx(-theta*cost, 1.0D+0)&
                       *u_hat(i + 1, j))
            uy_hat(1,j)=Amp*pi*(indexj - 1 - N)*(dcmplx(theta*cost, 1.0D+0)*&
                        gatherLu_hat(N-indexj+2) + dcmplx(-theta*cost, 1.0D+0)&
                         *u_hat(2, j))
            uy_hat(N/2+1,j)=Amp*pi*(indexj - 1 - N)*(dcmplx(theta*cost, 1.0D+0)*&
                        u_hat(N/2, j) + dcmplx(-theta*cost, 1.0D+0)&
                        *gatherRu_hat(N -indexj  + 2))
        end if
      end do
    end do
    if(myrank .eq. 0) then
        uy_hat(1, 1) = dcmplx(0.0D+0, 0.0D+0)
        uy_hat(N/2+1, 1) = dcmplx(0.0D+0, 0.0D+0)
    end if

    do j = 1, localN
      do i = 1, N/2
        rhs_hat(i, j) = ux_hat(i, j) + uy_hat(i, j) + &
                      dcmplx(C, 4*pi*lam*eps*(i-1))*u_hat(i, j)
      !  rhs_hat(i, j) = uy_hat(i, j) 
      end do
       rhs_hat(N/2+1, j) = ux_hat(N/2+1,j) + uy_hat(N/2+1, j) +&
                           dcmplx(C, 0.0D+0)*u_hat(N/2+1,j)
      ! rhs_hat(N/2+1, j) = uy_hat(N/2+1,j) 
    end do

    do j = 1,localN
        indexj = j + local_0_start
       do i = 1, N/2 + 1
        if(j .lt. N/2+2) then 
        u_hat(i, j) = 1.0D+0/(4*pi*pi*((i-1)*(i-1)+(indexj-1)*(indexj-1))&
                      *eps*dt+1)*(u_hat(i,j) + dt*rhs_hat(i,j))
        else 
        u_hat(i,j)= 1.0D+0/(4*pi*pi*((i-1)*(i-1)+(indexj-1-N)*(indexj-1-N))&
                    *eps*dt+1)*(u_hat(i,j) + dt*rhs_hat(i,j))
        end if
       end do
    end do

!    do j = 1, N/2
!      do i = 2, N/2
!        uy_hat(i,j)=Amp*pi*(j - 1)*(dcmplx(theta*cost, 1.0D+0)*&
!                    u_hat(i - 1, j) + dcmplx(-theta*cost, 1.0D+0)&
!                    *u_hat(i + 1, j))
!      end do
!    end do

!    do j = N/2 + 2, N
!      do i = 2, N/2
!        uy_hat(i, j)=Amp*pi*(j - 1 - N)*(dcmplx(theta*cost, 1.0D+0)*&
!                       u_hat(i - 1, j) + dcmplx(-theta*cost, 1.0D+0)&
!                       *u_hat(i + 1, j))
!      end do
!    end do

!    uy_hat(1, 1) = dcmplx(0.0D+0, 0.0D+0)
!    do j = 2, N/2
!      uy_hat(1, j) = Amp*pi*(j - 1)*(dcmplx(theta*cost, 1.0D+0)*&
!               conjg(u_hat(2,N-j+2)) + dcmplx(-theta*cost, 1.0D+0)&
!                     *u_hat(2, j))
!    end do
!
!    do j = N/2 + 2, N
!      uy_hat(1,j)=Amp*pi*(j - 1 - N)*(dcmplx(theta*cost, 1.0D+0)*&
!               conjg(u_hat(2,N-j+2)) + dcmplx(-theta*cost, 1.0D+0)&
!                     *u_hat(2, j))
!    end do
!
!    uy_hat(N/2+1, 1) = dcmplx(0.0D+0, 0.0D+0)
!    do j = 2, N/2
!      uy_hat(N/2+1,j)=Amp*pi*(j - 1)*(dcmplx(theta*cost, 1.0D+0)*&
!                        u_hat(N/2, j) + dcmplx(-theta*cost, 1.0D+0)&
!                        *conjg(u_hat(N/2, N - j + 2)))
!    end do
!
!    do j = N/2 + 2, N
!      uy_hat(N/2+1,j)=Amp*pi*(j - 1 - N)*(dcmplx(theta*cost, 1.0D+0)*&
!                        u_hat(N/2, j) + dcmplx(-theta*cost, 1.0D+0)&
!                        *conjg(u_hat(N/2, N - j + 2)))
!    end do
!
!    do i = 1, N/2 + 1
!      uy_hat(i, N/2 + 1) = dcmplx(0.0D+0, 0.0D+0)
!    end do
    deallocate(tmpv_up)
    deallocate(tmpv_dn)
    deallocate(tmpv)

    end 
      
