      subroutine sparse_CopyMatrix_Symmetrize(m , N, temp_hat1, BNlap)
      integer:: m, N, i, j, k
      complex*16 :: temp_hat1(1:N*N, 1:m), BNlap(1:(N/2+1)*N, 1:m)
      do k = 1, m
       do j = 1, N
         do i = 1, N/2 + 1
           temp_hat1(i+(j-1)*N,k) = BNlap(i+(j-1)*(N/2+1),k)
         end do
       end do
       ! Symmetrize
       do j = 2, N
         do i = N/2 + 2, N
           temp_hat1(i+(j-1)*N,k)=conjg(BNlap(N+2-i+(N+1-j)*(N/2+1),k))
         end do
       end do
        do i = N/2 + 2, N
             temp_hat1(i,k)=conjg(BNlap(N+2-i,k))
        end do
       end do
       end

       subroutine sparse_BuildBNx(m, N, pi, BNx,  BN)
       integer:: m , N, i, j, k
       real*8 :: pi
       complex*16::BNx(1:(N/2+1)*N, 1:m), BN(1:(N/2+1)*N, 1:m)
        do k = 1, m
           do j = 1, N
             do i = 1, N/2
                BNx(i+(j-1)*(N/2+1),k)= dcmplx(0.0D+0, 2*pi*(i-1))&
                    *BN(i+(j-1)*(N/2+1),k)
             end do
            end do

            do j = 1, N
                BNx(N/2+1+(j-1)*(N/2+1),k) = dcmplx(0.0D+0, 0.0D+0)
            end do
         end do
       end


       subroutine sparse_BuildBNy(m, N, pi, BNy,  BN)
       integer:: m , N, i, j, k
       real*8 :: pi
       complex*16::BNy(1:(N/2+1)*N, 1:m), BN(1:(N/2+1)*N, 1:m)
  
        do k = 1, m
          do j = 1, N/2
            do i = 1, N/2 + 1
              BNy(i+(j-1)*(N/2+1),k)= dcmplx(0.0D+0, 2*pi*(j-1))&
                *BN(i+(j-1)*(N/2+1),k)
            end do
          end do

          do j = N/2 + 2, N
            do i = 1, N/2 + 1
              BNy(i+(j-1)*(N/2+1),k)= dcmplx(0.0D+0, 2*pi*(j-1-N))&
                *BN(i+(j-1)*(N/2+1),k)
            end do
          end do

          do i = 1, N/2 + 1
            BNy(i+(N/2+1-1)*(N/2+1),k) = dcmplx(0.0D+0, 0.0D+0)
          end do
        end do      
       end


      subroutine sparse_BuildBNlap(m, N, pi, BNlap,  BN)
       integer :: m, N, i, j, k
       real*8::pi
       complex*16 :: BNlap(1:(N/2+1)*N, 1:m), BN(1:(N/2+1)*N, 1:m)
       do k = 1, m
         do j = 1, N/2 + 1
          do i = 1, N/2 + 1
            BNlap(i+(j-1)*(N/2+1),k)=-4*pi*pi*((i-1)*(i-1)+(j-1)*(j-1))&
              *BN(i+(j-1)*(N/2+1),k)
           end do
         end do

          do j = N/2 + 2, N
            do i = 1, N/2 + 1
              BNlap(i+(j-1)*(N/2+1),k)=-4*pi*pi*((i-1)*(i-1)+&
               (j-1-N)*(j-1-N))*BN(i+(j-1)*(N/2+1),k)
             end do
           end do
         end do
       end


       subroutine sparse_BuildRBBNxRBBNy1(m, N, cosx, Amp, Bx, By,&
                  RBNx, RBNy, RBBNx, RBBNy)
       integer :: m, N, i, j, k
       real*8 :: Bx(1:N, 1:N), By(1:N, 1:N), RBNx(1:N, 1:N, 1:m),&
                 RBNy(1:N,1:N,1:m),RBBNx(1:N, 1:N, 1:m),&
                 RBBNy(1:N, 1:N, 1:m),cosx(1:N), Amp
      forall(j = 1:N)
         Bx(1:N, j) = Amp*cosx(j)
      end forall

      ! The part of By without time t
      forall(j = 1:N)
         forall(i = 1:N)
             By(i, j) = Amp*cosx(i)
         end forall
      end forall

      !RBBNx(1:N, 1:N, 1:m) is Bx*RBNx in real space without time t
        do k = 1, m
          forall(j = 1:N)
              forall(i = 1:N)
                  RBBNx(i, j, k) = Bx(i, j)*RBNx(i, j, k)
               end forall
           end forall
        end do

        ! RBBNy(1:N, 1:N, 1:m) is By*RBNy in real space without time t
          do k = 1, m
            forall(j = 1:N)
                forall(i = 1:N)
                   RBBNy(i, j, k) = By(i, j)*RBNy(i, j, k)
                 end forall
             end forall
            end do
      end

       subroutine sparse_BuildRBBNxRBBNy2(m, N, sinx, Amp, theta,  Bx, By,&
                  RBNx, RBNy, RBBNx,RBBNy)
       integer :: m, N, i, j, k
       real*8 :: Bx(1:N, 1:N), By(1:N, 1:N), RBNx(1:N, 1:N, 1:m),&
                 RBNy(1:N,1:N,1:m), RBBNx(1:N, 1:N, 1:m), RBBNy(1:N, 1:N, 1:m),&
                 sinx(1:N), theta, Amp

      ! The part of Bx with time t
      forall(j = 1:N)
        Bx(1:N, j) = Amp*theta*sinx(j)
      end forall

     ! The part of By with time t
      forall(j = 1:N)
         forall(i = 1:N)
           By(i, j) = Amp*theta*sinx(i)
         end forall
        end forall

     ! RBBNx(1:N, 1:N, 1:m) is Bx*RBNx in real space with time t
       do k = 1, m
          forall(j = 1:N)
             forall(i = 1:N)
              RBBNx(i, j, k) = Bx(i, j)*RBNx(i, j, k)
             end forall
          end forall
        end do

     ! RBBNy(1:N, 1:N, 1:m) is By*RBNy in real space with time t
       do k = 1, m
       forall(j = 1:N)
         forall(i = 1:N)
             RBBNy(i, j, k) = By(i, j)*RBNy(i, j, k)
         end forall
        end forall
      end do
      end
