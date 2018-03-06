subroutine fftw_mpi_2d(N, ierr, u_hat, Nu_hat,  alloc_local, local_M, local_j_offset)
use, intrinsic :: iso_c_binding
include '/opt/fftw3/include/fftw3-mpi.f03'
integer :: ierr
!!integer*8 :: M 
integer*8 :: N  
!!integer(C_INTPTR_T), parameter :: N = 4  
!!integer(C_INTPTR_T), parameter :: L 
type(C_PTR) :: plan, cdata, c1data
complex(C_DOUBLE_COMPLEX), pointer :: data(:,:)
  
integer(C_INTPTR_T) :: i, j, alloc_local, local_M, local_j_offset
complex*16:: u_hat(1:N,1:local_M), Nu_hat(1:N/2+1, 1:local_M)
!!    call MPI_Init(ierr)
    call fftw_mpi_init()
! get local data size and allocate (note dimension reversal)
!    alloc_local = fftw_mpi_local_size_2d(N, N, ierr, &
!    local_M, local_j_offset)
!    localN = local_M
!    offset = local_j_offset
!    allocate(u_hat(1:N, 1:local_M))   
!    write(*,*)'localM:', local_M, local_j_offset, alloc_local 
    cdata = fftw_alloc_complex(alloc_local)
    call c_f_pointer(cdata, data, [N ,local_M])
!create MPI plan for in-place forward DFT (note dimension reversal)
    plan = fftw_mpi_plan_dft_2d(N, N, data, data, ierr, &
    FFTW_FORWARD, FFTW_MEASURE)
!initialize data to some function my_function(i,j)
    do j = 1, local_M
        do i = 1, N
            data(i, j) =  u_hat(i, j) 
        end do
    end do
! compute transform (as many times as desired)
   call fftw_mpi_execute_dft(plan, data, data)
   
   do j = 1, local_M
       do i = 1, N/2+1
         Nu_hat(i, j) = data(i, j)
!         write(*,*), data(i, j)
       end do
   end do

   call fftw_destroy_plan(plan)
   call fftw_free(cdata)
!   call MPI_FINALIZE(ierr)
end

subroutine getLocalInf(N, ierr, alloc_local,  localN, offset)
use, intrinsic :: iso_c_binding
include '/opt/fftw3/include/fftw3-mpi.f03'
integer :: ierr
integer*8 :: N, localN, offset 
integer(C_INTPTR_T) :: alloc_local, local_M, local_j_offset
    call fftw_mpi_init()
! get local data size and allocate (note dimension reversal)
    alloc_local = fftw_mpi_local_size_2d(N, N, ierr, &
    local_M, local_j_offset)
    localN = local_M
    offset = local_j_offset
end



