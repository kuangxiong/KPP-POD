/* Copyright (C) 
 * 2017 - kuangxiong
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 */
#include"../fun.h"
/* ============================================================================*/
/**
 * @Synopsis      this function is use to APOD algorithm to solve KPP　problem
 *
 * @Param         N        number of plane function
 * @Param         laml1    \lambda
 * @Param         eps      \eps
 * @Param         tau       paramter \tau
 * @Param         theta    \theta
 * @Param         dt       the step of time
 * @Param         t         
 * @Param         Amp       paramter 
 * @Param         y1        outer 
 * @Param         myrank    MPI_rank
 * @Param         nprocs    number of numprocs
 * @Param         MPI_Comm   Communication domain
 */
/* ============================================================================*/
void KPP_PW(int N, DOUBLE lam, DOUBLE eps, DOUBLE tau, DOUBLE theta, DOUBLE dt, DOUBLE t, 
			  DOUBLE Amp, DOUBLE cM, DOUBLE C, DOUBLE y1, int myrank, int nprocs, 
			  MPI_Comm comm)
{     
     int N1 = 1000, beta;
     long int alloc_local, localN, local_0_start;
     long int alloc_local1, localN1, local_0_start1;
     double *u, tmax, cost, cosx[N], sinx[N], h, *RtmpV;
     double *u1, cosx1[N1], sinx1[N1], *RtmpV1, localsum, localself, sum, selfsum;
	 
	 complex *u_hat, *CtmpV;
	 complex *u_hat1, *CtmpV1;
     fftw_plan p1, p2;
     fftw_plan Finep1, Finep2;
     int i, k, j, niter = 0, ONE=1, tmpN;

     if(N1%N !=0){
         printf("the Coarse Grid is wrong!\n");
         return ;
     }
     beta = N1/N;
/******************************************************************************
 *              computing for Coarse grid
 * ***************************************************************************/
     alloc_local = fftw_mpi_local_size_2d(N, N/2+1, comm, &localN, &local_0_start);
     u_hat = calloc(localN*(N/2+1), sizeof(*u_hat));
     u = calloc(localN * N, sizeof(*u));
     
	 CtmpV = calloc(alloc_local, sizeof(*CtmpV));
     RtmpV = calloc(2 * alloc_local, sizeof(*RtmpV));
     /************get initial u and u_hat**************/
     p1 = fftw_mpi_plan_dft_r2c_2d(N, N, RtmpV, CtmpV, comm, FFTW_MEASURE);
     p2 = fftw_mpi_plan_dft_c2r_2d(N, N, CtmpV, RtmpV, comm, FFTW_MEASURE);
  

	 for(i=0; i< localN; i++){
         for(j=0; j< N; j++){
            RtmpV[i*2*(N/2+1)+j] = 1;
         }
     }

     MPI_Barrier(comm);     
     fftw_execute(p1);
	 tmpN = localN*(N/2+1);
	 zcopy_(&tmpN, CtmpV, &ONE, u_hat, &ONE);
     h = 1.0/N;
	 for(i=0; i< N; i++){
		cosx[i] = cos(2*PI*i*h);
		sinx[i] = sin(2*PI*i*h);
	 }
/****************************************************************************
 *           computing for Fine Grid
 * **********************************************************************/
     alloc_local1 = fftw_mpi_local_size_2d(N1, N1/2+1, comm, &localN1, &local_0_start1);
     u_hat1 = calloc(localN1*(N1/2+1), sizeof(*u_hat1));
     u1 = calloc(localN1 * N1, sizeof(*u1));
     
	 CtmpV1 = calloc(alloc_local1, sizeof(*CtmpV1));
     RtmpV1 = calloc(2 * alloc_local1, sizeof(*RtmpV1));
     /************get initial u and u_hat**************/
     Finep1 = fftw_mpi_plan_dft_r2c_2d(N1, N1, RtmpV1, CtmpV1, comm, FFTW_MEASURE);
     Finep2 = fftw_mpi_plan_dft_c2r_2d(N1, N1, CtmpV1, RtmpV1, comm, FFTW_MEASURE);
  

	 for(i=0; i< localN1; i++){
         for(j=0; j< N1; j++){
            RtmpV1[i*2*(N1/2+1)+j] = 1;
         }
     }

     MPI_Barrier(comm);     
     fftw_execute(Finep1);
	 tmpN = localN1*(N1/2+1);
	 zcopy_(&tmpN, CtmpV1, &ONE, u_hat1, &ONE);
     h = 1.0/N1;
	 for(i=0; i< N1; i++){
		cosx1[i] = cos(2*PI*i*h);
		sinx1[i] = sin(2*PI*i*h);
     }
/******************************************************************************/
	 tmax = 0.01;
    
	 k = 1;
     t = 0;
     while(t < tmax)
//	 while(niter < 250)
	 {
        niter = niter + 1;
        cost = cos(t);
        t = t + dt;
        KPP_ComputePlane(N, t, dt, C, Amp, lam, theta, eps, cost, u_hat, myrank, nprocs,localN, local_0_start,  comm);
        KPP_ComputePlane(N1, t, dt, C, Amp, lam, theta, eps, cost, u_hat1, myrank, nprocs,localN1, local_0_start1, comm);
		k++;
     }
     
     tmpN = localN * (N/2+1);
	 zcopy_(&tmpN, u_hat, &ONE, CtmpV, &ONE);
     tmpN = localN1 * (N1/2+1);
	 zcopy_(&tmpN, u_hat1, &ONE, CtmpV1, &ONE);
//	 fftw_execute_dft_r2c(p2, CtmpV, RtmpV);
	 fftw_execute(p2);
	 fftw_execute(Finep2);
     for(i=0; i< localN; i++)
         for(j=0; j< N; j++)
             u[i*N+j]=RtmpV[i*2*(N/2+1)+j]/(N*N);
     for(i=0; i< localN1; i++)
         for(j=0; j< N1; j++)
             u1[i*N1+j]=RtmpV1[i*2*(N1/2+1)+j]/(N1*N1);

//     i = 4;
//     j=5;
     localsum = 0;
     localself = 0;
	 for(i=0; i< localN; i++){
		for(j=0; j< N; j++){
            localsum +=fabs(u[i*(N/2+1)+j]-u1[beta*i*(N1/2+1)+j])*fabs(u[i*(N/2+1)+j]-u1[beta*i*(N1/2+1)+j]);
            localself+=u[i*(N/2+1)+j]*u[i*(N/2+1)+j];
        }
     }
     tmpN = localN*N;
     MPI_Reduce(&localsum, &sum, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
     MPI_Reduce(&localself, &selfsum, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
     if(myrank==0)
        printf("%d\t%f\n", beta, sqrt(sum)/sqrt(selfsum));


	 free(CtmpV);
	 free(RtmpV);
	 free(u_hat);
	 free(u);

     fftw_destroy_plan(p1);
     fftw_destroy_plan(p2);
     fftw_destroy_plan(Finep1);
     fftw_destroy_plan(Finep2);

     free(CtmpV1);
	 free(RtmpV1);
	 free(u_hat1);
     free(u1);

}
