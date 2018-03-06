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
     long int alloc_local, localN, local_0_start;
     double *u, tmax, cost, cosx[N], sinx[N], h, *RtmpV;
	 
	 complex *u_hat, *CtmpV;
     fftw_plan p1, p2;
     int i, k, j, niter = 0, ONE=1, tmpN;
     FILE *fp, *fp1, *fp2;

     fp = fopen("nodes.text","w");
     fp1 = fopen("nodes1.text","w");
     fp2 = fopen("nodes2.text","w");
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
/******************************************************************************/
	 k = 1;
     t = 0;
     while(t < T)
	 {
        niter = niter + 1;
        cost = cos(t);
        t = t + dt;
        KPP_ComputePlane(N, t, dt, C, Amp, lam, theta, eps, cost, u_hat, myrank, nprocs,localN, local_0_start,  comm);
		k++;

        zcopy_(&tmpN, u_hat, &ONE, CtmpV, &ONE);

	    fftw_execute(p2);

        if(myrank==0){
            printf("11111111:%f\n", t);
            fprintf(fp, "%f\t%f\t%f\t%f\n", t, RtmpV[(localN-1)*2*(N/2+1)+256]/(N*N),
                                       RtmpV[(localN-1)*2*(N/2+1)+512]/(N*N), RtmpV[(localN-1)*2*(N/2+1)+768]/(N*N));
        }
        if(myrank==nprocs/2)
            fprintf(fp1, "%f\t%f\t%f\t%f\n", t, RtmpV[(localN-1)*2*(N/2+1)+256]/(N*N),
                                        RtmpV[(localN-1)*2*(N/2+1)+512]/(N*N), RtmpV[(localN-1)*2*(N/2+1)+768]/(N*N));
        if(myrank==nprocs-2)
            fprintf(fp2, "%f\t%f\t%f\t%f\n", t, RtmpV[(localN-1)*2*(N/2+1)+256]/(N*N),
                                        RtmpV[(localN-1)*2*(N/2+1)+512]/(N*N), RtmpV[(localN-1)*2*(N/2+1)+768]/(N*N)); 
     }

	 free(CtmpV);
	 free(RtmpV);
	 free(u_hat);
	 free(u);

     fftw_destroy_plan(p1);
     fftw_destroy_plan(p2);
     fclose(fp);
     fclose(fp1);
     fclose(fp2);

}
