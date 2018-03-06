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
 * 
 */
/**
 * @filename      KPP_ComputePlane.c
 * @Synopsis      
 * @author kuangxiong
 * @version 1
 * @date 2017-11-25
 */
#include"../fun.h"

/* ============================================================================*/
/**
 * @Synopsis      this function is used to compute plane solution of KPP problem
 *
 * @Param         N        input, number of basis function
 * @Param         t        input, time 
 * @Param         dt       input, the step of t
 * @Param         C        input, paramter
 * @Param         Amp      input, paramter 
 * @Param         lam      input , paramter \lambda
 * @Param         theta    input, \theta
 * @Param         eps      input, 
 * @Param         u_hat    input/output, input the curtime solution , output the next time solution
 * @Param         myrank   input   the id of process
 * @Param         nprocs   input   number of process
 * @Param         MPI_Comm input   communication domain
 */
/* ============================================================================*/
void KPP_ComputePlane_lz(int N, DOUBLE t, DOUBLE dt, DOUBLE C, DOUBLE Amp, DOUBLE lam, DOUBLE theta, 
         DOUBLE eps, DOUBLE cost, complex *u_hat, int myrank, int nprocs, long int local_n0, 
         long int local_0_start, MPI_Comm comm)
{   
    complex *ux_hat, *uy_hat, *rhs_hat, *localU, *localD, *tmpu_hat, *tmp_hat, *tmplocalL, *tmplocalR, *tmpL, *tmpR, a, b;
//    complex ux_hat[local_n0*(N/2+1)], uy_hat[local_n0*(N/2+1)], rhs_hat[local_n0*(N/2+1)], localU[local_n0*(N/2+1)], 
//			localD[N/2+1], tmpu_hat[local_n0*(N/2+1)], tmp_hat[local_n0*(N/2+1)], tmplocalL[local_n0], 
//			tmplocalR[local_n0], tmpL[N], tmpR[N], a, b;
    int i, j, k, n, tmpi, tmpj, befRank, nextRank, ONE, Index[local_n0];
    MPI_Status status;


//	printf("lamlam:%f\n", lam);

    ONE = 1;
    tmpu_hat = calloc(local_n0 * (N/2 + 1), sizeof(*tmpu_hat));
    tmp_hat = calloc(local_n0 * (N/2 + 1), sizeof(*tmp_hat));
    ux_hat = calloc(local_n0 * (N/2 + 1), sizeof(*ux_hat));
    uy_hat = calloc(local_n0 * (N/2 + 1), sizeof(*uy_hat));
    rhs_hat = calloc(local_n0 * (N/2 + 1), sizeof(*rhs_hat));
    localU = calloc(N/2 + 1, sizeof(*localU));
    localD = calloc(N/2 + 1, sizeof(*localD));
    tmplocalL = calloc(local_n0, sizeof(*tmplocalL));
    tmplocalR = calloc(local_n0, sizeof(*tmplocalR));
    tmpL = calloc(N, sizeof(*tmpL));
    tmpR = calloc(N, sizeof(*tmpR));
     
    for(i=0; i< local_n0; i++){
        tmplocalL[i] = conj(u_hat[i*(N/2+1)+1]);
        tmplocalR[i] = conj(u_hat[i*(N/2+1)+N/2-1]);
    }
    for(i=0; i< local_n0; i++){	
		Index[i] = local_0_start + i;
	}
    MPI_Allgather(tmplocalL, local_n0, MPI_C_DOUBLE_COMPLEX, tmpL, local_n0, MPI_C_DOUBLE_COMPLEX, comm);
    MPI_Allgather(tmplocalR, local_n0, MPI_C_DOUBLE_COMPLEX, tmpR, local_n0, MPI_C_DOUBLE_COMPLEX, comm);
    /****************exchange information between neighboring rank */
    n = N/2 + 1;
    befRank = (myrank - 1 + nprocs)%nprocs;
    nextRank = (myrank + 1)%nprocs;  
//    if(myrank==1)
//		for(i=0; i< N; i++)
//			printf("hello:%f\t%f\n", creal(tmpR[i]), cimag(tmpR[i]));

	MPI_Sendrecv(u_hat, n ,MPI_C_DOUBLE_COMPLEX, befRank, 990,
                 localD, n, MPI_C_DOUBLE_COMPLEX, nextRank, 990, comm, &status);

    MPI_Sendrecv(&u_hat[(local_n0-1)*n], n ,MPI_C_DOUBLE_COMPLEX, nextRank, 992,
                 localU, n, MPI_C_DOUBLE_COMPLEX, befRank, 992, comm, &status);
//	if(myrank==2)
//		for(i=0; i< N/2+1;i++)
//			printf("hahah111:%d\t%d\t%f\t%f\n",N/2+1, nextRank, creal(localU[i]), cimag(localU[i]));
   /**********************computing ux_hat******************************/ 
    for(i=1; i< local_n0-1; i++){
        for(j=0; j< N/2-1; j++){
            ux_hat[i*(N/2+1)+j] = Amp*((PI * j* (theta*cost + I) + lam/2*(1-theta*cost*I))
                    *u_hat[(i-1)*(N/2+1)+ j] + (PI*j*(-theta*cost + I) +lam/2
                    *(1 + theta*cost*I))*u_hat[(i+1)*(N/2+1) + j]);
        }
    }
//	if(myrank==0)
//		printf("ux_hat%f&%f\t%f&%f\n", ux_hat[1*(N/2+1)+N/2-1], u_hat[1*(N/2+1)+N/2-1]);
    for(j=0; j< N/2; j++){
        ux_hat[j] = Amp*((PI * j* (theta * cost + I) + lam/2*(1-theta*cost*I))
                            *localU[j] + (PI*j*(-theta*cost + I) +lam/2
                            *(1 + theta*cost*I))*u_hat[(N/2+1) + j]);
    }
        
    for(j=0; j< N/2; j++){
        ux_hat[(local_n0-1)*(N/2+1) + j] = Amp*((PI * j* (theta*cost + I) + lam/2*(1-theta*cost*I))
                            *u_hat[(local_n0-2)*(N/2+1)+j] + (PI*j*(-theta*cost + I) +lam/2
                            *(1 + theta*cost*I))*localD[j]);
    }
    for(i=1; i< local_n0-1; i++)
        ux_hat[i*(N/2 + 1) + N/2] = Amp*lam/2*((1-theta*cost*I)
                            *u_hat[(i-1)*(N/2+1)+ N/2] + (1 + theta*cost*I)*u_hat[(i+1)*(N/2+1)+N/2]);
       
    ux_hat[N/2] = Amp*lam/2*((1-theta*cost*I)*localU[N/2] + 
                   (1 + theta*cost*I)*u_hat[1*(N/2+1) + N/2]);
    ux_hat[(local_n0-1)*(N/2+1) + N/2] = Amp*lam/2*((1-theta*cost*I)
                  *u_hat[(local_n0-2)*(N/2+1)+N/2] + (1 + theta*cost*I)*localD[N/2]);
/**********************computing uy_hat*******************************/
    for(i=0; i< local_n0; i++){
        for(j=1; j< N/2; j++){
            if(Index[i] < N/2){
            uy_hat[i*(N/2+1) + j]=Amp*PI*i*((theta*cost + I)*u_hat[i * n +j-1]
                                 + (-theta*cost + I)*u_hat[i * n + j+1]);
            }
            else if(Index[i] == N/2){
            uy_hat[i*(N/2+1) + j] =0;
            }
            else{
            uy_hat[i*(N/2+1) + j]=Amp*PI*(Index[i]-N)*((theta*cost + I)*u_hat[i * n +j-1]
                                 + (-theta*cost + I)*u_hat[i * n + j+1]);
			
            }
        }
    }
//		if(myrank==0)
//			printf("hah1hahah:%e&%e\t%e&%e\n", u_hat[1*(N/2+1)+N/2-2], u_hat[1*(N/2+1)+N/2]);
    for(i=0; i< local_n0; i++){
        if(Index[i] == 0 || Index[i] == N/2){
            uy_hat[i * (N/2+1) + 0] = 0;
            uy_hat[i * (N/2+1) + N/2] = 0;
        }
        else if(Index[i] < N/2){
            uy_hat[i*(N/2+1) + N/2] = Amp*PI *i * ((theta*cost + I)* u_hat[i*(N/2+1) + N/2-1]
                                     +(-theta*cost + I) * tmpR[N-Index[i]]);
            uy_hat[i*(N/2+1)] =Amp*PI*i*((theta*cost+ I)*tmpL[N-Index[i]]
                                +(-theta*cost + I)*u_hat[i*(N/2+1) + 1]);
        }
        else{ 
            uy_hat[i*(N/2+1) + N/2] = Amp*PI *(Index[i]-N) * ((theta*cost + I)* u_hat[i*(N/2+1) + N/2-1]
                                     +(-theta*cost + I) * tmpR[N-Index[i]]);
            uy_hat[i*(N/2+1)] =Amp*PI*(Index[i]-N)*((theta*cost+ I)*tmpL[N-Index[i]]
                                +(-theta*cost + I)*u_hat[i*(N/2+1) + 1]); 
		   }
		}
//			if(myrank==0)
//				printf("hah1hahah:%d\t%f&%f\n", Index[1], u_hat[1*(N/2+1)+1]);
//		if(myrank==1)
//			for(i=0; i< N; i++)
//			printf("gg11:%f&%f\n", uy_hat[(local_n0-1)*(N/2+1)+N/2]);
/*******************************************************************/
//       if(myrank==0)
//		   printf("gg:%f\t%f\n", uy_hat[1*(N/2+1)+0]);
        for(i=0; i< local_n0; i++){
              for(j=0; j< N/2; j++){
                  rhs_hat[i*(N/2+1)+j] = ux_hat[i*(N/2+1)+j] + uy_hat[i*(N/2+1)+j] 
                                       +(C + I*(4*PI*lam*eps*j))*u_hat[i*(N/2+1)+j];
            }
        }
        for(i=0; i< local_n0; i++){
            rhs_hat[i*(N/2+1)+N/2] = ux_hat[i*(N/2+1)+j] + uy_hat[i*(N/2+1)+j]
                                   + C * u_hat[i*(N/2+1)+j];
        }
        
        for(i=0; i<local_n0; i++){
            for(j=0; j< N/2+1; j++){
                if(Index[i]< N/2+1){
                    u_hat[i*(N/2+1)+j] = 1.0/(4*PI*PI*(Index[i]*Index[i]+j*j)*eps*dt+1)*
                                          (u_hat[i*(N/2+1)+j]+ dt * rhs_hat[i*(N/2+1)+j]);
                }
                else{
                    u_hat[i*(N/2+1)+j] = 1.0/(4*PI*PI*((Index[i]-N)*(Index[i]-N)+j*j)*eps*dt+1)*
                                          (u_hat[i*(N/2+1)+j]+ dt * rhs_hat[i*(N/2+1)+j]);
                }
            }     
        }
		/***************************test*******************/
//        for(i=0; i< local_n0; i++)
//			for(j=0; j< (N/2+1); j++)
//				u_hat[i*(N/2+1)+j] = u_hat[i*(N/2+1)+j];
//	if(myrank ==0)
//		printf("googd:%e&%e\t%e&%e\n", uy_hat[1*(N/2+1)+N/2-1], u_hat[1*(N/2+1)+N/2-1]);
       free(tmpu_hat);
       free(tmp_hat);
       free(ux_hat);
	   free(uy_hat);
	   free(rhs_hat);
	   free(localU);
	   free(localD);
	   free(tmplocalL);
	   free(tmplocalR);
	   free(tmpL);
	   free(tmpR);
}
