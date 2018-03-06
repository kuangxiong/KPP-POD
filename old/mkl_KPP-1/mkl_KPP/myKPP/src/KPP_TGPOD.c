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
void KPP_TGPOD(int N, DOUBLE lam, DOUBLE eps, DOUBLE tau, DOUBLE theta, DOUBLE dt, DOUBLE t, 
			  DOUBLE Amp, DOUBLE cM, DOUBLE C, DOUBLE y1, int myrank, int nprocs, 
			  MPI_Comm comm1)
{     
     long int alloc_local, localN, local_0_start, ncup;
     double *u, cost, *S, *S1, cosx[N], sinx[N], h, *RtmpV, *uPOD, CoardT, 
            errind, tmpt, err, localerr, Norm, localNorm, zero, Cdt;
	 
	 complex *u_hat, *UP_hat, *gatherU, *gatherU1, testu_hat[N*(N/2+1)], *tmp4svd_hat,
			 *B, *B1, *BN1, *BN, *PODu0, *Firm, *Varm, *tmpFirm, *PODu, *CPODu, *tmp4svd_hat1,
			 *CtmpV, *BNlap, *TMat, *NotTMat, alpha, tmpvalue, beta, *mixMat, *err_hat, *Coaru_hat;
     fftw_plan p1, p2;
     int i, k, j, niter = 0, ONE=1, N1, N2, nPOD, nPOD1, nPOD2,  tmpN, tmpN1, CoarN, flag, CoarIntval;
	 FILE *fp;

	 fp = fopen("TG-APODerr.text","w");
     alloc_local = fftw_mpi_local_size_2d(N, N/2+1, comm1, &localN, &local_0_start);
     u_hat = calloc(localN*(N/2+1), sizeof(*u_hat));
     Coaru_hat = calloc(localN*(N/2+1), sizeof(*Coaru_hat));
     u = calloc(localN * N, sizeof(*u));
     uPOD = calloc(localN * N, sizeof(*uPOD));
     
	 CtmpV = calloc(alloc_local, sizeof(*CtmpV));
     RtmpV = calloc(2 * alloc_local, sizeof(*RtmpV));
     err_hat = calloc(localN*(N/2+1), sizeof(*err_hat));
	 CPODu = calloc(localN*(N/2+1), sizeof(*CPODu));
     N1 = localN * (N/2+1);
	 N2 = localN * N;
/************get initial u and u_hat**************/
     p1 = fftw_mpi_plan_dft_r2c_2d(N, N, RtmpV, CtmpV, comm1, FFTW_MEASURE);
     p2 = fftw_mpi_plan_dft_c2r_2d(N, N, CtmpV, RtmpV, comm1, FFTW_MEASURE);
  
	 for(i=0; i< localN; i++){
         for(j=0; j< N; j++){
            RtmpV[i*2*(N/2+1)+j] = 1;
         }
     }

     MPI_Barrier(comm1);     
     fftw_execute(p1);
	 tmpN = localN*(N/2+1);
	 zcopy_(&tmpN, CtmpV, &ONE, u_hat, &ONE);
     h = 1.0/N;
	 for(i=0; i< N; i++){
		cosx[i] = cos(2*PI*i*h);
		sinx[i] = sin(2*PI*i*h);
	 }
/*************************************************/
	 ncup = (int)(tmax/dt)+2;
     ncup = ncup/intval + 1; 
     zero = 0.0;
     CoarIntval = 10;
     Cdt = 10 * dt;	 
     CoardT = 0.06;
     gatherU = calloc(N1*ncup, sizeof(*gatherU));
     tmp4svd_hat = calloc(localN*N*ncup, sizeof(*tmp4svd_hat));
     S = calloc(ncup, sizeof(*S));
	 B = calloc(localN*N*N*N, sizeof(*B));
	 zcopy_(&N1,  u_hat, &ONE, gatherU, &ONE);

	 k = 1;
     t = 0;
     while(t <= tmax)
     {
        niter = niter + 1;
        cost = cos(t);
        t = t + dt;

        KPP_ComputePlane(N, t, dt, C, Amp, lam, theta, eps, cost, u_hat, myrank, nprocs,localN, local_0_start,  comm1);
        if(k%intval ==0)
		    zcopy_(&N1,  u_hat, &ONE, &gatherU[(k/intval)* N1], &ONE);
		k++;
        if(myrank==0)
            printf("111111111%f\n", t);
     }

	 KPP_Build_FFTMatrix1(N, ncup, localN, gatherU, tmp4svd_hat, comm1);
     KPP_mpi_svd(ncup, N, localN, tmp4svd_hat, S, B, nprocs, myrank);
     nPOD = KPP_GetPODNumber(S, ncup, gamma1);
     BN1 = calloc(nPOD*localN*N, sizeof(*BN1));
     BN = calloc(nPOD*localN*(N/2+1), sizeof(*BN));
	 for(i=0; i< nPOD; i++)
	 	zcopy_(&N2, &B[i*N2], &ONE, &BN1[i*N2], &ONE);
     for(k=0; k<nPOD; k++)
        for(i=0; i< localN; i++)
            for(j=0; j< N/2+1; j++)
                BN[k*localN*(N/2+1)+i*(N/2+1)+j] = BN1[k*localN*N+i*N+j];
/**********************************************************************************/   
     flag ==0;
	 zcopy_(&N1,  u_hat, &ONE, Coaru_hat, &ONE);
	 zcopy_(&N1,  u_hat, &ONE, CPODu,  &ONE);

//     printf("helloworld\n");

     MPI_Barrier(comm1);
     while(t<T)
     {
       tmpt = t;
/*****************************collect Coarse Grid PW solution ***********************/
       if(t + CoardT > T)
           CoarN = (T - t)/Cdt;
       else 
           CoarN = CoardT/Cdt;
//       printf("test:%d\n", CoarN);
       gatherU1 = calloc(N1*(CoarN/CoarIntval)+1, sizeof(*gatherU1));
       S1 = calloc(CoarN/CoarIntval + 1, sizeof(*S1));
	   B1 = calloc(localN*N*N*N, sizeof(*B1));
       tmp4svd_hat1 = calloc(localN*N*(CoarN/CoarIntval+1), sizeof(*tmp4svd_hat1));
       
       for(k=0; k< CoarN; k++){
          t = t + Cdt;
          KPP_ComputePlane(N, t, Cdt, C, Amp, lam, theta, eps, cost, Coaru_hat, myrank, nprocs,
                           localN, local_0_start,  comm1);
          if(k%CoarIntval ==0)
	         zcopy_(&N1,  Coaru_hat, &ONE, &gatherU1[(k/CoarIntval)* N1], &ONE);
       }
	   KPP_Build_FFTMatrix1(N, CoarN/CoarIntval + 1, localN, gatherU1, tmp4svd_hat1, comm1);
       KPP_mpi_svd(CoarN/CoarIntval + 1, N, localN, tmp4svd_hat1, S1, B1, nprocs, myrank);
//       if(myrank==0){
//       for(i=0; i< CoarN/CoarIntval+1; i++)
//           printf("S1:%f\n", S1[i]);
       
       nPOD1 = KPP_GetPODNumber(S1, CoarN/CoarIntval + 1, gamma2);
/******************************Build POD Matrix *************************************/
       free(S);   free(B);   free(gatherU1);  free(tmp4svd_hat1);       
       mixMat = calloc((nPOD+nPOD1)*localN*N, sizeof(*mixMat));
       S = calloc(nPOD+nPOD1, sizeof(*S));
	   B = calloc(localN*N*N*N, sizeof(*B)); 
       for(i=0; i< nPOD; i++){
          zcopy_(&N2, &BN1[i*N2], &ONE, &mixMat[i*N2], &ONE);
       }
       printf("number: %d\t%d\n", nPOD, nPOD1);
       for(i=0; i< nPOD1; i++){
          zcopy_(&N2, &B1[i*N2], &ONE, &mixMat[(i+nPOD)*N2], &ONE);
       }
       KPP_mpi_svd(nPOD+nPOD1, N, localN, mixMat, S, B, nprocs, myrank); 
       nPOD = KPP_GetPODNumber(S, nPOD+nPOD1, gamma1);
       nPOD = 10;
       if(myrank==0)
           for(i=0; i< 10; i++)
               printf("%f\n", S[i]);

       free(mixMat); free(B1), free(BN); free(BN1);

       Firm = calloc(nPOD*nPOD, sizeof(*Firm));
       tmpFirm = calloc(nPOD*nPOD, sizeof(*tmpFirm));
	   Varm = calloc(nPOD*nPOD, sizeof(*Varm));
       KPP_BuildPODMatrix(nPOD, N, localN, local_0_start, alloc_local, eps,
					   theta, lam, C, Amp, BN1, Firm, Varm, cosx, sinx, myrank, 
					   nprocs, comm1);
       BN1 = calloc(nPOD*localN*N, sizeof(*BN1));
       BN = calloc(nPOD*localN*(N/2+1), sizeof(*BN));
       PODu0 = calloc(nPOD, sizeof(*PODu0));
       PODu = calloc(nPOD, sizeof(*PODu));
	   for(i=0; i< nPOD; i++)
	      zcopy_(&N2, &B[i*N2], &ONE, &BN1[i*N2], &ONE);
       if(myrank==0){
           printf("%f&%f\t%f&%f\n", B[3*N2+3*N], B[3*N2+3*N+N/2+1 +3]);      
       }

       for(k=0; k<nPOD; k++)
         for(i=0; i< localN; i++)
            for(j=0; j< N/2+1; j++)
                BN[k*localN*(N/2+1)+i*(N/2+1)+j] = BN1[k*localN*N+i*N+j];
/**************************************************************************************/
/**********************get initial POD solution ****************************/
       KPP_GetInitialPOD(N, localN, nPOD, CPODu, BN1, PODu0, myrank, nprocs);   
       if(myrank==0)
       for(i=0; i< nPOD; i++)
            printf("hello:%d\t%f&%f\n",nPOD, PODu0[i]);

       t = tmpt; 
       while(t < tmpt + CoarN* Cdt)
       {
       for(i=0; i< nPOD; i++)
    		for(j=0; j< nPOD; j++)
    			tmpFirm[i*nPOD+j] = Firm[i*nPOD+j] + cos(t)*Varm[i*nPOD+j];
        cost = cos(t);
        t = t + dt;
        alpha = dt;
        beta = 0.0;
    	/*******************update PODu***************/
    	zgemv_("C", &nPOD, &nPOD, &alpha, tmpFirm, &nPOD, PODu0, &ONE, &beta, PODu, &ONE);    
    	for(i=0; i< nPOD; i++){
    		PODu[i] = PODu0[i] + PODu[i];
    	}       
     }
    /********************transform POD solution to PW solution *****************/
//       tmpvalue = 1.0;
//       zgemv_("N", &N1, &nPOD, &tmpvalue, BN, &N1, PODu, &ONE, &beta, Coaru_hat, &ONE);
     }
	 free(S);
	 free(B);
	 free(CtmpV);
	 free(RtmpV);
	 free(u_hat);
	 free(CPODu);
	 free(PODu);
	 free(PODu0);
	 free(gatherU);
	 free(tmp4svd_hat);
}
