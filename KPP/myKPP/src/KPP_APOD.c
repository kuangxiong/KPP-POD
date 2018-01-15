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
void KPP_APOD(int N, DOUBLE lam, DOUBLE eps, DOUBLE tau, DOUBLE theta, DOUBLE dt, DOUBLE t, 
			  DOUBLE Amp, DOUBLE cM, DOUBLE C, DOUBLE y1, int myrank, int nprocs, 
			  MPI_Comm MPI_COMM_WROLD)
{     
     long int alloc_local, localN, local_0_start, ncup;
     double *u, tmax, cost, *S, *S1, gamma, gamma1, gamma2, cosx[N], sinx[N], h, *RtmpV, *uPOD, 
            dT, errind, T;
	 
	 complex *u_hat, *UP_hat, *gatherU, *gatherU1, testu_hat[N*(N/2+1)], *tmp4svd_hat,
			 *B, *B1, *BN1, *BN, *PODu0, *Firm, *Varm, *tmpFirm, *PODu, *CPODu,
			 *CtmpV, *BNlap, *TMat, *NotTMat, alpha, tmpvalue, beta, *mixMat;
     fftw_plan p1, p2;
     int i, k, j, niter = 0, ONE=1, N1, N2, nPOD, nPOD1, nPOD2,  tmpN;
	 FILE *fp;

//	 fp = fopen("solution.text","w");
     alloc_local = fftw_mpi_local_size_2d(N, N/2+1, MPI_COMM_WORLD, &localN, &local_0_start);
     u_hat = calloc(localN*(N/2+1), sizeof(*u_hat));
     u = calloc(localN * N, sizeof(*u));
     uPOD = calloc(localN * N, sizeof(*uPOD));
     
	 CtmpV = calloc(alloc_local, sizeof(*CtmpV));
     RtmpV = calloc(2 * alloc_local, sizeof(*RtmpV));
//	 testu_hat = calloc(N*(N/2+1), sizeof(complex));
	 CPODu = calloc(localN*(N/2+1), sizeof(*CPODu));
     N1 = localN * (N/2+1);
	 N2 = localN * N;
	  
     /************get initial u and u_hat**************/
     p1 = fftw_mpi_plan_dft_r2c_2d(N, N, RtmpV, CtmpV, MPI_COMM_WORLD, FFTW_MEASURE);
     p2 = fftw_mpi_plan_dft_c2r_2d(N, N, CtmpV, RtmpV, MPI_COMM_WORLD, FFTW_MEASURE);
  

	 for(i=0; i< localN; i++){
         for(j=0; j< N; j++){
//           u[i*2*(N/2+1)+j] = 1;
            RtmpV[i*2*(N/2+1)+j] = 1;
         }
     }

     MPI_Barrier(MPI_COMM_WORLD);     
     fftw_execute(p1);
	 tmpN = localN*(N/2+1);
	 zcopy_(&tmpN, CtmpV, &ONE, u_hat, &ONE);
     h = 1.0/N;
	 for(i=0; i< N; i++){
		cosx[i] = cos(2*PI*i*h);
		sinx[i] = sin(2*PI*i*h);
	 }
/*************************************************/
	 T = 0.050;
     tmax = 0.002;
	 dT = 0.0001;
     gamma = 0.99999999999;
     gamma1 = 0.99999999999;
     gamma2 = 0.99999999999;
	 ncup = (int)(tmax/dt)+2;
     gatherU = calloc(N1*ncup, sizeof(*gatherU));
     tmp4svd_hat = calloc(localN*N*ncup, sizeof(*gatherU));
     S = calloc(ncup, sizeof(*S));
	 B = calloc(localN*N*N*N, sizeof(*B));
    
	 zcopy_(&N1,  u_hat, &ONE, gatherU, &ONE);

	 k = 1;
     t = 0;
//     printf("111111111111\n");
     while(t <= tmax)
     {
        niter = niter + 1;
        cost = cos(t);
        t = t + dt;
        KPP_ComputePlane(N, t, dt, C, Amp, lam, theta, eps, cost, u_hat, myrank, nprocs,localN, local_0_start,  MPI_COMM_WORLD);
		zcopy_(&N1,  u_hat, &ONE, &gatherU[k * N1], &ONE);
		k++;
     }
     MPI_Gather(u_hat, localN*(N/2+1), MPI_C_DOUBLE_COMPLEX, testu_hat , localN*(N/2+1), MPI_C_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);

//	printf("11111111111111\n");
//	 if(myrank==0)
//		 for(i=0; i< N; i++)
//			for(j=0; j< N/2+1; j++)
//		printf("%d\t%f\t%f\n",ncup, testu_hat[i*(N/2+1)+j]);
	 KPP_Build_FFTMatrix1(N, ncup, localN, gatherU, tmp4svd_hat);

     KPP_mpi_svd(ncup, N, localN, tmp4svd_hat, S, B, nprocs, myrank);
	 k=0;
	 if(myrank ==0)
		 for(i=0; i< ncup; i++)
			printf("S:%d\t%f\n",ncup,  S[i]);

	nPOD = KPP_GetPODNumber(S, ncup, gamma);
//	if(myrank==0)
//		printf("%d\t%d\n",ncup, nPOD);
/********************get initial POD solution****************/
	BN1 = calloc(nPOD*localN*N, sizeof(*BN1));
	BN = calloc(nPOD*localN*(N/2+1), sizeof(*BN));
	PODu0 = calloc(nPOD, sizeof(*PODu0));
	PODu = calloc(nPOD, sizeof(*PODu));
	for(i=0; i< nPOD; i++)
		zcopy_(&N2, &B[i*N2], &ONE, &BN1[i*N2], &ONE);

    for(k=0; k<nPOD; k++)
        for(i=0; i< localN; i++)
            for(j=0; j< N/2+1; j++)
                BN[k*localN*(N/2+1)+i*(N/2+1)+j] = BN1[k*localN*N+i*N+j];
//**************get initial POD solution********************/
    KPP_GetInitialPOD(N, localN, nPOD, u_hat, BN1, PODu0, myrank, nprocs);
	if(myrank==0)
		for(i=0; i< nPOD; i++)
			printf("PODu0:%f&%f\n", PODu0[i]);

/*************defined matrix for building POD matrix************/
	NotTMat=calloc(nPOD*localN*(N/2+1), sizeof(*NotTMat));
	TMat=calloc(nPOD*localN*(N/2+1), sizeof(*TMat));
    Firm = calloc(nPOD*nPOD, sizeof(*Firm));
    tmpFirm = calloc(nPOD*nPOD, sizeof(*tmpFirm));
	Varm = calloc(nPOD*nPOD, sizeof(*Varm));
//	KPP_BuildPODMatrix(nPOD, N, localN, local_0_start, alloc_local, eps,
//					   theta, lam, C, Amp, BN1, Firm, Varm, cosx, sinx, myrank, 
//					   nprocs);
	KPP_BuildPODMatrix1(nPOD, N, localN, local_0_start, alloc_local, eps,
					   theta, lam, C, Amp, BN1, NotTMat, TMat, Firm, Varm, cosx, sinx, myrank, 
					   nprocs);
//    if(myrank==0)
//        for(i=0; i< nPOD; i++)
//            for(j=0; j< nPOD; j++)
//                printf("%f&%f\n", Firm[i*nPOD+j]);
	tmpN = localN*(N/2+1);
    alpha = dt;
    beta = 0.0;
	while(t <= T) 
    {
    	for(i=0; i< nPOD; i++)
    		for(j=0; j< nPOD; j++)
    			tmpFirm[i*nPOD+j] = Firm[i*nPOD+j] + cos(t)*Varm[i*nPOD+j];
        t = t + dt;
    	/*******************update PODu***************/
    	zgemv_("C", &nPOD, &nPOD, &alpha, tmpFirm, &nPOD, PODu0, &ONE, &beta, PODu, &ONE);    
    	for(i=0; i< nPOD; i++){
    		PODu[i] = PODu0[i] + PODu[i];
    	}
//    	if(myrank==0){
//    		for(i=0; i< nPOD; i++)
//    			printf("%f&%f\n", PODu[i]);
//    	}
//		/*****************transform POD sulution to  plane solution************/
//		for(i=0; i< localN; i++){
//			for(j=0; j< N/2+1; j++){
//				CPODu[i*(N/2+1)+j] = 0;
//				for(k=0; k< nPOD; k++){
//					CPODu[i*(N/2+1)+j] += PODu[k]* BN1[k*localN*N+i*N+j]; 
//				}
//			}
//		}

        tmpvalue = 1.0;
        beta = 0;
        tmpN = localN*(N/2+1);
        zgemv_("N", &tmpN, &nPOD, &tmpvalue, BN, &tmpN, PODu, &ONE, &beta, CPODu, &ONE);
//		zcopy_(&tmpN, CPODu, &ONE, CtmpV, &ONE);
//		fftw_execute(p2);
//        for(i=0; i< localN; i++)
//			for(j=0; j< N; j++)
//				uPOD[i*N+j] = RtmpV[i*2*localN*(N/2+1)+j]/(N*N);  
		/***************************computing err indicator*********************/
        errind =KPP_GetErrIndicator(nPOD, localN, N, cost, dt, PODu, PODu0, BN1, NotTMat, TMat);
//      /***************************Updata POD Matrix**************************/
        tmpN = localN*(N/2+1);
        if(myrank==0)
            printf("errindicator:%f\n", errind);
        if(errind > 0.005){
            /*******************free Matrix***************************/
            free(gatherU); free(tmp4svd_hat); free(S); free(B); free(TMat);
            free(NotTMat); free(Firm); free(Varm); free(tmpFirm); free(PODu0);
            free(PODu);
            /*********************************************************/
            MPI_Barrier(MPI_COMM_WORLD);
            ncup = dT/dt +1;  
            tmp4svd_hat = calloc(localN*N*ncup, sizeof(*tmp4svd_hat));
            S1 = calloc(ncup, sizeof(*S1));
	        B1 = calloc(localN*N*N*N, sizeof(*B1));
            gatherU1 = calloc(N1*ncup, sizeof(*gatherU1));
            // transfer POD solution to PW solution 
            zgemv_("N", &tmpN, &nPOD, &tmpvalue, BN, &tmpN, PODu, &ONE, &beta, u_hat, &ONE);
	        zcopy_(&N1, u_hat, &ONE, gatherU, &ONE);
            k = 1;
            while(t< 0.00025 + dT){
            /***************************************************/
                KPP_ComputePlane(N, t, dt, C, Amp, lam, theta, eps, cost, u_hat, myrank, nprocs,
                                localN, local_0_start,  MPI_COMM_WORLD);
		        zcopy_(&N1, u_hat, &ONE, &gatherU[k * N1], &ONE);
		        k++;       
            }
	        KPP_Build_FFTMatrix1(N, ncup, localN, gatherU1, tmp4svd_hat);
            KPP_mpi_svd(ncup, N, localN, tmp4svd_hat, S1, B1, nprocs, myrank);
            nPOD1 = KPP_GetPODNumber(S, ncup, gamma);
            /**************assemble the mixed Matrix for SVD*****/
            mixMat = calloc((nPOD+nPOD1)*N1, sizeof(*mixMat));
            for(i=0; i< nPOD; i++){
                zcopy_(&N2, &B[i*N2], &ONE, &mixMat[i*N2], &ONE);
            }
            for(i=0; i< nPOD1; i++){
                zcopy_(&N2, &B1[i*N2], &ONE, &mixMat[(i+nPOD)*N2], &ONE);
            }

            S = calloc(nPOD+nPOD1, sizeof(*S));
	        B = calloc(localN*N*N*N, sizeof(*B)); 

	        KPP_Build_FFTMatrix1(N, nPOD+nPOD1, localN, mixMat, tmp4svd_hat);
            KPP_mpi_svd(nPOD+nPOD1, N, localN, tmp4svd_hat, S, B, nprocs, myrank); 
            nPOD = KPP_GetPODNumber(S, nPOD+nPOD1, gamma);
            BN1 = calloc(nPOD, sizeof(*BN1));
            for(i=0; i< nPOD; i++){
                zcopy_(&N2, &B[i*N2], &ONE, &BN1[i*N2], &ONE);
            }

            /****************update POD Matrix******************/
	        NotTMat=calloc(nPOD*localN*(N/2+1), sizeof(*NotTMat));
	        TMat=calloc(nPOD*localN*(N/2+1), sizeof(*TMat));
            Firm = calloc(nPOD*nPOD, sizeof(*Firm));
            tmpFirm = calloc(nPOD*nPOD, sizeof(*tmpFirm));
	        Varm = calloc(nPOD*nPOD, sizeof(*Varm));
            /****************************************************/
	        KPP_BuildPODMatrix1(nPOD, N, localN, local_0_start, alloc_local, eps, theta, lam, 
                                C, Amp, BN1, NotTMat, TMat, Firm, Varm, cosx, sinx, myrank, 
					            nprocs);
        }
//        /****************************************************************************/
        for(i=0; i< nPOD; i++){
            PODu0[i] = PODu[i];
        }
    }

	 free(gatherU);
	 free(tmp4svd_hat);
	 free(S);
	 free(B);
	 free(CtmpV);
	 free(RtmpV);
	 free(u_hat);
	 free(CPODu);
	 free(PODu);
	 free(PODu0);
}
