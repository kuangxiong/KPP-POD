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
void KPP_APOD_test(int N, DOUBLE lam, DOUBLE eps, DOUBLE tau, DOUBLE theta, DOUBLE dt, DOUBLE t, 
			  DOUBLE Amp, DOUBLE cM, DOUBLE C, DOUBLE y1, int myrank, int nprocs, 
			  MPI_Comm comm1)
{     
     long int alloc_local, localN, local_0_start, ncup;
     double *u, cost, *S, *S1, *SS1,  cosx[N], sinx[N], h, *RtmpV, *uPOD, 
            errind, tmpt, err, localerr, Norm, localNorm, zero;
	 
	 complex *u_hat, *UP_hat, *gatherU, *gatherU1, testu_hat[N*(N/2+1)], *tmp4svd_hat, *tmp4svd_hat1,
			 *B, *B1,*BB1,  *BN1, *BN, *PODu0, *Firm, *Varm, *tmpFirm, *PODu, *CPODu, *preu_hat,
			 *CtmpV, *BNlap, *TMat, *NotTMat, alpha, tmpvalue, beta, *mixMat, *err_hat;
     fftw_plan p1, p2;
     int i, k, j, niter = 0, ONE=1, N1, N2, nPOD, nPOD1, nPOD2,  tmpN, tmpN1, flag, tmpnPOD;
	 FILE *fp, *fp1;

	 fp = fopen("APODerr.text","w");
	 fp1 = fopen("APODerrind.text","w");
     alloc_local = fftw_mpi_local_size_2d(N, N/2+1, comm1, &localN, &local_0_start);
     u_hat = calloc(localN*(N/2+1), sizeof(*u_hat));
     preu_hat = calloc(localN*(N/2+1), sizeof(*preu_hat));
     u = calloc(localN * N, sizeof(*u));
     uPOD = calloc(localN * N, sizeof(*uPOD));
     
	 CtmpV = calloc(alloc_local, sizeof(*CtmpV));
     RtmpV = calloc(2 * alloc_local, sizeof(*RtmpV));
//	 testu_hat = calloc(N*(N/2+1), sizeof(complex));
     err_hat = calloc(localN*(N/2+1), sizeof(*err_hat));
	 CPODu = calloc(localN*(N/2+1), sizeof(*CPODu));
     N1 = localN * (N/2+1);
	 N2 = localN * N;
     zero = 0.0; 
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
     if(ncup%nprocs!=0)
     {
         if(myrank==0){
             printf("ncup:%d\n", ncup);
             printf("error: ncup should divided by nprocs!!!!\n");
         }
         return;
     
     }
     gatherU = calloc(N1*ncup, sizeof(*gatherU));
     tmp4svd_hat = calloc(localN*N*ncup, sizeof(*tmp4svd_hat));
     tmp4svd_hat1 = calloc(localN*N*ncup, sizeof(*tmp4svd_hat1));
     S = calloc(ncup, sizeof(*S));
	 B = calloc(localN*N*N*N, sizeof(*B));
     SS1 = calloc(ncup, sizeof(*SS1));
	 BB1 = calloc(localN*N*N*N, sizeof(*BB1));
    
	 zcopy_(&N1, u_hat, &ONE, gatherU, &ONE);

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
     }
//   MPI_Gather(u_hat, localN*(N/2+1), MPI_C_DOUBLE_COMPLEX, testu_hat , localN*(N/2+1), MPI_C_DOUBLE_COMPLEX, 0, comm1);
	 KPP_Build_FFTMatrix1(N, ncup, localN, gatherU, tmp4svd_hat, comm1);
	 KPP_Build_FFTMatrix1(N, ncup, localN, gatherU, tmp4svd_hat1, comm1);
//     for(i=0; i< localN*N; i++)
//     if(myrank==0)
//        for(i=0; i< localN; i++)
//            for(j=1; j< N/2+1; j++)
//             printf("%e&%e\t%e&%e\n", tmp4svd_hat[localN*N + i*N+j], tmp4svd_hat[localN*N+i*N+N-j]);

     KPP_mpi_svd2(ncup, N, localN, tmp4svd_hat, S, B, nprocs, myrank);
     KPP_mpi_ev1(ncup, N, localN, tmp4svd_hat1, SS1, BB1, nprocs, myrank);
	 k=0;
	 if(myrank ==0){
         printf("\nmethod:  \tSVD \t\t\t EV\n");
		 for(i=0; i< ncup; i++)
			printf("singular value:%f\t%f\n", S[i], SS1[i]);
     }
	nPOD = KPP_GetPODNumber(S, ncup, gamma1);
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
//    KPP_GetInitialPOD(N, localN, nPOD, u_hat, BN1, PODu0, myrank, nprocs);
//	if(myrank==0)
//		for(i=0; i< nPOD; i++)
//			printf("PODu0:%f&%f\n", PODu0[i]);
//
//    if(myrank==0)
//        for(i=0; i< localN; i++)
//            for(j=1; j< N/2+1; j++)
//                printf("%e&%e\t%e&%e \t%e&%e\t%e&%e\n", B[i*N+j], B[i*N+N-j], BB1[i*N+j], BB1[i*N+N-j]);
///*************defined matrix for building POD matrix************/
//	NotTMat=calloc(nPOD*localN*(N/2+1), sizeof(*NotTMat));
//	TMat=calloc(nPOD*localN*(N/2+1), sizeof(*TMat));
//    Firm = calloc(nPOD*nPOD, sizeof(*Firm));
//    tmpFirm = calloc(nPOD*nPOD, sizeof(*tmpFirm));
//	Varm = calloc(nPOD*nPOD, sizeof(*Varm));
////	KPP_BuildPODMatrix(nPOD, N, localN, local_0_start, alloc_local, eps,
////					   theta, lam, C, Amp, BN1, Firm, Varm, cosx, sinx, myrank, 
////					   nprocs);
//	KPP_BuildPODMatrix1(nPOD, N, localN, local_0_start, alloc_local, eps,
//					   theta, lam, C, Amp, BN1, NotTMat, TMat, Firm, Varm, cosx, sinx, myrank, 
//					   nprocs, comm1);
//	tmpN = localN*(N/2+1);
//	tmpN1 = 2*localN*(N/2+1);
//    alpha = dt;
//    beta = 0.0;
#if 0
//	while(t <= T) 	
    {
    	for(i=0; i< nPOD; i++)
    		for(j=0; j< nPOD; j++)
    			tmpFirm[i*nPOD+j] = Firm[i*nPOD+j] + cos(t)*Varm[i*nPOD+j];
        cost = cos(t);
        t = t + dt;
    	/*******************update PODu***************/
    	zgemv_("C", &nPOD, &nPOD, &alpha, tmpFirm, &nPOD, PODu0, &ONE, &beta, PODu, &ONE);    
    	for(i=0; i< nPOD; i++){
    		PODu[i] = PODu0[i] + PODu[i];
    	}
        /********************computing error of POD solution*******************/
		/***transform POD sulution to  plane solution******/
        tmpvalue = 1.0;
//      tmpN = localN*(N/2+1);
        
        zgemv_("N", &tmpN, &nPOD, &tmpvalue, BN, &tmpN, PODu, &ONE, &beta, CPODu, &ONE);
        /**************************************************/
    	zcopy_(&N1, u_hat, &ONE, preu_hat, &ONE);        
        KPP_ComputePlane(N, t, dt, C, Amp, lam, theta, eps, cost, u_hat, myrank, nprocs,localN, local_0_start,  comm1);
        for(i=0; i< localN; i++){
            for(j=0; j< N/2+1; j++){
                err_hat[i*(N/2+1)+j] = CPODu[i*(N/2+1)+j] - u_hat[i*(N/2+1)+j];   
            }
        }
        localerr = dnrm2_(&tmpN1, err_hat, &ONE);
        localerr = localerr * localerr;
        localNorm = dnrm2_(&tmpN1, u_hat, &ONE);
        localNorm = localNorm * localNorm;

        MPI_Allreduce(&localerr, &err, ONE, MPI_DOUBLE, MPI_SUM, comm1);
        MPI_Allreduce(&localNorm, &Norm, ONE, MPI_DOUBLE, MPI_SUM, comm1);
        if(myrank==0)
            fprintf(fp, "%f\t%f\n", t, sqrt(err/Norm));
        /*********************************************************************/
//		fftw_execute(p2);
//      for(i=0; i< localN; i++)
//	        for(j=0; j< N; j++)
//				uPOD[i*N+j] = RtmpV[i*2*localN*(N/2+1)+j]/(N*N);  
		/***************************computing err indicator*********************/
        errind =KPP_GetErrIndicator(nPOD, localN, N, cost, dt, PODu, PODu0, BN1, NotTMat, TMat, comm1);
//      /***************************Update POD Matrix**************************/
//        if(myrank==0)
//            printf("111111111111:%f\n", t);
//        tmpN = localN*(N/2+1);
//        tmpN1 = 2*localN*(N/2+1);
//            printf("%f\t%f\n", t, errind);

//        if(errind > errflag)
        if(1)
        { 
            
            tmpvalue = 1.0;
            // transfer POD solution to PW solution 
            t = t -dt;
            zgemv_("N", &tmpN, &nPOD, &tmpvalue, BN, &tmpN, PODu0, &ONE, &beta, CPODu, &ONE);
            /*******************free Matrix***************************/
            free(gatherU); free(tmp4svd_hat); free(S); free(B); free(TMat);
            free(NotTMat); free(Firm); free(Varm); free(tmpFirm); free(PODu0);
            free(PODu); free(BN); 
            /*********************************************************/
    		zcopy_(&N1, preu_hat, &ONE, u_hat, &ONE);
            MPI_Barrier(comm1);

            ncup = initer/intval +1;

            tmp4svd_hat = calloc(localN*N*ncup, sizeof(*tmp4svd_hat));
            tmp4svd_hat1 = calloc(localN*N*ncup, sizeof(*tmp4svd_hat1));
            gatherU = calloc(ncup*localN*(N/2+1), sizeof(*gatherU));
	        zcopy_(&N1, CPODu, &ONE, gatherU, &ONE);
            k = 1;

            /***************collecting PW solution fot updating POD Matrix*****************/
//            if(myrank==0)
//                printf("11111111:%f\t%f\n", t, errind);
            while(k <= initer ){
                cost = cos(t);
                t = t + dt;
                if(t > T){
                    printf("can i get here\n");
                    goto stop;
                }
//                if(myrank==0)
//                printf("111111111111:%f\n", t);
                                KPP_ComputePlane(N, t, dt, C, Amp, lam, theta, eps, cost, CPODu, myrank, nprocs,
                                localN, local_0_start,  comm1);
                if(k%intval==0)
    		        zcopy_(&N1, CPODu, &ONE, &gatherU[(k/intval) * N1], &ONE);
                /*******************computing error*******************************/

                KPP_ComputePlane(N, t, dt, C, Amp, lam, theta, eps, cost, u_hat, myrank, nprocs,
                                localN, local_0_start,  comm1);
                for(i=0; i< localN; i++){
                    for(j=0; j< N/2+1; j++){
                        err_hat[i*(N/2+1)+j] = CPODu[i*(N/2+1)+j] - u_hat[i*(N/2+1)+j];   
                    }
                }
                localerr = dnrm2_(&tmpN1, err_hat, &ONE);
                localerr = localerr * localerr;
                localNorm = dnrm2_(&tmpN1, u_hat, &ONE);
                localNorm = localNorm * localNorm;
                MPI_Allreduce(&localerr, &err, ONE, MPI_DOUBLE, MPI_SUM, comm1);
                MPI_Allreduce(&localNorm, &Norm, ONE, MPI_DOUBLE, MPI_SUM, comm1);
                if(myrank==0){
                    fprintf(fp, "%f\t%e\n", t, sqrt(err/Norm));
                    fprintf(fp1, "%f\t%e\n", t, zero); 
                }
              /*********************************************************************/
		        k++;                   
            }
//          printf("hello , i am here:%d\t%d\n", ncup, k/intval);
            S1 = calloc(ncup, sizeof(*S1));
	        B1 = calloc(localN*N*N*N, sizeof(*B1));
            SS1 = calloc(ncup, sizeof(*SS1));
	        BB1 = calloc(localN*N*N*N, sizeof(*BB1));
	        KPP_Build_FFTMatrix1(N, ncup, localN, gatherU, tmp4svd_hat, comm1);
	        KPP_Build_FFTMatrix1(N, ncup, localN, gatherU, tmp4svd_hat1, comm1);
            
            MPI_Barrier(comm1);
            
            KPP_mpi_svd(ncup, N, localN, tmp4svd_hat, S1, B1, nprocs, myrank);
            KPP_mpi_ev(ncup, N, localN, tmp4svd_hat1, SS1, BB1, nprocs, myrank);
            if(myrank==0)
                for(k=0; k< 6; k++)
                    for(i=0; i< localN; i++)
                        for(j=1; j< N/2+1; j++)
//                            if(fabs(B1[k*localN*N+i*N+j]-BB1[(ncup-1-k)*localN*N+i*N+j])>0.00000001)
//                                printf("error:%d\t %f&%f\t%f&%f\t%e&%e\n",k, BB1[(ncup-1-k)*localN*N+i*N+j], B1[k*localN*N+i*N+j], BB1[(ncup-1-k)*localN*N+i*N+j]-B1[k*localN*N+i*N+j]);
                            if(cimag(B1[k*localN*N+i*N+j]+B1[k*localN*N+i*N+N-j])>0.00000001)
                                printf("error:%d\t %e&%e\t%e&%e\t%e&%e\n",k, B1[k*localN*N+i*N+j], B1[k*localN*N + i*N + N-j], B1[k*localN*N+i*N+j]-B1[k*localN*N+i*N + N-j]);

            nPOD1 = KPP_GetPODNumber(S1, ncup, gamma2);
            if(myrank==0)
                for(i=0; i< ncup; i++)
                    printf("%d\t%f\t%f\n",ncup, S1[i], SS1[i]);
            /**************assemble the mixed Matrix for SVD*****/
            mixMat = calloc((nPOD+nPOD1)*localN*N, sizeof(*mixMat));

            for(i=0; i< nPOD; i++){
                zcopy_(&N2, &BN1[i*N2], &ONE, &mixMat[i*N2], &ONE);
//              zcopy_(&N2, &BN1[i*N2], &ONE, &tmp4svd_hat1[i*N2], &ONE);
            }
            for(i=0; i< nPOD1; i++){
                zcopy_(&N2, &B1[i*N2], &ONE, &mixMat[(i+nPOD)*N2], &ONE);
//              zcopy_(&N2, &B1[i*N2], &ONE, &tmp4svd_hat1[(i+nPOD)*N2], &ONE);
            }

            S = calloc(nPOD+nPOD1, sizeof(*S));
//          S = calloc(nPOD1, sizeof(*S));
	        B = calloc(localN*N*(nPOD+nPOD1), sizeof(*B)); 
//	        B = calloc(localN*N*N*N, sizeof(*B)); 
            if(myrank==0){
                printf("number: %d\t%d\n", nPOD, nPOD1);
     //           for(i=0; i< localN; i++)
     //               for(j=0; j< N/2-1; j++)
     //                   printf("11111:%f&%f\t%f&%f\n", mixMat[1*N2+i*N+N/2-j], mixMat[1*N2+i*N+N/2+j]);
            }

            tmpnPOD = nPOD+nPOD1;
//          KPP_mpi_svd(nPOD+nPOD1, N, localN, tmp4svd_hat1, SS1, BB1, nprocs, myrank); 
            KPP_mpi_ev(nPOD+nPOD1, N, localN, mixMat, S, B, nprocs, myrank); 
//          KPP_mpi_ev(nPOD1, N, localN, B1, S, B, nprocs, myrank); 
//          KPP_mpi_svd(nPOD1, N, localN, B1, S, B, nprocs, myrank); 
            printf("nPOD:%d\n", nPOD+nPOD1);
            nPOD = KPP_GetPODNumber(S, nPOD+nPOD1, gamma3);
//          nPOD = KPP_GetPODNumber(S, nPOD1, gamma1);
//            if(myrank==0){
//                printf("number: %d\t%d\n", nPOD, nPOD1);
//                for(i=0; i< localN; i++)
//                    for(j=0; j< N/2-1; j++)
//                        printf("hello:%f&%f\t%f&%f\n", B[1*N2+i*N+N/2-j], B[1*N2+i*N+N/2+j]);
//            }
//          
//            if(myrank==0)
//              for(i=0; i< 9; i++)
//                  printf("Ssss:%d\t%f\t%f\n",nPOD,  S[i], SS1[i]);
            if(myrank==0)
            {
                for(i=0; i< localN; i++)
                    for(j=1; j< N/2; j++)
                        if(fabs(B[i*N+j]-B[i*N+N-j])>0.00001)
                            printf("error:%f&%f\t%f&%f, %e\n", B[i*N+j], B[i*N+N-j], fabs(B[i*N+j]-B[i*N+N-j]));
            }
//           if(myrank==0)
//              for(i=0; i< nPOD; i++)
//                  printf("Ssss:%d\t%f\n",nPOD,  S[i]);
           free(S1);  free(B1); free(BN1); free(mixMat);
////
           BN1 = calloc(nPOD*localN*N, sizeof(*BN1));
           BN = calloc(nPOD*localN*(N/2+1), sizeof(*BN));
           PODu0 = calloc(nPOD, sizeof(*PODu0));
           PODu = calloc(nPOD, sizeof(*PODu));
           for(i=0; i< nPOD; i++){
               zcopy_(&N2, &B[(tmpnPOD-1-i)*N2], &ONE, &BN1[i*N2], &ONE);
//             zcopy_(&N2, &B[i*N2], &ONE, &BN1[i*N2], &ONE);
           }
           for(k=0; k< nPOD; k++){
               for(i=0; i< localN; i++){
                   for(j=0; j< N/2+1; j++){
                       BN[k*localN*(N/2+1)+i*(N/2+1)+j] = BN1[k*localN*N+i*N+j];
                   }
               }
           }
////         /****************update POD Matrix******************/
	        NotTMat=calloc(nPOD*localN*(N/2+1), sizeof(*NotTMat));
	        TMat=calloc(nPOD*localN*(N/2+1), sizeof(*TMat));
            Firm = calloc(nPOD*nPOD, sizeof(*Firm));
            tmpFirm = calloc(nPOD*nPOD, sizeof(*tmpFirm));
	        Varm = calloc(nPOD*nPOD, sizeof(*Varm));
            /****************************************************/
	        KPP_BuildPODMatrix1(nPOD, N, localN, local_0_start, alloc_local, eps, theta, lam, 
                                C, Amp, BN1, NotTMat, TMat, Firm, Varm, cosx, sinx, myrank, 
					            nprocs, comm1);
            KPP_GetInitialPOD(N, localN, nPOD, CPODu, BN1, PODu, myrank, nprocs);  
            /************************************************************************************
             * Attention!!! , PODu is double vector. however, the imaginary part  of PODu isn't 0, 
             * it's a small number, so we set imag of PODu is 0
             * **********************************************************************************/
            for(i=0; i< nPOD; i++)
                PODu[i]=creal(PODu[i]);
///***********************************test*****************************************/
        /******************test initial err*****************/
            if(myrank==0)
                for(i=0; i< nPOD; i++)
                    printf("PODu:%f&%f\n", PODu[i]);
            tmpvalue = 1.0;
//          transfer POD solution to PW solution 
            zgemv_("N", &tmpN, &nPOD, &tmpvalue, BN, &tmpN, PODu, &ONE, &beta, err_hat, &ONE);
                for(i=0; i< localN; i++){
                    for(j=0; j< N/2+1; j++){
                        err_hat[i*(N/2+1)+j] = CPODu[i*(N/2+1)+j] - err_hat[i*(N/2+1)+j];   
                    }
                }
                localerr = dnrm2_(&tmpN1, err_hat, &ONE);
                localerr = localerr * localerr;
                localNorm = dnrm2_(&tmpN1, CPODu, &ONE);
                localNorm = localNorm * localNorm;
                MPI_Allreduce(&localerr, &err, ONE, MPI_DOUBLE, MPI_SUM, comm1);
                MPI_Allreduce(&localNorm, &Norm, ONE, MPI_DOUBLE, MPI_SUM, comm1);
                if(myrank==0)
                    printf("initial err:%f\t%e\n", t, sqrt(err));
////  /*********************************************************************/
//        for(i=0; i< nPOD; i++){
//            PODu0[i] = PODu[i];
//        }
//
//    	for(i=0; i< nPOD; i++)
//    		for(j=0; j< nPOD; j++)
//    			tmpFirm[i*nPOD+j] = Firm[i*nPOD+j] + cos(t)*Varm[i*nPOD+j];
//        cost = cos(t);
//        t = t + dt;
//    	zgemv_("C", &nPOD, &nPOD, &alpha, tmpFirm, &nPOD, PODu0, &ONE, &beta, PODu, &ONE);    
//    	for(i=0; i< nPOD; i++){
//    		PODu[i] = PODu0[i] + PODu[i];
//    	}
//        errind =KPP_GetErrIndicator(nPOD, localN, N, cost, dt, PODu, PODu0, BN1, NotTMat, TMat, comm1);
//        if(myrank==0)
//        printf("errind:%f\t%e\n",t, errind);
//
//        goto stop1;
//#if 0
//#endif
        }
        
        else{   
        if(myrank==0)
            fprintf(fp1, "%f\t%f\n", t, errind);
        }
      /****************************************************************************/

        for(i=0; i< nPOD; i++){
            PODu0[i] = PODu[i];
//            if(myrank==0)
//                printf("%d\t%f&%f\n",nPOD,  PODu0[i]);
        }
     }
//     stop1:
#endif
	 free(S);
	 free(B);
	 free(CtmpV);
	 free(RtmpV);
	 free(u_hat);
	 free(CPODu);
	 free(PODu);
	 free(PODu0);
     stop:
	 free(gatherU);
	 free(tmp4svd_hat);
}
