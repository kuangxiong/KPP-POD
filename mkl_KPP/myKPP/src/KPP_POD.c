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
void KPP_POD(int N, DOUBLE lam, DOUBLE eps, DOUBLE tau, DOUBLE theta, DOUBLE dt, DOUBLE t, 
			  DOUBLE Amp, DOUBLE cM, DOUBLE C, DOUBLE y1, int myrank, int nprocs, MPI_Comm comm  
			  )
{     
     long int alloc_local, localN, local_0_start, ncup;
     double *u, cost, *S, *S1, cosx[N], sinx[N], h, *RtmpV, *uPOD, 
             errind, localerr, err, localNorm, Norm;
	 
	 complex *u_hat, *gatherU, testu_hat[N*(N/2+1)], *tmp4svd_hat,
			 *B, *B1, *BN1, *BN, *PODu0, *Firm, *Varm, *tmpFirm, *PODu, *CPODu,
			 *CtmpV, alpha, tmpvalue, beta, *err_hat, *NotTMat, *TMat;
     fftw_plan p1, p2;
     int i, k, j, niter = 0, ONE=1, N1, N2, nPOD, nPOD1, nPOD2,  tmpN, tmpN1;
	 FILE *fp, *fp1;

	 fp = fopen("PODerr.text","w");
	 fp1 = fopen("PODerrind.text","w");
     alloc_local = fftw_mpi_local_size_2d(N, N/2+1, comm, &localN, &local_0_start);
     u_hat = calloc(localN*(N/2+1), sizeof(*u_hat));
     err_hat = calloc(localN*(N/2+1), sizeof(*err_hat));
     u = calloc(localN * N, sizeof(*u));
     uPOD = calloc(localN * N, sizeof(*uPOD));
	 CPODu = calloc(localN*(N/2+1), sizeof(*CPODu));
     
	 CtmpV = calloc(alloc_local, sizeof(*CtmpV));
     RtmpV = calloc(2 * alloc_local, sizeof(*RtmpV));
//	 testu_hat = calloc(N*(N/2+1), sizeof(complex));
     N1 = localN * (N/2+1);
	 N2 = localN * N;
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
/*************************************************/
	 ncup = (int)(tmax/dt)+2;
     ncup = ncup/intval + 1;
     gatherU = calloc(N1*ncup, sizeof(*gatherU));
     tmp4svd_hat = calloc(localN*N*ncup, sizeof(*gatherU));
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
        
        KPP_ComputePlane(N, t, dt, C, Amp, lam, theta, eps, cost, u_hat, myrank, nprocs,localN, local_0_start,  comm);
		if(k%intval ==0){
            zcopy_(&N1,  u_hat, &ONE, &gatherU[(k/intval) * N1], &ONE);
        }
		k++;
     }
//     MPI_Gather(u_hat, localN*(N/2+1), MPI_C_DOUBLE_COMPLEX, testu_hat , localN*(N/2+1), MPI_C_DOUBLE_COMPLEX, 0, comm);

	 KPP_Build_FFTMatrix1(N, ncup, localN, gatherU, tmp4svd_hat, comm);
     
     KPP_mpi_svd(ncup, N, localN, tmp4svd_hat, S, B, nprocs, myrank);
	 k=0;
//	 if(myrank ==0)
//		 for(i=0; i< ncup; i++)
//			printf("S:%d\t%f\n",ncup,  S[i]);

	nPOD = KPP_GetPODNumber(S, ncup, gamma1);
//	if(myrank==0)
//		printf("%d\t%d\n",ncup, nPOD);
/********************get initial POD solution****************/
	NotTMat = calloc(nPOD*localN*(N/2+1), sizeof(*NotTMat));
	TMat = calloc(nPOD*localN*(N/2+1), sizeof(*TMat));

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
//	if(myrank==0)
//		for(i=0; i< nPOD; i++)
//			printf("PODu0:%f&%f\n", PODu0[i]);
//***************computing initial POD solution err*********/
//    tmpvalue = 1.0;
//    beta = 0;
//    zgemv_("N", &tmpN, &nPOD, &tmpvalue, BN, &tmpN, PODu0, &ONE, &beta, CPODu, &ONE);
//    for(i=0; i< localN; i++) 
//        for(j=0;  j< N/2+1; j++)
//           err_hat[i*(N/2+1)+j] = u_hat[i*(N/2+1)+j]-CPODu[i*(N/2+1)+j];
//    localerr = dnrm2_(&tmpN1, err_hat, &ONE);
//    localerr = localerr*localerr;
//    localNorm = dnrm2_(&tmpN1, u_hat, &ONE);
//    localNorm = localNorm * localNorm;
//    MPI_Allreduce(&localerr, &err, ONE, MPI_DOUBLE, MPI_SUM, comm);
//    MPI_Allreduce(&localNorm, &Norm, ONE, MPI_DOUBLE, MPI_SUM, comm);
//    if(myrank==0)
//       printf("initial POD err:%f\t%e\n", t, sqrt(err/Norm));
/*************defined matrix for building POD matrix************/
    Firm = calloc(nPOD*nPOD, sizeof(*Firm));
    tmpFirm = calloc(nPOD*nPOD, sizeof(*tmpFirm));
	Varm = calloc(nPOD*nPOD, sizeof(*Varm));

//    KPP_BuildPODMatrix(nPOD, N, localN, local_0_start, alloc_local, eps,
//					   theta, lam, C, Amp, BN1, Firm, Varm, cosx, sinx, myrank, 
//					   nprocs, comm);
	KPP_BuildPODMatrix1(nPOD, N, localN, local_0_start, alloc_local, eps,
					   theta, lam, C, Amp, BN1, NotTMat, TMat, Firm, Varm, cosx, sinx, myrank, 
					   nprocs, comm);
	tmpN = localN*(N/2+1);
	tmpN1 = 2*localN*(N/2+1);
    alpha = dt;
    beta = 0.0;
	while(t <= T) 
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
        /***********computing error indicator*********************/
        errind = KPP_GetErrIndicator(nPOD, localN, N, cost, dt, PODu, PODu0, BN1, NotTMat, 
                                    TMat, comm);
        if(myrank==0)
            fprintf(fp1, "%f\t%f\n", t, errind);
	    /*********transform POD sulution to  plane solution********/
        tmpvalue = 1.0;
        beta = 0;
        zgemv_("N", &tmpN, &nPOD, &tmpvalue, BN, &tmpN, PODu, &ONE, &beta, CPODu, &ONE);
        /***************************PW solution********************************/
        KPP_ComputePlane(N, t, dt, C, Amp, lam, theta, eps, cost, u_hat, myrank, nprocs,
                localN, local_0_start,  comm);
        /***************************computing err*****************************/
        for(i=0; i< localN; i++) 
            for(j=0;  j< N/2+1; j++)
               err_hat[i*(N/2+1)+j] = u_hat[i*(N/2+1)+j]-CPODu[i*(N/2+1)+j];
        MPI_Barrier(comm);
        localerr = dnrm2_(&tmpN1, err_hat, &ONE);
        localerr = localerr*localerr;
        localNorm = dnrm2_(&tmpN1, u_hat, &ONE);
        localNorm = localNorm * localNorm;
        MPI_Allreduce(&localerr, &err, ONE, MPI_DOUBLE, MPI_SUM, comm);
        MPI_Allreduce(&localNorm, &Norm, ONE, MPI_DOUBLE, MPI_SUM, comm);
        if(myrank==0)
            fprintf(fp, "%f\t%f\n", t, sqrt(err/Norm));
//		zcopy_(&tmpN, CPODu, &ONE, CtmpV, &ONE);
//		fftw_execute(p2);
//        for(i=0; i< localN; i++)
//			for(j=0; j< N; j++)
//				uPOD[i*N+j] = RtmpV[i*2*localN*(N/2+1)+j]/(N*N);  
//        /****************************************************************************/
        for(i=0; i< nPOD; i++){
            PODu0[i] = PODu[i];
        }
    }
    free(NotTMat);
    free(TMat);
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
