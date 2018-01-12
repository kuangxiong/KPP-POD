/*************************************************************************
	> File Name: KPP_BuildPODMatrix.c
	> Author:kuangxiong 
	> Mail:kuangxiong@lsec.cc.ac.cn 
	> Created Time: 2017年12月31日 星期日 22时33分52秒
 ************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include"../fun.h"

void
KPP_BuildPODMatrix(int nPOD, int N, int localN, int local_0_start, int alloc_local, double eps, double theta, 
				   double lam, double C, double Amp, complex *BN1,  complex *Firmx, complex *Varm, double *cosx, 
				   double *sinx, int myrank, int nprocs)
{
	int i, j, k, index[localN], N1, alpha, descA[9], ONE, ZERO, tmpN, tmpN1;
	complex *BN, *BNx, *BNy, *BNlap, *CtmpV, *tmphat1, value;
	double *RBNx, *RBNy, *RBN, *RtmpV, *Bx, *By, *RBBN, *RBBNx, *RBBNy;
	fftw_plan pford, pback;

	ONE = 1;
	ZERO = 0;
//	tmpN = nPOD;
//	tmpN1 = localN*(N/2+1);

	descA[0] = 1; descA[1] = 0; descA[2] = N*N; descA[3] = 1; descA[4]=1;
	descA[5] = 1; descA[6] = 0; descA[7] =0; descA[8] = localN*N; 

	for(i=0; i< localN; i++){
		index[i] = local_0_start + i;
	}

	for(i=0; i<nPOD*nPOD; i++){
		Firmx[i] = 0;
		Varm[i] = 0;
	}
	alpha = 1.0/(N*N);
	
	N1 = N*N;

	CtmpV = calloc(alloc_local, sizeof(*CtmpV));
	RtmpV = calloc(2 * alloc_local, sizeof(*RtmpV));
	pford = fftw_mpi_plan_dft_r2c_2d(N, N, RtmpV, CtmpV, MPI_COMM_WORLD, FFTW_MEASURE);
	pback = fftw_mpi_plan_dft_c2r_2d(N, N, CtmpV, RtmpV, MPI_COMM_WORLD, FFTW_MEASURE);

	Bx = calloc(localN *N, sizeof(*Bx));
	By = calloc(localN *N, sizeof(*By));
	BN = calloc(localN*(N/2+1)*nPOD, sizeof(*BN));
	BNx = calloc(localN*(N/2+1)*nPOD, sizeof(*BN));
	BNy = calloc(localN*(N/2+1)*nPOD, sizeof(*BN));
	BNlap = calloc(localN*(N/2+1)*nPOD, sizeof(*BNlap));
	RBN = calloc(localN *N * nPOD, sizeof(*RBN));
	RBBN = calloc(localN *N * nPOD, sizeof(*RBBN));
	RBNx = calloc(localN *N * nPOD, sizeof(*RBNx));
	RBNy = calloc(localN *N * nPOD, sizeof(*RBNy));
	RBBNx = calloc(localN *N * nPOD, sizeof(*RBBNx));
	RBBNy = calloc(localN *N * nPOD, sizeof(*RBBNy));
	tmphat1 = calloc(localN*N*nPOD, sizeof(*tmphat1));

	for(k=0; k< nPOD; k++)
		for(i=0; i< localN; i++)
			for(j=0; j< N/2+1; j++)
				BN[k*(localN*(N/2+1))+i*(N/2+1)+j]=BN1[k*(localN*N)+i*N+j]; 
/******************BNlap is the laplacian of BN************/
	for(k=0; k< nPOD; k++){
		for(i=0; i< localN; i++){
			for(j=0; j< N/2+1; j++){
				if(index[i] < N/2+1){
					BNlap[k*(localN*(N/2+1))+i*(N/2+1)+j] = -4*PI*PI*(i*i+j*j)
						*BN[k*(localN*(N/2+1))+i*(N/2+1)+j];
				}
				else{
					BNlap[k*(localN*(N/2+1))+i*(N/2+1)+j] = -4*PI*PI*((index[i]-N)*(index[i]-N)+j*j) 
						*BN[k*(localN*(N/2+1))+i*(N/2+1)+j];
				}
			}
		}
	}
//    tmpN = N*N;
//	pzdotc_(&tmpN, &value, &BN1[2*(localN*N)], &ONE, &ONE, descA,  &ONE, &BN1[1*(localN*N)], 
//			&ONE, &ONE, descA, &ONE);	
//	printf("hello:hello:%f&%f\n", value);
/******************copy and symmetrize(different with zliu)*******************/
    KPP_Build_FFTMatrix1(N, nPOD, localN, BNlap, tmphat1);
//	for(k=0;k<nPOD; k++)
//		for(i = 0; i< localN; i++)
//			for(j=0; j< N; j++)
//				tmphat1[k*(localN*N)+i*N+j] = conj(tmphat1[k*(localN*N)+i*N+j]);
/********************compute BN1*BNlap*******************************/
	tmpN = N*N;
	for(i=0; i< nPOD; i++){
		for(j = 0; j< nPOD; j++){
			pzdotc_(&tmpN, &value, &BN1[i*(localN*N)], &ONE, &ONE, descA,  &ONE, &tmphat1[j*(localN*N)], 
					&ONE, &ONE, descA, &ONE);	
			Firmx[i*nPOD+j] = Firmx[i*nPOD+j] + eps*value;
		}
	}
//	if(myrank==0){
//		for(i=0; i< localN*N; i++){
//				printf("hehe:%e&%e\t %e&%e\n", BN1[i], tmphat1[i]);
//			}
//			printf("\n");
//		}

/**************BNx is the derivatives of BN to x*********/	
	for(k=0; k< nPOD; k++){
		for(i=0; i< localN; i++){
			for(j=0; j< N/2; j++){
					BNx[k*(localN*(N/2+1))+i*(N/2+1)+j] = 2*PI*j*I
						*BN[k*(localN*(N/2+1))+i*(N/2+1)+j];
				}
					BNx[k*(localN*(N/2+1))+i*(N/2+1)+N/2] = 0; 
			}
		}
/******************copy and symmetrize(different with zliu)*******************/
    KPP_Build_FFTMatrix1(N, nPOD, localN, BNx, tmphat1);
//    if(myrank==0){
//	for(i=0; i< localN; i++)
//		for(j=0; j< N; j++)
//			printf("%d\t%f&%f\n",myrank, tmphat1[localN*N+i*N+j]);
//	}
//	MPI_Barrier(MPI_COMM_WORLD);
//    if(myrank==1){
//	for(i=0; i< localN; i++)
//		for(j=0; j< N; j++)
//			printf("%d\t%f&%f\n",myrank, tmphat1[localN*N+i*N+j]);
//	}
//	MPI_Barrier(MPI_COMM_WORLD);
/*******************compute eps*BN*BNlap+2*eps*lam*BN*BNx************************/
	tmpN = N*N;
	for(i=0; i< nPOD; i++){
		for(j = 0; j< nPOD; j++){
			pzdotc_(&tmpN, &value, &BN1[i*(localN*N)], &ONE, &ONE, descA,  &ONE, &tmphat1[j*(localN*N)], 
					&ONE, &ONE, descA, &ONE);	
			Firmx[i*nPOD+j] = Firmx[i*nPOD+j] + 2*lam*eps*value;
		}
	}
/*******************compute eps*BN*BNlap+2*eps*lam*BN*BNx + CI*************************/
	for(i=0; i<nPOD; i++){
		Firmx[i*nPOD+i] = Firmx[i*nPOD+i] + C;
	}
	for(i=0; i< localN; i++){
		for(j=0; j< N; j++){
			Bx[i*N+j] = cosx[j];
			By[i*N+j] = cosx[index[i]];
		}
	}
/*******************compute lam (\psi, cos(2pi y)\psi)****************************/
///****************BNy is the derivatives of BN to y*********/	
	for(k=0; k< nPOD; k++){
		for(i=0; i< localN; i++){
			for(j=0; j< N/2+1; j++){
				if(index[i]< N/2)
					BNy[k*(localN*(N/2+1))+i*(N/2+1)+j] = 2*PI*index[i]*I
						*BN[k*(localN*(N/2+1))+i*(N/2+1)+j];
				else if(index[i] == N/2)
					BNy[k*(localN*(N/2+1))+i*(N/2+1)+j] = 0;

				else
					BNy[k*(localN*(N/2+1))+i*(N/2+1)+j] = 2*PI*(index[i]-N)*I
						*BN[k*(localN*(N/2+1))+i*(N/2+1)+j];
				}
			}
		}
/**********************RBBNx is the derivatives of BN to x in real space*************/
	tmpN = localN*(N/2+1);
	for(k=0; k< nPOD; k++){
	/****************RBN is BN in real space*************/
		zcopy_(&tmpN, &BN[k*localN*(N/2+1)], &ONE, CtmpV, &ONE);
		fftw_execute(pback);
		for(i=0; i< localN; i++)
			for(j=0; j< N; j++)
				RBN[k*(localN*N)+i*N+j] = RtmpV[i*2*(N/2+1)+j]/(N*N);		
    /*************compute RBNx*******************/		
		zcopy_(&tmpN, &BNx[k*localN*(N/2+1)], &ONE, CtmpV, &ONE);
//		fftw_execute_dft_c2r(pback, CtmpV, RtmpV);
		fftw_execute(pback);
		for(i=0; i< localN; i++)
			for(j=0; j< N; j++)
				RBNx[k*(localN*N)+i*N+j] = RtmpV[i*2*(N/2+1)+j]/(N*N);
	/************compute RBNy****************/
		zcopy_(&tmpN, &BNy[k*localN*(N/2+1)], &ONE, CtmpV, &ONE);
//		fftw_execute_dft_c2r(pback, CtmpV, RtmpV);
		fftw_execute(pback);
		for(i=0; i< localN; i++)
			for(j=0; j< N; j++)
				RBNy[k*(localN*N)+i*N+j] = RtmpV[i*2*(N/2+1)+j]/(N*N);
	}

	for(k=0;k<nPOD; k++){
		for(i=0; i< localN; i++){
			for(j=0; j< N; j++){
				RBBNx[k*(localN*N)+i*N+j] = By[i*N+j]*RBNx[k*(localN*N)+i*N+j];
				RBBNy[k*(localN*N)+i*N+j] = Bx[i*N+j]*RBNy[k*(localN*N)+i*N+j];
			    RBBN[k*(localN*N)+i*N+j] = RBBNx[k*(localN*N)+i*N+j] + RBBNy[k*(localN*N)+i*N+j]
										   + lam*By[i*N+j]*RBN[k*(localN*N)+i*N+j];
			}
		}
	}
/****************transform real space RBBN to K space **************************/
	tmpN = localN * N;
	for(k=0;k<nPOD; k++){
		for(i=0; i< localN; i++)
			for(j=0; j<N; j++)
				RtmpV[i*2*(N/2+1)+j] = RBBN[k*localN*N+i*N+j];
//		fftw_execute_dft_r2c(pback, RtmpV, CtmpV);
		fftw_execute(pford);
		for(i=0; i< localN; i++){
			for(j=0; j< N/2+1; j++){
				tmphat1[k*(localN*N)+i*N+j] = CtmpV[i*(N/2+1)+j];
			}
			for(j=N/2+1; j< N; j++){
				tmphat1[k*(localN*N)+i*N+j] = conj(CtmpV[i*(N/2+1)+N-j]);
			}
		}
	}	
	tmpN = N*N;
	for(i=0; i< nPOD; i++){
		for(j = 0; j< nPOD; j++){
			pzdotc_(&tmpN, &value, &BN1[i*(localN*N)], &ONE, &ONE, descA,  &ONE, &tmphat1[j*(localN*N)], 
					&ONE, &ONE, descA, &ONE);	
			Firmx[i*nPOD+j] = Firmx[i*nPOD+j] + value;
		}
	}
//	if(myrank==0){
//		printf("test:%e&%e\n", value);
//		for(i=0; i< nPOD; i++){
//			for(j=0; j< nPOD; j++){
//				printf("%f&%f\t", Firmx[i*nPOD+j]);
//			}
//			printf("\n");
//		}
//	}
/***************************the part of t*******************************/
	MPI_Barrier(MPI_COMM_WORLD);
	for(i=0; i< localN; i++){
		for(j=0; j< N; j++){
			Bx[i*N+j] = Amp * theta * sinx[j];
			By[i*N+j] = Amp * theta * sinx[index[i]];
		}
	}
	for(k=0;k<nPOD; k++){
		for(i=0; i< localN; i++){
			for(j=0; j< N; j++){
				RBBNx[k*(localN*N)+i*N+j] = By[i*N+j] * RBNx[k*(localN*N)+i*N+j];
				RBBNy[k*(localN*N)+i*N+j] = Bx[i*N+j] * RBNy[k*(localN*N)+i*N+j];
			    RBBN[k*(localN*N)+i*N+j] = RBBNx[k*(localN*N)+i*N+j] + RBBNy[k*(localN*N)+i*N+j]
										 +lam * By[i*N+j] * RBN[k*(localN*N)+i*N+j];
			}
		}
	}
//	if(myrank==0){
//		printf("111111111alpha:%f&%f\n", lam, 1.0/(N*N));
//		for(i=0; i< localN; i++){
//			for(j=0; j< N; j++)
//				printf("%f\t", RBBN[i*N+j]);
//			printf("\n");	
//		}
//	}
    
	for(k=0;k<nPOD; k++){
		for(i=0; i< localN; i++)
			for(j=0; j<N; j++)
				RtmpV[i*2*(N/2+1)+j] = RBBN[k*localN*N+i*N+j];

		fftw_execute(pford);

		for(i=0; i< localN; i++){
			for(j=0; j< N/2+1; j++){
				tmphat1[k*(localN*N)+i*N+j] = CtmpV[i*(N/2+1)+j];
			}
			for(j=N/2+1; j< N; j++){
				tmphat1[k*(localN*N)+i*N+j] = conj(CtmpV[i*(N/2+1)+N-j]);
			}
		}
	}
	
	
	tmpN = N * N;
	for(i=0; i< nPOD; i++){
		for(j = 0; j< nPOD; j++){
			pzdotc_(&tmpN, &value, &BN1[i*(localN*N)], &ONE, &ONE, descA,  &ONE, &tmphat1[j*(localN*N)], 
					&ONE, &ONE, descA, &ONE);	
			Varm[i*nPOD+j] = Varm[i*nPOD+j] + value;
		}
	}
//    if(myrank==0){
//	for(i=0; i< localN; i++)
//		for(j=0; j< N; j++)
//			printf("%d\t%f&%f\n",myrank, RBBN[i*N+j]);
//	}
//	MPI_Barrier(MPI_COMM_WORLD);
//    if(myrank==1){
//	for(i=0; i< localN; i++)
//		for(j=0; j< N; j++)
//			printf("%d\t%f&%f\n",myrank, RBBN[i*N+j]);
//	}
//	MPI_Barrier(MPI_COMM_WORLD);
//	if(myrank==0){
//		for(i=0; i< nPOD; i++){
//			for(j=0; j< nPOD; j++){
//				printf("%f&%f\t", Firmx[i*nPOD+j]);
//			}
//			printf("\n");
//		}
//	}

#if 0
#endif
/****************transform real space RBBN to K space **************************/
	free(Bx);
	free(By);
	free(RBBN);
	free(RBBNx);
	free(RBBNy);
	free(BN);
	free(BNx);
	free(BNy);
	free(RBNx);
	free(RBNy);
	free(BNlap);
	free(tmphat1);
}
