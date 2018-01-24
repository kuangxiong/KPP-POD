/*************************************************************************
	> File Name: KPP_GetInitialPOD.c
	> Author:kuangxiong 
	> Mail:kuangxiong@lsec.cc.ac.cn 
	> Created Time: 2018年01月03日 星期三 09时18分51秒
 ************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include"../fun.h"

void
KPP_GetInitialPOD(int N, int localN, int nPOD, complex *u_hat,  complex *BN1, complex *PODu0, int myrank, int numprocs)
{
	int descA[9], k, i, j, tmpN, ONE=1, ictxt, nprow, npcol;
	complex tmpV[localN*N], tmpu_hat[localN*N], value;
	

	tmpN = N*N;
	descA[0]=1; descA[1]=0; descA[2] = N*N; descA[3] = 1; descA[4] = 1;
	descA[5] = 1; descA[6] = 0; descA[7] = 0; descA[8] = localN*N;

	nprow = numprocs;
	npcol = 1;

	Cblacs_pinfo(&myrank, &numprocs);
	Cblacs_get(-1, 0, &ictxt);
	Cblacs_gridinit(&ictxt, "Row", nprow, npcol);

	for(i=0; i< localN; i++){
		for(j=0; j< N; j++){
			if(j<N/2+1){
				tmpu_hat[i*N+j] = u_hat[i*(N/2+1)+j]; 
			}
			else{
				tmpu_hat[i*N+j] = conj(u_hat[i*(N/2+1)+N-j]);		
			}
		}
	}
	
/***********************************************/	
	k = 0;
	for(k=0; k< nPOD; k++)
	{
		for(i=0; i< localN; i++){
			for(j=0; j< N; j++){
				tmpV[i*N+j] = BN1[k*(localN*N)+i*N+j];
			}
		}
		pzdotc_(&tmpN, &value, tmpV, &ONE, &ONE, descA, &ONE, tmpu_hat, 
						   &ONE, &ONE, descA, &ONE);

//		printf("%d\t%f&%f\n", myrank, value);
		PODu0[k] = value;
	}

//	if(myrank==0)
//		for(i=0; i< localN; i++)
//			for(j=0; j< N; j++)
//			printf("%d\t%f&%f\t%f&%f\n", myrank, tmpu_hat[i*(N)+j], tmpV[i*(N)+j]);

}

