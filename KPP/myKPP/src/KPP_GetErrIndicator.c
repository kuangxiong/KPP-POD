/*************************************************************************
	> File Name: KPP_GetErrIndiCator.c
	> Author:kuangxiong 
	> Mail: kuangxiong@lsec.cc.ac.cn
	> Created Time: 2018年01月14日 星期日 16时29分20秒
 ************************************************************************/
#include"stdio.h"
#include<stdlib.h>
#include"math.h"
#include"../fun.h"
#include"mpi.h"

double 
KPP_GetErrIndicator(int nPOD, int localN, int N, double cost, double dt,complex *curPOD, 
        complex *prePOD, complex *BN1, complex *NotTMat, complex *TMat)
{
    double Errind, ErrlocalNorm, ErrNorm, localNorm, Norm ;
    int i, j, k, tmpN, tmpM, ONE = 1;
    complex *curu, *preu, *lgNotTMat, *BN, *lgTMat, *tmpu, alpha, beta;
    curu = calloc(localN*(N/2+1), sizeof(*curu));
    preu = calloc(localN*(N/2+1), sizeof(*preu));
    tmpu = calloc(localN*(N/2+1), sizeof(*tmpu));
    BN = calloc(nPOD*localN*(N/2+1), sizeof(*BN));

    lgNotTMat = calloc(localN*(N/2+1), sizeof(*lgNotTMat));
    lgTMat = calloc(localN*(N/2+1), sizeof(*lgTMat));

    for(k=0; k< nPOD; k++)
        for(i=0; i< localN; i++)
            for(j=0; j< (N/2+1); j++)
                BN[k*localN*(N/2+1)+i*(N/2+1)+j] = BN1[k*localN*N+i*N+j];

//    for(i=0; i< nPOD; i++)
//        printf("%f&%f\n", curPOD[i]);



//    for(k=0; k< nPOD; k++){
//        for(i=0; i< localN; i++){
//            for(j=0; j< N/2+1; j++){
//                lgNotTMat[i*(N/2+1)+j] += NotTMat[k*localN*(N/2+1)+i*(N/2+1)+j]* prePOD[k];
//                lgTMat[i*(N/2+1)+j] += TMat[k*localN*(N/2+1)+i*(N/2+1)+j]*prePOD[k];
////              curu[i*(N/2+1)+j] += BN[k*localN*(N/2+1)+i*(N/2+1)+j]*curPOD[k];
////              preu[i*(N/2+1)+j] += BN[k*localN*(N/2+1)+i*(N/2+1)+j]*prePOD[k];
//            }
//        }
//    }
    alpha = 1.0;
    beta = 0.0;
    tmpM = localN*(N/2+1);
    tmpN = nPOD;
    zgemv_("N", &tmpM, &tmpN, &alpha, NotTMat, &tmpM, prePOD, &ONE, &beta, lgNotTMat, &ONE);
    zgemv_("N", &tmpM, &tmpN, &alpha, TMat, &tmpM, prePOD, &ONE, &beta, lgTMat, &ONE);
    zgemv_("N", &tmpM, &tmpN, &alpha, BN, &tmpM, curPOD, &ONE, &beta, curu, &ONE);
    zgemv_("N", &tmpM, &tmpN, &alpha, BN, &tmpM, prePOD, &ONE, &beta, preu, &ONE);



    for(i=0; i< localN; i++){
        for(j=0; j< (N/2+1); j++){
            tmpu[i*(N/2+1)+j] = lgNotTMat[i*(N/2+1)+j] + cost * lgTMat[i*(N/2+1)+j];
            tmpu[i*(N/2+1)+j] = (curu[i*(N/2+1)+j]-preu[i*(N/2+1)+j])/dt - tmpu[i*(N/2+1)+j];
        }
    }
    
    tmpN = 2*localN*(N/2+1);
    ErrlocalNorm = dnrm2_(&tmpN, tmpu, &ONE);
    ErrlocalNorm = ErrlocalNorm*ErrlocalNorm;
    localNorm = dnrm2_(&tmpN, preu, &ONE);
    localNorm = localNorm*localNorm;
    MPI_Allreduce(&ErrlocalNorm, &ErrNorm, ONE, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&localNorm, &Norm, ONE, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    Errind = sqrt(ErrNorm/Norm);
//  Errind = sqrt(Norm);
    free(curu);
    free(preu);
    free(BN);
    free(tmpu);
    free(lgNotTMat);
    free(lgTMat);
//    printf("%e\n", Errind);
    return Errind;
}

