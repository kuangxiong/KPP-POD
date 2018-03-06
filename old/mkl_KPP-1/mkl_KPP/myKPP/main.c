/**********************************************
 * @filename      main.c
 * @Synopsis      
 * @author kuangxiong
 * @version 1
 * @date 2017-11-25
 *****************************************************/
#include"fun.h"

int 
main(int argc, char* argv[])
{
    int myrank, nprocs, N;
    DOUBLE eps, tau, theta, dt, laml, lamr, laml1, vaeps, y1, y2, Amp, t, cM, C;
    
    MPI_Init(&argc, &argv);
//    MPI_Comm *comm;
//    MPI_Comm_dup(MPI_COMM_WORLD, comm);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
/***************initialize variables****************/
    N = 128;
    laml = 0;
    lamr = 10;
    eps = 0.05;
    tau = 1.0;
    theta = 1.0;
    dt = 0.000001;
    Amp = 1.0;
    vaeps = 5e-3;
/***************************************************/
    if(N%nprocs!=0){
        printf("err:nprocs must be division by N!!!\n");
        MPI_Finalize();
        return 0;
    }
    if(myrank == 0){
		printf("step size:%f,\t theta: %f,\t Amp:%f\n", dt, theta, Amp);
		printf("Parameter:%f\n", tau);
		printf("The convergence criteria for the outer loop:%f\n", vaeps);
	}
/*************compute c_M****************/
    laml1 = lamr - (sqrt(5.0)-1)/2*(lamr-laml);

    cM = eps * laml1 * laml1 + laml1*Amp*sqrt(theta*theta + 1.0)+1.0/tau;
    C = -laml1*Amp*sqrt(theta*theta + 1.0);
/***************several adaptive ways to solver KPP problem*************/
//	KPP_PW(N, laml1, eps, tau, theta, dt, t, Amp, cM, C, y1, myrank, nprocs, MPI_COMM_WORLD);
//	KPP_TGPW(N, laml1, eps, tau, theta, dt, t, Amp, cM, C, y1, myrank, nprocs, MPI_COMM_WORLD);
//	testKPP_PW(N, laml1, eps, tau, theta, dt, t, Amp, cM, C, y1, myrank, nprocs, MPI_COMM_WORLD);
//	KPP_POD(N, laml1, eps, tau, theta, dt, t, Amp, cM, C, y1, myrank, nprocs, MPI_COMM_WORLD);
//	KPP_TGPOD(N, laml1, eps, tau, theta, dt, t, Amp, cM, C, y1, myrank, nprocs, MPI_COMM_WORLD);
	KPP_TGAPOD(N, laml1, eps, tau, theta, dt, t, Amp, cM, C, y1, myrank, nprocs, MPI_COMM_WORLD);
//	KPP_APOD(N, laml1, eps, tau, theta, dt, t, Amp, cM, C, y1, myrank, nprocs, MPI_COMM_WORLD);

    if(myrank ==0){
        printf("lambda:%f\t, eps:%f\t, cM:%f\n", laml1, eps, cM);
		printf("stop time:%f\n", t);
		printf("left point C*:%f\n", y1);
	}
    MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();
    return 0;
}
