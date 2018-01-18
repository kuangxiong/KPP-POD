/*************************************************************************
	> File Name: KPP_GetPODNumber.c
	> Author:kuangxiong 
	> Mail:kuangxiong@lsec.cc.ac.cn 
	> Created Time: 2017年12月29日 星期五 15时31分28秒
 ************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int
KPP_GetPODNumber(double *S, int ncup, double gamma)
{
	int i;
	double sum1, partsum;
	sum1 = 0;
	partsum = 0;
	for(i=0; i< ncup; i++)
		sum1 = sum1 + S[i];
	for(i=0; i<ncup; i++){
		partsum = partsum + S[i];
		if(partsum > gamma*sum1)
			return i+1;
	}
}

