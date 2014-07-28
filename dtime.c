/*
 * =====================================================================================
 *
 *       Filename:  dtime.c
 *
 *    Description:  call gettimeofday()
 *
 *        Version:  1.0
 *        Created:  07/25/2014 04:20:28 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Cao Zongyan (), zycao@sccas.cn
 *        Company:  Supercomputing Center, CNIC, CAS, China
 *
 * =====================================================================================
 */
#include "dtime.h"

double dtime()
{
	double tseconds = 0.0;
	struct timeval mytime;	
	gettimeofday(&mytime, (struct timezone*)0);
	tseconds = (double)(mytime.tv_sec + mytime.tv_usec*1.0e-6);
	return (tseconds);
}

