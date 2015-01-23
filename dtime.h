/*
 * =====================================================================================
 *
 *       Filename:  dtime.h
 *
 *    Description:  Dtime function, call gettimeofday()
 *
 *        Version:  1.0
 *        Created:  07/25/2014 04:23:05 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Cao Zongyan (), zycao@sccas.cn
 *        Company:  Supercomputing Center, CNIC, CAS, China
 *
 * =====================================================================================
 */

#ifndef _DTIME_H_
#define _DTIME_H_

#include "offload.h"
#include <sys/time.h>

__OffloadFunc_Macro__
inline double dtime() // Return double style timestamp using gettimeofday().
{
	double tseconds;
	struct timeval mytime;

	gettimeofday(&mytime, NULL);
	tseconds = (double)(mytime.tv_sec + mytime.tv_usec*1.0e-6);
	return (tseconds);
} /* dtime() */

#endif /* _DTIME_H_ */
