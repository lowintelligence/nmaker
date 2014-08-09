/*
 * =====================================================================================
 *
 *       Filename:  offload.h
 *
 *    Description:  Offload controller head file
 *
 *        Version:  1.0
 *        Created:  07/24/2014 02:23:34 PM
 *       Revision:  none
 *       Compiler:  icc
 *
 *         Author:  Cao Zongyan (), zycao@sccas.cn
 *        Company:  Supercomputing Center, CNIC, CAS, China
 *
 * =====================================================================================
 */

#ifndef _OFFLOAD_H_
#define _OFFLOAD_H_

#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <pthread.h>

#ifdef __INTEL_OFFLOAD 
	#define __OffloadVar_Macro__ __declspec(target(mic)) 
	#define __OffloadFunc_Macro__ __attribute__ ((target(mic)))
#else
	#define __OffloadVar_Macro__
	#define __OffloadFunc_Macro__
#endif

#endif
