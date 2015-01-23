/*
 * =====================================================================================
 *
 *       Filename:  dbginfo.h
 *
 *    Description:  Debugging info macro defines.
 *
 *        Version:  1.0
 *        Created:  12/11/2014 06:33:29 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Cao Zongyan (), zycao@sccas.cn
 *        Company:  Supercomputing Center, CNIC, CAS, China
 *
 * =====================================================================================
 */

#ifndef _DBGINFO_H_
#define _DBGINFO_H_

#include "offload.h"
#include <stdio.h>

__OffloadVar_Macro__
extern struct _IO_FILE *stdin;      /*  Standard input stream.  */
__OffloadVar_Macro__
extern struct _IO_FILE *stdout;     /*  Standard output stream.  */
__OffloadVar_Macro__
extern struct _IO_FILE *stderr;     /*  Standard error output stream.  */

__OffloadFunc_Macro__
extern int fprintf (FILE *__restrict __stream, __const char *__restrict __format, ...);

#ifdef NMK_DEBUG

#ifdef NMK_DBG_STDOUT
#define DBG_OUT stdout
#else
#define DBG_OUT stderr
#endif /* NMK_DBG_STDOUT */

#define DBG_MEMORY	1
#define DBG_DATA	2
#define DBG_CALC	4
#define DBG_LOGIC	8
#define DBG_TIME	16
#define DBG_TMP		32768
#define DBG_ALL		65535

__OffloadVar_Macro__
extern int dbg_level;
#define DBG_INFOL(dbgcode, fmt, args...) {\
	if (dbgcode & dbg_level) \
		fprintf(DBG_OUT, "[DBG %d]%s(%d): "fmt" ", dbgcode, __FILE__, __LINE__, ##args);\
}
#define DBG_INFOFL(dbgcode, dbgfile, fmt, args...) {\
	if (dbgcode & dbg_level) \
		fprintf(dbgfile, "[DBG %d]%s(%d): "fmt" ", dbgcode, __FILE__, __LINE__, ##args);\
}
#define DBG_INFO(dbgcode, fmt, args...) {\
	if (dbgcode & dbg_level) \
		fprintf(DBG_OUT, fmt, ##args);\
}
#define DBG_INFOF(dbgcode, dbgfile, fmt, args...) {\
	if (dbgcode & dbg_level) \
		fprintf(dbgfile, fmt, ##args);\
}

#else /* NMK_DEBUG */

#define DBG_INFOL(...)
#define DBG_INFOFL(...)
#define DBG_INFO(...)
#define DBG_INFOF(...)
  
#endif /* NMK_DEBUG */

#define LOG_OUT stdout

#define LOG_COMMON	1
#define LOG_FATAL	2
#define LOG_VB1		4
#define LOG_VB2		8
#define LOG_TMP		32768
#define LOG_ALL		65535

__OffloadVar_Macro__
extern int log_level;
#define LOG_INFOL(logcode, fmt, args...) {\
	if (logcode & log_level) \
		fprintf(LOG_OUT, "[INFO %d]%s(%d): "fmt" ", logcode, __FILE__, __LINE__, ##args);\
}
#define LOG_INFOFL(logcode, logfile, fmt, args...) {\
	if (logcode & log_level) \
		fprintf(logfile, "[INFO %d]%s(%d): "fmt" ", logcode, __FILE__, __LINE__, ##args);\
}
#define LOG_INFO(logcode, fmt, args...) {\
	if (logcode & log_level) \
		fprintf(LOG_OUT, fmt, ##args);\
}
#define LOG_INFOF(logcode, logfile, fmt, args...) {\
	if (logcode & log_level) \
		fprintf(logfile, fmt, ##args);\
}
#define CLOG(fmt, args...) {\
	if (LOG_COMMON & log_level) \
		fprintf(LOG_OUT, fmt, ##args);\
}
#define CLOGF(logfile, fmt, args...) {\
	if (LOG_COMMON & log_level) \
		fprintf(logfile, fmt, ##args);\
}

#endif /* _DBGINFO_H_ */
