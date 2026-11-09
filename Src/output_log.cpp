////////////////////////////////////////////////////////////////////////
/*!
  \file
  \brief Output log file driver.

  The integration log file is divided into a "pre-step"
  and a "post-step" output.

  \authors A. Mignone(andrea.mignone@unito.it)\n
           B. Vaidya

  \date    Apr 16, 2026
*/
////////////////////////////////////////////////////////////////////////
  
#include "gpluto.hpp"
#include <sys/stat.h>
#include <sys/types.h>

#ifdef PARALLEL
static FILE *g_flog;
static char log_file_name[512];
#endif

////////////////////////////////////////////////////////////////////////
/*    The next set of functions provides basic functionalities to

     - set the log file
     - formatted output to the log file through the printLog() and 
       printLog() functions                                               
*/              
////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////
/*!
  Close a previously opened log file
*/
////////////////////////////////////////////////////////////////////////
void LogFileClose(void)
{
#ifdef PARALLEL
  if (g_flog != NULL) fclose(g_flog);
#endif
}

////////////////////////////////////////////////////////////////////////
/*!
  Flushes the output buffer of a log file stream.
*/
////////////////////////////////////////////////////////////////////////
void LogFileFlush(void)
{
#ifdef PARALLEL
  // print ("\n> LogFileFlush 00\n");
  if (g_flog != NULL){
    // print ("\n> LogFileFlush 01\n");
    fflush (g_flog);
    // print ("\n> LogFileFlush 02\n");
  }
#endif
}

////////////////////////////////////////////////////////////////////////
/*!
  Open log file name in parallel mode.
  Each processor has its own log file "pluto.g_rank.log", unless the
  MULTIPLE_LOG_FILES flag has been set to FALSE.
  
  \param [in] output_dir  the name of the output directory
  \param [in] cmd         pointer to cmd line option structure.
*/
////////////////////////////////////////////////////////////////////////
void LogFileOpen(char *log_dir, const char *mode)
{
#ifdef PARALLEL

  snprintf(log_file_name, sizeof(log_file_name), "%s/pluto.%d.log", log_dir, g_rank);

  #if MULTIPLE_LOG_FILES == NO
  if (g_rank != 0) return;
  #endif

  g_flog = fopen(log_file_name, mode);

  // -- check that we have a valid directory name -- 

  if (g_flog == NULL){
    printf ("! SetLogFile(): %s cannot be written.\n",log_file_name);
    printf ("  Using current directory instead.\n");
    snprintf (log_file_name, sizeof(log_file_name),"./pluto.%d.log",g_rank);
    g_flog = fopen(log_file_name, mode);
  }
#endif
} 

////////////////////////////////////////////////////////////////////////
/*!
  Define print function to be used in serial/parallel mode
  
  Note: when MULTIPLE_LOG_FILES == NO, only proc #0 does the
  actual printing.
*/
////////////////////////////////////////////////////////////////////////
void print (const char *fmt, ...)
{
#if MULTIPLE_LOG_FILES == NO
  if (g_rank != 0) return;
#endif

  va_list args;
  va_start(args, fmt);

#ifdef PARALLEL
  vfprintf(g_flog, fmt, args);
#else
  vprintf(fmt, args);
#endif

  va_end(args);
}

////////////////////////////////////////////////////////////////////////
/*!
  Define printLog function [All processors will write]
*/
////////////////////////////////////////////////////////////////////////
void printLog (const char *fmt, ...)
{
  va_list args;
  va_start(args, fmt);

#ifdef PARALLEL
  // --------------------------------------------
  // File may not be opened if
  // MULTIPLE_LOG_FILES is set to NO.
  // --------------------------------------------

  #if MULTIPLE_LOG_FILES == NO
  if (g_flog == NULL)  g_flog = fopen(log_file_name, "a");
  #endif

  vfprintf(g_flog, fmt, args);
#else
  vprintf(fmt, args);
#endif

  va_end(args);
}

////////////////////////////////////////////////////////////////////////
/*!
*/
////////////////////////////////////////////////////////////////////////
char *IndentString()
{
  static char str[64];

#if 00
  if      (g_stepNumber < 10)     sprintf (str,"%7s"," ");
  else if (g_stepNumber < 100)    sprintf (str,"%8s"," ");
  else if (g_stepNumber < 1000)   sprintf (str,"%9s"," ");
  else if (g_stepNumber < 10000)  sprintf (str,"%10s"," ");
  else if (g_stepNumber < 100000) sprintf (str,"%11s"," ");
  else                            sprintf (str,"%12s"," ");
#endif
//  sprintf (str,"%10s"," ");
  snprintf (str,sizeof(str), "%10s"," ");
  return str;
}
