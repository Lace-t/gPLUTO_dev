////////////////////////////////////////////////////////////////////////
/*!
  \file
  \brief Provide a simple interface for measuring function execution time.

  The ClockBeg() and ClockEnd() function may be called to time a
  given function:

  void Myfun(...){
    ClockBeg("FuncName");
   
    < function body >;


    ClockEnd("FuncName");
  }

  ClockEnd() will also automatically print the time elapsed since ClockBeg()

  \authors A. Mignone (andrea.mignone@unito.it)
  \date    June 25, 2024
*/
////////////////////////////////////////////////////////////////////////
#include "gpluto.hpp"

#define CLOCK_MAX_FUNCTIONS 64
#ifdef PARALLEL
static double tbeg[CLOCK_MAX_FUNCTIONS];
static double tend[CLOCK_MAX_FUNCTIONS];
#else
static clock_t tbeg[CLOCK_MAX_FUNCTIONS];
static clock_t tend[CLOCK_MAX_FUNCTIONS];
#endif
static string func_name_stack[CLOCK_MAX_FUNCTIONS];
static string func_null = "Func NULL";

////////////////////////////////////////////////////////////////////////
void ClockBeg(string func_name)
/*! 
  Start clock for function 'func_name'
*/
////////////////////////////////////////////////////////////////////////
{
  int id;
  static int first_call = YES;
  
  // 0. Initialize all stack members to "NULL" 

  if (first_call) {
    for (id = 0; id < CLOCK_MAX_FUNCTIONS; id++) {
      func_name_stack[id] = func_null;
    }
    first_call = NO;
  }

  // 1. Loop on CLCOK_MAX_FUNCTIONS

  for (id = 0; id < CLOCK_MAX_FUNCTIONS; id++){
    if (func_name_stack[id] == func_name) { 
      string msg =   "! ClockBeg(): function " + func_name
                   + " already started but never ended";
      QuitPLUTO(msg, 1);
    }

    if (func_name_stack[id] == func_null) {
      func_name_stack[id] = func_name;
      break;
    }
  }

  #ifdef PARALLEL
  tbeg[id] = MPI_Wtime();
  #else
  tbeg[id] = clock();
  #endif
}

////////////////////////////////////////////////////////////////////////
void ClockEnd(string func_name)
/*! 
 * Stop clock for function 'func_name'
*/
////////////////////////////////////////////////////////////////////////
{
  int    id;
  double wtime;

  // 0. Search function id 

  for (id = 0; id < CLOCK_MAX_FUNCTIONS; id++){
    if (func_name_stack[id] == func_name) {
      break;
    }
  }
  if (id == CLOCK_MAX_FUNCTIONS){
    string msg = "! ClockEnd(): function "+func_name+" not known";
    QuitPLUTO(msg, 1);
  }

  // 1. Clock end & print

  #ifdef PARALLEL
  tend[id] = MPI_Wtime();
  wtime    = (tend[id] - tbeg[id]);
  #else
  tend[id] = clock();
  wtime    = (double)(tend[id] - tbeg[id]) / CLOCKS_PER_SEC;
  #endif

  print ("> ClockEnd(); ** function %s [id = %d] took = %12.6e (sec) **\n", 
         func_name.c_str(), id, wtime);

  // 2. Reset stack member
  
  func_name_stack[id] = func_null;
}
