////////////////////////////////////////////////////////////////////////
/*! 
  \file  
  \brief Parse command line options.

  Parse command line options at runtime and set values in
  the Cmd_Line structure.

  \authors A. Mignone (andrea.mignone@unito.it)
  \date    Nov 22, 2025
*/
////////////////////////////////////////////////////////////////////////
#include "gpluto.hpp"

static void PrintUsage();

////////////////////////////////////////////////////////////////////////
/*!
  Parse command line options. Error messages will be output to
  stdout.
  
  \param[in]  argc       argument count
  \param[in]  argv       argument vector
  \param[out] cmd        the command-line structure.
*/
////////////////////////////////////////////////////////////////////////
void ParseCmdLineArgs (int argc, char *argv[], cmdLine *cmd)
{
  int i;
  int dir;
  int maxsize = 32;

  // 1. Set default values first

  snprintf(cmd->input_file, maxsize, "pluto.ini");  
  cmd->restart   = NO;
  cmd->h5restart = NO;
  cmd->prestart  = NO;   // If restart is enabled, restart also particles 
  cmd->maxsteps  = -1;
  cmd->write     = YES;
  cmd->makegrid  = NO; 
  cmd->jet       = -1;   // means option is not used
  cmd->xres      = -1;   // means no grid resizing  

  cmd->nproc[IDIR] = -1; // means autodecomp will be used 
  cmd->nproc[JDIR] = -1;
  cmd->nproc[KDIR] = -1;

  // 2. Parallelize all (active) directions

  #ifdef PARALLEL  
  cmd->parallel_dim[IDIR] = (INCLUDE_IDIR ? YES:NO); 
  cmd->parallel_dim[JDIR] = (INCLUDE_JDIR ? YES:NO); 
  cmd->parallel_dim[KDIR] = (INCLUDE_KDIR ? YES:NO); 
  #endif

  // 3. Parse Command Line Options.
  // 
  //    Note: Since at this time the output directory
  //          is not known and "pluto.log" has not
  //          been opened yet, we use printf to issue
  //          error messages.

  for (i = 1; i < argc ; i++){

    if (!strcmp(argv[i],"-dec")) {

      // Start reading integers at i+1

      for (dir = 0; dir < DIMENSIONS; dir++){

        if ((++i) >= argc){
          if (g_rank == 0){
            DIM_SELECT(printf ("! You must specify -dec n1\n");  ,
                       printf ("! You must specify -dec n1  n2\n");  ,
                       printf ("! You must specify -dec n1  n2  n3\n");)
          }
          QuitPLUTO(" ", 1);
        }
        cmd->nproc[dir] = atoi(argv[i]);
        if (cmd->nproc[dir] == 0){
          if (g_rank == 0) {
            printf ("! Incorrect number of processor for dir = %d \n", dir);
          }
          QuitPLUTO(" ", 1);
        }
      }

    }else if (!strcmp(argv[i],"-i")) {

      snprintf (cmd->input_file, maxsize, "%s",argv[++i]);
    
    }else if (!strcmp(argv[i],"-makegrid")) {

      cmd->makegrid = YES;

    }else if (!strcmp(argv[i],"-maxsteps")){

      if ((++i) >= argc){
        if (g_rank == 0) printf ("! You must specify -maxsteps nn\n");
        QuitPLUTO(" ", 1);
      }else{
        cmd->maxsteps = atoi(argv[i]);
        if (cmd->maxsteps < 0) {
          if (g_rank == 0)
            printf ("! You must specify -maxsteps nn, with nn >= 0 \n");
            QuitPLUTO(" ", 1);
        }
      }

    }else if (!strcmp(argv[i],"-no-write")) {

      cmd->write = NO;

    }else if (!strcmp(argv[i],"-no-x1par")) {

      cmd->parallel_dim[IDIR] = NO;

    }else if (!strcmp(argv[i],"-no-x2par")) {

      cmd->parallel_dim[JDIR] = NO;

    }else if (!strcmp(argv[i],"-no-x3par")) {

      cmd->parallel_dim[KDIR] = NO;

    }else  if (   !strcmp(argv[i],"-restart")
               || !strcmp(argv[i],"-frestart")
               || !strcmp(argv[i],"-h5restart")) {

     // Default restart is last written file (-1)

      if (!strcmp(argv[i], "-restart")) {
        cmd->restart  = YES;  // can only take YES/NO values 
        cmd->prestart = YES;  // Restart particles as well   
      }else if (!strcmp(argv[i], "-frestart")) {
        cmd->restart  = YES;  // can only take YES/NO values 
        cmd->prestart = NO;   // Do not restart particles    
      }else{
        cmd->h5restart = YES;
      }
      cmd->nrestart = -1;   // the file number to restart from

      if ((++i) < argc){
        char *endptr;
        cmd->nrestart = (int)strtol(argv[i], &endptr, 10);

        // if a non-numerical character is encountered, 
        //  cmd->nrestart should reset to -1

        if (endptr == argv[i]){
          i--;
          cmd->nrestart = -1;
        }
      }

    } else if (!strcmp(argv[i],"-x1jet")) {

      cmd->jet = IDIR;

    } else if (!strcmp(argv[i],"-x2jet")) {

      cmd->jet = JDIR;

    } else if (!strcmp(argv[i],"-x3jet")) {

      cmd->jet = KDIR;

    }else if (!strcmp(argv[i],"-xres")){

      if ((++i) >= argc){
        if (g_rank == 0) printf ("! You must specify -xres nn\n");
        QuitPLUTO(" ", 1);
      }else{
        cmd->xres = atoi(argv[i]);
        if (cmd->xres <= 1) {
          if (g_rank == 0) printf ("! You must specify -xres nn, with nn > 1 \n");
          QuitPLUTO(" ", 1);
        }
      }

    }else if (!strcmp(argv[i],"--help")){

      PrintUsage();
      QuitPLUTO(" ", 1);

    }else{
      if (g_rank == 0) printf ("! Unknown option '%s'\n",argv[i]);
      QuitPLUTO(" ", 1);
    }
  }

  // 4. Disable domain decomposition in the 
  //    direction specified by -xnjet

  if      (cmd->jet == IDIR) cmd->parallel_dim[IDIR] = NO;
  else if (cmd->jet == JDIR) cmd->parallel_dim[JDIR] = NO;
  else if (cmd->jet == KDIR) cmd->parallel_dim[KDIR] = NO;

}
////////////////////////////////////////////////////////////////////////
/*!
  
*/
////////////////////////////////////////////////////////////////////////
void PrintUsage()
{

  if (g_rank != 0) return;
  
  printf ("Usage: pluto [options]\n\n");
  printf ("           or \n\n");
  printf ("       mpirun -np NP ./pluto [options]\n\n");
  printf ("[options] are:\n\n");
  printf (" -dec n1 [n2] [n3]\n");  
  printf ("    Enable user-defined parallel decomposition mode. The integers\n");
  printf ("    n1, n2 and n3 specify the number of processors along the x1,\n");
  printf ("    x2, and x3 directions. There must be as many integers as the\n");
  printf ("    number of dimensions and their product must equal the total\n");
  printf ("    number of processors used by mpirun or an error will occurr.\n\n"); 

  printf (" -frestart n\n");
  printf ("    Restart computations for the fluid (no particles) from the\n");
  printf ("    n-th output file in double in precision format (.dbl).\n\n");

  printf (" -i <name>\n");
  printf ("    Use <name> as initialization file instead of pluto.ini.\n\n");

  printf (" --help\n");
  printf ("    Show this option summary.\n\n");

  printf (" -h5restart n\n");
  printf ("    Restart computations from the n-th output file in HDF5\n");
  printf ("    double precision format (.dbl.h5).\n\n");

  printf (" -makegrid\n");
  printf ("    Generate grid only, do not start computations.\n\n");

  printf (" -maxsteps n\n");
  printf ("    Stop computations after n steps.\n\n");

  printf (" -no-write\n");
  printf ("    Do not write data to disk.\n\n");
  
  printf (" -no-x1par, -no-x2par, -no-x3par\n");
  printf ("    Do not perform parallel domain decomposition along the x1, x2\n");
  printf ("    or x3 direction, respectively.\n\n");

  printf (" -restart n\n");
  printf ("    Restart computations from the n-th output file in double in\n");
  printf ("    precision format (.dbl).\n\n");
  
  printf (" -x1jet, -x2jet, -x3jet\n");
  printf ("    Exclude from integration regions of zero pressure gradient\n");
  printf ("    that extends up to the end of the domain in x1, x2 or x3\n");
  printf ("    direction, respectively. This option is specifically\n");
  printf ("    designed for jets propagating along one of the coordinate\n");
  printf ("    axis. In parallel mode, parallel decomposition is not\n");
  printf ("    performed along the selected direction.\n\n");

  printf (" -xres n1\n");
  printf ("    Set the grid resolution in the x1 direction to n1 zones\n");
  printf ("    by overriding pluto.ini. Cell aspect ratio is preserved by\n");
  printf ("    modifying the grid resolution in the other coordinate\n"); 
  printf ("    directions accordingly.\n");

}
