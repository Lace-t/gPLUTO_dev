////////////////////////////////////////////////////////////////////////
/*!
  \file
  \brief Output grid to "grid.out"

  \author A. Mignone (andrea.mignone@unito.it)
  \date   March 17, 2026
*/
////////////////////////////////////////////////////////////////////////
#include "gpluto.hpp"

void OutputGrid (Grid *grid, Runtime *rtime)
{
  int dir;

  // ------------------------------------------------------
  // Write grid file (grid.out)
  // ------------------------------------------------------

  if (g_rank == 0){
    int i;
    char fname[MAX_STRING_SIZE+16];
    time_t time_now;
    time(&time_now);
    FILE *fg;

    snprintf (fname, sizeof(fname), "%s/grid.out", rtime->output_dir);
    fg = fopen(fname,"w");
    fprintf (fg, "# ******************************************************\n");
    fprintf (fg, "# gPLUTO %s Grid File\n",gPLUTO_VERSION);

    fprintf (fg, "# Generated on  %s",asctime(localtime(&time_now)));

    fprintf (fg, "# \n");
    fprintf (fg,"# DIMENSIONS: %d\n",DIMENSIONS);
    #if GEOMETRY == CARTESIAN
    fprintf (fg, "# GEOMETRY:   CARTESIAN\n");
    #elif GEOMETRY == CYLINDRICAL
    fprintf (fg, "# GEOMETRY:   CYLINDRICAL\n");
    #elif GEOMETRY == SPHERICAL
    fprintf (fg, "# GEOMETRY:   SPHERICAL\n");
    #endif
    DIM_LOOP(dir){
      fprintf (fg, "# X%d: [% f, % f], %d point(s), %d ghosts\n", dir+1,
               grid->xbeg_glob[dir], grid->xend_glob[dir],
               grid->nint_glob[dir], grid->nghost[dir]);
    }
    fprintf (fg, "# ******************************************************\n");
    for (dir = 0; dir < 3; dir++) {
      int ngh  = grid->nghost[dir];
      int gbeg = grid->gbeg[dir];
      int gend = grid->gend[dir];
      fprintf (fg, "%d \n", gend - gbeg + 1);
      for (i = gbeg; i <= gend; i++) {
       fprintf (fg, " %d   %18.12e    %18.12e\n",
                i-ngh+1, grid->xL_glob[dir][i], grid->xR_glob[dir][i]);
      }
    }
    fclose(fg);
  }  // if (g_rank == 0) 
}