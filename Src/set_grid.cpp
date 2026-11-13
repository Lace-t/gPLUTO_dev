////////////////////////////////////////////////////////////////////////
/*!
  \file
  \brief Create global and local grid.

  \author A. Mignone (andrea.mignone@unito.it)
  \date   March 24, 2026
*/
////////////////////////////////////////////////////////////////////////
#include "gpluto.hpp"
#include "grid_input.hpp"

////////////////////////////////////////////////////////////////////////
void SetGrid (Runtime *rtime, int *procs, Grid *grid, int cf)
/*!
  
  \param [in]  rtime   pointer to Runtime structure, used for defining
  the global grid.
  \param [in]  procs   an array containing the number of processors
  along the three diretions (parallel only).
  \param [out] grid    a pointer to the Grid structure
  \param [in]  cf      coarsening factor: a de-refinement ratio specifying 
  how much the grid should be de-refined.
*/
////////////////////////////////////////////////////////////////////////
{
  int i,j,k, dir;
  int ngh, gbeg;
  int lbeg, lend, offset;
  int include_dir[] = {INCLUDE_IDIR, INCLUDE_JDIR, INCLUDE_KDIR};
  int gsize[3];
  int lsize[3];
  static int first_call = 1;

  if (first_call) print ("\n> Generating grid...\n");

  // ------------------------------------------------------
  // 1. Set global grid properties
  // ------------------------------------------------------

  for (dir = 0; dir < 3; dir++) {

    gsize[dir] = rtime->npoint[dir];

    // 1a. De-refine grid if required

    if (cf > 1 && include_dir[dir]) {
      if ((gsize[dir] % cf) != 0){
        QuitPLUTO ("SetGrid(): coarsening factor incompatible",1);
      }
      gsize[dir] /= cf;  
    }

    // 1b. If a direction is not included, one zone will
    //     be set and no ghost zone.

    ngh = GetNghost()*include_dir[dir];
    grid->nghost[dir] = ngh;

    // 1c. Set sizes (interior and total number of zones)

    grid->gbeg[dir]   = ngh;
    grid->gend[dir]   = grid->gbeg[dir] + gsize[dir] - 1;
    grid->nint_glob[dir] = gsize[dir];
    grid->ntot_glob[dir] = gsize[dir] + 2*ngh;

    grid->nproc[dir]  = procs[dir];

    grid->xC_glob[dir] = ARRAY_1D(grid->ntot_glob[dir], double);
    grid->xL_glob[dir] = ARRAY_1D(grid->ntot_glob[dir]+1, double);
//    grid->xR_glob[dir] = grid->xL_glob[dir] + 1;
    grid->xR_glob[dir] = ARRAY_1D(grid->ntot_glob[dir], double);
    grid->dx_glob[dir] = ARRAY_1D(grid->ntot_glob[dir], double);

    #if MAPPED_GRID == YES
    grid->dxi_dx_glob[dir] = ARRAY_1D(grid->ntot_glob[dir], double);
    #endif
  }

  // ------------------------------------------------------
  // 2. Create global grid
  // ------------------------------------------------------

  for (dir = 0; dir < 3; dir++){

    double *xC   = grid->xC_glob[dir];
    double *xL   = grid->xL_glob[dir];
    double *xR   = grid->xR_glob[dir];
    double *dx   = grid->dx_glob[dir];
    double *xbeg = grid->xbeg_glob;
    double *xend = grid->xend_glob;

    ngh  = grid->nghost[dir];
    gbeg = grid->gbeg[dir];

    // 2a. Create interior grid

    xbeg[dir] = rtime->xbeg[dir];
    xend[dir] = rtime->xend[dir];

    if (rtime->grid_spec[dir].count) {
      std::vector<double> edges;
      try { edges=GridInput::Edges(rtime->grid_spec[dir],cf); }
      catch (const std::exception &error) {
        QuitPLUTO(string("SetGrid(): ")+error.what(),1);
      }
      int count=grid->nint_glob[dir];
      if (edges.size()!=size_t(count+1)) QuitPLUTO("SetGrid(): patch cell count mismatch",1);
      double first=edges[1]-edges[0],last=edges[count]-edges[count-1];
#if MAPPED_GRID == NO
      for(int j=1;j<count;++j)
        if(std::abs((edges[j+1]-edges[j])-first)>1.e-10*first)
          QuitPLUTO("Nonuniform ini grid requires MAPPED_GRID=YES (no user MappedGrid function needed)",1);
#endif
      auto edge=[&](int j) {
        if(j<0) return edges[0]+j*first;
        if(j>count) return edges[count]+(j-count)*last;
        return edges[j];
      };
      for(i=0;i<grid->ntot_glob[dir];++i) {
        xL[i]=edge(i-gbeg); xR[i]=edge(i-gbeg+1);
        dx[i]=xR[i]-xL[i]; xC[i]=.5*(xL[i]+xR[i]);
      }
#if MAPPED_GRID == YES
      grid->dxi[dir]=1.0/count;
      // Same metric discretization as the existing MappedGrid path.
      for(i=0;i<grid->ntot_glob[dir];++i) {
        double span=edge(i-gbeg+2)-edge(i-gbeg-1);
        grid->dxi_dx_glob[dir][i]=24.0*grid->dxi[dir]/(27.0*dx[i]-span);
      }
#endif
      print("> PLUTO-style grid [dir = %d, patches = %d]\n",dir,rtime->grid_spec[dir].count);
      continue;
    }

    // 2b. Mapped grid definition

    #if MAPPED_GRID == YES
    print ("> Mapped grid generation [dir = %d]\n",dir);

    grid->dxi[dir] = 1.0/grid->nint_glob[dir];
    double *dxi_dx = grid->dxi_dx_glob[dir];
    for (i = 0; i < grid->ntot_glob[dir]; i++){
      double xiL = (i - gbeg)*grid->dxi[dir];
      double xiR = (i - gbeg + 1)*grid->dxi[dir];

      xL[i] = MappedGrid (grid, dir, xiL);
      xR[i] = MappedGrid (grid, dir, xiR);

      #if MAPPED_GRID_FORCE_PERIODS == YES
      if (i < gbeg) {
        double L = grid->xend_glob[dir] - grid->xbeg_glob[dir];
        xL[i] = MappedGrid (grid, dir, xiL+1) - L;
        xR[i] = MappedGrid (grid, dir, xiR+1) - L;
      }
      if (i > grid->gend[dir]) {
        double L = grid->xend_glob[dir] - grid->xbeg_glob[dir];
        xL[i] = MappedGrid (grid, dir, xiL-1) + L;
        xR[i] = MappedGrid (grid, dir, xiR-1) + L;
      }
      #endif

      // 2b2. Verify that grid remains monotonic

      if (xL[i] >= xR[i]){
        printLog ("! SetGrid(): grid does not appear to be monotonic\n");
        printLog ("             dir = %d, xL[%d] = %12.6e, xR[%d] = %12.6e\n",
                   dir, i, xL[i], i, xR[i]);
        QuitPLUTO("",1);
      }

      // 2b3. Check that no Nan is produced

      if (isnan(xL[i]) || isnan(xR[i])){
        printLog ("! SetGrid(): nan found at zone %d\n",i);
        printLog ("             dir = %d, xL[%d] = %12.6e, xR[%d] = %12.6e\n",
                   dir, i, xL[i], i, xR[i]);
        QuitPLUTO("",1);
      }

      dx[i] = xR[i] - xL[i];
      xC[i] = 0.5*(xR[i] + xL[i]);
    }

    // 2b4. Compute 4th-order approximation to 
    //      dxi/dx = 1/(dx/dxi) from edge values.

    for (i = 1; i < grid->ntot_glob[dir]-1; i++){
      double h = grid->dxi[dir];
      dxi_dx[i] = 24.0*h/(27.0*(xR[i] - xL[i]) - (xR[i+1] - xL[i-1]));
    }

    #else

    for (i = 0; i < grid->ntot_glob[dir]; i++){
      dx[i] = (xend[dir] - xbeg[dir])/grid->nint_glob[dir];
      xC[i] = xbeg[dir] + (i - gbeg + 0.5)*dx[i];

      xL[i] = xbeg[dir] + (i - gbeg)*dx[i];
      xR[i] = xbeg[dir] + (i - gbeg + 1)*dx[i];
    }
    #endif
  }

  // ------------------------------------------------------
  // 3. Create local grid.
  //    In parallel mode, when the coarsening factor is > 1, 
  //    the local processor grid indices must be recalculated
  //    since AL_Get_bounds references to the original grid:
  //
  //
  //    beg -= ngh, end -= ngh-1
  //    beg /= cf , end /= cf
  //    beg += ngh, end += ngh-1
  // ------------------------------------------------------

  int  beg[3] = {0,0,0};
  int  end[3] = {0,0,0};
  #ifdef PARALLEL
  
#if USE_AL == NO
for (int dir = 0; dir < DIMENSIONS; dir++) {
  beg[dir] = new_SZ.beg[dir];
  end[dir] = new_SZ.end[dir];
}

if (cf > 1){
  DIM_LOOP(dir){
    beg[dir] -= grid->nghost[dir]; 
    end[dir] -= grid->nghost[dir]-1;
    beg[dir] /= cf;
    end[dir] /= cf;
    beg[dir] += grid->nghost[dir]; 
    end[dir] += grid->nghost[dir]-1;
  }
}
#else
  AL_Get_bounds  (SZ, beg, end, grid->nghost, AL_C_INDEXES);

  if (cf > 1){
    DIM_LOOP(dir){
      beg[dir] -= grid->nghost[dir]; 
      end[dir] -= grid->nghost[dir]-1;
      beg[dir] /= cf;
      end[dir] /= cf;
      beg[dir] += grid->nghost[dir]; 
      end[dir] += grid->nghost[dir]-1;
    }
  }
#endif
  #else
  for (dir = 0; dir < 3; dir++) {
    beg[dir] = grid->gbeg[dir];
    end[dir] = grid->gend[dir];
  }
  #endif

// printLog("SetGrid(): beg = {%d, %d, %d}; end = {%d, %d, %d}\n", 
//            beg[IDIR], beg[JDIR], beg[KDIR], 
//            end[IDIR], end[JDIR], end[KDIR]);

  for (dir = 0; dir < 3; dir++){
    ngh = grid->nghost[dir];

    lsize[dir] = end[dir] - beg[dir] + 1;

    grid->beg[dir] = beg[dir];
    grid->end[dir] = end[dir];

    grid->lbeg[dir] = lbeg = ngh;
    grid->lend[dir] = lend = lbeg + lsize[dir] - 1;

    grid->nint[dir] = lsize[dir];
    grid->ntot[dir] = lsize[dir] + 2*ngh;

    // 3a. Here ntot_stag[JDIR][dir] gives the total number of zones
    //     for the k-th vector component in the dir direction, e.g.
    //
    //     ntot_stag[JDIR][IDIR] = nx 
    //     ntot_stag[JDIR][JDIR] = ny+1
    //     ntot_stag[JDIR][KDIR] = nz 

    grid->ntot_stag[IDIR][dir] = grid->ntot[dir] + (dir == IDIR ? 1:0);
    grid->ntot_stag[JDIR][dir] = grid->ntot[dir] + (dir == JDIR ? 1:0);
    grid->ntot_stag[KDIR][dir] = grid->ntot[dir] + (dir == KDIR ? 1:0);

    offset = grid->beg[dir] - ngh;
    grid->xC[dir] = grid->xC_glob[dir] + offset;
    grid->xR[dir] = grid->xR_glob[dir] + offset;
    grid->xL[dir] = grid->xL_glob[dir] + offset;
    grid->dx[dir] = grid->dx_glob[dir] + offset;
    #if MAPPED_GRID == YES
    grid->dxi_dx[dir] = grid->dxi_dx_glob[dir] + offset;
    #endif

    grid->xbeg[dir] = grid->xL[dir][lbeg];
    grid->xend[dir] = grid->xR[dir][lend];

    grid->lbound[dir] = grid->rbound[dir] = 0;
    if (grid->beg[dir] == grid->gbeg[dir]) {
      grid->lbound[dir] = rtime->left_bound[dir];
    }

    if (grid->end[dir] == grid->gend[dir]) {
      grid->rbound[dir] = rtime->right_bound[dir];
    }
  }

// printLog("SetGrid(): lsize = {%d, %d, %d}\n", 
//            lsize[IDIR], lsize[JDIR], lsize[KDIR]);

// printLog("SetGrid(): nint = {%d, %d, %d}\n", 
//            grid->nint[IDIR], grid->nint[JDIR], grid->nint[KDIR]);

  // 4. Set geometry

  SetGeometry(grid);

  // 5. Get the maximum number of points and
  //    minimum cell length among all direction.

#if PHYSICS != GRMHD
  double dl_min[3];
  for (dir = 0; dir < 3; dir++) dl_min[dir] = 1.e30;

  for (k = grid->lbeg[KDIR]; k <= grid->lend[KDIR]; k++) {
  for (j = grid->lbeg[JDIR]; j <= grid->lend[JDIR]; j++) {
  for (i = grid->lbeg[IDIR]; i <= grid->lend[IDIR]; i++) {
    double dl;
    Metric metric;
    #if GEOMETRY != CARTESIAN
    MetricSet (metric, grid->rC[i], grid->sC[j]);
    #endif

    dl           = grid->dx[IDIR][i]*metric.h1; 
    dl_min[IDIR] = MIN (dl_min[IDIR], dl);

    dl           = grid->dx[JDIR][j]*metric.h2; 
    dl_min[JDIR] = MIN (dl_min[JDIR], dl);
 
    dl           = grid->dx[KDIR][k]*metric.h3; 
    dl_min[KDIR] = MIN (dl_min[KDIR], dl);
  }}}

  #ifdef PARALLEL
  MPI_Allreduce (dl_min, grid->dl_min, 3, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
  #else
  for (dir = 0; dir < 3; dir++) grid->dl_min[dir] = dl_min[dir];
  #endif
#endif

  // 6. Print domain specifications

  if (first_call){
    print ("\n  Global grid:\n");
    DIM_LOOP(dir){
      if (fabs(grid->xend_glob[dir] - grid->xbeg_glob[dir]) >= 1.e4) { // Use scientific notation 
        print ("  X%d: [ %12.6e, %12.6e], %d cells (%d ghosts);", dir+1,
                grid->xbeg_glob[dir], grid->xend_glob[dir],
                grid->nint_glob[dir], grid->nghost[dir]);
      }else{                // Use floating point notation 
        print ("  X%d: [ %8.4f, %8.4f], %d cells (%d ghosts);", dir+1,
                grid->xbeg_glob[dir], grid->xend_glob[dir],
                grid->nint_glob[dir], grid->nghost[dir]);
      }
      print (" gbeg, gend = [%d, %d]\n", grid->gbeg[dir], grid->gend[dir]);
    }

    print ("\n");
    print ("  Local grid:\n");
    DIM_LOOP(dir){
      if (fabs(grid->xend[dir] - grid->xbeg[dir]) >= 1.e4) { // Use scientific notation 
        print ("  X%d: [ %-13.6e, %-13.6e], %d cells;", dir+1,
                grid->xbeg[dir], grid->xend[dir], grid->nint[dir]);
      }else{
        print ("  X%d: [ %8.4f, %8.4f], %d cells;", dir+1,
                grid->xbeg[dir], grid->xend[dir], grid->nint[dir]);
      }
      print (" lbeg, lend = [%d, %d]; beg, end = [%d, %d]\n",
                grid->lbeg[dir], grid->lend[dir],
                grid->beg[dir], grid->end[dir]);
    }
  }
  first_call = 0;

  grid->nmax = MAX(grid->ntot[IDIR], grid->ntot[JDIR]);
  grid->nmax = MAX(grid->nmax, grid->ntot[KDIR]);
}
