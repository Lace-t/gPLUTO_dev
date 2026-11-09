////////////////////////////////////////////////////////////////////////
/*!
  \file
  \brief Set boundary conditions.

  The Boundary() function sets both internal and physical boundary
  conditions on all sides of the computational domain.
  It is used to fill ghost zones of both cell-centered and face-centered
  d arrays.\n
  The type of boundary conditions at the leftmost or rightmost side of a
  given grid is specified by the integers <c> grid.lbound[dir] </c> or
  <c> grid.rbound[dir] </c>, respectively.
  When this value is different from zero, the local processor borders
  the physical boundary and the admissible values for \c lbound or \c
  rbound are OUTFLOW, REFLECTIVE, AXISYMMETRIC, EQTSYMMETRIC, PERIODIC,
  SHEARING, ROTADTED or USERDEF.
  Conversely, when this value is zero (internal boundary), two neighboring
  processors that share the same side need to fill ghost zones by exchanging
  d values.
  This step is done here only for parallel computations on static grids.

  Predefined physical boundary conditions are handled by the
  following functions:

  - OutflowBoundary()
  - PeriodicBoundary()
  - ReflectiveBoundary()


  \note For staggered mesh+periodic b.c, the leftmost value of \e each
        processor is overwritten with a value NX1 zones to the right:


      -->Before:

         proc #0     proc #1     proc #2     proc #3
        |.......|   |.......|   |.......|   |.......|
        A       B   C       D   E       F   G       H

      --> PeriodicBoundary()
      --> After:

         proc #0     proc #1     proc #2     proc #3
        |.......|   |.......|   |.......|   |.......|
        H       B   B       D   D       F   F       H

  \author A. Mignone  (andrea.mignone@unito.it)
          V. Berta    (vittoria.berta@unito.it)
          M. Rossazza (marco.rossazza@unito.it)
          
  \date   Jan 27, 2026
*/
////////////////////////////////////////////////////////////////////////
#include "gpluto.hpp"
#define CT_DEBUG 0

////////////////////////////////////////////////////////////////////////
void BoundaryConservative (Data *data, Grid *grid)
/*!
  Set boundary conditions on one or more sides of the computational
  domain.
  
  \param [in,out] d  pointer to PLUTO Data structure containing the
  solution array (including centered and staggered
  fields)
  \param [in]     grid   pointer to grid structure.
*/
////////////////////////////////////////////////////////////////////////
{ 
  int  is;
  int  side[6] = {X1_BEG, X1_END, X2_BEG, X2_END, X3_BEG, X3_END};
  int  type[6], vsign[NVAR];
  int  *lbound = grid->lbound;
  int  *rbound = grid->rbound;
  int  par_dim[3] = {0, 0, 0};
  int  ib,ie,jb,je,kb,ke;

  int  ibeg = grid->lbeg[IDIR], iend = grid->lend[IDIR];
  int  jbeg = grid->lbeg[JDIR], jend = grid->lend[JDIR];
  int  kbeg = grid->lbeg[KDIR], kend = grid->lend[KDIR];
  int  nx1_tot = grid->ntot[IDIR];
  int  nx2_tot = grid->ntot[JDIR];
  int  nx3_tot = grid->ntot[KDIR];

  RBox bbox;
  #pragma acc enter data copyin(bbox)
  #pragma omp target enter data map(to: bbox)

  #pragma acc enter data create (vsign[:NVAR])
  #pragma omp target enter data map(alloc: vsign[:NVAR])

  // ------------------------------------------------------
  // 0. Assign b.c. using conservative variables 
  // ------------------------------------------------------

  RBox boxt;
  RBoxDefine(0, nx1_tot - 1, 0, nx2_tot - 1, 0, nx3_tot - 1, &boxt);
  #pragma acc enter data copyin(boxt)
  #pragma omp target enter data map(to: boxt)

  // 0b. Check the number of processors in each direction +
  //     Call userdef internal boundary with side == 0 

  DIM_EXPAND(par_dim[0] = grid->nproc[IDIR] > 1;  ,
             par_dim[1] = grid->nproc[JDIR] > 1;  ,
             par_dim[2] = grid->nproc[KDIR] > 1;)

// Do not allow internal boundary on conservative variables at the moment.
//  #if INTERNAL_BOUNDARY == YES
//  UserDefBoundary (data, NULL, 0, grid);
//  #endif

  #if (PHYSICS == ResRMHD) && (DIVE_CONTROL == CONSTRAINED_TRANSPORT)
  // RBoxDefine (0, nx1_tot-1, 0, nx2_tot-1, 0, nx3_tot-1, &bbox);
  CT_ComputeCharge(data, &boxt, grid);
  #endif

// -- 0c. Exchange data between processors -- 

  #ifdef PARALLEL
  ParallelExchange (data->Uc, data->Vs, grid);
  #endif

  // ------------------------------------------------------
  // 2. Main loop on computational domain sides
  // ------------------------------------------------------

  type[0] = lbound[IDIR]*INCLUDE_IDIR; type[1] = rbound[IDIR]*INCLUDE_IDIR;
  type[2] = lbound[JDIR]*INCLUDE_JDIR; type[3] = rbound[JDIR]*INCLUDE_JDIR;
  type[4] = lbound[KDIR]*INCLUDE_KDIR; type[5] = rbound[KDIR]*INCLUDE_KDIR;

mynvtxstart_("Boundary",SILVER);

  for (is = 0; is <= (2*DIMENSIONS - 1); is++){

    if (type[is] == 0) continue;  // No physical boundary or non-active  *
                                   // * dimension: skip                     

    // 2a. Define boundary boxes, sweeping direction.

    ib = 0; ie = grid->ntot[IDIR]-1;
    jb = 0; je = grid->ntot[JDIR]-1;
    kb = 0; ke = grid->ntot[KDIR]-1;

    if (side[is] == X1_BEG) {
      ib  = 0;    ie  = ibeg-1;
    } else if (side[is] == X1_END) {
      ib  = iend+1; ie  = nx1_tot-1;
    } else if (side[is] == X2_BEG) {
      jb  = 0;    je  = jbeg-1;
    } else if (side[is] == X2_END) {
      jb  = jend+1; je  = nx2_tot-1;
    } else if (side[is] == X3_BEG) {
      kb  = 0; ke  = kbeg-1;
    } else if (side[is] == X3_END) {
      kb  = kend+1; ke  = nx3_tot-1;
    }

    RBoxDefine (ib, ie, jb, je, kb, ke, &bbox);
    #pragma acc update device(bbox)
    #pragma omp target update to(bbox)

    // 2b. Apply boundary conditions.

    if (type[is] == OUTFLOW) {
   
      // ----------------------------------------
      // 2b1: [OUTFLOW] Boundary Conditions.
      // ----------------------------------------

      mynvtxstart_("OUTFLOW",YELLOW);
      OutflowBoundary (data->Uc, data->Vs, data->q, &bbox, side[is]);
      mynvtxstop_();

      // Recover normal magnetic field and average it to cell center. 
      // Transverse components are assigned consistently with
      // cell-centered quantities.

      #ifdef STAGGERED_MHD
      CT_FillStaggeredFields (data->Vs, data->q, side[is], grid); 
      #endif  // STAGGERED_MHD 

      #pragma acc wait

      // Correction step for total conserved energy (is this really needed ?)

      #if (PHYSICS == MHD) && (defined STAGGERED_MHD) && (HAVE_ENERGY)
      int i,j,k;
      if (side[is] == X1_BEG) {
        #pragma acc parallel loop collapse(3) seq present(bbox, data)
        BOX_LOOP(bbox,i,j,k){
          double dBy = data->Vs[BX2s](ibeg,j+1,k) - data->Vs[BX2s](ibeg,j,k);
          data->Uc(i,j,k,ENG) = data->Uc(i+1,j,k,ENG) + data->Vs[BX1s](i+1,j,k)*dBy;
        }
      }
      if (side[is] == X1_END) {
        #pragma acc parallel loop collapse(3) seq present(bbox, data)
        BOX_LOOP(bbox,i,j,k){
          double dBy = data->Vs[BX2s](iend,j+1,k) - data->Vs[BX2s](iend,j,k);
          data->Uc(i,j,k,ENG) = data->Uc(i-1,j,k,ENG) - data->Vs[BX1s](i,j,k)*dBy;
        }
      }
      if (side[is] == X2_BEG) {
        #pragma acc parallel loop collapse(3) seq present(bbox, data)
        BOX_LOOP(bbox,i,j,k){
          double dBx = data->Vs[BX1s](i+1,jbeg,k) - data->Vs[BX1s](i,jbeg,k);
          data->Uc(i,j,k,ENG) = data->Uc(i,j+1,k,ENG) + data->Vs[BX2s](i,j+1,k)*dBx;
        }
      }
      if (side[is] == X2_END) {
        #pragma acc parallel loop collapse(3) seq present(bbox, data)
        BOX_LOOP(bbox,i,j,k){
          double dBx = data->Vs[BX1s](i+1,jend,k) - data->Vs[BX1s](i,jend,k);
          data->Uc(i,j,k,ENG) = data->Uc(i,j-1,k,ENG) - data->Vs[BX2s](i,j,k)*dBx;
        }
      }
      #endif

    }else if (  (type[is] == REFLECTIVE)
             || (type[is] == AXISYMMETRIC)
             || (type[is] == EQTSYMMETRIC)){

      // ----------------------------------------
      // 2b2. [REFLECTIVE/AXISYMMETRIC/EQTSYMMETRIC]
      //       Boundary Conditions.
      // ----------------------------------------

      FlipSign (side[is], type[is], vsign);
      #pragma acc update device (vsign[:NVAR])
      #pragma omp target update to(vsign[:NVAR])

      ReflectiveBoundary(data->Uc, data->Vs, vsign, &bbox, side[is]);

    }else if (type[is] == PERIODIC){  // -- Periodic B.C. (serial or 1 proc) -- 

      // ----------------------------------------
      // 2b3. [PERIODIC] Boundary Conditions.
      //      Assigned only if the direction is not 
      //      parallel (par_dim[is/2] == NO).
      //      NOTE: for staggered meshes we overwrite 
      //      the leftmost face value with the 
      //      rightmost one, as done in parallel.
      // ----------------------------------------

       if (!par_dim[is/2]) {
         mynvtxstart_("PERIODIC",YELLOW);
         PeriodicBoundary(data->Uc, data->Vs, &bbox, side[is], grid);
         mynvtxstop_();
       }

    }else if (type[is] == ROTATED) { 

      // --------------------------------------------------
      // 6d. [ROTATED] Boundary Conditions.
      // --------------------------------------------------

      RotateBoundary (data->Uc, data->Vs, &bbox, side[is], grid);

    }else if (type[is] == USERDEF) {

      // --------------------------------------------------
      // 6f. [USERDEF] Boundary Conditions.
      //      NOTE: this section should not be called since 
      //      a userdef b.c. does not typically define the 
      //      conservative variables in the ghost zones.
      //      Since U is only used to compute the Laplacian,
      //      we replace a USERDEF boundary with an
      //      OUTFLOW b.c. (see section 6a)
      // --------------------------------------------------

      #ifdef GLM_MHD
      int i,j,k;
      BOX_LOOP(bbox, i,j,k) {
        data->Uc(i,j,k,PSI_GLM) = 0.0;
        #ifdef PHI_GLM
        data->Uc(i,j,k,PHI_GLM) = 0.0;
        #endif
      }
      #endif

      UserDefBoundary (data, &bbox, side[is], grid);

      #ifdef STAGGERED_MHD
      CT_FillStaggeredFields (data->Vs, data->q, side[is], grid); 
      #endif

/*
      double vp[NVAR], up[NVAR];

      BOX_LOOP(bbox, i,j,k){
        NVAR_LOOP(nv) vp[nv] = data->Vpc(i,j,k,nv);
        PrimToCons(vp, up);
        NVAR_LOOP(nv) data->Upc(i,j,k,nv) = up[nv];
      }
       
      for (k = bbox.kbeg; k <= bbox.kend; k++){
      for (j = bbox.jbeg; j <= bbox.jend; j++){
      for (i = bbox.jbeg; j <= bbox.jend; j++){
        DIM_EXPAND(d2ux = HO_Average4D<IDIR>(data->Upc,i,j,k,nv);   ,
                   d2uy = HO_Average4D<JDIR>(data->Upc,i,j,k,nv);   ,
                   d2uz = HO_Average4D<KDIR>(data->Upc,i,j,k,nv);)

        data->Uc(i,j,k,nv) = data->Upc(i,j,k,nv) + (d2ux + d2uy + d2uz)/24.0;
      }}}
*/

    } // end if (type[is]) 

    // ------------------------------------------------------
    // 7. Redefine cell-center field in ghost zones from
    // staggered components.
    // Note that this defines cell-centered data->Uc.
    // ------------------------------------------------------

    #ifdef STAGGERED_MHD
    mynvtxstart_("CT_AverageStaggeredFields",BLUE);
    CT_AverageStaggeredFields (data, 0, &bbox, grid);
    mynvtxstop_();
    #endif

  } // end loop on sides 

mynvtxstop_();

  // ------------------------------------------------------
  // 8. Compute charge for the next time level
  // ------------------------------------------------------

  #if (PHYSICS == ResRMHD) && (DIVE_CONTROL == CONSTRAINED_TRANSPORT)
  RBoxDefine (0, nx1_tot-1, 0, nx2_tot-1, 0, nx3_tot-1, &bbox);

  CT_ComputeCharge (data, &bbox, grid);
  #endif

  #pragma acc exit data delete(vsign[:NVAR], bbox, boxt)
  #pragma omp target exit data map(delete: vsign[:NVAR], bbox, boxt)
}

