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
  SHEARING or USERDEF.
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

  \author A. Mignone (andrea.mignone@unito.it)\n
          M. Rossazza
  \date   Feb 27, 2026
*/
////////////////////////////////////////////////////////////////////////
#include "gpluto.hpp"

////////////////////////////////////////////////////////////////////////
/*!
  Set boundary conditions on one or more sides of the computational
  domain.
  
  \param [in,out] data  pointer to PLUTO Data structure containing the
  solution array (including centered and staggered
  fields)
  \param [in]  iarr     an integer determining on which array b.c. are
  imposed
  \param [in]  grid     pointer to grid structure.
*/
////////////////////////////////////////////////////////////////////////
void Boundary (Data *data, int iarr, Grid *grid)
{
  int  is;
  int  side[6] = {X1_BEG, X1_END, X2_BEG, X2_END, X3_BEG, X3_END};
  int  type[6], sbeg, send, vsign[NVAR];
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
  
  const Ary4D& Vc      = (iarr == 0 ? data->Vc: data->Vpc);  
  const Ary3D (&Vs)[6] = (iarr == 0 ? data->Vs: data->Vps);

  RBox bbox;
  #pragma acc enter data copyin(bbox) 

  #pragma acc enter data create (vsign[:NVAR])
  #pragma omp target enter data map(to: vsign)

  // 0. Check the number of processors in each direction

  DIM_EXPAND(par_dim[0] = grid->nproc[IDIR] > 1;  ,
             par_dim[1] = grid->nproc[JDIR] > 1;  ,
             par_dim[2] = grid->nproc[KDIR] > 1;)

  RBoxDefine (grid->lbeg, grid->lend, &bbox);
  #pragma acc update device(bbox)
  #pragma omp target enter data map(to: bbox)

  //#pragma omp target update to(bbox)

  // 1. Call userdef internal boundary with side == 0, 
  //    compute charge if necessary.

  #if INTERNAL_BOUNDARY == YES
  UserDefBoundary (data, &bbox, iarr == 0 ? 0:-1, grid);
  #endif

  #if (PHYSICS == ResRMHD) && (DIVE_CONTROL == CONSTRAINED_TRANSPORT)
  CT_ComputeCharge(data, &bbox, grid);
  #endif

  // 2. Exchange data between processors

  mynvtxstart_("ParallelExcahnge",GREEN);

  #ifdef PARALLEL
  ParallelExchange (Vc, Vs, grid);
  #if (PHYSICS == ResRMHD) && (DIVE_CONTROL == CONSTRAINED_TRANSPORT)
  ParallelExchangeDouble(data->q, grid);
  #endif
  #endif

  mynvtxstop_();

  // 3. Boundaries are imposed on ALL sides: a loop from 
  //    sbeg = 0 to send = 2*3 - 1 is performed. 

  sbeg = 0;
  send = 2*3 - 1;

  // -------------------------------------------
  // 4. Main loop on computational domain sides
  // -------------------------------------------

  type[0] = lbound[IDIR]*INCLUDE_IDIR; type[1] = rbound[IDIR]*INCLUDE_IDIR;
  type[2] = lbound[JDIR]*INCLUDE_JDIR; type[3] = rbound[JDIR]*INCLUDE_JDIR;
  type[4] = lbound[KDIR]*INCLUDE_KDIR; type[5] = rbound[KDIR]*INCLUDE_KDIR;

mynvtxstart_("Boundary",SILVER);

  for (is = sbeg; is <= send; is++){

    if (type[is] == 0) continue;  // No physical boundary or non-active  
                                  // dimension: skip                     

    // 4a. Define boundary boxes, sweeping direction.

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

    // 4b. Apply boundary conditions.

    if (type[is] == SINGULAR_AXIS){  // -- Singular axis condition -- 

      // 4c. [SINGULAR_AXIS] Boundary Conditions.
      
      #if GEOMETRY == CYLINDRICAL || (GEOMETRY == SPHERICAL && DIMENSIONS == 3)
      SingularAxisBoundary (Vc, Vs, &bbox, side[is], grid);
      #endif

    }else if (type[is] == OUTFLOW) {

      // 4d1. [OUTFLOW] Boundary Conditions.

      mynvtxstart_("OUTFLOW",YELLOW);
      OutflowBoundary (Vc, Vs, data->q, &bbox, side[is]);
      mynvtxstop_();

      // 4d2. Recover normal magnetic field and average
      //      it to cell center. Transverse components
      //      are assigned consistently with
      //      cell-centered quantities.

      #ifdef STAGGERED_MHD
      CT_FillStaggeredFields (Vs, data->q, side[is], grid); 
      #endif  // STAGGERED_MHD 

      #pragma acc wait

    }else if (  (type[is] == REFLECTIVE)
             || (type[is] == AXISYMMETRIC)
             || (type[is] == EQTSYMMETRIC)){

      // 4e. [REFLECTIVE/AXISYMMETRIC/EQTSYMMETRIC]
      //     Boundary Conditions.

      FlipSign (side[is], type[is], vsign);
      #pragma acc update device (vsign[:NVAR])
      #pragma omp target update to(vsign[:NVAR])

      mynvtxstart_("REFLECTIVE",YELLOW);
      ReflectiveBoundary(Vc , Vs , vsign, &bbox, side[is]);
      mynvtxstop_();

    }else if (type[is] == PERIODIC){  // -- Periodic B.C. (serial or 1 proc) -- 

      // 4f. [PERIODIC] Boundary Conditions.
      //     Assigned only if the direction is not parallel
      //     (par_dim[is/2] == NO).
      //     NOTE: for staggered meshes we overwrite the
      //     leftmost face value with the rightmost one, as
      //     done in parallel.

      if (!par_dim[is/2]) {
        mynvtxstart_("PERIODIC",YELLOW);
        PeriodicBoundary(Vc , Vs, &bbox, side[is], grid);
        mynvtxstop_();
      }
    }else if (type[is] == ROTATED) { 

      // 4g. [ROTATED] Boundary Conditions.

      RotateBoundary (Vc, Vs, &bbox, side[is], grid);

    }else if (type[is] == USERDEF) {

      // 4h. [USERDEF] Boundary Conditions.

      #ifdef GLM_MHD  
      int i,j,k;    
      BOX_LOOP(bbox,i,j,k) {
        Vc(i,j,k,PSI_GLM) = 0.0;
        #ifdef PHI_GLM
        Vc(i,j,k,PHI_GLM) = 0.0;
        #endif
      }
      #endif

      UserDefBoundary (data, &bbox, side[is], grid);

      // Assign normal component of staggered B 
      // using the div.B = 0 condition         

      #ifdef STAGGERED_MHD
      CT_FillStaggeredFields (Vs, data->q, side[is], grid); 
      #endif

    } // end if (type[is]) 

    // 5. Redefine cell-center field in ghost zones from
    //    staggered components.
    //    Note that this defines cell-centered data->Uc
    //    and then we need to further copy on data->Vc.

    #ifdef STAGGERED_MHD
    #ifdef HIGH_ORDER
    if (type[is] != OUTFLOW) continue;
    #endif
    int i,j,k,nv;

    mynvtxstart_("CT_AverageStaggeredFields",BLUE);
    CT_AverageStaggeredFields (data, 0, &bbox, grid);
    mynvtxstop_();

    // have to move this here from inside the kernel for a bug on Leonardo
    double u[NVAR]; 

    #pragma acc parallel loop collapse(3) present(data, grid) private(u[:NVAR])
    #pragma omp target teams distribute parallel for collapse(3) private(u) //map(present: data[:1])
    for(k = kb; k <= ke; k++){
    for(j = jb; j <= je; j++){
    for(i = ib; i <= ie; i++){
      #pragma acc loop seq
      NVAR_LOOP(nv) u[nv] = data->Uc(i,j,k,nv);
      #if GEOMETRY != CARTESIAN
      Metric metric;
      MetricSet(metric, grid->rC[i], grid->sC[j]);
      ConsToOrtho (u, metric);
      #endif

      DIM_EXPAND(Vc(i,j,k,BX1) = u[BX1];  ,
                 Vc(i,j,k,BX2) = u[BX2];  ,
                 Vc(i,j,k,BX3) = u[BX3];)
      #if (PHYSICS == ResRMHD) && (DIVE_CONTROL == CONSTRAINED_TRANSPORT)
      DIM_EXPAND(Vc(i,j,k,EX1) = u[EX1];  ,
                 Vc(i,j,k,EX2) = u[EX2];  ,
                 Vc(i,j,k,EX3) = u[EX3];)
      #endif
    }}}

    #endif

  } // end loop on sides 
mynvtxstop_();

  // ------------------------------------------------------
  // 6. Compute charge for the next time level
  // ------------------------------------------------------

  #if (PHYSICS == ResRMHD) && (DIVE_CONTROL == CONSTRAINED_TRANSPORT)
  RBoxDefine (0, nx1_tot-1, 0, nx2_tot-1, 0, nx3_tot-1, &bbox);
  #pragma acc update device(bbox)
  #pragma omp target update to(bbox)
  CT_ComputeCharge (data, &bbox, grid);
  #endif

  #pragma acc exit data delete(vsign[:NVAR], bbox)
  #pragma omp target exit data map(delete: vsign[:NVAR], bbox)
}


#if GEOMETRY == CYLINDRICAL || (GEOMETRY == SPHERICAL && DIMENSIONS == 3)
////////////////////////////////////////////////////////////////////////
void SingularAxisBoundary (const Ary4D &Vc, const Ary3D (&Vs)[6], 
                              RBox *box, int side, Grid *grid) 
/*! 
  
*/
////////////////////////////////////////////////////////////////////////
{
  int   i, j, k, nv;
  int   jm;
  int jend = grid->lend[JDIR];

  #if GEOMETRY == CYLINDRICAL
  int im;
  int ibeg = grid->lbeg[IDIR];
  int nx2  = grid->nint[JDIR];
  if (side == X1_BEG){  // X1_BEG boundary 

    // 1a. Ghost zones are filled by copying from the opposite azimuth 
    //    (shift by pi) and flipping the sign of the polar-basis vector 
    //    components (vr --> -vr, vphi --> -vphi and likewise for B)

    #pragma acc parallel loop collapse(3) present(Vc, box)
    BOX_LOOP(*box,i,j,k){

      jm = j + nx2/2;
      if (jm > jend) jm -= nx2;
      im = 2*ibeg - i - 1;

      NVAR_LOOP(nv) Vc(i,j,k,nv) = Vc(im,jm,k,nv); 

      // Odd-symmetry across the origin for in-plane vectors.

      Vc(i,j,k,VX1) *= -1.0;
      Vc(i,j,k,VX2) *= -1.0;

      #if (PHYSICS == MHD) || (PHYSICS == RMHD)
      Vc(i,j,k,BX1) *= -1.0;
      Vc(i,j,k,BX2) *= -1.0;
      #endif

      #if RADIATION
      Vc(i,j,k,FR1) *= -1.0;
      Vc(i,j,k,FR2) *= -1.0;
      #endif

      #ifdef STAGGERED_MHD
      Vs[BX1s](i,j,k) = -Vs[BX1s](im+1,jm, k);  
      #endif
    }

    // 1b. Repeat for staggered fields

    #if (defined STAGGERED_MHD) && INCLUDE_JDIR
    #pragma acc parallel loop collapse(3) present(Vs, box)
    BOXJ_LOOP(*box,i,j,k){
      jm = j + nx2/2;
      if (jm > jend+1) jm -= nx2;   // "+1" since BX2s is staggered along JDIR
      im = 2*ibeg - i - 1;

      Vs[BX2s](i,j,k) = -Vs[BX2s](im, jm, k);
    }
    #endif

    #if (defined STAGGERED_MHD) && INCLUDE_KDIR
    #pragma acc parallel loop collapse(3) present(Vs, box)
    BOXK_LOOP(*box,i,j,k){
      jm = j + nx2/2;
      if (jm > jend) jm -= nx2;
      im = 2*ibeg - i - 1;

      Vs[BX3s](i,j,k) = Vs[BX3s](im, jm, k);
    }
    #endif

  } 
  #endif  // GEOMETRY == CYLINDRICAL

  #if GEOMETRY == SPHERICAL && DIMENSIONS == 3
  int km;
  int nx3  = grid->nint[KDIR];
  int jbeg = grid->lbeg[JDIR];
  int kend = grid->lend[KDIR];
  if (side == X2_BEG || side == X2_END){  // X2_BEG/END boundary 

    // 2a. Ghost zones are filled by copying from the opposite azimuth 
    //    (shift by pi) and flipping the sign of the polar-basis vector 
    //     components (vth --> -vth, vphi --> -vphi and likewise for B)

    #pragma acc parallel loop collapse(3) present(Vc, box)
    BOX_LOOP(*box,i,j,k){
      km = k + nx3/2;
      if (km > kend) km -= nx3;

      if (side == X2_BEG) jm = 2*jbeg - j - 1; // Mirror point (North Pole)
      if (side == X2_END) jm = 2*jend - j + 1; // Mirror point (South Pole)

      NVAR_LOOP(nv) Vc(i,j,k,nv) = Vc(i,jm,km,nv); 

      // Odd-symmetry across the origin for in-plane vectors.

      Vc(i,j,k,VX2) *= -1.0;
      Vc(i,j,k,VX3) *= -1.0;
      #if PHYSICS == MHD
      Vc(i,j,k,BX2) *= -1.0;
      Vc(i,j,k,BX3) *= -1.0;
      #endif

      #if RADIATION
      Vc(i,j,k,FR2) *= -1.0;
      Vc(i,j,k,FR3) *= -1.0;
      #endif

      #ifdef STAGGERED_MHD
      if (side == X2_BEG) Vs[BX2s](i,j  ,k) = -Vs[BX2s](i,jm+1,km); 
      if (side == X2_END) Vs[BX2s](i,j+1,k) = -Vs[BX2s](i,jm  ,km); 
      #endif
    }

    // 2b. Repeat for staggered fields

    #ifdef STAGGERED_MHD
    #pragma acc parallel loop collapse(3) present(Vs, box)
    BOXK_LOOP(*box,i,j,k){
      km = k + nx3/2;
      if (km > kend+1) km -= nx3;   // "+1" since BX3s is staggered along KDIR

      if (side == X2_BEG) jm = 2*jbeg - j - 1; // Mirror point (North Pole)
      if (side == X2_END) jm = 2*jend - j + 1; // Mirror point (South Pole)

      Vs[BX3s](i,j,k) = -Vs[BX3s](i,jm,km);
    }
    #endif

    #ifdef STAGGERED_MHD
    #pragma acc parallel loop collapse(3) present(Vs, box)
    BOXI_LOOP(*box,i,j,k){
      km = k + nx3/2;
      if (km > kend) km -= nx3;

      if (side == X2_BEG) jm = 2*jbeg - j - 1; // Mirror point (North Pole)
      if (side == X2_END) jm = 2*jend - j + 1; // Mirror point (South Pole)

      Vs[BX1s](i,j,k) = Vs[BX1s](i,jm,km);
    }
    #endif
  }
  #endif // GEOMETRY == SPHERICAL

  #ifndef STAGGERED_MHD
  #endif
}
#endif // GEOMETRY == CYLINDRICAL || (GEOMETRY == SPHERICAL && DIMENSIONS == 3)

////////////////////////////////////////////////////////////////////////
/*!
  Implements periodic boundary conditions in serial mode or when
  one processor only handle the periodic direction.
  Note: for staggered fields, we follow the convention of "left" 
  ownership and the rightmost interface must be replaced with its 
  opposite value:
  
  \verbatim
  |...|...|---|---|---|---|...|...|
  C   D   A   B   C   D   A   B   C
  *-------------->*
  *-------------->*                           
  *-------------->*                           
  \endverbatim
  
  \param [in,out] q     a 3D array requiring ghost zone filling
  \param [in]     box   pointer to a RBox structure defining the
  extent of the boundary region
  \param [in]     side  the side of the computational domain.
*/
////////////////////////////////////////////////////////////////////////
void PeriodicBoundary (const Ary4D &Vc,const Ary3D (&Vs)[6], RBox *box,
                       int side, Grid *grid)
{
  int  i,j,k,nv;
  #pragma omp target enter data map(to: box[:1]) //, grid[:1])

  if (side == X1_BEG || side == X1_END){

    int di = (side == X1_BEG ? grid->nint[IDIR]:-grid->nint[IDIR]);

    #pragma acc parallel loop collapse(4) present(box, grid) async(1)
    #pragma omp target teams distribute parallel for collapse(4)
    for (nv = 0; nv <  NVAR; nv++){
      BOX_LOOP(*box, i,j,k) Vc(i,j,k,nv) = Vc((i+di),j,k,nv);
    }

    #ifdef STAGGERED_MHD 
    #if INCLUDE_IDIR
    if (side == X1_BEG){
      #pragma acc parallel loop collapse(3) present(box, grid) async(2)
      #pragma omp target teams distribute parallel for collapse(3)
      BOX_LOOP(*box, i,j,k) {
        CT_ASSIGN(Vs[BX1s](i,j,k) = Vs[BX1s](i+di,j,k);  ,
                  Vs[EX1s](i,j,k) = Vs[EX1s](i+di,j,k);)
      }
    }else{
      #pragma acc parallel loop collapse(3) present(box, grid) async(2)
      #pragma omp target teams distribute parallel for collapse(3)
      BOXI_LOOP(*box, i,j,k) {
        CT_ASSIGN(Vs[BX1s](i,j,k) = Vs[BX1s](i+di,j,k);  ,
                  Vs[EX1s](i,j,k) = Vs[EX1s](i+di,j,k);)
      }
    }
    #endif

    #if INCLUDE_JDIR
    #pragma acc parallel loop collapse(3) present(box, grid) async(3)
    #pragma omp target teams distribute parallel for collapse(3) 
    BOXJ_LOOP(*box, i,j,k) {
      CT_ASSIGN(Vs[BX2s](i,j,k) = Vs[BX2s]((i+di),j,k); ,
                Vs[EX2s](i,j,k) = Vs[EX2s]((i+di),j,k);)
    } 
    #endif

    #if INCLUDE_KDIR
    #pragma acc parallel loop collapse(3) present(box, grid) async(4)
    #pragma omp target teams distribute parallel for collapse(3) 
    BOXK_LOOP(*box, i,j,k) {
      CT_ASSIGN(Vs[BX3s](i,j,k) = Vs[BX3s]((i+di),j,k); ,
                Vs[EX3s](i,j,k) = Vs[EX3s]((i+di),j,k);)
    } 
    #endif
    #endif // STAGGERED_MHD 

    #pragma acc wait

  }else if (side == X2_BEG || side == X2_END){

    int dj = (side == X2_BEG ? grid->nint[JDIR]:-grid->nint[JDIR]);

    #pragma acc parallel loop collapse(4) present(box, grid) async(1)
    #pragma omp target teams distribute parallel for collapse(4) 
    for (nv = 0; nv <  NVAR; nv++){
      BOX_LOOP(*box, i,j,k) Vc(i,j,k,nv) = Vc(i,(j+dj),k,nv);
    }

    #ifdef STAGGERED_MHD
    #if INCLUDE_IDIR
    #pragma acc parallel loop collapse(3) present(box, grid)  async(2)
    #pragma omp target teams distribute parallel for collapse(3) 
    BOXI_LOOP(*box,i,j,k) {
      CT_ASSIGN(Vs[BX1s](i,j,k) = Vs[BX1s](i,(j+dj),k);  ,
                Vs[EX1s](i,j,k) = Vs[EX1s](i,(j+dj),k);)
    }
    #endif

    #if INCLUDE_JDIR
    if (side == X2_BEG){
      #pragma acc parallel loop collapse(3) present(box, grid) async(3)
      #pragma omp target teams distribute parallel for collapse(3) 
      BOX_LOOP(*box, i,j,k) {
        CT_ASSIGN(Vs[BX2s](i,j,k) = Vs[BX2s](i,j+dj,k);  ,
                  Vs[EX2s](i,j,k) = Vs[EX2s](i,j+dj,k);)
      }
    }else{
      #pragma acc parallel loop collapse(3) present(box, grid) async(3)
      #pragma omp target teams distribute parallel for collapse(3) 
      BOXJ_LOOP(*box, i,j,k) {
        CT_ASSIGN(Vs[BX2s](i,j,k) = Vs[BX2s](i,j+dj,k);  ,
                  Vs[EX2s](i,j,k) = Vs[EX2s](i,j+dj,k);)
      }
    }
    #endif

    #if INCLUDE_KDIR
    #pragma acc parallel loop collapse(3) present(box, grid) async(4)
    #pragma omp target teams distribute parallel for collapse(3) 
    BOXK_LOOP(*box,i,j,k) {
      CT_ASSIGN(Vs[BX3s](i,j,k) = Vs[BX3s](i,(j+dj),k); ,
                Vs[EX3s](i,j,k) = Vs[EX3s](i,(j+dj),k);)
    }
    #endif
    #endif // STAGGERED_MHD 

    #pragma acc wait

  }else if (side == X3_BEG || side == X3_END){

    int dk = (side == X3_BEG ? grid->nint[KDIR]:-grid->nint[KDIR]);

    #pragma acc parallel loop collapse(4) present(box, grid) async(1)
    #pragma omp target teams distribute parallel for collapse(4) 
    for (nv = 0; nv < NVAR; nv++){
      BOX_LOOP(*box, i,j,k) Vc(i,j,k,nv) = Vc(i,j,(k+dk),nv);
    }

    #ifdef STAGGERED_MHD
    #if INCLUDE_IDIR
    #pragma acc parallel loop collapse(3) present(box, grid) async(2)
    #pragma omp target teams distribute parallel for collapse(3) 
    BOXI_LOOP(*box,i,j,k) {
      CT_ASSIGN(Vs[BX1s](i,j,k) = Vs[BX1s](i,j,(k+dk)); ,
                Vs[EX1s](i,j,k) = Vs[EX1s](i,j,(k+dk));)
    }
    #endif

    #if INCLUDE_JDIR
    #pragma acc parallel loop collapse(3) present(box, grid) async(3)
    #pragma omp target teams distribute parallel for collapse(3) 
    BOXJ_LOOP(*box, i,j,k) {
      CT_ASSIGN(Vs[BX2s](i,j,k) = Vs[BX2s](i,j,(k+dk)); ,
                Vs[EX2s](i,j,k) = Vs[EX2s](i,j,(k+dk));)
    }
    #endif

    #if INCLUDE_KDIR
    if (side == X3_BEG){
      #pragma acc parallel loop collapse(3) present(box, grid) async(4)
      #pragma omp target teams distribute parallel for collapse(3) 
      BOX_LOOP(*box, i,j,k) {
        CT_ASSIGN(Vs[BX3s](i,j,k) = Vs[BX3s](i,j,k+dk);  ,
                  Vs[EX3s](i,j,k) = Vs[EX3s](i,j,k+dk);)
      }
    }else{
      #pragma acc parallel loop collapse(3) present(box, grid) async(4)
      #pragma omp target teams distribute parallel for collapse(3) 
      BOXK_LOOP(*box, i,j,k) {
        CT_ASSIGN(Vs[BX3s](i,j,k) = Vs[BX3s](i,j,k+dk);  ,
                  Vs[EX3s](i,j,k) = Vs[EX3s](i,j,k+dk);)
      }
    }
    #endif
    #endif // STAGGERED_MHD 

    #pragma acc wait
  }

  #ifndef STAGGERED_MHD
  #endif

}

////////////////////////////////////////////////////////////////////////
/*!
  Impose zero-gradient boundary conditions on 'q' on
  the boundary side specified by 'side'.
  The input array 'q' must not represent the normal component
  of a staggered magnetic fied.
  
  \param [in,out] q     a 3D array requiring ghost zone filling
  \param [in]     box   pointer to a RBox structure defining the
  extent of the boundary region
  \param [in]     side  the side of the computational domain.
*/
////////////////////////////////////////////////////////////////////////
void OutflowBoundary (const Ary4D &Vc, const Ary3D (&Vs)[6], const Ary3D &q,  
                      RBox *box, int side)
{
  int  i, j, k, nv;

  if (side == X1_BEG || side == X1_END){
    int i0 = (side == X1_BEG ? box->iend+1:box->ibeg-1);

    #pragma acc parallel loop collapse(4) present(box) async(1)
    #pragma omp target teams distribute parallel for collapse(4) 
    for (nv = 0; nv <  NVAR; nv++){
      BOX_LOOP(*box, i,j,k) Vc(i,j,k,nv) = Vc(i0,j,k,nv);
    }

    #if (defined STAGGERED_MHD) && INCLUDE_JDIR
    #pragma acc parallel loop collapse(3) present(box) async(2)
    #pragma omp target teams distribute parallel for collapse(3) 
    BOXJ_LOOP(*box, i,j,k) {
      CT_ASSIGN(Vs[BX2s](i,j,k) = Vs[BX2s](i0,j,k);  ,
                Vs[EX2s](i,j,k) = Vs[EX2s](i0,j,k);)
    }
    #endif

    #if (defined STAGGERED_MHD) && INCLUDE_KDIR
    #pragma acc parallel loop collapse(3) present(box) async(3)
    #pragma omp target teams distribute parallel for collapse(3) 
    BOXK_LOOP(*box, i,j,k) {
      CT_ASSIGN(Vs[BX3s](i,j,k) = Vs[BX3s](i0,j,k);  ,
                Vs[EX3s](i,j,k) = Vs[EX3s](i0,j,k);)
    }
    #endif

    #if (PHYSICS == ResRMHD) && (DIVE_CONTROL == CONSTRAINED_TRANSPORT)
    #pragma acc parallel loop collapse(3) present(box) async(4)
    #pragma omp target teams distribute parallel for collapse(3) 
    BOX_LOOP(*box, i,j,k) q(i,j,k) = q(i0,j,k);
    #endif

    #pragma acc wait

  }else if (side == X2_BEG || side == X2_END){

    int j0 = (side == X2_BEG ? box->jend+1:box->jbeg-1);

    #pragma acc parallel loop collapse(4) present(box) async(1)
    #pragma omp target teams distribute parallel for collapse(4) 
    for (nv = 0; nv <  NVAR; nv++){
      BOX_LOOP(*box, i,j,k) Vc(i,j,k,nv) = Vc(i,j0,k,nv);
    }

    #if (defined STAGGERED_MHD) && INCLUDE_IDIR
    #pragma acc parallel loop collapse(3) present(box) async(2)
    #pragma omp target teams distribute parallel for collapse(3) 
    BOXI_LOOP(*box,i,j,k){
      CT_ASSIGN(Vs[BX1s](i,j,k) = Vs[BX1s](i,j0,k);  ,
                Vs[EX1s](i,j,k) = Vs[EX1s](i,j0,k);)
    }
    #endif

    #if (defined STAGGERED_MHD) && INCLUDE_KDIR
    #pragma acc parallel loop collapse(3) present(box) async(3)
    #pragma omp target teams distribute parallel for collapse(3) 
    BOXK_LOOP(*box,i,j,k) {
      CT_ASSIGN(Vs[BX3s](i,j,k) = Vs[BX3s](i,j0,k);  ,
                Vs[EX3s](i,j,k) = Vs[EX3s](i,j0,k);)
    }
    #endif

    #if (PHYSICS == ResRMHD) && (DIVE_CONTROL == CONSTRAINED_TRANSPORT)
    #pragma acc parallel loop collapse(3) present(box) async(4)
    #pragma omp target teams distribute parallel for collapse(3) 
    BOX_LOOP(*box, i,j,k) q(i,j,k) = q(i,j0,k);
    #endif

    #pragma acc wait

  }else if (side == X3_BEG || side == X3_END){

    int k0 = (side == X3_BEG ? box->kend+1:box->kbeg-1);

    #pragma acc parallel loop collapse(4) present(box) async(1)
    #pragma omp target teams distribute parallel for collapse(4) 
    for (nv = 0; nv < NVAR; nv++){
      BOX_LOOP(*box, i,j,k) Vc(i,j,k,nv) = Vc(i,j,k0,nv);
    }

    #if (defined STAGGERED_MHD) && (INCLUDE_IDIR)
    #pragma acc parallel loop collapse(3) present(box) async(2)
    #pragma omp target teams distribute parallel for collapse(3) 
    BOXI_LOOP(*box,i,j,k) {
      CT_ASSIGN(Vs[BX1s](i,j,k) = Vs[BX1s](i,j,k0);  ,
                Vs[EX1s](i,j,k) = Vs[EX1s](i,j,k0);)
    }
    #endif

    #if (defined STAGGERED_MHD) && (INCLUDE_JDIR)
    #pragma acc parallel loop collapse(3) present(box) async(4)
    #pragma omp target teams distribute parallel for collapse(3) 
    BOXJ_LOOP(*box, i,j,k) {
      CT_ASSIGN(Vs[BX2s](i,j,k) = Vs[BX2s](i,j,k0);  ,
                Vs[EX2s](i,j,k) = Vs[EX2s](i,j,k0);)
    }
    #endif

    #if (PHYSICS == ResRMHD) && (DIVE_CONTROL == CONSTRAINED_TRANSPORT)
    #pragma acc parallel loop collapse(3) present(box) async(6)
    #pragma omp target teams distribute parallel for collapse(3) 
    BOX_LOOP(*box, i,j,k) q(i,j,k) = q(i,j,k0);
    #endif

    #pragma acc wait
  }

  #ifndef STAGGERED_MHD
  #endif
  #if DIVE_CONTROL != CONSTRAINED_TRANSPORT
  #endif
}

////////////////////////////////////////////////////////////////////////
/*!
  Reverse the sign of vector components with respect to axis side.
  Depending on type, one needs to symmetrize or anti-symmetrize:
  
  - REFLECTIVE:    \n
  o   Vn -> -Vn,  Bn -> -Bn,  En ->  En  \n
  o   Vp ->  Vp,  Bp ->  Bp,  Ep -> -Ep  \n
  o   Vt ->  Vt,  Bt ->  Bt,  Et -> -Et
  
  - AXISYMMETRIC:  \n
  o   Vn -> -Vn,  Bn -> -Bn,  En -> -En \n
  o   Vp ->  Vp,  Bp ->  Bp,  Ep ->  Ep \n
  o   Vt -> -Vt,  Bt -> -Bt,  Et -> -Et
  
  - EQTSYMMETRIC:  \n
  o   Vn -> -Vn,  Bn ->  Bn,  En -> -En  \n
  o   Vp ->  Vp,  Bp -> -Bp,  Ep ->  Ep  \n
  o   Vt ->  Vt,  Bt -> -Bt,  Et ->  Et
  
  where (n) is the normal components, (p) and (t)
  are the transverse (or poloidal and toroidal for
  cylindrical and spherical coordinates) components.
  
  If the radiation module is used, the transformation law
  for the radiation flux is the same as for the velocities.
  
  \param [in]  side  boundary side
  \param [in]  type boundary condition type
  \param [out] vsign an array of values (+1 or -1) giving the sign
*/
////////////////////////////////////////////////////////////////////////
void FlipSign (int side, int type, int *vsign)
{
  int  nv, dir;
  #ifdef BX1
  int  BXn = BX1, BXt = BX2, BXb = BX3;
  #endif
  #ifdef EX1
  int  EXn = EX1, EXt = EX2, EXb = EX3;
  #endif

  if (side == X1_BEG || side == X1_END) dir = IDIR;
  if (side == X2_BEG || side == X2_END) dir = JDIR;
  if (side == X3_BEG || side == X3_END) dir = KDIR;

  int  VXn = VX1+dir;
  #if RADIATION
  int  FRn = FR1+dir; 
  #endif

  // Set which variable flips sign

  NVAR_LOOP(nv) vsign[nv] = 1;

  #ifdef BX1
  INDEX_CYCLE (dir, BXn, BXt, BXb);
  #endif
  #ifdef EX1
  INDEX_CYCLE (dir, EXn, EXt, EXb);
  #endif

  vsign[VXn] = -1.0;
  #if (PHYSICS == MHD) || (PHYSICS == RMHD)  || (PHYSICS == ResRMHD) || (PHYSICS == GRMHD)
  vsign[BXn] = -1.0;
  #endif
  #if DUST_FLUID == YES
  vsign[VXn_D] = -1.0;
  #endif
  #if RADIATION
  vsign[FRn] = -1.0;
  #endif

  if (type == REFLECTIVE){

    // vsign[VXn] = -1;
    #ifdef BX1
    vsign[BXn] = -1;
    #endif
    #ifdef EX1
    vsign[EXt] = -1;
    vsign[EXb] = -1;
    #endif

  }else if (type == EQTSYMMETRIC){

    // vsign[VXn] = -1;
    #ifdef BX1
    vsign[BXn] =  1;
    vsign[BXt] = -1;
    vsign[BXb] = -1;
    #ifdef GLM_MHD
    vsign[PSI_GLM] = -1;
    #endif
    #ifdef EX1
    vsign[EXn] = -1;
    #endif
    #endif

  }else if (type == AXISYMMETRIC){

    // vsign[VXn] = -1;
    #if GEOMETRY != CARTESIAN
    vsign[iVPHI] = -1;
    #if (PHYSICS == MHD) || (PHYSICS == RMHD)  || (PHYSICS == ResRMHD) || (PHYSICS == GRMHD)
    vsign[iBPHI] = -1;
    #endif
    #endif

    #if DUST_FLUID == YES
    vsign[VXt_D] = -1;
    #endif

    #ifdef EX1
    vsign[EXn] = -1;
    vsign[EXt] = -1;
    #endif

    #if RADIATION && GEOMETRY != CARTESIAN
    vsign[iFRPHI] = -1;
    #endif
  }
}

////////////////////////////////////////////////////////////////////////
/*!
  Implements periodic boundary conditions in serial mode or when
  one processor only handle the periodic direction.
  
  \param [in,out] q     a 3D array requiring ghost zone filling
  \param [in]     box   pointer to a RBox structure defining the
  extent of the boundary region
  \param [in]     side  the side of the computational domain.
*/
////////////////////////////////////////////////////////////////////////
void ReflectiveBoundary (const Ary4D &Vc,const Ary3D (&Vs)[6], int *vsign, 
                         RBox *box, int side)
{
  int  i, j, k, nv;

  if (side == X1_BEG || side == X1_END){
    int ibeg, iend, i0;

    if (side == X1_BEG) {
      ibeg = box->iend+1;
      i0   = 2*ibeg-1;
    }else{
      iend = box->ibeg-1;
      i0   = 2*iend+1;
    }

    #pragma acc parallel loop collapse(4) present(box, vsign[:NVAR]) async(1)
    #pragma omp target teams distribute parallel for collapse(4) 
    for (nv = 0; nv < NVAR; nv++){
      BOX_LOOP(*box, i,j,k) {
        Vc(i,j,k,nv) = vsign[nv]*Vc((i0-i),j,k,nv);
      }
    }

    #if (defined STAGGERED_MHD) && INCLUDE_IDIR
    if (side == X1_BEG){
      #pragma acc parallel loop collapse(3) present(box, vsign[:NVAR]) async(2)
      #pragma omp target teams distribute parallel for collapse(3) 
      BOX_LOOP(*box, i,j,k) {
        Vs[BX1s](i,j,k) = vsign[BX1]*Vs[BX1s](i0-i+1,j,k);
      }
    }else{
      #pragma acc parallel loop collapse(3) present(box, vsign[:NVAR]) async(2)
      #pragma omp target teams distribute parallel for collapse(3) 
      BOX_LOOP(*box, i,j,k){
         Vs[BX1s](i+1,j,k) = vsign[BX1]*Vs[BX1s](i0-i,j,k);
      }
    }
    #endif

    #if (defined STAGGERED_MHD) && INCLUDE_JDIR
    #pragma acc parallel loop collapse(3) present(box, vsign[:NVAR]) async(3)
    #pragma omp target teams distribute parallel for collapse(3) 
    BOXJ_LOOP(*box, i,j,k) {
      Vs[BX2s](i,j,k) = vsign[BX2]*Vs[BX2s](i0-i,j,k);
    }
    #endif

    #if (defined STAGGERED_MHD) && INCLUDE_KDIR
    #pragma acc parallel loop collapse(3) present(box, vsign[:NVAR]) async(4)
    #pragma omp target teams distribute parallel for collapse(3) 
    BOXK_LOOP(*box, i,j,k) {
      Vs[BX3s](i,j,k) = vsign[BX3]*Vs[BX3s](i0-i,j,k);
    }
    #endif

    #pragma acc wait

  }else if (side == X2_BEG || side == X2_END){

    int jbeg, jend, j0;

    if (side == X2_BEG) {
      jbeg = box->jend+1;
      j0   = 2*jbeg-1;
    }else{
      jend = box->jbeg-1;
      j0   = 2*jend+1;
    }

    #pragma acc parallel loop collapse(4) present(box, vsign[:NVAR]) async(1)
    #pragma omp target teams distribute parallel for collapse(4) 
    for (nv = 0; nv <  NVAR; nv++){
      BOX_LOOP(*box, i,j,k) {
        Vc(i,j,k,nv) = vsign[nv]*Vc(i,j0-j,k,nv);
      }
    }

    #if (defined STAGGERED_MHD) && INCLUDE_IDIR
    #pragma acc parallel loop collapse(3) present(box, vsign[:NVAR]) async(2)
    #pragma omp target teams distribute parallel for collapse(3) 
    BOXI_LOOP(*box, i,j,k) {
      Vs[BX1s](i,j,k) = vsign[BX1]*Vs[BX1s](i,j0-j,k);
    }
    #endif

    #if (defined STAGGERED_MHD) && INCLUDE_JDIR
    if (side == X2_BEG){
      #pragma acc parallel loop collapse(3) present(box, vsign[:NVAR]) async(3)
      #pragma omp target teams distribute parallel for collapse(3) 
      BOX_LOOP(*box, i,j,k) {
        Vs[BX2s](i,j,k) = vsign[BX2]*Vs[BX2s](i,j0-j+1,k);
      }
    }else{
      #pragma acc parallel loop collapse(3) present(box, vsign[:NVAR]) async(3)
      #pragma omp target teams distribute parallel for collapse(3) 
      BOX_LOOP(*box, i,j,k){
         Vs[BX2s](i,j+1,k) = vsign[BX2]*Vs[BX2s](i,j0-j,k);
      }
    }
    #endif

    #if (defined STAGGERED_MHD) && INCLUDE_KDIR
    #pragma acc parallel loop collapse(3) present(box, vsign[:NVAR]) async(4)
    #pragma omp target teams distribute parallel for collapse(3) 
    BOXK_LOOP(*box, i,j,k) {
      Vs[BX3s](i,j,k) = vsign[BX3]*Vs[BX3s](i,j0-j,k);
    }
    #endif

    #pragma acc wait

  }else if (side == X3_BEG || side == X3_END){

    int kbeg, kend, k0;

    if (side == X3_BEG) {
      kbeg = box->kend+1;
      k0   = 2*kbeg-1;
    }else{
      kend = box->kbeg-1;
      k0   = 2*kend+1;
    }

    #pragma acc parallel loop collapse(4) present(box, vsign[:NVAR]) async(1)
    #pragma omp target teams distribute parallel for collapse(4) 
    for (nv = 0; nv <  NVAR; nv++){
      BOX_LOOP(*box, i,j,k) {
        Vc(i,j,k,nv) = vsign[nv]*Vc(i,j,k0-k,nv);
      }
    }

    #if (defined STAGGERED_MHD) && INCLUDE_IDIR
    #pragma acc parallel loop collapse(3) present(box, vsign[:NVAR]) async(2)
    #pragma omp target teams distribute parallel for collapse(3) 
    BOXI_LOOP(*box, i,j,k) {
      Vs[BX1s](i,j,k) = vsign[BX1]*Vs[BX1s](i,j,k0-k);
    }
    #endif

    #if (defined STAGGERED_MHD) && INCLUDE_JDIR
    #pragma acc parallel loop collapse(3) present(box, vsign[:NVAR]) async(3)
    #pragma omp target teams distribute parallel for collapse(3) 
    BOXJ_LOOP(*box, i,j,k) {
      Vs[BX2s](i,j,k) = vsign[BX2]*Vs[BX2s](i,j,k0-k);
    }
    #endif

    #if (defined STAGGERED_MHD) && INCLUDE_KDIR
    if (side == X3_BEG){
      #pragma acc parallel loop collapse(3) present(box, vsign[:NVAR]) async(4)
      #pragma omp target teams distribute parallel for collapse(3) 
      BOX_LOOP(*box, i,j,k) {
        Vs[BX3s](i,j,k) = vsign[BX3]*Vs[BX3s](i,j,k0-k+1);
      }
    }else{
      #pragma acc parallel loop collapse(3) present(box, vsign[:NVAR]) async(4)
      #pragma omp target teams distribute parallel for collapse(3) 
      BOX_LOOP(*box, i,j,k){
        Vs[BX3s](i,j,k+1) = vsign[BX3]*Vs[BX3s](i,j,k0-k);
      }
    }
    #endif

    #pragma acc wait
  }

  #ifndef STAGGERED_MHD
  #endif
}

