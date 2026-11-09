////////////////////////////////////////////////////////////////////////
/*!
  \file
  \brief HLL Riemann solver for HD.

  Solve the Riemann problem for the HD equations using the
  single-state HLL solver by Toro.

  On input, this function takes left and right primitive state vectors
  \c data->sweep->vL and \c data->sweep->vR at zone edge i-1/2.
  On output, return flux and pressure vectors at the same interface
  \c i-1/2 (note that the \c i refers to \c i-1/2).

  Also during this step, compute maximum wave propagation speed (cmax)
  for  explicit time step computation.

  \b Reference:
   - "Riemann Solver and Numerical Methods for Fluid Dynamics"
      by E.F. Toro (Chapter 10)
   - "On Godunov-Type Method near Low Densities"
     by B. Einfeldt, C.D. Munz, P.L. Roe, JCP 92, 273-295 (1991)

  \authors A. Mignone  (andrea.mignone@unito.it)  \n
           M. Rossazza (marco.rossazza@unito.it)  \n
           V. Berta    (vittoria.berta@unito.it)  \n
           A. Suriano  (alessio.suriano@unito.it)
  \date    Feb 12, 2026
*/
////////////////////////////////////////////////////////////////////////
#include "gpluto.hpp"

template<int dir>
void HLL_Solver (Data *data, timeStep *Dts, Grid *grid, RBox *box);

////////////////////////////////////////////////////////////////////////
/*!
*/
////////////////////////////////////////////////////////////////////////
void HLL_Solver_T (Data *data, timeStep *Dts, Grid *grid, RBox *box, int dir)
{
  if      (dir == IDIR) HLL_Solver<IDIR>(data, Dts, grid, box);
  else if (dir == JDIR) HLL_Solver<JDIR>(data, Dts, grid, box);
  else if (dir == KDIR) HLL_Solver<KDIR>(data, Dts, grid, box);
}

////////////////////////////////////////////////////////////////////////
/*!
  Solve Riemann problem for the HD equations
  using the HLL Riemann solver.
*/
////////////////////////////////////////////////////////////////////////
template<int dir>
void HLL_Solver (Data *data, timeStep *Dts, Grid *grid, RBox *box)
{
  int      i, j, k, nv;
  int      ibeg, jbeg, kbeg;
  int      iend, jend, kend;
  int      ni, nj;
  uint32_t err_kernel = 0;  // Cumulative kernel error flag 

  double *x1 = (dir == IDIR ? grid->xL[IDIR]:grid->xC[IDIR]); 
  double *x2 = (dir == JDIR ? grid->xL[JDIR]:grid->xC[JDIR]); 
  double *x3 = (dir == KDIR ? grid->xL[KDIR]:grid->xC[KDIR]); 

  RunConfig &run_config = *(data->run_config);

  // ------------------------------------------------------
  // 1. Compute strides and box permutation
  // ------------------------------------------------------

  ARRAY_STRIDES (grid, dir, ni, nj);

  // Cycle on directional indexes to preserve sweep-> ordering

  RBoxPermuteIndices<dir>(box, ibeg, iend, jbeg, jend, kbeg, kend);
  iend++;  // Riemann solver need one more interface value

  // ------------------------------------------------------
  // 2. Loop over domain box and solve Riemann problem
  //    in the direction given by dir.
  // ------------------------------------------------------

  #pragma omp target teams distribute parallel for collapse(3)
  #pragma acc parallel loop collapse(3) present(data, run_config, Dts, grid) \
                                        copy(err_kernel)
  for (k = kbeg; k <= kend; k++){
  for (j = jbeg; j <= jend; j++){
  for (i = ibeg; i <= iend; i++){
      
    // Compute offsets and assign local pointers

    long int offset  = ni*(j + nj*k);
    long int offset1 = NVAR*ni*(j + nj*k);

    Ary1D cmax(&Dts->cmax[offset],ni);
    Ary1D cs2(&data->sweep->cs2[offset],ni);
    Ary1D press(&data->sweep->press[offset],ni);
    Ary2D vL(&data->sweep->vL[offset1],ni,NVAR);
    Ary2D vR(&data->sweep->vR[offset1],ni,NVAR);
    Ary2D flux(&data->sweep->flux[offset1],ni,NVAR);

    // Private variable declaration

    int ti, tj, tk;
    double ul[NVAR], ur[NVAR];
    double vl[NVAR], vr[NVAR];
    double fL[NVAR], fR[NVAR];
    double Sl, Sr;
    double pl, pr;

    State stateL, stateR;

    // Set state and "true" indices

    if constexpr (dir == IDIR)  {ti = i; tj = j; tk = k;}
    if constexpr (dir == JDIR)  {ti = j; tj = i; tk = k;}
    if constexpr (dir == KDIR)  {ti = j; tj = k; tk = i;}

    NFLUID_LOOP(nv) {
      vl[nv] = vL(i,nv);
      vr[nv] = vR(i,nv);
    }

    data->err(ti,tj,tk) = 0;   // reset error 

    // 2a. Compute sound speed, cons variables and fluxes;
    //     compute HLL speeds.

    PrimToCons (vl, ul, run_config);
    PrimToCons (vr, ur, run_config);

    STATE_SET(stateL, vl, ul, NULL, x1[ti], x2[tj], x3[tk]);
    STATE_SET(stateR, vr, ur, NULL, x1[ti], x2[tj], x3[tk]);

    SoundSpeed2 (stateL, run_config);
    SoundSpeed2 (stateR, run_config);
    double a2L = stateL.a2;
    double a2R = stateR.a2;

    Flux<dir> (stateL, fL, pl);
    Flux<dir> (stateR, fR, pr);

    HLL_Speed<dir> (stateL, stateR, Sl, Sr);

    // 2b. Get max propagation speed for dt computation
  
    double bL = MIN(Sl, 0.0);
    double bR = MAX(0.0, Sr);

    cmax(i) = MAX(fabs(Sl), fabs(Sr));    
    cs2(i)  = (bR*a2R - bL*a2L)/(bR - bL);

    // 2c. Compute Flux if SL > 0 or SR < 0

    if (Sl > 0.0){
      #pragma acc loop seq
      for (nv = 0; nv < NFLX; nv++) flux(i,nv) = fL[nv];
      press(i) = pl;

    }else if (Sr < 0.0){
      #pragma acc loop seq
      for (nv = 0; nv < NFLX; nv++) flux(i,nv) = fR[nv];
      press(i) = pr;

    }else{

    // 2d.Compute Flux when  SL < 0 < SR

      double scrh = 1.0 / (Sr - Sl);

      #pragma acc loop seq
      for (nv = 0; nv < NFLX; nv++) {
        flux(i,nv) = Sl*Sr*(ur[nv] - ul[nv]) + Sr*fL[nv] - Sl*fR[nv];
        flux(i,nv) *= scrh;
      }
      press(i) = (Sr*pl - Sl*pr)*scrh;
    }

    // 2e. Compute tracer fluxes

    #if NSCL > 0
    if (flux(i,RHO) > 0.0){
      #pragma acc loop seq
      NSCL_LOOP(nv) flux(i,nv) = flux(i,RHO)*vl[nv];
    }else{
      #pragma acc loop seq
      NSCL_LOOP(nv) flux(i,nv) = flux(i,RHO)*vr[nv];
    }
    #endif

  }}} // End Domain loop (i,j,k)

  // ------------------------------------------------------
  // 3. Error Check / Quit condition
  // ------------------------------------------------------

  #if WARNING_MESSAGES > 0 
  Where ("HLL_Solver()", data, err_kernel, box, grid);
  #endif

  #if WARNING_MESSAGES >= 0 
  if (err_kernel & ERR_QUIT) QuitPLUTO("from HLL_Solver()",err_kernel);
  #endif
}
