////////////////////////////////////////////////////////////////////////
/*!
  \file
  \brief HLL Riemann solver for RHD.

  \authors A. Mignone  (andrea.mignone@unito.it)
           M. Rossazza (marco.rossazza@unito.it)
           A. Suriano  (alessio.suriano@unito.it)
  \date    Dec 17, 2025
*/
////////////////////////////////////////////////////////////////////////
#include "gpluto.hpp"

template<int dir>
void HLL_Solver (Data *, timeStep *, Grid *, RBox *);

////////////////////////////////////////////////////////////////////////
void HLL_Solver_T (Data *data, timeStep *Dts, Grid *grid, RBox *box, int dir)
/*!
*/
////////////////////////////////////////////////////////////////////////
{
  if      (dir == IDIR) HLL_Solver<IDIR>(data, Dts, grid, box);
  else if (dir == JDIR) HLL_Solver<JDIR>(data, Dts, grid, box);
  else if (dir == KDIR) HLL_Solver<KDIR>(data, Dts, grid, box);
}

////////////////////////////////////////////////////////////////////////
template<int dir>
void HLL_Solver (Data *data, timeStep *Dts, Grid *grid, RBox *box)
/*!
  Solve Riemann problem for the adiabatic/isothermal RHD equations
  using the HLL Riemann solver.
*/
////////////////////////////////////////////////////////////////////////
{
  int      i, j, k, nv;
  int      ibeg, jbeg, kbeg;
  int      iend, jend, kend;
  int  ni, nj;
  uint32_t err_kernel = 0;  // Cumulative kernel error flag 

  double *x1 = (dir == IDIR ? grid->xL[IDIR]:grid->xC[IDIR]); 
  double *x2 = (dir == JDIR ? grid->xL[JDIR]:grid->xC[JDIR]); 
  double *x3 = (dir == KDIR ? grid->xL[KDIR]:grid->xC[KDIR]); 

  RunConfig &run_config = *(data->run_config);

  // ------------------------------------------------------
  // 1. Compute strides and box permutation
  // ------------------------------------------------------

  ARRAY_STRIDES (grid, dir, ni, nj);

  // Cycle on directional indexes to preserve sweep ordering

  RBoxPermuteIndices<dir>(box, ibeg, iend, jbeg, jend, kbeg, kend);
  iend++;  // Riemann solver need one more interface value

  // ------------------------------------------------------
  // 2. Loop over domain box and solve Riemann problem
  //    in the direction given by dir.
  // ------------------------------------------------------

  #pragma omp target teams distribute parallel for collapse(3) 
  #pragma acc parallel loop collapse(3) present(data, run_config, Dts, grid)\
                                        copy(err_kernel)
  for (k = kbeg; k <= kend; k++){
  for (j = jbeg; j <= jend; j++){
  for (i = ibeg; i <= iend; i++){

    // Compute offsets and assign local pointers
      
    long int offset = ni*(j + nj*k);
    long int offset1 = NVAR*ni*(j + nj*k);

    Ary1D press(&data->sweep->press[offset],ni);
    Ary1D cmax(&Dts->cmax[offset],ni);
    Ary2D vL(&data->sweep->vL[offset1],ni,NVAR);
    Ary2D vR(&data->sweep->vR[offset1],ni,NVAR);
    Ary2D flux(&data->sweep->flux[offset1],ni,NVAR);
    Ary2D fL(&data->sweep->fL[offset1],ni,NVAR);
    Ary2D fR(&data->sweep->fR[offset1],ni,NVAR);

    // Private variable declaration

    int ti, tj, tk;
    double ul[NVAR], ur[NVAR];
    double vl[NVAR], vr[NVAR];
    double fl[NVAR], fr[NVAR];
    double Sl, Sr, pl, pr;
    State stateL, stateR;

    if constexpr (dir == IDIR) {ti = i; tj = j; tk = k;}
    if constexpr (dir == JDIR) {ti = j; tj = i; tk = k;}
    if constexpr (dir == KDIR) {ti = j; tj = k; tk = i;}

    for (nv = 0; nv < NVAR; nv++) {
      vl[nv] = vL(i,nv);
      vr[nv] = vR(i,nv);
    }
    PrimToCons (vl, ul, run_config);
    PrimToCons (vr, ur, run_config);

    data->err(ti,tj,tk) = 0;   // reset error 

    // 2a. - compute sound speed, cons variables and fluxes
    //     - compute HLL speeds.
      
    STATE_SET(stateL, vl, ul, NULL, x1[ti], x2[tj], x3[tk]);
    STATE_SET(stateR, vr, ur, NULL, x1[ti], x2[tj], x3[tk]);

    SoundSpeed2 (stateL, run_config);
    SoundSpeed2 (stateR, run_config);

    Flux<dir> (stateL, fl, pl);
    Flux<dir> (stateR, fr, pr);

    HLL_Speed<dir> (stateL, stateR, Sl, Sr);

    // 2b. Get max propagation speed for dt comp.

    double scrh = MAX(fabs(Sl), fabs(Sr));
    cmax(i)  = scrh;

    // 2c. Compute Flux if SL > 0 or SR < 0

    if (Sl > 0.0){

      #pragma acc loop seq
      NFLX_LOOP(nv) flux(i,nv) = fl[nv];
      press(i) = pl;

    }else if (Sr < 0.0){
      #pragma acc loop seq
      NFLX_LOOP(nv) flux(i,nv) = fr[nv];
      press(i) = pr;

    }else{

      // 2d. Compute Flux when  SL < 0 < SR

      scrh = 1.0/(Sr - Sl);
      #pragma acc loop seq
      NFLX_LOOP(nv) {
        flux(i,nv)  = Sl*Sr*(ur[nv] - ul[nv]) + Sr*fl[nv] - Sl*fr[nv];
        flux(i,nv) *= scrh;
      }
      press(i) = (Sr*pl - Sl*pr)*scrh;
    } // end if (SL < 0, SR > 0) 

    // 2e. Check for NaNs in flux

    #pragma acc loop seq
    for (nv = 0; nv < NFLX; nv++){
      if (isnan(flux(i,nv))) {
        ATOMIC_ERR_UPDATE(err_kernel, ERR_NAN, ERR_QUIT);
        KERNEL_ERR_UPDATE(data->err(ti,tj,tk),  ERR_NAN, ERR_QUIT); 
      }
    }

    // 2f. Compute tracer fluxes

    #if NSCL > 0
    if (flux(i,RHO) > 0.0){
      #pragma acc loop seq
      NSCL_LOOP(nv) flux(i,nv) = flux(i,RHO)*vl[nv];
    }else{
      #pragma acc loop seq
      NSCL_LOOP(nv) flux(i,nv) = flux(i,RHO)*vr[nv];
    }
    #endif
    
  }}} // End loop (i,j,k)

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
