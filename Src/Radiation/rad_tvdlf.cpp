// /////////////////////////////////////////////////////////////////////
/*!
  \file
  \brief Lax-Friedrichs (Rusanov) Riemann solver for the radiation fields.

  \authors J. D. Melon Fuksman (fuksman@mpia.de)
  \date    May 8, 2026
*/
// /////////////////////////////////////////////////////////////////////

#include "gpluto.hpp"

template<int dir>
void Rad_LF_Solver (Data *, timeStep *, Grid *, RBox *);

void Rad_LF_Solver_T (Data *data, timeStep *Dts, Grid *grid, RBox *box, int dir){
  if      (dir == IDIR) Rad_LF_Solver<IDIR>(data, Dts, grid, box);
  else if (dir == JDIR) Rad_LF_Solver<JDIR>(data, Dts, grid, box);
  else if (dir == KDIR) Rad_LF_Solver<KDIR>(data, Dts, grid, box);
}

template<int dir>
////////////////////////////////////////////////////////////////////////
/*!
  Solve the radiation Riemann problem using a Rusanov flux.

  \tparam        dir   Sweep direction (IDIR, JDIR, or KDIR).
  \param[in,out] data  Pointer to data struct.
  \param[in,out] Dts   Pointer to time step struct.
  \param[in]     grid  Pointer to grid struct.
  \param[in]     box   Pointer to sweep box.
*/
////////////////////////////////////////////////////////////////////////
void Rad_LF_Solver (Data *data, timeStep *Dts, Grid *grid, RBox *box)
 {
  int i, j, k, nv;
  int ibeg, jbeg, kbeg;
  int iend, jend, kend;
  int ni, nj;
  uint32_t err_kernel = 0;

  double *x1 = (dir == IDIR ? grid->xL[IDIR]:grid->xC[IDIR]);
  double *x2 = (dir == JDIR ? grid->xL[JDIR]:grid->xC[JDIR]);
  double *x3 = (dir == KDIR ? grid->xL[KDIR]:grid->xC[KDIR]);
  RunConfig &run_config = *(data->run_config);

  ARRAY_STRIDES (grid, dir, ni, nj);
  RBoxPermuteIndices<dir>(box, ibeg, iend, jbeg, jend, kbeg, kend);
  iend++;

  // --------------------------------------------------------
  //  1. Loop over domain box and solve Riemann problem
  //     in the direction given by dir.
  // --------------------------------------------------------

  #pragma omp target teams distribute parallel for collapse(3)
  #pragma acc parallel loop collapse(3) present(data, run_config, Dts, grid) \
                                        copy(err_kernel)
  for (k = kbeg; k <= kend; k++){
  for (j = jbeg; j <= jend; j++){
  for (i = ibeg; i <= iend; i++){

    int ti, tj, tk;
    double ul[NVAR], ur[NVAR];
    double vl[NVAR], vr[NVAR];
    double fl[NVAR], fr[NVAR];
    State stateL, stateR;
    #if RADIATION_NR
    double cr = run_config.rad_reduced_c;
    #endif

    if constexpr (dir == IDIR) {ti = i; tj = j; tk = k;}
    if constexpr (dir == JDIR) {ti = j; tj = i; tk = k;}
    if constexpr (dir == KDIR) {ti = j; tj = k; tk = i;}

    data->err(ti,tj,tk) = 0; // reset error

    long int offset = ni*(j + nj*k);
    long int offset1 = NVAR*ni*(j + nj*k);

    Ary1D cmax(&Dts->cmax[offset],ni);
    Ary2D vL(&data->sweep->vL[offset1],ni,NVAR);
    Ary2D vR(&data->sweep->vR[offset1],ni,NVAR);
    Ary2D flux(&data->sweep->flux[offset1],ni,NVAR);

      #pragma acc loop seq
      NVAR_LOOP(nv) {
        vl[nv] = vL(i,nv);
        vr[nv] = vR(i,nv);
      }

      // ----------------------------------------------------
      //  2a. Compute sound speed, cons variables and fluxes.
      // ----------------------------------------------------

      // Flux limiting (|F|<E)
      LimitRadFlux(vl) ;
      LimitRadFlux(vr) ;
      #pragma acc loop seq
      NRAD_LOOP(nv) {
        ul[nv] = vl[nv];
        ur[nv] = vr[nv];
      }

      STATE_SET(stateL, vl, ul, NULL, x1[ti], x2[tj], x3[tk]);
      STATE_SET(stateR, vr, ur, NULL, x1[ti], x2[tj], x3[tk]);
      stateL.i = stateR.i = ti;
      stateL.j = stateR.j = tj;
      stateL.k = stateR.k = tk;
      
      // L/R fluxes
      RadFlux<dir> (stateL, fl);
      RadFlux<dir> (stateR, fr);

      // --------------------------------------------------------
      //  2b. Compute max eigenvalue and Lax-Friedrichs fluxes
      // --------------------------------------------------------

      {
        double Sl;
        double Sr;
        RadSpeed<dir> (stateL, stateR, Sl, Sr, grid, run_config);

        // Set maximum speeds considering both radiation and fluid
        double scrh = MAX(fabs(Sl), fabs(Sr));
        #if RADIATION_NR
        cmax(i) = cr*scrh;
        #else
        cmax(i) = MAX(cmax(i),scrh);
        #endif

        #pragma acc loop seq
        NRAD_LOOP(nv) {
          flux(i,nv) =  0.5*(fl[nv] + fr[nv]) - 0.5*scrh*(ur[nv] - ul[nv]);
        }
      }

      #if RADIATION_NR
      #pragma acc loop seq
      NRAD_LOOP(nv) {
        flux(i,nv) *= cr;
      }
      #endif

      #pragma acc loop seq
      NRAD_LOOP(nv) {
        if (isnan(flux(i,nv))) {
          ATOMIC_ERR_UPDATE(err_kernel, ERR_NAN, ERR_QUIT);
          KERNEL_ERR_UPDATE(data->err(ti,tj,tk), ERR_NAN, ERR_QUIT);
        }
      }
  }}}

  #if WARNING_MESSAGES > 0
  Where ("Rad_LF_Solver()", data, err_kernel, box, grid);
  #endif

  #if WARNING_MESSAGES >= 0
  if (err_kernel & ERR_QUIT) QuitPLUTO("from Rad_LF_Solver()", err_kernel);
  #endif

}