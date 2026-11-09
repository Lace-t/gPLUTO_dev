////////////////////////////////////////////////////////////////////////
/*!
  \file
  \brief Lax-Friedrechs (Rusanov) Riemann solver for HD.

  Solve the Riemann problem for the adiabatic and isothermal HD
  equations using the Lax-Friedrichs Rusanov Riemann solver with local
  maximum characteristic speed:
  \f[
      \hat{F}_{i+\HALF} = \frac{F^L_{i+\HALF} + F^R_{i+\HALF}}{2}
   - c^{\rm max}_{i+\HALF}\frac{U^R_{i+\HALF} - U^L_{i-\HALF}}{2}
     \qquad{\rm where}\quad c^{\rm max} = \lambda_{\rm max}\left(
          \frac{V^R_{i+\HALF} + V^L_{i+\HALF}}{2}\right)
  \f]
  where \f$\lambda_{\rm \max}\f$ is a function of the arithmetic
  average between the left and the right states.

  On input, this function takes left and right primitive state vectors
  \c data->sweep->vL and \c data->sweep->vR at zone edge i-1/2.
  On output, return flux and pressure vectors at the same interface
  \c i-1/2 (note that the \c i refers to \c i-1/2).

  Also during this step, compute maximum wave propagation speed (cmax)
  for  explicit time step computation.

  \b Reference:
    - "Comparison of Some Flux Correced Transport and Total Variation
       Diminishing Numerical Schemes for Hydrodynamics and Magnetohydrodynamics
       Problems", Toth and Odstrcil, JCP (1996), 128,82

  \authors A. Mignone  (andrea.mignone@unito.it)\n
           M. Rossazza (marco.rossazza@unito.it)\n
           V. Berta    (vittoria.berta@unito.it)\n
           A. Suriano  (alessio.suriano@unito.it)
  \date    Feb 12, 2026
*/
////////////////////////////////////////////////////////////////////////
#include "gpluto.hpp"

template<int dir>
void LF_Solver (Data *, timeStep *, Grid *, RBox *);

////////////////////////////////////////////////////////////////////////
void LF_Solver_T (Data *data, timeStep *Dts, Grid *grid, RBox *box, int dir)
/*!
*/
////////////////////////////////////////////////////////////////////////
{
  if      (dir == IDIR) LF_Solver<IDIR>(data, Dts, grid, box);
  else if (dir == JDIR) LF_Solver<JDIR>(data, Dts, grid, box);
  else if (dir == KDIR) LF_Solver<KDIR>(data, Dts, grid, box);
}

////////////////////////////////////////////////////////////////////////
/*!
  Solve the Riemann problem between L/R states using a
  Rusanov-Lax Friedrichs flux.
*/
////////////////////////////////////////////////////////////////////////
template<int dir>
void LF_Solver (Data *data, timeStep *Dts, Grid *grid, RBox *box)
{
  int  i, j, k, nv;
  constexpr int  VXn = VX1 + dir;
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

    long int offset = ni*(j + nj*k);
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
    double vLR[NVAR];
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

    // 2a. Compute conservative variables

    PrimToCons (vl, ul, run_config);
    PrimToCons (vr, ur, run_config);

    STATE_SET(stateL, vl, ul, NULL, x1[ti], x2[tj], x3[tk]);
    STATE_SET(stateR, vr, ur, NULL, x1[ti], x2[tj], x3[tk]);

    // 2b. Compute average state in order to get the local max 
    //     characteristic velocity.
    //     In order to avoid underestimating this speed at strong 
    //     reflecting  boundaries (or when vxL = -vxR in general), 
    //     we average the absolute value of the normal velocity.
      
    #pragma acc loop seq
    for (nv = 0; nv < NFLX; nv++) vLR[nv] = 0.5*(vl[nv] + vr[nv]);

    vLR[VXn] = 0.5*(fabs(vl[VXn]) + fabs(vr[VXn]));

    // 2c. Compute sound speed, max and min signal speeds, pressure 
    //     and fluxes

    {
    State stateLR;
    STATE_SET(stateLR, vLR, NULL, NULL, x1[ti], x2[tj], x3[tk]);
    SoundSpeed2 (stateLR, run_config);
    double a = sqrt(stateLR.a2);
    double cmaxLR = vLR[VXn] + a;
    double cminLR = vLR[VXn] - a;
    double cLR = MAX(fabs(cminLR), fabs(cmaxLR));
    cmax(i) = cLR;
    cs2(i)  = stateLR.a2;
    }

    Flux<dir> (stateL, fL, pl);
    Flux<dir> (stateR, fR, pr);

    press(i) = 0.5*(pl + pr);

    // 2d. Compute fluxes

    #pragma acc loop seq
    for (nv = 0; nv < NFLX; nv++) {
      flux(i,nv) = 0.5*(fL[nv] + fR[nv]) - 0.5*cmax(i)*(ur[nv] - ul[nv]);
    }

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

  }}} // End Domain loop (i,j,k)

  // ------------------------------------------------------
  // 3. Error Check / Quit condition
  // ------------------------------------------------------

  #if WARNING_MESSAGES > 0 
  Where ("LF_Solver()", data, err_kernel, box, grid);
  #endif

  #if WARNING_MESSAGES >= 0 
  if (err_kernel & ERR_QUIT) QuitPLUTO("from LF_Solver()",err_kernel);
  #endif
}
