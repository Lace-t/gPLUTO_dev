// /////////////////////////////////////////////////////////////////////
/*!
  \file
  \brief HLL Riemann solver for the (M1) radiation transport equations.

  Solve the Riemann problem for the radiation transport equations using
  a single-sweep HLLC solver.

  On input, this function takes left and right primitive sweep vectors
  \c sweep->vL and \c sweep->vR at zone edge \c i+1/2;
  On output, return fluxes of the radiation fields at the same interface
  \c i+1/2 (note that the \c i refers to \c i+1/2).

  Also during this step, compute maximum wave propagation speed (cmax)
  for  explicit time step computation, considering both radiation and
  MHD fields.

  \authors J. D. Melon Fuksman (fuksman@mpia.de)

  \date    May 9, 2026
  
  \b References
 -  Melon Fuksman, J. D., and Mignone, A. 2019, "A Radiative
    Transfer Module for Relativistic Magnetohydrodynamics in the
    PLUTO Code" ApJS, 242, 20.
*/
// /////////////////////////////////////////////////////////////////////
#include "gpluto.hpp"

template<int dir>
void Rad_HLLC_Solver (Data *, timeStep *, Grid *, RBox *);

void Rad_HLLC_Solver_T (Data *data, timeStep *Dts, Grid *grid, RBox *box, int idir){
  if      (idir == IDIR) Rad_HLLC_Solver<IDIR>(data, Dts, grid, box);
  else if (idir == JDIR) Rad_HLLC_Solver<JDIR>(data, Dts, grid, box);
  else if (idir == KDIR) Rad_HLLC_Solver<KDIR>(data, Dts, grid, box);
}

template<int dir>
////////////////////////////////////////////////////////////////////////
/*!
  Solve the RHD Riemann problem using the HLLC Riemann solver.

  \tparam        dir   Sweep direction (IDIR, JDIR, or KDIR).
  \param[in,out] data  Pointer to data struct.
  \param[in,out] Dts   Pointer to time step struct.
  \param[in]     grid  Pointer to grid struct.
  \param[in]     box   Pointer to sweep box.
*/
////////////////////////////////////////////////////////////////////////
void Rad_HLLC_Solver (Data *data, timeStep *Dts, Grid *grid, RBox *box)
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
  GEN_INDEX_CYCLE (dir, FR1, FRn, FRt, FRb);
  RBoxPermuteIndices<dir>(box, ibeg, iend, jbeg, jend, kbeg, kend);
  iend++;

  // --------------------------------------------------------
  //  1. Loop over domain box and solve Riemann problem
  //    in the direction given by dir.
  // --------------------------------------------------------

  #pragma omp target teams distribute parallel for collapse(3)
  #pragma acc parallel loop collapse(3) present(data, run_config, Dts, grid) \
                                        copy(err_kernel)
  for (k = kbeg; k <= kend; k++){
  for (j = jbeg; j <= jend; j++){
  for (i = ibeg; i <= iend; i++){

    int ti, tj, tk;
    double usl[NVAR], usr[NVAR];
    double ul[NVAR], ur[NVAR];
    double vl[NVAR], vr[NVAR];
    double fl[NVAR], fr[NVAR];
    double Sl, Sr, ps;
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
    Ary1Duint32_t flag(&data->sweep->flag[offset],ni);
    Ary2D vL(&data->sweep->vL[offset1],ni,NVAR);
    Ary2D vR(&data->sweep->vR[offset1],ni,NVAR);
    Ary2D flux(&data->sweep->flux[offset1],ni,NVAR);
    Ary2D fL(&data->sweep->fL[offset1],ni,NVAR);
    Ary2D fR(&data->sweep->fR[offset1],ni,NVAR);

      #pragma acc loop seq
      NVAR_LOOP(nv) {
        vl[nv] = vL(i,nv);
        vr[nv] = vR(i,nv);
      }

      // ----------------------------------------------------
      //  2a. - Compute sound speed, cons variables and fluxes,
      //      - compute HLL speeds.
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

      // Speed calculation, flag for HLL if diffusion limiter is applied
      int diff_cond = RadSpeed<dir> (stateL, stateR, Sl, Sr, grid, run_config);
      if (diff_cond) flag(i) |= FLAG_HLL;

      // L/R fluxes
      RadFlux<dir> (stateL, fl);
      RadFlux<dir> (stateR, fr);

      // Set maximum speeds considering both radiation and fluid
      double scrh = MAX(fabs(Sl), fabs(Sr));
      #if RADIATION_NR
      cmax(i) = cr*scrh;
      #else
      cmax(i) = MAX(cmax(i),scrh);
      #endif

      // --------------------------------------------------------
      //  2b. Compute Flux if SrL > 0 or SrR < 0
      // --------------------------------------------------------

      // Switch to tvdlf flux if speeds or enr are small
      if ( (fabs(Sl) < 1e-300 && fabs(Sr) < 1e-300) ||
           (ul[ENR] < 1e-300 && ur[ENR] < 1e-300) ){

        #pragma acc loop seq
        NRAD_LOOP(nv) flux(i,nv) = 0.5*(fl[nv] + fr[nv] - scrh*(ur[nv] - ul[nv]));

      } else if (Sl > 0.0){

        #pragma acc loop seq
        NRAD_LOOP(nv) flux(i,nv) = fl[nv];

      }else if (Sr < 0.0){
        #pragma acc loop seq
        NRAD_LOOP(nv) flux(i,nv) = fr[nv];

      }else{

        // ----------------------------------------------------
        //  3. Compute Flux when  SrL < 0 < SrR
        // ----------------------------------------------------

        bool switch_to_hll = false;
        double us;

        #if RADIATION_DIFF_LIMITING || SHOCK_FLATTENING == MULTID
        if ((flag(i) & FLAG_HLL) || (flag(i-1) & FLAG_HLL)){
          switch_to_hll = true;
        }
        #endif

        if (!switch_to_hll) {

          // ---------------------------------------
          //  3a. Get *U
          // ---------------------------------------
        
          double vxl, vxr;
          {
            double fr2L = ul[FR1]*ul[FR1] + ul[FR2]*ul[FR2] + ul[FR3]*ul[FR3] ;
            double fr2R = ur[FR1]*ur[FR1] + ur[FR2]*ur[FR2] + ur[FR3]*ur[FR3] ;

            double cosL = (fr2L > 1e-200) ? ul[FRn]/sqrt(fr2L) : ul[FRn]/1e-100 ;
            double cosR = (fr2R > 1e-200) ? ur[FRn]/sqrt(fr2R) : ur[FRn]/1e-100 ;
          
            double chiL = ul[ENR]*ul[ENR];
            double chiR = ur[ENR]*ur[ENR];
          
            fr2L = fr2L/chiL;
            fr2R = fr2R/chiR;

            chiL = (3.0+4.0*fr2L)/(5.0+2.0*sqrt(4.0-3.0*fr2L));
            chiR = (3.0+4.0*fr2R)/(5.0+2.0*sqrt(4.0-3.0*fr2R));
          
            vxl = (fr2L > 1e-100) ? (1.5*chiL-0.5)*cosL/sqrt(fr2L) : 0.0 ;
            vxr = (fr2R > 1e-100) ? (1.5*chiR-0.5)*cosR/sqrt(fr2R) : 0.0 ;
          }

          double AL = Sl*ul[ENR] - fl[ENR];
          if (AL > 0.0) AL = 0.0;
          double AR = Sr*ur[ENR] - fr[ENR];
          if (AR < 0.0) AR = 0.0;
        
          double BL = Sl*ul[FRn] - fl[FRn];
          double BR = Sr*ur[FRn] - fr[FRn];
        
          double fpL = ul[FRt]*ul[FRt] + ul[FRb]*ul[FRb] ;
          double fpR = ur[FRt]*ur[FRt] + ur[FRb]*ur[FRb] ;
        
          double eeL = 1e-10*ul[ENR] ;
          double eeR = 1e-10*ur[ENR] ;

          // Return HLL flux in special cases
          if( (fabs(AL)<eeL && fabs(AR)<eeR) || (fpL<eeL*eeL && fpR<eeR*eeR) ){
            switch_to_hll = true;
          }

          if (!switch_to_hll) {
            double a = AR*Sl - AL*Sr;
            double b = AL + BL*Sr - AR - BR*Sl;
            double c = BR - BL;

            scrh = b*b - 4.0*a*c;
            if (scrh < 0.0) scrh = 0.0;
            scrh = -0.5*(b + DSIGN(b)*sqrt(scrh));
            us = c/scrh;

            if (fabs(Sl - us) < 1e-100 || fabs(Sr - us) < 1e-100) {
              switch_to_hll = true;
            }

            if (!switch_to_hll) {
              ps = (AL*us - BL)/(1.0 - us*Sl);
              usl[FRn] = (Sl*(ul[ENR] + ps) - ul[FRn])*us/(Sl - us);
              usr[FRn] = (Sr*(ur[ENR] + ps) - ur[FRn])*us/(Sr - us);

              if (isnan(usl[FRn]) || isnan(usr[FRn])) {
                switch_to_hll = true;
              }
            }
          }

          if (!switch_to_hll) {

            usl[FRt] = ul[FRt]*(Sl - vxl)/(Sl - us);
            usr[FRt] = ur[FRt]*(Sr - vxr)/(Sr - us);
            usl[FRb] = ul[FRb]*(Sl - vxl)/(Sl - us);
            usr[FRb] = ur[FRb]*(Sr - vxr)/(Sr - us);

            usl[ENR] = ul[ENR] + (usl[FRn]-ul[FRn])/Sl ;
            usr[ENR] = ur[ENR] + (usr[FRn]-ur[FRn])/Sr ;

            if (isnan(usl[ENR]) || isnan(usr[ENR]) ||
                isnan(usl[FRt]) || isnan(usr[FRt]) ||
                isnan(usl[FRb]) || isnan(usr[FRb])) {
              switch_to_hll = true;
            }

          }
        }

        if (switch_to_hll) {

          scrh = 1.0/(Sr - Sl);
          #pragma acc loop seq
          NRAD_LOOP(nv) {
            flux(i,nv)  = Sl*Sr*(ur[nv] - ul[nv]) + Sr*fl[nv] - Sl*fr[nv];
            flux(i,nv) *= scrh;
          }

        } else {

          // ---------------------------------------
          //  3b. Compute HLLC flux
          // ---------------------------------------

          if (us >= 0.0){
            #pragma acc loop seq
            NRAD_LOOP(nv) {
              flux(i,nv) = fl[nv] + Sl*(usl[nv] - ul[nv]);
            }
          } else {
            #pragma acc loop seq
            NRAD_LOOP(nv) {
              flux(i,nv) = fr[nv] + Sr*(usr[nv] - ur[nv]);
            }
          }

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
  Where ("Rad_HLLC_Solver()", data, err_kernel, box, grid);
  #endif

  #if WARNING_MESSAGES >= 0
  if (err_kernel & ERR_QUIT) QuitPLUTO("from Rad_HLLC_Solver()", err_kernel);
  #endif

}