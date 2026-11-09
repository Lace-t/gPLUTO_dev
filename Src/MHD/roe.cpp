////////////////////////////////////////////////////////////////////////
/*!
  \file
  \brief Implementation of the Roe Riemann solver for the MHD equations.

  Solve the Riemann problem for the adiabatic and isothermal MHD equations
  using the linearized Riemann solver of Roe.
  The implementation follows the approach outlined by Cargo & Gallice (1997).

  The isothermal version is recovered by taking the limit of
  \f$ \bar{a}^2 \f$ for \f$ \gamma \to 1\f$ which gives
  (page 451) \f$\bar{a}^2 \to {\rm g\_isoSoundSpeed2} + X\f$
  where \c X is defined as in the adiabatic case.
  This follows by imposing zero jump across the entropy wave
  (first of Eq. 4.20, page 452) giving
  \f$ \Delta p = (\bar{a}^2 - X)\Delta\rho\f$.
  Then all the terms like
  \f$ [X\Delta\rho + \Delta p] \to \bar{a}^2\Delta\rho \f$.

  When the macro HLL_HYBRIDIZATION is set to YES, we revert to the HLL
  solver whenever an unphysical state appear in the solution.

  The macro CHECK_ROE_MATRIX can be used to verify that the
  characteristic decomposition reproduces the Roe matrix.
  This can be done also when BACKGROUND_FIELD is set to YES.

  On input, this function takes left and right primitive state vectors
  \c data->sweep->vL and \c data->sweep->vR at zone edge i-1/2;
  On output, return flux and pressure vectors at the same interface
  \c i-1/2 (note that the \c i refers to \c i-1/2).

  Also during this step, compute maximum wave propagation speed (cmax)
  for  explicit time step computation.

  \b References:
    - "Roe Matrices for Ideal MHD and Systematic Construction
       of Roe Matrices for Systems of Conservation Laws",
       P. Cargo, G. Gallice, JCP (1997), 136, 446.

  \authors A. Mignone  (andrea.mignone@unito.it)  \n
           M. Rossazza (marco.rossazza@unito.it)  \n
           C. Zanni   (zanni@oato.inaf.it)        \n
           S. Truzzi  (stefano.truzzi@unito.it)   \n
           A. Suriano  (alessio.suriano@unito.it)
  \date    Feb 12, 2026
*/
////////////////////////////////////////////////////////////////////////
#include "gpluto.hpp"

#define sqrt_1_2  (0.70710678118654752440)

#ifndef HLL_HYBRIDIZATION
  #define HLL_HYBRIDIZATION    NO
#endif
#ifndef CHECK_ROE_MATRIX
  #define CHECK_ROE_MATRIX     NO
#endif

#ifndef ROE_FIX_DEGENERACY
  #define ROE_FIX_DEGENERACY   NO
#endif

// Function prototype
template<int dir>
void Roe_Solver (Data *, timeStep *, Grid *, RBox *);

////////////////////////////////////////////////////////////////////////
/*!
*/
////////////////////////////////////////////////////////////////////////
void Roe_Solver_T (Data *data, timeStep *Dts, Grid *grid, RBox *box, int dir)
{
  if      (dir == IDIR) Roe_Solver<IDIR>(data, Dts, grid, box);
  else if (dir == JDIR) Roe_Solver<JDIR>(data, Dts, grid, box);
  else if (dir == KDIR) Roe_Solver<KDIR>(data, Dts, grid, box);
}

////////////////////////////////////////////////////////////////////////
/*!
*/
////////////////////////////////////////////////////////////////////////
template<int dir>
void Roe_Solver (Data *data, timeStep *Dts, Grid *grid, RBox *box)
{
  int      i, j, k, g, nv;
  int      ibeg, jbeg, kbeg;
  int      iend, jend, kend;
  int      ni, nj;
  uint32_t err_kernel = 0;  // Cumulative kernel error flag 

  RunConfig &run_config = *(data->run_config);

  #ifdef GLM_MHD
  double glm_ch = run_config.glm_ch;
  double ch2    = glm_ch*glm_ch;
  #endif

  #if HLL_HYBRIDIZATION == YES
  int ifail=0;
  double dvel;
  double prr;
  #endif
  double *x1 = (dir == IDIR ? grid->xL[IDIR]:grid->xC[IDIR]); 
  double *x2 = (dir == JDIR ? grid->xL[JDIR]:grid->xC[JDIR]); 
  double *x3 = (dir == KDIR ? grid->xL[KDIR]:grid->xC[KDIR]); 

  #if EOS == IDEAL
  double gmm = run_config.gamma;
  #endif

  // ------------------------------------------------------
  // 0. Compatibility checks
  // ------------------------------------------------------

  #if (defined STAGGERED_MHD) && (CT_EMF_AVERAGE == UCT_GFORCE)
  QuitPLUTO("! Roe_Solver(): not compatible with UCT_GFORCE", 1);
  #endif

  // ------------------------------------------------------
  // 1. Compute strides and box permutation
  // ------------------------------------------------------

  ARRAY_STRIDES (grid, dir, ni, nj);

  GEN_INDEX_CYCLE (dir, VX1, VXn, VXt, VXb);
  GEN_INDEX_CYCLE (dir, MX1, MXn, MXt, MXb);
  GEN_INDEX_CYCLE (dir, BX1, BXn, BXt, BXb);
  GEN_INDEX_CYCLE (dir, IDIR, DIRn, DIRt, DIRb);

  // Cycle on directional indexes to preserve sweep ordering
  RBoxPermuteIndices<dir>(box, ibeg, iend, jbeg, jend, kbeg, kend);
  iend++;  // Riemann solver need one more interface value

  // ------------------------------------------------------
  // 2. Loop over domain box and solve Riemann problem
  //    in the direction given by dir.
  // ------------------------------------------------------

  #pragma omp target teams distribute parallel for collapse(3) 
  #pragma acc parallel loop collapse(3) present(data, run_config, grid, Dts) \
                                        copy(err_kernel)
  for (k = kbeg; k <= kend; k++){
  for (j = jbeg; j <= jend; j++){
  for (i = ibeg; i <= iend; i++){
  
    // Compute offsets and assign local pointers

    long int offset = ni*(j + nj*k);
    long int offset1 = NVAR*ni*(j + nj*k);

    Ary1D Bn(&data->sweep->Bn[offset],ni);
    Ary1D PSIn(&data->sweep->PSIn[offset],ni);
    Ary1D cmax(&Dts->cmax[offset],ni);
    Ary1D s_cs2(&data->sweep->cs2[offset],ni);
    Ary1D press(&data->sweep->press[offset],ni);
    Ary1Duint32_t flag(&data->sweep->flag[offset],ni);

    Ary2D vL(&data->sweep->vL[offset1],ni,NVAR);
    Ary2D vR(&data->sweep->vR[offset1],ni,NVAR);
    Ary2D flux(&data->sweep->flux[offset1],ni,NVAR);
    Ary2D fL(&data->sweep->fL[offset1],ni,NVAR);
    Ary2D fR(&data->sweep->fR[offset1],ni,NVAR);

    // Private variable declaration

    double Rc[NVAR][NVAR];
    double eta[NVAR], lambda[NVAR], alambda[NVAR];

    int switch_to_hll = 0;
    int ti, tj, tk;

    double vl[NVAR], ul[NVAR], fl[NVAR];
    double vr[NVAR], ur[NVAR], fr[NVAR];
    double Sl, pl, a2L;
    double Sr, pr, a2R;
    double bgf[3];
    double a2;
    double alpha_f,alpha_s, beta_y, beta_z;

    State stateL, stateR;

    // Set state and "true" indices

    if constexpr (dir == IDIR) {ti = i; tj = j; tk = k;}
    if constexpr (dir == JDIR) {ti = j; tj = i; tk = k;}
    if constexpr (dir == KDIR) {ti = j; tj = k; tk = i;}

    NFLUID_LOOP(nv) {
      vl[nv] = vL(i,nv);
      vr[nv] = vR(i,nv);
    }

    STATE_SET(stateL, vl, ul, bgf, x1[ti], x2[tj], x3[tk]);
    STATE_SET(stateR, vr, ur, bgf, x1[ti], x2[tj], x3[tk]);

    data->err(ti,tj,tk) = 0;   // reset error 

    // 2a. GLM pre-conditioner

    #ifdef GLM_MHD
    GLM_Solve (vl, vr, run_config, dir);
    #endif
    PrimToCons (vl, ul, run_config);
    PrimToCons (vr, ur, run_config);

    // 2b. Compute background field, sound speeds & fluxes

    #if BACKGROUND_FIELD == YES
    BackgroundField (stateL, run_config); // stateL->bgf = stateR->bgf
    #endif

    SoundSpeed2 (stateL, run_config);
    SoundSpeed2 (stateR, run_config);
    a2L = stateL.a2;
    a2R = stateR.a2;

    Flux<dir> (stateL, fl, pl, run_config);
    Flux<dir> (stateR, fr, pr, run_config);

    // 2c. Keep Roe / Switch to HLL ?

    #if HLL_SWITCH == YES
    if ((flag(i) & FLAG_HLL) || (flag(i-1) & FLAG_HLL)){
      switch_to_hll = 1;
    }
    #endif

    double dV[NFLX], dU[NFLX]; // has to be outside of if(!switch_to_hll)
                                  // otherwise GPU run does not work
                                  // (maybe compiler bug)                    

    if (!switch_to_hll){ // Use ROE 

      #pragma acc loop seq
      for (nv = 0; nv < NFLX; nv++) {
        dV[nv] = vr[nv] - vl[nv];
        dU[nv] = ur[nv] - ul[nv];
        #pragma acc loop seq
        for (g = 0; g < NFLX; g++) Rc[nv][g] = 0.0;
      }

      // 2d. Compute Roe averages

      double X  = dV[BXn]*dV[BXn] + dV[BXt]*dV[BXt] + dV[BXb]*dV[BXb];
      double sl, sr;
      {
        double sqr_rho_L = sqrt(vl[RHO]);
        double sqr_rho_R = sqrt(vr[RHO]);

        X /= (sqr_rho_L + sqr_rho_R)*(sqr_rho_L + sqr_rho_R)*2.0;

        sl = sqr_rho_L/(sqr_rho_L + sqr_rho_R);
        sr = sqr_rho_R/(sqr_rho_L + sqr_rho_R);
      }

      double rho = sr*vl[RHO] + sl*vr[RHO];
      double sqrt_rho = sqrt(rho);

      double Bx = sr*vl[BXn] + sl*vr[BXn];
      double By = sr*vl[BXt] + sl*vr[BXt];
      double Bz = sr*vl[BXb] + sl*vr[BXb];

      #if BACKGROUND_FIELD == YES
      // -- Define field B0 and total B. B1 is the deviation -- 
      double B0x = bgf[BXn-BX1];
      double B1x = sr*vl[BXn] + sl*vr[BXn];
      Bx = B0x + B1x;
      double B0y = bgf[BXt-BX1];
      double B1y = sr*vl[BXt] + sl*vr[BXt];
      By = B0y + B1y;
      double B0z = bgf[BXb-BX1];
      double B1z = sr*vl[BXb] + sl*vr[BXb];
      Bz = B0z + B1z;
      #endif

      double sBx = (Bx >= 0.0 ? 1.0 : -1.0);

      double bx = Bx/sqrt_rho;
      double by = By/sqrt_rho;
      double bz = Bz/sqrt_rho;

      double bt2   = by*by + bz*bz;
      double b2    = bx*bx + bt2;
      double Btmag = sqrt(bt2*rho);

      double u = sl*vl[VXn] + sr*vr[VXn];
      double v = sl*vl[VXt] + sr*vr[VXt];
      double w = sl*vl[VXb] + sr*vr[VXb];

      // 2e. Compute enthalpy and sound speed

      #if EOS == ISOTHERMAL
      a2 = 0.5*(a2L + a2R) + X;  // in most cases a2L = a2R
                                    // for isothermal MHD 
      #elif EOS == IDEAL
      double Hgas; // gas enthalpy 
      double vel2    = u*u + v*v + w*w;
      double g1 = gmm - 1.0;
      {
        double vdm = u*dU[MXn] + v*dU[MXt] + w*dU[MXb];
        #if BACKGROUND_FIELD == YES // BdB = B1.dB1 (deviation only) 
        double BdB = B1x*dU[BXn] + B1y*dU[BXt] + B1z*dU[BXb];
        #else
        double BdB = Bx*dU[BXn] + By*dU[BXt] + Bz*dU[BXb];
        #endif
        dV[PRS] = g1*((0.5*vel2 - X)*dV[RHO] - vdm + dU[ENG] - BdB);
      }
      {
        double HL = (ul[ENG] + pl)/vl[RHO];
        double HR = (ur[ENG] + pr)/vr[RHO];
        double H  = sl*HL + sr*HR;   // total enthalpy 
        
        #if BACKGROUND_FIELD == YES
        double scrh = B1x*Bx + B1y*By + B1z*Bz;;
        Hgas = H - scrh/rho; 
        #else
        Hgas = H - b2;       
        #endif
      }

      a2 = (2.0 - gmm)*X + g1*(Hgas - 0.5*vel2);
      if (a2 < 0.0) {
        ATOMIC_ERR_UPDATE(err_kernel, ERR_RIEMANN_SOLVE_FAIL);
        KERNEL_ERR_UPDATE(data->err(ti,tj,tk), ERR_RIEMANN_SOLVE_FAIL);
      }

      #endif // EOS == IDEAL 

      // 2f. Compute fast and slow magnetosonic speeds.
      //     The following expression appearing in the definitions
      //     of the fast magnetosonic speed
      //
      //     (a^2 - b^2)^2 + 4*a^2*bt^2 = (a^2 + b^2)^2 - 4*a^2*bx^2
      //
      //     is always positive and avoids round-off errors.
      //
      //     Note that we always use the total field to compute the
      //     characteristic speeds.

      double a,cs,ca,cf,scrh, sqDelta;
      double dab2 = a2 - b2;
      {
        double ca2  = bx*bx;
        scrh = dab2*dab2 + 4.0*bt2*a2;
        sqDelta = sqrt(scrh);

        double cf2 = 0.5*(a2 + b2 + sqDelta);
        double cs2 = a2*ca2/cf2;   // -- same as 0.5*(a2 + b2 - scrh) -- 

        cf = sqrt(cf2);
        cs = sqrt(cs2);
        ca = sqrt(ca2);
        a  = sqrt(a2);
      } // ca2,cs2,cf2 

      if (cf == cs) {
        alpha_f = 1.0;
        alpha_s = 0.0;
      }else{
        alpha_f = 0.5*( dab2/sqDelta + 1.0);
        alpha_s = 0.5*(-dab2/sqDelta + 1.0);

        alpha_f = MAX(0.0, alpha_f);
        alpha_s = MAX(0.0, alpha_s);
        alpha_f = sqrt(alpha_f);
        alpha_s = sqrt(alpha_s);
      }

      if (Btmag > 0.0) {
        beta_y = By/Btmag;
        beta_z = Bz/Btmag;
      } else {
        beta_z = sqrt_1_2;
        beta_y = sqrt_1_2;
      }

      // 2f1. Fix degeneracy in the Btmag = 0 limit, when less than 
      //      3 dimensions are included.
      //      This is useful to remove "noise" in the z-component
      //      when the problem is strictly planar.

      #if ROE_FIX_DEGENERACY == YES
        #if INCLUDE_KDIR == NO
        beta_y = (dir == IDIR); beta_z = 1.0-beta_y;
        #elif INCLUDE_JDIR == NO
        beta_y = (dir == IDIR); beta_z = 1.0-beta_y;
        #endif
      #endif

      // 2g.  Compute non-zero entries of conservative eigenvectors (Rc),
      //      wave strength L*dU (=eta) for all 8 (or 7) waves using the
      //      expressions given by Eq. [4.18]--[4.21].
      //      Fast and slow eigenvectors are multiplied by a^2 while
      //      jumps are divided by a^2.
      //
      //      Notes:
      //      - the expression on the paper has a typo in the very last term
      //        of the energy component: it should be + and not - !
      //      - with background field splitting: additional terms must be
      //        added to the energy component for fast, slow and Alfven waves.
      //        To obtain energy element, conservative eigenvector (with
      //        total field) must be multiplied by | 0 0 0 0 -B0y -B0z 1 |.
      //        Also, H - b2 does not give gas enthalpy. A term b0*btot must
      //        be added and eta (wave strength) should contain total field
      //        and deviation's delta.

      // Fast wave:  u - c_f 

      g = KFASTM;
      lambda[KFASTM] = u - cf;

      scrh  = alpha_s*cs*sBx;

      double beta_dv = beta_y*dV[VXt] + beta_z*dV[VXb];
      double beta_dB = beta_y*dV[BXt] + beta_z*dV[BXb];

      Rc[RHO][g] = alpha_f;

      Rc[MXn][g] = alpha_f*lambda[KFASTM];

      Rc[MXt][g] = alpha_f*v + scrh*beta_y;
      Rc[MXb][g] = alpha_f*w + scrh*beta_z;
      Rc[BXn][g] = 0.0;
      Rc[BXt][g] = alpha_s*a*beta_y/sqrt_rho;
      Rc[BXb][g] = alpha_s*a*beta_z/sqrt_rho;
      #if EOS == IDEAL
      double beta_v  = beta_y*v       + beta_z*w;
      Rc[ENG][g] =   alpha_f*(Hgas - u*cf) + scrh*beta_v
                  + alpha_s*a*Btmag/sqrt_rho;
      #if BACKGROUND_FIELD == YES
      Rc[ENG][g] -= B0y*Rc[BXt][g] + B0z*Rc[BXb][g];
      #endif
      eta[g] =   alpha_f*(X*dV[RHO] + dV[PRS]) + rho*scrh*beta_dv
              - rho*alpha_f*cf*dV[VXn]        + sqrt_rho*alpha_s*a*beta_dB;
      #elif EOS == ISOTHERMAL
      eta[g] =   alpha_f*(0.0*X + a2)*dV[RHO] + rho*scrh*beta_dv
              - rho*alpha_f*cf*dV[VXn] + sqrt_rho*alpha_s*a*beta_dB;
      #endif
      eta[g] *= 0.5/a2;

      // Fast wave:  u + c_f

      g = KFASTP;
      lambda[KFASTP] = u + cf;

      Rc[RHO][g] = alpha_f;
      Rc[MXn][g] = alpha_f*lambda[KFASTP];
      Rc[MXt][g] = alpha_f*v - scrh*beta_y;
      Rc[MXb][g] = alpha_f*w - scrh*beta_z;
      Rc[BXn][g] = 0.0;
      Rc[BXt][g] = Rc[BXt][KFASTM];
      Rc[BXb][g] = Rc[BXb][KFASTM];
      #if EOS == IDEAL
      Rc[ENG][g] =   alpha_f*(Hgas + u*cf) - scrh*beta_v
                  + alpha_s*a*Btmag/sqrt_rho;
      #if BACKGROUND_FIELD == YES
      Rc[ENG][g] -= B0y*Rc[BXt][g] + B0z*Rc[BXb][g];
      #endif
      eta[g] =   alpha_f*(X*dV[RHO] + dV[PRS]) - rho*scrh*beta_dv
              + rho*alpha_f*cf*dV[VXn]        + sqrt_rho*alpha_s*a*beta_dB;
      #elif EOS == ISOTHERMAL
      eta[g] =   alpha_f*(0.0*X + a2)*dV[RHO] - rho*scrh*beta_dv
              + rho*alpha_f*cf*dV[VXn]      + sqrt_rho*alpha_s*a*beta_dB;
      #endif
      eta[g] *= 0.5/a2;

      // Entropy wave:  u

      #if EOS == IDEAL
      g = KENTRP;
      lambda[KENTRP] = u;

      Rc[RHO][g] = 1.0;
      Rc[MXn][g] = u;
      Rc[MXt][g] = v;
      Rc[MXb][g] = w;
      Rc[BXn][g] = 0.0;
      Rc[BXt][g] = 0.0;
      Rc[BXb][g] = 0.0;
      Rc[ENG][g] = 0.5*vel2 + (gmm - 2.0)/g1*X;

      eta[g] = ((a2 - X)*dV[RHO] - dV[PRS])/a2;
      #endif

      // div.B wave (u): this wave exists when:
      //
      //   1) 8 wave formulation
      //   2) CT, since we always have 8 components, but it
      //      carries zero jump.
      //
      //   With GLM, KDIVB is replaced by KPSI_GLMM, KPSI_GLMP and these
      //   two waves should not enter in the Riemann solver (eta = 0.0)
      //   since the 2x2 linear system formed by (B,psi) has already
      //   been solved.
      //

      #ifdef GLM_MHD
      lambda[KPSI_GLMP] =  glm_ch;
      lambda[KPSI_GLMM] = -glm_ch;
      eta[KPSI_GLMP] = eta[KPSI_GLMM] = 0.0;
      #else
      g = KDIVB;
      lambda[KDIVB] = u;
      Rc[RHO][g] = 0.0;
      Rc[MXn][g] = 0.0;
      Rc[MXt][g] = 0.0;
      Rc[MXb][g] = 0.0;
      Rc[BXn][g] = 0.0;
      Rc[BXt][g] = 0.0;
      Rc[BXb][g] = 0.0;
      #if EOS == IDEAL
      Rc[ENG][g] = 0.0;
      #endif
      #if DIVB_CONTROL == EIGHT_WAVES
      Rc[BXn][g] = 1.0;
      eta[g]    = dU[BXn];
      #else
      Rc[BXn][g] = eta[g] = 0.0;
      #endif
      #endif

      // Slow wave:  u - c_s

      scrh = alpha_f*cf*sBx;

      g = KSLOWM;
      lambda[KSLOWM] = u - cs;

      Rc[RHO][g] = alpha_s;
      Rc[MXn][g] = alpha_s*lambda[KSLOWM];
      Rc[MXt][g] = alpha_s*v - scrh*beta_y;
      Rc[MXb][g] = alpha_s*w - scrh*beta_z;
      Rc[BXn][g] = 0.0;
      Rc[BXt][g] = - alpha_f*a*beta_y/sqrt_rho;
      Rc[BXb][g] = - alpha_f*a*beta_z/sqrt_rho;
      #if EOS == IDEAL
      Rc[ENG][g] =   alpha_s*(Hgas - u*cs) - scrh*beta_v
                   - alpha_f*a*Btmag/sqrt_rho;
      #if BACKGROUND_FIELD == YES
      Rc[ENG][g] -= B0y*Rc[BXt][g] + B0z*Rc[BXb][g];
      #endif
      eta[g] =   alpha_s*(X*dV[RHO] + dV[PRS]) - rho*scrh*beta_dv
              - rho*alpha_s*cs*dV[VXn]        - sqrt_rho*alpha_f*a*beta_dB;
      #elif EOS == ISOTHERMAL
      eta[g] =   alpha_s*(0.*X + a2)*dV[RHO] - rho*scrh*beta_dv
               - rho*alpha_s*cs*dV[VXn]      - sqrt_rho*alpha_f*a*beta_dB;
      #endif
      eta[g] *= 0.5/a2;

      // Slow wave:  u + c_s

      g = KSLOWP;
      lambda[KSLOWP] = u + cs;

      Rc[RHO][g] = alpha_s;
      Rc[MXn][g] = alpha_s*lambda[KSLOWP];
      Rc[MXt][g] = alpha_s*v + scrh*beta_y;
      Rc[MXb][g] = alpha_s*w + scrh*beta_z;
      Rc[BXn][g] = 0.0;
      Rc[BXt][g] = Rc[BXt][KSLOWM];
      Rc[BXb][g] = Rc[BXb][KSLOWM];
      #if EOS == IDEAL
      Rc[ENG][g] =   alpha_s*(Hgas + u*cs) + scrh*beta_v
                  - alpha_f*a*Btmag/sqrt_rho;
      #if BACKGROUND_FIELD == YES
      Rc[ENG][g] -= B0y*Rc[BXt][g] + B0z*Rc[BXb][g];
      #endif
      eta[g] =   alpha_s*(X*dV[RHO] + dV[PRS]) + rho*scrh*beta_dv
              + rho*alpha_s*cs*dV[VXn]        - sqrt_rho*alpha_f*a*beta_dB;
      #elif EOS == ISOTHERMAL
      eta[g] =   alpha_s*(0.*X + a2)*dV[RHO] + rho*scrh*beta_dv
              + rho*alpha_s*cs*dV[VXn]      - sqrt_rho*alpha_f*a*beta_dB;
      #endif

      eta[g] *= 0.5/a2;

      // Alfven wave:  u - c_a

      g = KALFVM;
      lambda[KALFVM] = u - ca;

      Rc[RHO][g] = 0.0;
      Rc[MXn][g] = 0.0;
      Rc[MXt][g] = - rho*beta_z;
      Rc[MXb][g] = + rho*beta_y;
      Rc[BXn][g] = 0.0;
      Rc[BXt][g] = - sBx*sqrt_rho*beta_z;
      Rc[BXb][g] =   sBx*sqrt_rho*beta_y;
      #if EOS == IDEAL
      Rc[ENG][g] = - rho*(v*beta_z - w*beta_y);
      #if BACKGROUND_FIELD == YES
      Rc[ENG][g] -= B0y*Rc[BXt][g] + B0z*Rc[BXb][g];
      #endif
      #endif
      eta[g] = + beta_y*dV[VXb]               - beta_z*dV[VXt]
               + sBx/sqrt_rho*(beta_y*dV[BXb] - beta_z*dV[BXt]);

      eta[g] *= 0.5;

      // Alfven wave:  u + c_a

      g = KALFVP;
      lambda[KALFVP] = u + ca;

      Rc[RHO][g] = 0.0;
      Rc[MXn][g] = 0.0;
      Rc[MXt][g] = - Rc[MXt][KALFVM];
      Rc[MXb][g] = - Rc[MXb][KALFVM];
      Rc[BXn][g] = 0.0;
      Rc[BXt][g] =   Rc[BXt][KALFVM];
      Rc[BXb][g] =   Rc[BXb][KALFVM];
      #if EOS == IDEAL
      Rc[ENG][g] = - Rc[ENG][KALFVM];
      #if BACKGROUND_FIELD == YES
      Rc[ENG][g] -= B0y*Rc[BXt][g] + B0z*Rc[BXb][g];
      #endif
      #endif
      eta[g] = - beta_y*dV[VXb]               + beta_z*dV[VXt]
               + sBx/sqrt_rho*(beta_y*dV[BXb] - beta_z*dV[BXt]);

      eta[g] *= 0.5;

      // 2h. Compute maximum signal velocity

      cmax(i) = fabs(u) + cf;
      s_cs2(i)  = a2;

      #pragma acc loop seq
      for (nv = 0; nv < NFLX; nv++) alambda[nv] = fabs(lambda[nv]);

      // 2l. Entropy Fix

      double delta = 1.e-6;

      if (alambda[KFASTM] < 0.5*delta) {
        alambda[KFASTM] = lambda[KFASTM]*lambda[KFASTM]/delta + 0.25*delta;
      }
      if (alambda[KFASTP] < 0.5*delta) {
        alambda[KFASTP] = lambda[KFASTP]*lambda[KFASTP]/delta + 0.25*delta;
      }

      if (alambda[KSLOWM] < 0.5*delta) {
        alambda[KSLOWM] = lambda[KSLOWM]*lambda[KSLOWM]/delta + 0.25*delta;
      }
      if (alambda[KSLOWP] < 0.5*delta) {
        alambda[KSLOWP] = lambda[KSLOWP]*lambda[KSLOWP]/delta + 0.25*delta;
      }

      // 2i.  Compute Roe numerical flux

      #pragma acc loop seq
      for (nv = 0; nv < NFLX; nv++) {
        double dflux = 0.0;

        dflux += alambda[KFASTM]*eta[KFASTM]*Rc[nv][KFASTM];
        dflux += alambda[KFASTP]*eta[KFASTP]*Rc[nv][KFASTP];
        #if HAVE_ENERGY
        dflux += alambda[KENTRP]*eta[KENTRP]*Rc[nv][KENTRP];
        #endif
        #if DIVB_CONTROL != DIV_CLEANING
        dflux += alambda[KDIVB]*eta[KDIVB]*Rc[nv][KDIVB];
        #endif
        dflux += alambda[KSLOWM]*eta[KSLOWM]*Rc[nv][KSLOWM];
        dflux += alambda[KSLOWP]*eta[KSLOWP]*Rc[nv][KSLOWP];
        dflux += alambda[KALFVM]*eta[KALFVM]*Rc[nv][KALFVM];
        dflux += alambda[KALFVP]*eta[KALFVP]*Rc[nv][KALFVP];
        #if DIVB_CONTROL == DIV_CLEANING
        dflux += alambda[KPSI_GLMM]*eta[KPSI_GLMM]*Rc[nv][KPSI_GLMM];
        dflux += alambda[KPSI_GLMP]*eta[KPSI_GLMP]*Rc[nv][KPSI_GLMP];
        #endif

        flux(i,nv) = 0.5*(fl[nv] + fr[nv] - dflux);
      }

      press(i) = 0.5*(pl + pr);

      // 2j. Check the Roe matrix condition, FR - FL = A*(UR - UL)
      //     where A*(UR - UL) = R*lambda*eta.

      #if CHECK_ROE_MATRIX == YES
      #pragma acc loop seq
      for (nv = 0; nv < NFLX; nv++){
        dV[nv] = fr[nv] - fl[nv];
        dV[nv] -= lambda[KFASTM]*eta[KFASTM]*Rc[nv][KFASTM];
        dV[nv] -= lambda[KFASTP]*eta[KFASTP]*Rc[nv][KFASTP];
        #if HAVE_ENERGY
        dV[nv] -= lambda[KENTRP]*eta[KENTRP]*Rc[nv][KENTRP];
        #endif
        #if DIVB_CONTROL != DIV_CLEANING
        dV[nv] -= lambda[KDIVB]*eta[KDIVB]*Rc[nv][KDIVB];
        #endif
        dV[nv] -= lambda[KSLOWM]*eta[KSLOWM]*Rc[nv][KSLOWM];
        dV[nv] -= lambda[KSLOWP]*eta[KSLOWP]*Rc[nv][KSLOWP];
        dV[nv] -= lambda[KALFVM]*eta[KALFVM]*Rc[nv][KALFVM];
        dV[nv] -= lambda[KALFVP]*eta[KALFVP]*Rc[nv][KALFVP];
        #if DIVB_CONTROL == DIV_CLEANING
        dV[nv] -= lambda[KPSI_GLMM]*eta[KPSI_GLMM]*Rc[nv][KPSI_GLMM];
        dV[nv] -= lambda[KPSI_GLMP]*eta[KPSI_GLMP]*Rc[nv][KPSI_GLMP];
        #endif
        // if (fabs(dV[nv]) > 1.e-4){
        //     local_error          |= ERR_MATRIX_NO_COND;
        //     data->err(ti,tj,tk)  |= ERR_NEGATIVE_SOUND_SPEED;
        // }
      }
      #endif

      // 2k. Hybridize with HLL solver: replace occurences of unphysical
      //      states (p < 0, rho < 0) with HLL Flux. Reference:
      //
      //     "A Positive Conservative Method for MHD based based on HLL
      //      and Roe methods", P. Janhunen, JCP (2000), 80, 649

      Sl  = lambda[KFASTM];
      Sr  = lambda[KFASTP];

      #if HLL_HYBRIDIZATION == YES
      double Uv[NFLX];
      if (Sl < 0.0 && Sr > 0.0){
        ifail = 0;

        // Check left state

        #if EOS == ISOTHERMAL
        Uv[RHO] = ul[RHO] + (flux(i,RHO) - fl[RHO])/Sl;
        ifail  = (Uv[RHO] < 0.0);
        #else
        NFLX_LOOP(nv) {
          Uv[nv] = ul[nv] + (flux(i,nv) - fl[nv])/Sl;
        }
        Uv[MXn] += (press(i) - pl)/Sl;

        vel2  = Uv[MX1]*Uv[MX1] + Uv[MX2]*Uv[MX2] + Uv[MX3]*Uv[MX3];
        b2    = Uv[BX1]*Uv[BX1] + Uv[BX2]*Uv[BX2] + Uv[BX3]*Uv[BX3];
        prr    = Uv[ENG] - 0.5*vel2/Uv[RHO] - 0.5*b2;
        ifail = (prr < 0.0) || (Uv[RHO] < 0.0);
        #endif

        // Check right state

        #if EOS == ISOTHERMAL
        Uv[RHO] = ur[RHO] + (flux(i,RHO) - fr[RHO])/Sr;
        ifail  = (Uv[RHO] < 0.0);
        #else
        for (nv = 0; nv < NFLX;  nv++){
          Uv[nv] = ur[nv] + (flux(i,nv) - fr[nv])/Sr;
        }
        Uv[MXn] += (press(i) - pr)/Sr;

        vel2  = Uv[MX1]*Uv[MX1] + Uv[MX2]*Uv[MX2] + Uv[MX3]*Uv[MX3];
        b2    = Uv[BX1]*Uv[BX1] + Uv[BX2]*Uv[BX2] + Uv[BX3]*Uv[BX3];
        prr    = Uv[ENG] - 0.5*vel2/Uv[RHO] - 0.5*b2;
        ifail = (prr < 0.0) || (Uv[RHO] < 0.0);
        #endif

        #if DIMENSIONS > 1
        // Use the HLL flux function if the interface
        //  lies within a strong shock.
        //  The effect of this switch is visible
        //  in the Mach reflection test.

        dvel = vr[VXn] - vl[VXn];
        #if EOS == ISOTHERMAL
        scrh  = fabs(vl[RHO] - vr[RHO]);
        scrh /= MIN(vl[RHO], vr[RHO]);
        #else
        scrh  = fabs(vl[PRS] - vr[PRS]);
        scrh /= MIN(vl[PRS], vr[PRS]);
        #endif
        if (scrh > 1.0 && (vr[VXn] < vl[VXn])){
          scrh = 1.0/(Sr - Sl);
          #pragma acc loop seq
          for (nv = 0; nv < NFLX; nv++) {
            flux(i,nv) = Sl*Sr*(ur[nv] - ul[nv]) +
                      Sr*fl[nv] - Sl*fr[nv];
            flux(i,nv) *= scrh;
          }
          press(i) = (Sr*pl - Sl*pr)*scrh;
        }
        #endif // DIMENSIONS > 1
      }
      #endif  // HLL_HYBRIDIZATION == YES 
    } else { // Use HLL 

      HLL_Speed<dir> (stateL, stateR, Sl, Sr);
      double scrh = MAX(fabs(Sl), fabs(Sr));
      cmax(i) = scrh;

      double bL = MIN(0.0, Sl);
      double bR = MAX(0.0, Sr);

      scrh = 1.0 / (bR - bL);
      #pragma acc loop seq
      for (nv = 0; nv < NFLX; nv++) {
        flux(i,nv)  = bL*bR*(ur[nv] - ul[nv]) + bR*fl[nv] - bL*fr[nv];
        flux(i,nv) *= scrh;
      }
      press(i) = (bR*pl - bL*pr)*scrh;

    }  // end if (switch_to_hll) 

    // 2l. Compute Powell's / GLM interface states

    #if DIVB_CONTROL == EIGHT_WAVES
    Bn(i) = 0.5*(vl[BXn] + vr[BXn]);
    #elif DIVB_CONTROL == DIV_CLEANING
    Bn(i)   = flux(i,PSI_GLM)/ch2;
    PSIn(i) = flux(i,BXn);
    #endif

    // 2m. Compute tracer fluxes

    #if NSCL > 0
    if (flux(i,RHO) > 0.0){
      #pragma acc loop seq
      NSCL_LOOP(nv) flux(i,nv) = flux(i,RHO)*vl[nv];
    }else{
      #pragma acc loop seq
      NSCL_LOOP(nv) flux(i,nv) = flux(i,RHO)*vr[nv];
    }
    #endif

    // 2n. Check for NaNs in flux

    #pragma acc loop seq
    for (nv = 0; nv < NFLX; nv++){
      if (isnan(flux(i,nv))) {
        ATOMIC_ERR_UPDATE(err_kernel, ERR_NAN, ERR_QUIT);
        KERNEL_ERR_UPDATE(data->err(ti,tj,tk),  ERR_NAN, ERR_QUIT); 
      }
    }

    // 2o. Store emf-flux components for CT

    #if (DIVB_CONTROL == CONSTRAINED_TRANSPORT)
    #if (CT_EMF_AVERAGE == ARITHMETIC) || (CT_EMF_AVERAGE == CT_CONTACT)
    CT_Store_Arith<dir>(data, flux, i, ti, tj, tk);
    #elif CT_EMF_AVERAGE == CT_FLUX
    CT_Store_Flux<dir>(data, fl, fr, flux, i, ti, tj, tk);
    #elif CT_EMF_AVERAGE == UCT_HLL
    CT_Store_UCTHll<dir>(data, vl, vr, Sl, Sr, ti, tj, tk);
    #elif CT_EMF_AVERAGE == UCT_HLLD
    CT_Store_UCTHlld<dir> (data, ul, ur, vl, vr, Sl, Sr, ti, tj, tk);
    #endif
    #endif

    // 2p. Conservative Energy correction step
    //        [Experimental]
  
    //   #if (CT_ENERGY_CORRECTION == 2) || (CT_ENERGY_CORRECTION == 3)
    //   double EgL, EgR;
    //   double fEgL   = HD_EnergyFlux(vl, &EgL, dir, run_config);
    //   double fEgR   = HD_EnergyFlux(vr, &EgR, dir, run_config);
    //   double fEg_lf = 0.5*(fEgL + fEgR) - 0.5*cmax(i)*(EgR - EgL); 

    //   flux(i,ENG) = fEg_lf;
    //   #endif

  }}}  // End for (i,j,k)

  // ------------------------------------------------------
  // 3. Error Check / Quit condition
  // ------------------------------------------------------

  #if WARNING_MESSAGES > 0 
  Where ("Roe_Solver()", data, err_kernel, box, grid);
  #endif

  int quit_condition = (err_kernel & ERR_NAN || err_kernel & ERR_RIEMANN_SOLVE_FAIL);

  #if FALLBACK_MODE != NO
  quit_condition *= (data->integration_mode != DEFAULT);
  Fallback_Trigger(data->err, data->flag, box);
  #endif

  if (quit_condition) {
    QuitPLUTO("from Roe_Solver(): failure (proabably a2 < 0 or NaN in flux)",
               err_kernel);
  }

}
#undef CHECK_ROE_MATRIX
