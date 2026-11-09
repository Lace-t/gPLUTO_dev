#include "gpluto.hpp" 
 
//////////////////////////////////////////////////////////////////////// 
void Init (double *us, double x1, double x2, double x3, RunConfig &run_config) 
/*!
*/
//////////////////////////////////////////////////////////////////////// 
{  
  double gmm = 5./3.;	
  double *param = run_config.inputParam;
  run_config.gamma = gmm;
  
  // ------------------------------------------------------
  // 0. Set problem units 
  // ------------------------------------------------------

  run_config.unit_length   = 100.0e3*CONST_pc;
  run_config.unit_velocity = CONST_c;
  run_config.unit_density  = 0.05*CONST_amu;

  // ------------------------------------------------------
  // 1. Define problem parameters
  // ------------------------------------------------------
  
  double rx2 = param[RADX]*param[RADX];
  double ry2 = param[RADY]*param[RADY];
  double rad = sqrt(rx2*(x1*x1) + ry2*(x2*x2));  
   
  double h   = 2.5*param[TEMP]+sqrt(9.*param[TEMP]*param[TEMP]/4. + 1);
  double bt  = sqrt(param[NU]*h*param[GAMMA_UP]*param[GAMMA_UP]*param[SIGMA]);

  us[VX1] = 0.0; 
  us[VX2] = 0.0; 
  us[VX3] = 0.0;
  us[BX1] = 0.5*bt;
  us[BX2] = 0.0;
  us[BX3] = 0.0;
  us[TRC] = 0.0; 
  us[TRC+1] = 1.0;
   
  us[RHO] = pow(1+pow(rad,2),-1.5*param[BETA]); 
  us[PRS] = param[TEMP]*param[NU]; 

  #ifdef PSI_GLM
  us[PSI_GLM] = 0.0;
  #endif

} 

////////////////////////////////////////////////////////////////////////
void InitDomain (Data *data, Grid *grid)
/*!
*/
////////////////////////////////////////////////////////////////////////
{
}

////////////////////////////////////////////////////////////////////////
void Analysis(Data *data, Grid *grid)
/*!
 */ 
//////////////////////////////////////////////////////////////////////// 
{ 

} 
 
//////////////////////////////////////////////////////////////////////// 
void UserDefBoundary (const Data *data, RBox *box, int side, Grid *grid) 
/*!  
  restart;
  eqBn   := Bn   = (BL   + BR)/2   - (psiR - psiL)/2:
  eqpsin := psin = (psiL + psiR)/2 - (BR   - BL)/2:
  solve ({eqBn, eqpsin}, {psiL, psin});
*/
//////////////////////////////////////////////////////////////////////// 
{ 
  int    i, j, k; 
  int    ibeg = grid->lbeg[IDIR];
  #ifdef PSI_GLM
  iend = grid->lend[IDIR];
  #endif
  int    nv;
  double x1, x2, prof; 
  double rjc;
  double vjup;
  double h, bt; 
  RunConfig &run_config = *(data->run_config);
  const double *param = run_config.inputParam;

  vjup = sqrt(1.0 - 1.0/(param[GAMMA_UP]*param[GAMMA_UP])); 
  
  rjc	= 0.01; 
  h   = 2.5*param[TEMP]+sqrt(9.*param[TEMP]*param[TEMP]/4. + 1);
  bt  = sqrt(param[NU]*h*param[GAMMA_UP]*param[GAMMA_UP]*param[SIGMA]);
 
  // ------------------------------------------------------
  // 0. Assign jet value from internal boundary 
  // ------------------------------------------------------

  if (side == 0) { 
    //nst = g_stepNumber;
    #pragma acc parallel loop collapse(3) present(box,grid,data,param) 
    BOX_LOOP(*box,i,j,k) {

      x1 = grid->xC[IDIR][i];
      x2 = grid->xC[JDIR][j];
      // rr = sqrt(x2*x2)/rjc;

      // double bp = (bt*0.5)/(1. + ((rr/0.5)*(rr/0.5)));
	
      if ( (sqrt(x2*x2) < rjc) && (x1 >= 0.0) && (x1 <= 0.01) ){
        data->Vc(i,j,k,RHO) =  param[NU];

        data->Vc(i,j,k,VX1) =  vjup;  
        data->Vc(i,j,k,VX2) =  0.0;
        data->Vc(i,j,k,VX3) =  0.0;

        data->Vc(i,j,k,BX1) =  0.5*bt; // bp;
        data->Vc(i,j,k,BX2) =  0.0; // (0.5*bt)*rr; 
        data->Vc(i,j,k,BX3) =  0.0;

        #ifdef PSI_GLM
        if (grid->xC[IDIR][i+1] > 0.01){
          double Bn = 0.5*bt;
          double BL = Bn, BR = data->Vc(i+1,j,k,BX1);
          double psiR = data->Vc(i+1,j,k,PSI_GLM);
          data->Vc(i,j,k,PSI_GLM) = 2.0*Bn + psiR - (BL + BR);
        }
        if (grid->xC[IDIR][i-1] < 0.0){
          double Bn = 0.5*bt;
          double BR = data->Vc(i-1,j,k,BX1), BL = Bn;
          double psiL = data->Vc(i-1,j,k,PSI_GLM);
          data->Vc(i,j,k,PSI_GLM) = psiL -2.0*Bn + (BL + BR);
        }
        #endif

        data->Vc(i,j,k,PRS) =  param[TEMP]*param[NU]; 
        data->Vc(i,j,k,TRC) =  1.;
        data->Vc(i,j,k,TRC+1) = 0.0;
        data->flag(i,j,k)   |= FLAG_INTERNAL_BOUNDARY;  
        // data->flag(i,j,k)   |= FLAG_MINMOD;
        // data->flag(i,j,k)   |= FLAG_HLL;
      }
    }
    // Convert to conservative
    PrimToCons3D(data, 0, grid, box);
  }

  // ------------------------------------------------------
  // 1. Assign jet value from regular boundary 
  // ------------------------------------------------------

  if (side == X1_BEG){  // -- X2_BEG boundary -- 
    #pragma acc parallel loop collapse(3) present(box, data)
    BOX_LOOP(*box,i,j,k){
      double vjet[NVAR], vout[NVAR];

      // 1a. Set jet values 

      x1 = grid->xC[IDIR][i];
      x2 = grid->xC[JDIR][j];
      // rr = sqrt(x2*x2)/rjc;
      // double bp = (bt*0.5)/(1. + ((rr/0.5)*(rr/0.5)));
      vjet[RHO] =  param[NU];

      vjet[VX1] =  vjup;  
      vjet[VX2] =  0.0;
      vjet[VX3] =  0.0;

      vjet[BX1] =  0.5*bt; // bp;
      vjet[BX2] =  0.0;    // (0.5*bt)*rr; 
      vjet[BX3] =  0.0;

      #ifdef PSI_GLM
      vjet[PSI_GLM] =  0.0;

      double Bn = 0.5*bt;
      double BL = Bn, BR = data->Vc(box->iend+1,j,k,BX1);
      double psiR = data->Vc(box->iend+1,j,k,PSI_GLM);
      data->Vc(i,j,k,PSI_GLM) = 2.0*Bn + psiR - (BL + BR);
      #endif

      vjet[PRS] =  param[TEMP]*param[NU]; 
      vjet[TRC] =  1.;
      vjet[TRC+1] = 0.0;

      // 1b. Copy and reflect ambient medium 

      NVAR_LOOP(nv) vout[nv] = data->Vc(i,j,k,nv);
      vout[VX1] *= -1.0;
      vout[BX2] *= -1.0;
      vout[BX3] *= -1.0;

      // 1c. Assign b.c. on primitive variables.

      prof = (sqrt(x2*x2) < rjc);
      NVAR_LOOP(nv) data->Vc(i,j,k,nv) = vout[nv] - (vout[nv] - vjet[nv])*prof;

      // 1d. Lower accuracy close to boundary

      if (prof > -1.e-3){
        data->flag(ibeg,j,k) |= FLAG_MINMOD;
        data->flag(ibeg,j,k) |= FLAG_HLL;
      }
    }

    #ifdef STAGGERED_MHD
    #pragma acc parallel loop collapse(3) present(box, data)
    BOXJ_LOOP(*box,i,j,k){
      double vout_BX2 = -data->Vs[BX2s](2*ibeg-i-1,j,k);
      double vjet_BX2 = 0.0;
      x2 = grid->xL[JDIR][j];
      prof = (fabs(x2/rjc) <= 1.0);
      data->Vs[BX2s](i,j,k) = vout_BX2 - (vout_BX2 - vjet_BX2)*prof;
    }
    #endif

  }
}

////////////////////////////////////////////////////////////////////////
#pragma acc routine seq
void BodyForceVector(const State& state, const RunConfig& run_config, double *g)
/*!
*/
////////////////////////////////////////////////////////////////////////
{
}

#if CT_EMF_BOUNDARY == YES
void CT_EMF_Boundary(const Data *data, Grid *grid) 
{
  int i,j,k; 
  int ibeg = grid->lbeg[IDIR], iend = grid->lend[IDIR] + INCLUDE_IDIR;
  int jbeg = grid->lbeg[JDIR], jend = grid->lend[JDIR] + INCLUDE_JDIR;
  int kbeg = grid->lbeg[KDIR], kend = grid->lend[KDIR] + INCLUDE_KDIR;
  EMF *emf = data->emf;

  for (k = kbeg; k <= kend; k++){ 
  for (j = jbeg; j <= jend; j++){
  for (i = ibeg; i <= iend; i++){

    //interior
    if ((data->flag(i,j,k) & FLAG_INTERNAL_BOUNDARY)) {
      emf->Ee[KDIR](i,j,k) = 0.0; 
    }

    // Left
    if (  !(data->flag(i-1,j,k) & FLAG_INTERNAL_BOUNDARY)
        && (data->flag(i,j,k)   & FLAG_INTERNAL_BOUNDARY)) {
      emf->Ee[KDIR](i,j,k) = 0.0; 
    }

    // Right
    if (  !(data->flag(i,j,k)   & FLAG_INTERNAL_BOUNDARY)
        && (data->flag(i-1,j,k) & FLAG_INTERNAL_BOUNDARY)) {
      emf->Ee[KDIR](i,j,k) = 0.0; 
    }

    // Bottom
    if (  !(data->flag(i,j-1,k) & FLAG_INTERNAL_BOUNDARY)
        && (data->flag(i,j,k)   & FLAG_INTERNAL_BOUNDARY)) {
      emf->Ee[KDIR](i,j,k) = 0.0; 
    }

    // Top
    if (  !(data->flag(i,j,k)   & FLAG_INTERNAL_BOUNDARY)
        && (data->flag(i,j-1,k) & FLAG_INTERNAL_BOUNDARY)) {
      emf->Ee[KDIR](i,j,k) = 0.0; 
    }
  }}}
}
#endif