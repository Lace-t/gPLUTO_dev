////////////////////////////////////////////////////////////////////////
/*!
  \file
  \brief GLM module implementation.

  Contains functions for the GLM module.

  \authors A. Mignone    (andrea.mignone@unito.it)
  \date    Apr 1, 2026

  \b Reference
     - "Hyperbolic Divergence Cleaning for the MHD equations"
        Dedner et al., JCP (2002) 175, 645
     - "A second-order unsplit Godunov scheme for cell-centered MHD:
        the CTU-GLM scheme" Mignone \& Tzeferacos, JCP (2010) 229, 2117
*/
////////////////////////////////////////////////////////////////////////
#include "gpluto.hpp"

#ifdef GLM_MHD

////////////////////////////////////////////////////////////////////////
void GLM_Init (const Data *data, Grid *grid)
/*!
  Initialize the maximum propagation speed glm_ch at the beginning
  of integration cycle.
*/
////////////////////////////////////////////////////////////////////////
{
  double maxc;
  double glm_ch = data->run_config->glm_ch;

  if (glm_ch > 0.0) return;

#if PHYSICS == MHD
  int    i,j,k;
  double lambda, tau, A2, scrh0, scrh1;
  double a2, b2, bn2, bt2, At2, cf;

  RBox box;
  RBoxDefine (grid->lbeg, grid->lend, &box);

  maxc = 0.0;

#pragma acc data copy(maxc)
{
  #pragma omp target teams distribute parallel for collapse(3) \
              reduction(max:maxc,lambda) map(to:data->run_config)// map(to: box , grid)
  #pragma acc parallel loop collapse(3) present(data, data->run_config) copyin(box) \
                                        reduction(max:maxc) 
  BOX_LOOP(box,i,j,k){

    #if EOS != ISOTHERMAL
    a2 = data->run_config->gamma*data->Vc(i,j,k,PRS)/data->Vc(i,j,k,RHO);
    #elif EOS == ISOTHERMAL
    double a = data->run_config->isoSoundSpeed;
    a2 = a*a;
    #endif

    tau = 1.0/data->Vc(i,j,k,RHO);
    // sqrt_rho = sqrt(data->Vc(i,j,k,RHO));

    // 1. Compute the fast waves in the three directions
    //    simultanously

    // X1-dir 

    bn2 =   data->Vc(i,j,k,BX1)*data->Vc(i,j,k,BX1);                   // > 0 
    bt2 =   data->Vc(i,j,k,BX2)*data->Vc(i,j,k,BX2) 
          + data->Vc(i,j,k,BX3)*data->Vc(i,j,k,BX3);   // > 0 

    b2  = bn2 + bt2;     //  >0     
    A2  = b2*tau;        //  >0     
    At2 = bt2*tau;

    scrh1 = a2 - A2;
    scrh0 = sqrt(scrh1*scrh1 + 4.0*a2*At2);      //   >0   

    cf = sqrt(0.5*(a2 + A2 + scrh0));

    lambda = fabs(data->Vc(i,j,k,VX1)) + cf;
    maxc = MAX(maxc,lambda);

    // X2-dir 

    #if DIMENSIONS > 1
    bn2 =   data->Vc(i,j,k,BX2)*data->Vc(i,j,k,BX2);   // > 0 
    bt2 =   data->Vc(i,j,k,BX1)*data->Vc(i,j,k,BX1) 
          + data->Vc(i,j,k,BX3)*data->Vc(i,j,k,BX3);   // > 0 

    b2  = bn2 + bt2;     //  >0  
    A2  = b2*tau;        //  >0  
    At2 = bt2*tau;

    scrh1 = a2 - A2;
    scrh0 = sqrt(scrh1*scrh1 + 4.0*a2*At2);      //   >0   

    cf = sqrt(0.5*(a2 + A2 + scrh0));

    lambda = fabs(data->Vc(i,j,k,VX2)) + cf;
    maxc = MAX(maxc,lambda);
    #endif

    // X3-dir 

    #if DIMENSIONS == 3
    bn2 = data->Vc(i,j,k,BX3)*data->Vc(i,j,k,BX3);                   // > 0 
    bt2 = data->Vc(i,j,k,BX1)*data->Vc(i,j,k,BX1) + data->Vc(i,j,k,BX2)*data->Vc(i,j,k,BX2);   // > 0 

    b2  = bn2 + bt2;     //  >0    
    A2  = b2*tau;        //  >0    
    At2 = bt2*tau;

    scrh1 = a2 - A2;
    scrh0 = sqrt(scrh1*scrh1 + 4.0*a2*At2);      //   >0   

    cf = sqrt(0.5*(a2 + A2 + scrh0));

    lambda = fabs(data->Vc(i,j,k,VX3)) + cf;
    maxc = MAX(maxc,lambda);
    #endif

  }
}
  glm_ch = maxc;

#elif (PHYSICS == RMHD) || (PHYSICS == ResRMHD)

  if (glm_ch < 0.0) glm_ch = maxc = 1.0;

#endif

#ifdef PARALLEL
  double gmaxc;
  MPI_Allreduce (&maxc, &gmaxc, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  maxc = gmaxc;
  glm_ch = maxc;
#endif

  data->run_config->glm_ch = glm_ch;
  #pragma acc update device(data->run_config)
  #pragma omp target update to(data->run_config)

}

////////////////////////////////////////////////////////////////////////
/*!
  Include the damping source term of the Lagrangian multiplier
  equation in a split fashion for the mixed GLM formulation.
  Ref. Mignone & Tzeferacos, JCP (2010) 229, 2117, Equation (27).
*/
////////////////////////////////////////////////////////////////////////
void GLM_Source (Data *data, double dt, Grid *grid)
{
  int    i,j,k;
  double scrh;
  double glm_ch = data->run_config->glm_ch;
  DIM_EXPAND(double *dx1 = grid->dx[IDIR];,
             double *dx2 = grid->dx[JDIR];,
             double *dx3 = grid->dx[KDIR];)
  RBox box;
  box.ibeg = grid->lbeg[IDIR], box.iend = grid->lend[IDIR];
  box.jbeg = grid->lbeg[JDIR], box.jend = grid->lend[JDIR];
  box.kbeg = grid->lbeg[KDIR], box.kend = grid->lend[KDIR];

  #ifdef _USE_NVTX
  mynvtxstart_("GLM_Source",RED);
  #endif

  // ------------------------------------------------------
  // 1. Apply source step to Lagrangian multiplier:
  //    psi = psi(0)*exp(-GLM_ALPHA*(dt/dh)*ch) 
  // ------------------------------------------------------

  #if (PHYSICS == MHD) || (PHYSICS == RMHD) || (PHYSICS == ResRMHD)
  #pragma omp target teams distribute parallel for collapse(3)
  #pragma acc parallel loop collapse(3) present(data) copyin(box)
  BOX_LOOP(box,i,j,k){

    double dl[3], dtdh, dh;
    Metric metric;
    #if GEOMETRY != CARTESIAN
    MetricSet (metric, grid->rC[i], grid->sC[j]);
    #endif

    // 1a. Compute cell length (dh).
    //     The choice to take the maximum (instead of the minimum or other 
    //     norms) is justified only by validation or test robustness 
    //     (see, e.g., conf #14 in RMHD/Jet.) 

    DIM_EXPAND(dl[IDIR] = metric.h1*dx1[i];,
               dl[JDIR] = metric.h2*dx2[j];,
               dl[KDIR] = metric.h3*dx3[k];)

    DIM_EXPAND(dh = dl[IDIR];         ,
               dh = MAX(dh, dl[JDIR]);,
               dh = MAX(dh, dl[KDIR]);)
 
    // cp   = sqrt(dh*glm_ch/GLM_ALPHA);

    // 1b. Integrate exactly: 
    //     psi = psi(0)*exp(-alpha*dtdh*ch)

    dtdh = dt/dh;
    scrh = dtdh*glm_ch*GLM_ALPHA;
    scrh = exp(-scrh);

    data->Uc(i,j,k,PSI_GLM) *= scrh;
    data->Vc(i,j,k,PSI_GLM) *= scrh;
    #if (PHYSICS == ResRMHD) && (DIVE_CONTROL == DIV_CLEANING)
    data->Vc(i,j,k,PHI_GLM) *= scrh;
    data->Uc(i,j,k,PHI_GLM) *= scrh;
    #endif
  }
  #endif

  // ------------------------------------------------------
  // 2. For high order, recover point values
  // ------------------------------------------------------

  #ifdef HIGH_ORDER
  BoundaryConservative (data, grid);  // Set b.c. on data->Uc (1 ghost)
  PointValue(data, grid, &box);       // Obtain point values Upc, Vpc
  #endif

#ifdef _USE_NVTX
mynvtxstop_();
#endif

  return;
}
#endif // GLM_MHD 
