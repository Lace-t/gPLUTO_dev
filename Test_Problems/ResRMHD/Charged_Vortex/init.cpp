////////////////////////////////////////////////////////////////////////
/*!
  \file
  \brief Vortex problem

*/
////////////////////////////////////////////////////////////////////////
#include "gpluto.hpp"

////////////////////////////////////////////////////////////////////////
void Init (double *v, double x1, double x2, double x3, RunConfig &run_config)
/*!
  
*/
////////////////////////////////////////////////////////////////////////
{
  double *inputParam = run_config.inputParam;
  double q0  = inputParam[CRG0];
  double p0  = inputParam[P0];
  double r2  = x1*x1 + x2*x2;
  double r   = sqrt(r2);
  double r_1 = 1.0/r;
  double Gamma1;
  double scrh = r2 + 1;
  double c = x1*r_1, s = x2*r_1, Er, Bz, vphi;

  run_config.gamma = (4.0/3.0);

  Gamma1  = run_config.gamma/(run_config.gamma - 1); 
  
  v[RHO] = 1.0;

  Er     = 0.5*q0*r/scrh;
  v[EX1] = Er*c;
  v[EX2] = Er*s;
  v[EX3] = 0.0;

  v[BX1] = 0.0;
  v[BX2] = 0.0;
  v[BX3] = Bz = sqrt(scrh*scrh - q0*q0*0.25)/scrh;

  vphi = -Er/Bz;

  v[VX1] = -vphi*s;
  v[VX2] =  vphi*c;
  v[VX3] = 0.0;

  #ifdef CRG
  v[CRG] = q0/(scrh*scrh);
  #endif
 
  scrh    = (4.0*scrh - q0*q0)/(scrh*(4.0 - q0*q0));
  v[PRS]  = pow(scrh, 0.5*Gamma1)*(p0 + v[RHO]/Gamma1);
  v[PRS] -= v[RHO]/Gamma1;

  #ifdef PHI_GLM
  v[PHI_GLM] = 0.0;
  #endif
  #ifdef PSI_GLM
  v[PSI_GLM] = 0.0;
  #endif

  #if ENTROPY_SWITCH != NO
  v[ENTR] = Entropy(v, run_config);
  #endif

}

////////////////////////////////////////////////////////////////////////
void InitDomain (Data *data, Grid *grid)
/*! 
  Assign initial condition by looping over the computational domain.
  Called after the usual Init() function to assign initial conditions
  on primitive variables.
  Value assigned here will overwrite those prescribed during Init().
*/
////////////////////////////////////////////////////////////////////////
{
}

////////////////////////////////////////////////////////////////////////
void Analysis (Data *data, Grid *grid)
/* 
 *
*/
////////////////////////////////////////////////////////////////////////
{
}
////////////////////////////////////////////////////////////////////////
void UserDefBoundary (const Data *data, RBox *box, int side, Grid *grid) 
/*
 *
*/
////////////////////////////////////////////////////////////////////////
{
  int    i,j,k, nv;
  double *x = grid->xC[IDIR], *xl = grid->xL[IDIR];
  double *y = grid->xC[JDIR], *yl = grid->xL[JDIR];
  double *z = grid->xC[KDIR]; 
  #ifdef HIGH_ORDER
  double *xr = grid->xR[IDIR]; 
  double *yr = grid->xR[JDIR];
  double *zr = grid->xR[KDIR], *zl = grid->xL[KDIR];
  #endif
  double *dx = grid->dx[IDIR];
  RunConfig &run_config = *(data->run_config);
  #if DIVE_CONTROL == CONSTRAINED_TRANSPORT
  double r2, scrh;
  double *inputParam = run_config.inputParam;
  #endif

  if (side == X1_BEG || side == X1_END){ 
    #pragma acc parallel loop collapse(3) present (data, run_config)
    BOX_LOOP(*box,i,j,k){
      double v[16];
      Init (v, x[i], y[j], z[k], run_config);
      NVAR_LOOP(nv) data->Vc(i,j,k,nv) = v[nv];

      #ifdef HIGH_ORDER
      int Ng = HO_STARTUP_NGAUSS;
      double u_av[16], v_av[16];
      State state;
      STATE_SET(state, v_av, u_av, NULL, x[i], y[j], z[k]);
      state.flag = &(data->flag(i,j,k));

      InitGaussAverage (NULL, xl[i], xr[i], yl[j], yr[j], zl[k], zr[k], 
                        u_av, Ng, Ng, Ng, run_config);
      NVAR_LOOP(nv) data->Vpc(i,j,k,nv) = v[nv];
      NVAR_LOOP(nv) data->Vc(i,j,k,nv)  = v_av[nv];
      NVAR_LOOP(nv) data->Uc(i,j,k,nv)  = u_av[nv];
      #endif

      #if DIVE_CONTROL == CONSTRAINED_TRANSPORT
      r2   = x[i]*x[i] + y[j]*y[j];
      scrh = r2 + 1.0;
      data->q(i,j,k) = inputParam[CRG0]/(scrh*scrh);
      #endif
    }

    #pragma acc parallel loop collapse(3) present (data, run_config)
    BOXJ_LOOP(*box, i,j,k){
      double v[16];

      Init (v, x[i], yl[j], z[k], run_config);
      #ifdef STAGGERED_MHD
      data->Vs[BX2s](i,j,k) = v[BX2];
      #if DIVE_CONTROL == CONSTRAINED_TRANSPORT
      data->Vs[EX2s](i,j,k) = v[EX2];
      #endif
      #endif
    }
  }

  if (side == X2_BEG || side == X2_END){ 
    #pragma acc parallel loop collapse(3) present (data, run_config)
    BOX_LOOP(*box,i,j,k){
      double v[16];
      Init (v, x[i], y[j], z[k], run_config);
      NVAR_LOOP(nv) data->Vc(i,j,k,nv) = v[nv];

      #ifdef HIGH_ORDER
      int Ng = HO_STARTUP_NGAUSS;
      double u_av[16], v_av[16];
      State state;
      STATE_SET(state, v_av, u_av, NULL, x[i], y[j], z[k]);
      state.flag = &(data->flag(i,j,k));

      InitGaussAverage (NULL, xl[i], xr[i], yl[j], yr[j], zl[k], zr[k], 
                        u_av, Ng, Ng, Ng, run_config);

      NVAR_LOOP(nv) data->Vpc(i,j,k,nv) = v[nv];
      NVAR_LOOP(nv) data->Vc(i,j,k,nv)  = v_av[nv];
      NVAR_LOOP(nv) data->Uc(i,j,k,nv)  = u_av[nv];
      #endif

      #if DIVE_CONTROL == CONSTRAINED_TRANSPORT
      r2   = x[i]*x[i] + y[j]*y[j];
      scrh = r2 + 1.0;
      data->q(i,j,k) = inputParam[CRG0]/(scrh*scrh);
      #endif
    }

    #pragma acc parallel loop collapse(3) present (data, run_config)
    BOXI_LOOP(*box, i,j,k){
      double v[16];
      if (i < box->iend) Init (v, xl[i], y[j], z[k], run_config);
      else               Init (v, x[i-1] + dx[i-1]*0.5, y[j], z[k], run_config);
      #ifdef STAGGERED_MHD
      data->Vs[BX1s](i,j,k) = v[BX1];
      #if DIVE_CONTROL == CONSTRAINED_TRANSPORT
      data->Vs[EX1s](i,j,k) = v[EX1];
      #endif
      #endif 
    }
  }
}

