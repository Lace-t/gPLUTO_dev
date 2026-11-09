////////////////////////////////////////////////////////////////////////
/*!
  \file
  \brief Stellar wind test problem.

  Sets initial condition for a spherically symmetric radial wind blowing 
  from the origin of coordinates, see [Mig14].
  The initial condition consists of a constant-density ambient medium with 
  \f[
     \rho = \rho_a \,,\qquad  p = p_a\,,\qquad \mathbf{v} = v_{csm}\mathbf{\hat{j}} 
  \f]
  where \c v_csm is the velocity of the star with respect to the background.
  The wind is injected using the \c INTERNAL_BOUNDARY where flow quantities are
  kept constant in time and equal to 
  \f[
     r^2v_r\rho = \rho_0 V_0^2r_0^2 \,,\qquad
     \mathbf{v} = V_0\mathbf{\hat{r}} \,,\qquad
     p          = \frac{c_s^2}{\Gamma}\rho^\Gamma
  \f]
  These value are defined through the  UserDefBoundary() function when
  \c side is equal to 0. 

  Dimensions are chosen so that the spherical wind shell has radius 1,
  density 1 and velocity 1 (\f$ r_0 = V_0 = \rho_0 = 1\f$).
  
  The input parameters that control the problem dynamics are
  
  - <tt>param[CS_WIND]</tt>: sets the sound speed in the wind region;
  - <tt>param[RHO_AMB]</tt>: sets the ambient density
  - <tt>param[CS_AMB]</tt>:  sets the ambient sound speed;
  - <tt>param[V_CSM]</tt>:   sets the velocity of the star with 
                             respect to the background;

  Configurations #01-05 and #07-08 work in 2D cylindrical axisymmetric
  coordinates while conf. #06 is 3D Cartesian.

  <CENTER>
  Conf.|PHYS|GEOMETRY   |DIM|REC        |TIME_STEP  |Solver
  -----|----|-----------|---|-----------|-----------|-------
   #01 |HD  |CYLINDRICAL|3  |LINEAR     |RK2        |hll
   #03 |HD  |CYLINDRICAL|3  |LINEAR     |RK2        |hllc
   #05 |HD  |CYLINDRICAL|3  |LINEAR     |RK2        |hll
   #06 |HD  |CARTESIAN  |3  |LINEAR     |RK2        |roe
   #07 |HD  |CYLINDRICAL|3  |PARABOLIC  |RK3        |roe
  <\CENTER>
  
  \author A. Mignone (andrea.mignone@unito.it)
  \date   Apr 25, 2026

  \b References
     - [Mig14]: "High-order conservative reconstruction schemes for finite
       volume methods in cylindrical and spherical coordinates",
       Mignone, JCP (2014) 270, 784 
*/
////////////////////////////////////////////////////////////////////////
#include "gpluto.hpp"

////////////////////////////////////////////////////////////////////////
/*!
*/
////////////////////////////////////////////////////////////////////////
void Init (double *us, double x1, double x2, double x3, RunConfig &run_config)
{
  double rho, cs;
  double R;
  double *param = run_config.inputParam;

  rho = param[RHO_AMB];
  cs  = param[CS_AMB];
  us[RHO] = rho;
  us[PRS] = cs*cs*rho/run_config.gamma;

  #if GEOMETRY == CARTESIAN
  us[VX1] = 0.0;         
  us[VX2] = 0.0;         
  us[VX3] = -param[V_CSM];       

  R = sqrt(x1*x1 + x2*x2 + x3*x3);
  #elif GEOMETRY == CYLINDRICAL
  us[VX1] = 0.0;
  us[VX2] = 0.0;
  us[VX3] = -param[V_CSM];

  R = sqrt(x1*x1 + x3*x3);
  #elif GEOMETRY == SPHERICAL
  us[VX1] =  param[V_CSM]*cos(x2);
  us[VX2] = -param[V_CSM]*sin(x2);
  us[VX3] = 0.0;         

  R = x1;
  #endif

  #if NTRACER > 0
  us[TRC] = (R <= 1.0 ? 1.0:0.0);
  #endif

  #if ENTROPY_SWITCH 
  us[ENTR] = Entropy(us, run_config);
  #endif

  run_config.small_prs = 1.e-5;
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
/*!
*/
////////////////////////////////////////////////////////////////////////
{
}

////////////////////////////////////////////////////////////////////////
/*!
  Sets inflow boundary condition at the top boundary (side == X2_END)
  and the stellar wind region when side == 0.
*/
////////////////////////////////////////////////////////////////////////
void UserDefBoundary (const Data *data, RBox * box, int side, Grid *grid) 
{
  int   i, j, k;
  double  r, r0, cs;
  double  Vwind  = 1.0, rho, vr;
  double *x1 = grid->xC[IDIR];
  #if GEOMETRY == CARTESIAN
  double *x2 = grid->xC[JDIR];
  #endif
  double *x3 = grid->xC[KDIR];
  double *param = data->run_config->inputParam;
  double gamma  = data->run_config->gamma;
  #if ENTROPY_SWITCH 
  int nv;
  RunConfig &run_config = *(data->run_config);
  #endif

  if (side == 0){
    r0 = 1.0;
    cs = param[CS_WIND];
    #pragma acc parallel loop collapse(3) present(data, box)
    BOX_LOOP(*box,i,j,k){ 
      #if GEOMETRY == CARTESIAN
      r  = sqrt(x1[i]*x1[i] + x2[j]*x2[j] + x3[k]*x3[k]);
      if (r <= r0){  
        vr    = tanh(r/r0/0.1)*Vwind;
        rho   = Vwind*r0*r0/(vr*r*r);
        data->Vc(i,j,k,RHO) = rho;
        data->Vc(i,j,k,VX1) = Vwind*x1[i]/r;
        data->Vc(i,j,k,VX2) = Vwind*x2[j]/r;
        data->Vc(i,j,k,VX3) = Vwind*x3[k]/r;
        data->Vc(i,j,k,PRS) = cs*cs/gamma*pow(rho,gamma);
        data->flag(i,j,k)  |= FLAG_INTERNAL_BOUNDARY;

        #if ENTROPY_SWITCH  // Need to recompute entropy ? 
        double v[NVAR];
        NVAR_LOOP(nv) v[nv] = data->Vc(i,j,k,nv);
        data->Vc(i,j,k,ENTR) = Entropy(v,run_config);
        #endif 
      }    
      #elif GEOMETRY == CYLINDRICAL
      r  = sqrt(x1[i]*x1[i] + x3[k]*x3[k]);
      if (r <= r0){  
        vr    = tanh(r/r0/0.1)*Vwind;
        rho   = Vwind*r0*r0/(vr*r*r);
        data->Vc(i,j,k,RHO) = rho;
        data->Vc(i,j,k,VX1) = Vwind*x1[i]/r;
        data->Vc(i,j,k,VX2) = 0.0;
        data->Vc(i,j,k,VX3) = Vwind*x3[k]/r;
        data->Vc(i,j,k,PRS) = cs*cs/gamma*pow(rho,gamma);
        data->flag(i,j,k)  |= FLAG_INTERNAL_BOUNDARY;

        #if ENTROPY_SWITCH  // Need to recompute entropy ? 
        double v[NVAR];
        NVAR_LOOP(nv) v[nv] = data->Vc(i,j,k,nv);
        data->Vc(i,j,k,ENTR) = Entropy(v,run_config);
        #endif 
     }    
     #endif

    }
    
    // Convert to conservative
    PrimToCons3D (data, 0, grid, box);
  }

  if (side == X2_END){  // -- X2_END boundary -- 
    cs  = param[CS_AMB];
    rho = param[RHO_AMB];
    #pragma acc parallel loop collapse(3) present(data, box)
    BOX_LOOP(*box,i,j,k){
      #if GEOMETRY == CYLINDRICAL
      data->Vc(i,j,k,VX1) = 0.0;
      data->Vc(i,j,k,VX2) = -param[V_CSM];  
      data->Vc(i,j,k,VX3) = 0.0;  
      data->Vc(i,j,k,RHO) =  rho;       
      data->Vc(i,j,k,PRS) =  cs*cs*rho/gamma; 
      #endif
      #if ENTROPY_SWITCH // Recompute entropy
      double v[NVAR];
      v[RHO] = data->Vc(i,j,k,RHO);
      v[PRS] = data->Vc(i,j,k,PRS);
      v[ENTR] = Entropy(v,run_config);
      #endif
    }
  }

  if (side == X3_BEG){  // -- X3_BEG boundary -- 
  }
 
  if (side == X3_END){  // -- X3_END boundary -- 
    cs  = param[CS_AMB];
    rho = param[RHO_AMB];
    #pragma acc parallel loop collapse(3) present(data, box)
    BOX_LOOP(*box,i,j,k){
      #if GEOMETRY == CARTESIAN
      data->Vc(i,j,k,VX1) = 0.0;
      data->Vc(i,j,k,VX2) = 0.0;
      data->Vc(i,j,k,VX3) = -param[V_CSM];  
      data->Vc(i,j,k,RHO) =  rho;       
      data->Vc(i,j,k,PRS) =  cs*cs*rho/gamma; 
      #elif GEOMETRY == CYLINDRICAL
      data->Vc(i,j,k,VX1) = 0.0;
      data->Vc(i,j,k,VX2) = 0.0;  
      data->Vc(i,j,k,VX3) = -param[V_CSM];  
      data->Vc(i,j,k,RHO) =  rho;       
      data->Vc(i,j,k,PRS) =  cs*cs*rho/gamma; 
      #endif

      #if ENTROPY_SWITCH // Recompute entropy
      double v[NVAR];
      v[RHO] = data->Vc(i,j,k,RHO);
      v[PRS] = data->Vc(i,j,k,PRS);
      v[ENTR] = Entropy(v,run_config);
      #endif
    }
  }
}

