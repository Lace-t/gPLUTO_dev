////////////////////////////////////////////////////////////////////////
/*!
  \file
  \brief Hydro 2D Jet problem
  <CENTER>
  Conf.|PHYS|GEOMETRY   |DIM|REC        |TIME_STEP  |Solver
  -----|----|-----------|---|-----------|-----------|-------
   #01 |HD  |CYLINDRICAL|3  |LINEAR     |RK2        |hll
   #02 |HD  |SPHERICAL  |2  |LINEAR     |RK2        |hll
   #04 |HD  |SPHERICAL  |2  |LINEAR     |RK2        |hll
  <\CENTER>

  \author A. Mignone (andrea.mignone@unito.it)

  \date   June 08, 2014

  \b References
     - Sod, G. A. (1978).
       "A Survey of Several Finite Difference Methods for Systems of
        Nonlinear Hyperbolic Conservation Laws".
        JCP (1978) 27:1-31

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
  double Ta = 50.0;  // Ambient temperature 

  v[RHO] = inputParam[ETA];
  v[VX1] = 0.0;
  v[VX2] = 0.0;
  v[VX3] = 0.0;

  #if EOS == IDEAL
  v[PRS]  = 1.0/run_config.gamma;
  #elif EOS == PVTE_LAW
  v[PRS] = Pressure(v,Ta);
  #endif

#if PHYSICS == MHD
  #if GEOMETRY == CYLINDRICAL
  v[BX1] = 0.0;
  v[BX2] = 0.0;
  v[BX3] = inputParam[BMAG];
  #elif GEOMETRY == SPHERICAL
  v[BX1] = inputParam[BMAG]/(x1*x1); // Monopole field 
  v[BX2] = 0.0;
  v[BX3] = 0.0;
  #endif
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
/*!
*/
////////////////////////////////////////////////////////////////////////
{
  int   i, j, k, nv;
  int ibeg = grid->lbeg[IDIR], iend = grid->lend[IDIR];
  int jbeg = grid->lbeg[JDIR], jend = grid->lend[JDIR];
  int kbeg = grid->lbeg[KDIR], kend = grid->lend[KDIR];

  double  *x1, *x2, *x3;
  double cs, Tj, vj[NVAR], vinj;
  RunConfig &run_config = *(data->run_config);
  double *inputParam = run_config.inputParam;
  
  x1 = grid->xC[IDIR];  // -- array pointer to x1 coordinate -- 
  x2 = grid->xC[JDIR];  // -- array pointer to x2 coordinate -- 
  x3 = grid->xC[KDIR];  // -- array pointer to x3 coordinate -- 

  // ----------------------------------------------
  // Define jet values with cylindrical coordinates 
  // ----------------------------------------------

  vj[RHO] = 1.0;
  #if EOS == IDEAL
  vj[PRS] = 1.0/run_config.gamma;         // -- Pressure-matched jet -- 
  vj[VX2] = inputParam[MACH];  // -- Sound speed is one in this case -- 
  #elif EOS == PVTE_LAW
  Tj      = inputParam[TJET];
  vj[PRS] = Pressure(vj,Tj);
  cs      = sqrt(vj[PRS]/vj[RHO]); // -- For simplicity, isothermal cs is used -- 
  vj[VX2] = inputParam[MACH]*cs;
  #endif
  
  vinj = vj[VX2];

#if GEOMETRY == CYLINDRICAL
  if (side == X3_BEG){     // -- select the boundary side -- 
    #pragma acc parallel loop collapse(3) present(data, box)
    BOX_LOOP(*box,i,j,k){   // -- Loop over boundary zones -- 
      if (fabs(x1[i]) <= 1.0){   // -- set jet values for r <= 1 -- 
        double vphi = vinj*x1[i]*inputParam[VPHI];

        data->Vc(i,j,k,RHO) = vj[RHO];
        data->Vc(i,j,k,VX1) = 0.0;
        data->Vc(i,j,k,VX2) = vphi;
        data->Vc(i,j,k,VX3) = vinj;
        data->Vc(i,j,k,PRS) = vj[PRS]; 
      }else{                // -- reflective boundary for r > 1 --
        NVAR_LOOP(nv) data->Vc(i,j,k,nv) = data->Vc(i,j,2*kbeg-k-1,nv);
        data->Vc(i,j,k,VX3) *= -1.0;
      }
    }
  }
#endif


#if GEOMETRY == SPHERICAL
  if (side == X1_BEG){     // -- select the boundary side -- 
    #pragma acc parallel loop collapse(3) present(data, box)
    BOX_LOOP(*box,i,j,k){   // -- Loop over boundary zones -- 
      if (fabs(x2[j]) <= 0.1){   // -- set jet values for r <= 1 -- 
        double vphi = vinj*x2[j]*inputParam[VPHI];

        data->Vc(i,j,k,RHO) = vj[RHO];
        data->Vc(i,j,k,VX1) = vinj;
        data->Vc(i,j,k,VX2) = 0.0;
        data->Vc(i,j,k,VX3) = vphi;
        data->Vc(i,j,k,PRS) = vj[PRS]; 
      }else{                // -- reflective boundary for r > 1 --
        NVAR_LOOP(nv) data->Vc(i,j,k,nv) = data->Vc(2*ibeg-i-1,j,k,nv);
        data->Vc(i,j,k,VX1) *= -1.0;
      }

      #if PHYSICS == MHD
      data->Vc(i,j,k,BX1) = inputParam[BMAG]/(x1[i]*x1[i]);
      data->Vc(i,j,k,BX2) = 0.0;
      data->Vc(i,j,k,BX3) = 0.0;
      #endif

    }
  }
#endif
}
