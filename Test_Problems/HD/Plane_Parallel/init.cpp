#include "gpluto.hpp"

////////////////////////////////////////////////////////////////////////
void Init (double *v, double x1, double x2, double x3, RunConfig &run_config)
/*!

<CENTER>
Conf.|PHYS|GEOMETRY   |DIM|REC        |TIME_STEP  |Solver
-----|----|-----------|---|-----------|-----------|-------
 #01 |HD  |SPHERICAL  |2  |LINEAR     |RK2        |hll
<\CENTER>

*/
////////////////////////////////////////////////////////////////////////
{
  double *inputParam = run_config.inputParam;

  v[RHO] = 1.0 + inputParam[DN_PROF]*pow(x1*sin(x2), 2.0);
  v[VX1] =  inputParam[MACH]*cos(x2);
  v[VX2] = -inputParam[MACH]*sin(x2);
  v[VX3] = 0.0;
  v[PRS] = 1.0/run_config.gamma;

  #if PHYSICS == MHD 
  v[IBR]     =  inputParam[MAG_FIELD]*cos(x2);
  v[IBTHETA] = -inputParam[MAG_FIELD]*sin(x2); 
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
  Assign user-defined boundary conditions.
  
  \param [in,out] d  pointer to the PLUTO data structure containing
  cell-centered primitive quantities (data->Vc) and 
  staggered magnetic fields (data->Vs, when used) to 
  be filled.
  \param [in] box    pointer to a RBox structure containing the lower
  and upper indices of the ghost zone-centers/nodes
  or edges at which data values should be assigned.
  \param [in] side   specifies the boundary side where ghost zones need
  to be filled. It can assume the following 
  pre-definite values: X1_BEG, X1_END,
  X2_BEG, X2_END, 
  X3_BEG, X3_END.
  The special value side == 0 is used to control
  a region inside the computational domain.
  \param [in] grid  pointer to an array of Grid structures.
*/
////////////////////////////////////////////////////////////////////////
{
  int     i, j, k, nv;
  int nx2 = grid->nint[JDIR];
  double v[256];
  double  *x1, *x2, *x3;
  RunConfig &run_config = *(data->run_config);

  x1 = grid->xC[IDIR];
  x2 = grid->xC[JDIR];
  x3 = grid->xC[KDIR];

  #ifdef PARALLEL
    #error This setup cannot be parallelized owing to the b.c.
  #endif

  if (side == X1_BEG){

    #pragma acc parallel loop collapse(3) present(data, run_config, box)
    BOX_LOOP(*box,i,j,k){
      Init (v, x1[i], x2[j], x3[k], run_config);

      #pragma acc loop seq
      NVAR_LOOP(nv) data->Vc(i,j,k,nv) = v[nv];
    }

  } else if (side == X1_END) {

    #pragma acc parallel loop collapse(3) present(data, run_config, box)
    BOX_LOOP(*box,i,j,k){
      Init (v, x1[i], x2[j], x3[k], run_config);

      #pragma acc loop seq
      NVAR_LOOP(nv) data->Vc(i,j,k,nv) = v[nv];
    }
    
  }

  if (side == X2_BEG){

    #pragma acc parallel loop collapse(3) present(data, run_config, box)
    BOX_LOOP(*box,i,j,k){
      Init (v, x1[i], x2[j], x3[k], run_config);

//      #pragma acc loop seq
//      NVAR_LOOP(nv) data->Vc(i,j,k,nv) = v[nv];
      
      data->Vc(i,j,k,RHO) =  1.0;
      data->Vc(i,j,k,VX1) = -data->Vc(i,j + nx2,k,VX2);
      data->Vc(i,j,k,VX2) =  data->Vc(i,j + nx2,k,VX1);
      data->Vc(i,j,k,PRS) =  1.0/run_config.gamma;
    }

  } else if (side == X2_END) {

    #pragma acc parallel loop collapse(3) present(data, run_config, box)
    BOX_LOOP(*box,i,j,k){
      Init (v, x1[i], x2[j], x3[k], run_config);

//      #pragma acc loop seq
//      NVAR_LOOP(nv) data->Vc(i,j,k,nv) = v[nv];

      data->Vc(i,j,k,RHO) =  1.0;
      data->Vc(i,j,k,VX1) =  data->Vc(i,j - nx2,k,VX2);
      data->Vc(i,j,k,VX2) = -data->Vc(i,j - nx2,k,VX1);
      data->Vc(i,j,k,PRS) =  1.0/run_config.gamma;
    }
    
  }
}

