#include "gpluto.hpp"

////////////////////////////////////////////////////////////////////////
/*!
*/
////////////////////////////////////////////////////////////////////////
void Init (double *us, double x1, double x2, double x3, RunConfig &run_config)
{
   run_config.gamma = 1.4;


   us[VX1] = (x1 > 0.6 && x2 < 0.2 ? 0.0:3.0);
   us[VX2] = 0.0;
   us[VX3] = 0.0;
   us[TRC] = 0.0;

   us[RHO] = 1.4;
   us[PRS] = 1.0;
us[VX1] = (x1 < 1.0 ? 3.0:0.0);
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
*/
////////////////////////////////////////////////////////////////////////
{
}

////////////////////////////////////////////////////////////////////////
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
void UserDefBoundary (const Data *data, RBox *box, int side, Grid *grid) 
{
  int   i, j, k;
  double *x = grid->xC[IDIR];
  double *y = grid->xC[JDIR];

  if (side == 0){
    #pragma acc parallel loop collapse(3) present(data, box)
    BOX_LOOP(*box,i,j,k){

double y1 = 0.5 + 1.2*(x[i]-1.5);
double y2 = 0.5 - 0.8*(x[i]-1.5);

if (y[j] < y1 && y[j] > y2 && x[i] < 1.7) {
        data->flag(i,j,k) |= FLAG_INTERNAL_BOUNDARY;   
        data->Vc(i,j,k,RHO) = 10.0;
        data->Vc(i,j,k,VX1) = 0.0;
        data->Vc(i,j,k,VX2) = 0.0;
        data->Vc(i,j,k,VX3) = 0.0;
        data->Vc(i,j,k,PRS) = 1.0;
}



    //   if (y[j] <= 0.2 && x[i] >= 0.6) {
    //     data->flag(i,j,k) |= FLAG_INTERNAL_BOUNDARY;   
    //     data->Vc(i,j,k,RHO) = 1.4;
    //     data->Vc(i,j,k,VX1) = 0.0;
    //     data->Vc(i,j,k,VX2) = 0.0;
    //     data->Vc(i,j,k,VX3) = 0.0;
    //     data->Vc(i,j,k,PRS) = 1.0;
    //   }
    }
  }
  
  if (side == X1_BEG){
    #pragma acc parallel loop collapse(3) present(data, box)
    BOX_LOOP(*box,i,j,k){
      data->Vc(i,j,k,RHO) = 1.4;
      data->Vc(i,j,k,VX1) = 3.0;
      data->Vc(i,j,k,VX2) = 0.0;
      data->Vc(i,j,k,VX3) = 0.0;
      data->Vc(i,j,k,PRS) = 1.0;
    }
  } 
}
