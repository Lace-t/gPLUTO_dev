////////////////////////////////////////////////////////////////////////
/*!
//to be commented by D.Crocco
  \file
  \brief Compute the rhs for parabolic terms


UpdateStage:

  <X1dir>:
    Compute viscous flux+sources in IDIR:
          data->fluidRiemannSolver(data, Dts, grid, &sweep_box);
          ParabolicFlux<IDIR>()
    Add it to RightHandSide()
 
  <X2dir>:
    Compute viscous flux+source in JDIR:
          data->fluidRiemannSolver(data, Dts, grid, &sweep_box);
          ParabolicFlux<JDIR>()
    Add it to RightHandSide()

  <X3dir>:
    Compute viscous flux+sources in KDIR:
          data->fluidRiemannSolver(data, Dts, grid, &sweep_box);
          ParabolicFlux<KDIR>()
    Add it to RightHandSide()

 
  \author Z. Ahamne
  \date   Apr 19, 2024
*/
////////////////////////////////////////////////////////////////////////
#include "gpluto.hpp"

/*void ParabolicFlux<dir> (Data *data, RBox *box, Grid *grid)
{
  #if INCLUDE_VISCOSITY
  ViscousFlux<dir> ()
  #endif

  #if INCLUDE_RESISTIVITY
  ResistiveFlux<dir> ()
  #endif

  .. 


  


}

XXXFlux<dir>()
// Compute Parabolic flux for XXX process (XXX = Viscosity, Resistivity, etc...
//   at faces in the dir direction)
{
  for (k = kbeg; k <= kend; k++){
  for (j = jbeg; j <= jend; j++){
    for (i = ibeg; i <= iend; i++){
      int ti, tj, tk;
      ti = 
      tj = 
      tk = 

      sweep->flux_par(i,nv) = 

// Remember to compute dt estimate (or related coefficients)
    }
  }}
}


int ParabolicRHS(Data *data, RBox *box, Grid *grid)
{
  int i,j,k;
  double *dx1 = grid->dx[IDIR];
  double *dx2 = grid->dx[JDIR];
  double *dx3 = grid->dx[KDIR];



  BOX_LOOP(*box, i,j,k){

    #if AMBICYLINDRICAL_DIFFUSION

    dBx = (data->Vc(i+1,j,k,BX1) - data->Vc(i-1,j,k,BX1))/dx;  

    Jx = dBy/dz - dBz/dx... 

    ... compute here ...

    Fx1(i,j,k,BX1) = ... 
    Fx1(i,j,k,BX2) = ... 
    Fx1(i,j,k,BX3) = ... 

    Fx2(i,j,k,BX1) = ... 
    Fx2(i,j,k,BX2) = ... 
    Fx3(i,j,k,BX3) = ... 

    Fx3(i,j,k,BX1) = ... 
    Fx3(i,j,k,BX2) = ... 
    Fx3(i,j,k,BX3) = ... 

    #endif

  }

  BOX_LOOP(*box, i,j,k){
    dU += dt/dx*(Fx1(i+1,j,k,nv) - Fx1(i,j,k,nv));
  }



}
*/