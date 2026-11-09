////////////////////////////////////////////////////////////////////////
/*! 
  \file  
  \brief Add source terms before the FARGO advection step.
  
  This function is called prior to the FARGO advection algorithm to add
  source terms.
  At present, we use it only for the energy equation in the
  ShearingBox module:
  \f[ 
     \pd{E'}{t} + w\pd{E'}{y} = (B_xB_y - \rho v_xv'_y)\pd{w}{x}
  \f]
  where \f$w = -q\Omega x\f$.
  The discretization follows the algorithm of [GS10], see
  Eq. (51) and (63) of that paper.

  \b Reference
     - [GS10] "Implementation of the shearing box approximation in Athena",
       Stone & Gardiner, ApJS (2010) 189, 142.
 
  \author A. Mignone (andrea.mignone@unito.it)
  \date   Apr 30, 2020
*/
////////////////////////////////////////////////////////////////////////
#include "gpluto.hpp"

////////////////////////////////////////////////////////////////////////
void FARGO_Source(Ary4D &UU, Ary3D (&Vs)[6], double dt, Grid *grid)
/*!
  
  \param [in,out] UU     array of conserved variables
  \param [in]     Vs     array of staggered magnetic fields
  \param [in]     dt     the current time increment
  \param [in]     grid   pointer to an array of Grid structures
*/
////////////////////////////////////////////////////////////////////////
{
#ifdef SHEARINGBOX
  #if HAVE_ENERGY || !INCLUDE_JDIR
  int    i,j,k;
  double scrh, rho, mx, my;
  double *x  = grid->x[IDIR];
  double *xr = grid->xr[IDIR];
  double *dx = grid->dx[IDIR];
  double *dz = grid->dx[KDIR];  
  #if PHYSICS == MHD
  double Bx, By;
  double w, wp, wm;
  #endif
  int ibeg    = grid->lbeg[IDIR], iend = grid->lend[IDIR];
  int jbeg    = grid->lbeg[JDIR], jend = grid->lend[JDIR];
  int kbeg    = grid->lbeg[KDIR], kend = grid->lend[KDIR];
  
  scrh = SB_Q*SB_OMEGA;
  for(k=kbeg; k <= kend; k++){
  for(j=jbeg; j <= jend; j++){
  for(i=ibeg; i <= iend; i++){
    rho = UU(i,j,k,RHO);
    mx  = UU(i,j,k,MX1);
    my  = UU(i,j,k,MX2);
    #if PHYSICS == MHD
    Bx  = UU(i,j,k,BX1);
    By  = UU(i,j,k,BX2);
    #endif

    // ------------------------------------------------------
    // In the axisymmetric case, By is not a staggered
    // component and not treated during FARGO_ShiftSolution().
    // Contributions from d(w*Bx)/dx and d(w*Bz)/dz are
    // added here.
    // ------------------------------------------------------
  
    #if (PHYSICS == MHD) && !INCLUDE_JDIR   // Axisymmetric ShearingBox 
    wp = -scrh*xr[i];
    wm = -scrh*xr[i-1];
    w  = -scrh*x[i];
    
    UU[k][j][i][BX2] += dt*( (wp*Vs[BX1s](i,j,k) - wm*Vs[BX1s](i-1,j,k))/dx[i]
                            + w*(Vs[BX1s](i,j,k) - Vs[BX1s](i,j,k-1))/dz[k] );
    #endif

    #if HAVE_ENERGY && (PHYSICS == MHD)
    UU(i,j,k,ENG) += - dt*scrh*Bx*(By - 0.5*dt*Bx*scrh)
                        + dt*scrh*mx*my/rho;
    #endif                        
    #if HAVE_ENERGY && (PHYSICS == HD)
    UU(i,j,k,ENG) += dt*scrh*mx*my/rho;
    #endif                        
    
  }}}

  #endif
#endif // SHEARINGBOX 
}
