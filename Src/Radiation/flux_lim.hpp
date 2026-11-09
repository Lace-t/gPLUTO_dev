// /////////////////////////////////////////////////////////////////////
/*!
  \file 
  \brief  Radiation flux limiting

  Function used to limit the radiation flux in case it exceeds the radiation
  energy density.
  
  \authors J. D. Melon Fuksman (fuksman@mpia.de)
  \date    Dec 12, 2019

*/
// /////////////////////////////////////////////////////////////////////

#pragma omp declare target
#pragma acc routine seq
////////////////////////////////////////////////////////////////////////
/*!
  Fix the condition |F|<E on the L/R states at cell interfaces.

  \param[in,out] u  Vector of conserved or primitive fields.
*/
////////////////////////////////////////////////////////////////////////
inline void LimitRadFlux (double *u)
{
  double f_abs;
  double n1, n2, n3;

  f_abs = sqrt (u[FR1]*u[FR1] + u[FR2]*u[FR2] + u[FR3]*u[FR3]) ;

  if( f_abs > u[ENR] ) {
    if( f_abs > 1e-50 ){
      n1 = u[FR1]/f_abs ;
      n2 = u[FR2]/f_abs ;
      n3 = u[FR3]/f_abs ;
      u[FR1] = u[ENR]*n1 ;
      u[FR2] = u[ENR]*n2 ;
      u[FR3] = u[ENR]*n3 ;
    } else {
      n1 = u[FR1]/1e-50 ;
      n2 = u[FR2]/1e-50 ;
      n3 = u[FR3]/1e-50 ;
      u[FR1] = u[ENR]*n1 ;
      u[FR2] = u[ENR]*n2 ;
      u[FR3] = u[ENR]*n3 ;
    }
  }

  return;
}
#pragma omp end declare target

#pragma omp declare target
#pragma acc routine seq
////////////////////////////////////////////////////////////////////////
/*!
  Check whether the total reconstructed radiation flux is > Er.
  If so, flatten distribution.

  \param[in]     v   Cell-centered state.
  \param[in,out] vL  Left reconstructed interface states.
  \param[in,out] vR  Right reconstructed interface states.
  \param[in]     i   Interface index.
*/
////////////////////////////////////////////////////////////////////////
inline void RadFluxLimFlatten (Ary2D &v, Ary2D &vL, Ary2D &vR, const int i)
{
// --  Flatten distribution if |Frad| > Erad  --
  double fm, fp;
  int k;

  fm = vR(i,FR1)*vR(i,FR1)  + vR(i,FR2)*vR(i,FR2) + vR(i,FR3)*vR(i,FR3);
  fm = sqrt(fm);
  fp = vL(i+1,FR1)*vL(i+1,FR1)  + vL(i+1,FR2)*vL(i+1,FR2) + vL(i+1,FR3)*vL(i+1,FR3);
  fp = sqrt(fp);

  if (fm >= vR(i,ENR) || fp >= vL(i+1,ENR)){
    #pragma acc loop seq
    for (k = FR3; k >= ENR; --k) vR(i,k) = vL(i+1,k) = v(i,k);
  }
}
#pragma omp end declare target
