////////////////////////////////////////////////////////////////////////
/*!
  \file
  \brief Store emf-related components obtained during 1D sweeps.

   Collections of inline function used to store emf-related components
   during Riemann solvers.

   Here \c emf->e (and \c emf->e_dff) are matrices containing the 
   components of electric field or the velocity (depending on the scheme).
   During the x1-sweep, for instance, <tt>{n,t,b} = {x,y,z}</tt> we have

       when dir = DIRn = IDIR:  (Bxn, Bxt, Bxb) = (Bx, By, Bz)
       emf->e[DIRb][DIRn] = emf->e[KDIR][IDIR] = -flux(BXt)  = Ez
       emf->e[DIRt][DIRn] = emf->e[JDIR][IDIR] =  flux(BXb)  = Ey 

       when dir = DIRn = JDIR:  (Bxn, Bxt, Bxb) = (By, Bz, Bx)
       emf->e[DIRb][DIRn] = emf->e[IDIR][JDIR] = -flux(BXt)  = Ex
       emf->e[DIRt][DIRn] = emf->e[KDIR][JDIR] =  flux(BXb)  = Ez 

       when dir = DIRn = KDIR:  (Bxn, Bxt, Bxb) = (Bz, Bx, By)
       emf->e[DIRb][DIRn] = emf->e[JDIR][KDIR] = -flux(BXt)  = Ey
       emf->e[DIRt][DIRn] = emf->e[IDIR][KDIR] =  flux(BXb)  = Ex 

  or, if we store velocities:

       when dir = DIRn = IDIR:  
       emf->e[DIRb][DIRn] = +vy
       emf->e[DIRt][DIRn] = -vz

       when dir = DIRn = JDIR:  
       emf->e[DIRb][DIRn] = +vz
       emf->e[DIRt][DIRn] = -vx

       when dir = DIRn = KDIR:  
       emf->e[DIRb][DIRn] = +vx
       emf->e[DIRt][DIRn] = -vy

   Signs are taken so that the final composition formula for a UCT-like 
   scheme are of the type

       Ez = (aL*vL*BL + aR*vR*BR) +- (dR*BR - dL*BL)


  \authors M. Rossazza  \n
           A. Mignone (andrea.mignone@unito.it)
  \date    May 15, 2025
*/
////////////////////////////////////////////////////////////////////////

#if (CT_EMF_AVERAGE == ARITHMETIC) || (CT_EMF_AVERAGE == CT_CONTACT)

// Function prototype
/**
 * @brief  Store full electric field from Riemann solver flux.
 * @tparam dir         Direction index.
 * @param  d           Pointer to Data structure.
 * @param  flux        Reference to Ary2D flux array.
 * @param  i           Index in the array.
 * @param  ti, tj, tk  Indices in 3D space.
 */
template<int dir>
inline void CT_Store_Arith(Data *, Ary2D &, int, int, int, int);

// Explicit template instantiations
/** \cond SuppressExplicitInstantiations */
template void CT_Store_Arith<IDIR>(Data *, Ary2D &, int, int, int, int);
template void CT_Store_Arith<JDIR>(Data *, Ary2D &, int, int, int, int);
template void CT_Store_Arith<KDIR>(Data *, Ary2D &, int, int, int, int);
/** \endcond*/
////////////////////////////////////////////////////////////////////////
#pragma acc routine seq
template<int dir>
inline void CT_Store_Arith(Data *data, Ary2D &flux, int i, int ti, int tj, int tk)
/*! Arithmetic average
*/
////////////////////////////////////////////////////////////////////////
{
  EMF  *emf = data->emf;
  
  GEN_INDEX_CYCLE (dir, BX1, BXn, BXt, BXb);
  GEN_INDEX_CYCLE (dir, IDIR, DIRn, DIRt, DIRb);

  emf->e[DIRb][DIRn](ti,tj,tk) = -flux(i,BXt); 
  emf->e[DIRt][DIRn](ti,tj,tk) =  flux(i,BXb);
  #if CT_EMF_AVERAGE == CT_CONTACT
  int s;
  if      (flux(i,RHO) >  eps_CT_CONTACT) s =  1;
  else if (flux(i,RHO) < -eps_CT_CONTACT) s = -1;
  else s = 0;
  emf->sv[DIRn](ti,tj,tk) = s;
  #endif
}
#endif

#if CT_EMF_AVERAGE == CT_FLUX
// Function prototype
template<int dir>
inline void CT_Store_Flux(Data *, double *, double *, Ary2D &, int, int, int, int);

// Explicit template instantiations
/** \cond SuppressExplicitInstantiations */
template void CT_Store_Flux<IDIR>(Data *, double *, double *, Ary2D &, 
                                  int, int, int, int);
template void CT_Store_Flux<JDIR>(Data *, double *, double *, Ary2D &, 
                                  int, int, int, int);
template void CT_Store_Flux<KDIR>(Data *, double *, double *, Ary2D &,  
                                  int, int, int, int);
/** \endcond */

////////////////////////////////////////////////////////////////////////
/*! Separately store point-flux (electric field) and diffusive 
  contributions.                                                    */
////////////////////////////////////////////////////////////////////////
#pragma acc routine seq
template<int dir>
inline void CT_Store_Flux(Data *data, double *fl, double *fr, Ary2D &flux, 
                          int i, int ti, int tj, int tk)
{
  EMF  *emf = data->emf;
  
  GEN_INDEX_CYCLE (dir, BX1, BXn, BXt, BXb);
  GEN_INDEX_CYCLE (dir, IDIR, DIRn, DIRt, DIRb);

  emf->e[DIRb][DIRn](ti,tj,tk) =  -0.5*(fl[BXt] + fr[BXt]);
  emf->e[DIRt][DIRn](ti,tj,tk) =   0.5*(fl[BXb] + fr[BXb]);
  emf->e_dff[DIRb][DIRn](ti,tj,tk) = -(flux(i,BXt) - 0.5*(fl[BXt] + fr[BXt]));
  emf->e_dff[DIRt][DIRn](ti,tj,tk) =   flux(i,BXb) - 0.5*(fl[BXb] + fr[BXb]);
}
#endif

// #if CT_EMF_AVERAGE == UCT_HLL
#if PHYSICS != GRMHD
// Function prototype
template<int dir>
inline void CT_Store_UCTHll(Data *, double *, double *, double, double, 
                            int, int, int);

// Explicit template instantiations
/** \cond SuppressExplicitInstantiations */
template void CT_Store_UCTHll<IDIR>(Data *, double *, double *, 
                                    double, double, int, int, int);
template void CT_Store_UCTHll<JDIR>(Data *, double *, double *, 
                                    double, double, int, int, int);
template void CT_Store_UCTHll<KDIR>(Data *, double *, double *,  
                                    double, double, int, int, int);
/** \endcond */

////////////////////////////////////////////////////////////////////////
/*! Store transverse velocities, weight and diffusion coefficients obtained 
  from Riemann solver.
*/
////////////////////////////////////////////////////////////////////////
#pragma acc routine seq
template<int dir>
inline void CT_Store_UCTHll(Data *data, double *vl, double *vr, 
                            double Sl, double Sr, 
                            int ti, int tj, int tk)
{
  EMF  *emf = data->emf;
  
  GEN_INDEX_CYCLE (dir, VX1, VXn, VXt, VXb);
  GEN_INDEX_CYCLE (dir, IDIR, DIRn, DIRt, DIRb);

  double aR_1 = MAX(0.0, Sr);
  double aL_1 = MIN(0.0, Sl);

  double scrh = 1.0/(aR_1 - aL_1);

  emf->e[DIRb][DIRn](ti,tj,tk) =  (aR_1*vl[VXt] - aL_1*vr[VXt])*scrh;
  emf->e[DIRt][DIRn](ti,tj,tk) = -(aR_1*vl[VXb] - aL_1*vr[VXb])*scrh;

  emf->aL[DIRn](ti,tj,tk) =  aR_1*scrh;       // Right weight
  emf->aR[DIRn](ti,tj,tk) = -aL_1*scrh;       // Left weight
  emf->dR[DIRn](ti,tj,tk) = -aL_1*aR_1*scrh;  // Right diffusion coeff
  emf->dL[DIRn](ti,tj,tk) = -aL_1*aR_1*scrh;  // Left  diffusion coeff
}

#else // PHYSICS == GRMHD

// Function prototype that takes also the local metric
template<int dir>
inline void CT_Store_UCTHll(Data *, double *, double *, double, double, 
                            int, int, int, Metric *);

// Explicit template instantiations
/** \cond SuppressExplicitInstantiations */
template void CT_Store_UCTHll<IDIR>(Data *, double *, double *, 
                                    double, double, int, int, int, Metric *);
template void CT_Store_UCTHll<JDIR>(Data *, double *, double *, 
                                    double, double, int, int, int, Metric *);
template void CT_Store_UCTHll<KDIR>(Data *, double *, double *,  
                                    double, double, int, int, int, Metric *);
/** \endcond */

/* ********************************************************************* */
/*! Store transverse velocities, weight and diffusion coefficients obtained 
 *  from Riemann solver.
 *
 *********************************************************************** */
#pragma acc routine seq
template<int dir>
inline void CT_Store_UCTHll(Data *d, double *vl, double *vr, 
                            double Sl, double Sr, 
                            int ti, int tj, int tk, Metric *G)
{
  EMF  *emf = d->emf;
  
  GEN_INDEX_CYCLE (dir, VX1, VXn, VXt, VXb);
  GEN_INDEX_CYCLE (dir, IDIR, DIRn, DIRt, DIRb);

  double aR_1 = MAX(0.0, Sr);
  double aL_1 = MIN(0.0, Sl);

  double scrh = 1.0/(aR_1 - aL_1);

  // Use the transport velocity
  vl[VXt] = vl[VXt]*G->lapse - G->shift[VXt-VX1];   
  vl[VXb] = vl[VXb]*G->lapse - G->shift[VXb-VX1];   
  vr[VXt] = vr[VXt]*G->lapse - G->shift[VXt-VX1];   
  vr[VXb] = vr[VXb]*G->lapse - G->shift[VXb-VX1];   

  emf->e[DIRb][DIRn](ti,tj,tk) =  (aR_1*vl[VXt] - aL_1*vr[VXt])*scrh/1.;
  emf->e[DIRt][DIRn](ti,tj,tk) = -(aR_1*vl[VXb] - aL_1*vr[VXb])*scrh/1.;

  emf->aL[DIRn](ti,tj,tk) =  aR_1*scrh;       // Right weight
  emf->aR[DIRn](ti,tj,tk) = -aL_1*scrh;       // Left weight
  emf->dR[DIRn](ti,tj,tk) = -aL_1*aR_1*scrh;  // Right diffusion coeff
  emf->dL[DIRn](ti,tj,tk) = -aL_1*aR_1*scrh;  // Left  diffusion coeff
}
#endif
// #endif

#if CT_EMF_AVERAGE == UCT_HLLD
// Function prototype
template<int dir>
inline void CT_Store_UCTHlld(Data *, double *, double *, 
                                     double *, double *, 
                                     double, double, int, int, int);

// Explicit template instantiations
/** \cond SuppressExplicitInstantiations */
template void CT_Store_UCTHlld<IDIR>(Data *, double *, double *, double *, double *, 
                                     double, double, int, int, int);
template void CT_Store_UCTHlld<JDIR>(Data *, double *, double *, double *, double *, 
                                     double, double, int, int, int);
template void CT_Store_UCTHlld<KDIR>(Data *, double *, double *, double *, double *, 
                                     double, double, int, int, int);
/** \endcond */

////////////////////////////////////////////////////////////////////////
#pragma acc routine seq
template<int dir>
inline void CT_Store_UCTHlld(Data *data, double *ul, double *ur, 
                                      double *vl, double *vr,
                                      double Sl, double Sr, 
                                      int ti, int tj, int tk)
/*! 
  Store transverse velocities, weight and diffusion coefficients obtained 
  from Riemann solver.
*/
////////////////////////////////////////////////////////////////////////
{
  EMF  *emf = data->emf;
  
  GEN_INDEX_CYCLE (dir, VX1, VXn, VXt, VXb);
  GEN_INDEX_CYCLE (dir, MX1, MXn, MXt, MXb);
  GEN_INDEX_CYCLE (dir, BX1, BXn, BXt, BXb);
  GEN_INDEX_CYCLE (dir, BX1s, BXns, BXts, BXbs);
  GEN_INDEX_CYCLE (dir, IDIR, DIRn, DIRt, DIRb);

  int    switch_to_hll = 0;
  double scrh;
  double aR_1 = MAX(0.0, Sr);
  double aL_1 = MIN(0.0, Sl);

  #ifdef HIGH_ORDER
  double Bn  = data->Vps[BXns](ti,tj,tk);
  #else
  double Bn  = data->Vs[BXns](ti,tj,tk);
  #endif

  double dul = Sl - vl[VXn];
  double dur = Sr - vr[VXn];

  #if EOS == IDEAL
  double pL = vl[PRS] + 0.5*(vl[BX1]*vl[BX1] + vl[BX2]*vl[BX2] + vl[BX3]*vl[BX3]);
  double pR = vr[PRS] + 0.5*(vr[BX1]*vr[BX1] + vr[BX2]*vr[BX2] + vr[BX3]*vr[BX3]);
  scrh  = 1.0/(dur*ur[RHO] - dul*ul[RHO]);
  double Sc = (dur*ur[MXn] - dul*ul[MXn] - pR + pL)*scrh;
  double usl_RHO = ul[RHO]*dul/(Sl - Sc);
  double usr_RHO = ur[RHO]*dur/(Sr - Sc);

  if (usl_RHO < 0.0 || usr_RHO < 0.0) switch_to_hll = 1;
 
  double sqrL = sqrt(usl_RHO);
  double sqrR = sqrt(usr_RHO);

  double Sal = Sc - fabs(Bn)/sqrL;
  double Sar = Sc + fabs(Bn)/sqrR;

  double chiL  = (vl[VXn] - Sc)*(Sl - Sc)/(Sal + Sl - 2.0*Sc);
  double chiR  = (vr[VXn] - Sc)*(Sr - Sc)/(Sar + Sr - 2.0*Sc);
  #elif EOS == ISOTHERMAL
  scrh    = 1.0/(Sr - Sl);
  double rho_h   = (ur[RHO]*dur - ul[RHO]*dul)*scrh;
  //double Sc = 0.5*(Sal + Sar);
  double Sc = (Sl*ur[RHO]*dur - Sr*ul[RHO]*dul*scrh)/rho_h;

  if (rho_h < 0.0) switch_to_hll = 1;

  double sqrho_h = sqrt(rho_h);
  double Sal = Sc - fabs(Bn)/sqrho_h;
  double Sar = Sc + fabs(Bn)/sqrho_h;
  double chiL  = (vl[VXn] - Sc)*(Sl - Sc)/(Sal + Sl - 2.0*Sc);
  double chiR  = (vr[VXn] - Sc)*(Sr - Sc)/(Sar + Sr - 2.0*Sc);
  #endif

  double nuL  = (Sal + Sl)/(fabs(Sal) + fabs(Sl));
  double nuR  = (Sar + Sr)/(fabs(Sar) + fabs(Sr));
  double nuLR = (Sal + Sar)/(fabs(Sal) + fabs(Sar));

  if ( fabs(Sar - Sal) > 1.e-9*fabs(Sr-Sl)){
    emf->dL[dir](ti,tj,tk) =   0.5*(chiL*nuL - chiL*nuLR)
                             + 0.5*(fabs(Sal) - nuLR*Sal);

    emf->dR[dir](ti,tj,tk) =   0.5*(chiR*nuR - chiR*nuLR)
                             + 0.5*(fabs(Sar) - nuLR*Sar);
    emf->aL[dir](ti,tj,tk) = 0.5*(1.0 + nuLR);
    emf->aR[dir](ti,tj,tk) = 0.5*(1.0 - nuLR);
  }else{   // HLLC, degenerate limit Bx -> 0 
    emf->dL[dir](ti,tj,tk) = 0.5*chiL*nuL + 0.5*fabs(Sal);
    emf->dR[dir](ti,tj,tk) = 0.5*chiR*nuR + 0.5*fabs(Sar);

    emf->aL[dir](ti,tj,tk) = 0.5;
    emf->aR[dir](ti,tj,tk) = 0.5;
  }

  scrh = 1.0/(aR_1 - aL_1);

  // HLL Diffusion coefficients 

  if (switch_to_hll){
    emf->aL[DIRn](ti,tj,tk) =  aR_1*scrh;       // Right weight
    emf->aR[DIRn](ti,tj,tk) = -aL_1*scrh;       // Left weight
    emf->dR[DIRn](ti,tj,tk) = -aL_1*aR_1*scrh;  // Right diffusion coeff
    emf->dL[DIRn](ti,tj,tk) = -aL_1*aR_1*scrh;  // Left  diffusion coeff
  }

  emf->e[DIRb][DIRn](ti,tj,tk) =  (aR_1*vl[VXt] - aL_1*vr[VXt])*scrh;
  emf->e[DIRt][DIRn](ti,tj,tk) = -(aR_1*vl[VXb] - aL_1*vr[VXb])*scrh;

  emf->e_dff[DIRb][DIRn](ti,tj,tk) = 0.0;
  emf->e_dff[DIRt][DIRn](ti,tj,tk) = 0.0;
}
#endif

#if CT_EMF_AVERAGE == UCT_GFORCE
#if PHYSICS != GRMHD
// Function prototype
template<int dir>
inline void CT_Store_GForce(Data *, double *, double *, 
                                    double *, double *, 
                                    double, double, 
                                    int, int, int);

// Explicit template instantiations
/** \cond SuppressExplicitInstantiations */
template void CT_Store_GForce<IDIR>(Data *, double *, double *, double *, double *, 
                                    double, double, int, int, int);
template void CT_Store_GForce<JDIR>(Data *, double *, double *, double *, double *, 
                                    double, double, int, int, int);
template void CT_Store_GForce<KDIR>(Data *, double *, double *, double *, double *, 
                                    double, double, int, int, int);
/** \endcond */

////////////////////////////////////////////////////////////////////////
#pragma acc routine seq
template<int dir>
inline void CT_Store_GForce(Data *data, double *vl, double *vr, 
                                     double *uLW, double *vLW, 
                                     double om, double lambda,
                                     int ti, int tj, int tk)
/*!
  Store transverse velocities and diffusion coefficients obtained 
  from Riemann solver.
*/
////////////////////////////////////////////////////////////////////////
{
  EMF  *emf = data->emf;
  
  GEN_INDEX_CYCLE (dir, VX1, VXn, VXt, VXb);
  GEN_INDEX_CYCLE (dir, MX1, MXn, MXt, MXb);
  GEN_INDEX_CYCLE (dir, IDIR, DIRn, DIRt, DIRb);

  double dtdx = 1.0/lambda;
  // double vxn = 0.5*(vl[VXn] + vr[VXn]);
  double vxt = 0.5*(vl[VXt] + vr[VXt]);
  double vxb = 0.5*(vl[VXb] + vr[VXb]);

  double vxn_lw = vLW[VXn];
  double vxt_lw = vLW[VXt];
  double vxb_lw = vLW[VXb];
  double fn = dtdx*vxn_lw;
 
  emf->e[DIRb][DIRn](ti,tj,tk) = -(-om*(vxt_lw - 0.5*fn*(vr[VXt] - vl[VXt])) 
                                   - (1.0 - om)*vxt);
  emf->e[DIRt][DIRn](ti,tj,tk) =   -om*(vxb_lw - 0.5*fn*(vr[VXb] - vl[VXb])) 
                                   - (1.0 - om)*vxb;

  //emf->dL[DIRn](ti,tj,tk) = 0.5*(om*dtdx*vl[VXn]*vxn_lw + (1.0 - om)*lambda);
  //emf->dR[DIRn](ti,tj,tk) = 0.5*(om*dtdx*vr[VXn]*vxn_lw + (1.0 - om)*lambda);

  emf->dL[DIRn](ti,tj,tk) = 0.5*(om*fn*vl[VXn] + (1.0 - om)*lambda);
  emf->dR[DIRn](ti,tj,tk) = 0.5*(om*fn*vr[VXn] + (1.0 - om)*lambda);

}
#else // PHYSICS == GRMHD
      
// Function prototype that takes also the local metric
template<int dir>
inline void CT_Store_GForce(Data *, double *, double *, 
                                    double *, double *, 
                                    double, double, 
                                    int, int, int, Metric *);

// Explicit template instantiations
/** \cond SuppressExplicitInstantiations */
template void CT_Store_GForce<IDIR>(Data *, double *, double *, double *, double *, 
                                    double, double, int, int, int, Metric *);
template void CT_Store_GForce<JDIR>(Data *, double *, double *, double *, double *, 
                                    double, double, int, int, int, Metric *);
template void CT_Store_GForce<KDIR>(Data *, double *, double *, double *, double *, 
                                    double, double, int, int, int, Metric *);
/** \endcond */

/* ********************************************************************* */
#pragma acc routine seq
template<int dir>
inline void CT_Store_GForce(Data *d, double *vl, double *vr, 
                                     double *uLW, double *vLW, 
                                     double om, double lambda,
                                     int ti, int tj, int tk, Metric *G)
/*!
 *  Store transverse velocities and diffusion coefficients obtained 
 *  from Riemann solver.
 *
 *********************************************************************** */
{
  EMF  *emf = d->emf;
  
  GEN_INDEX_CYCLE (dir, VX1, VXn, VXt, VXb);
  GEN_INDEX_CYCLE (dir, MX1, MXn, MXt, MXb);
  GEN_INDEX_CYCLE (dir, IDIR, DIRn, DIRt, DIRb);

  double dtdx = 1.0/lambda;

  // Use the transport velocity
  for (int nv = 0; nv < 3; ++nv) {
    vl[VX1+nv]  =  vl[VX1+nv]*G->lapse - G->shift[nv];   
    vr[VX1+nv]  =  vr[VX1+nv]*G->lapse - G->shift[nv];   
    vLW[VX1+nv] = vLW[VX1+nv]*G->lapse - G->shift[nv];   
  }

  double vxn = 0.5*(vl[VXn] + vr[VXn]);
  double vxt = 0.5*(vl[VXt] + vr[VXt]);
  double vxb = 0.5*(vl[VXb] + vr[VXb]);

  double vxn_lw = vLW[VXn];
  double vxt_lw = vLW[VXt];
  double vxb_lw = vLW[VXb];
  double fn = dtdx*vxn_lw;
 
  emf->e[DIRb][DIRn](ti,tj,tk) = -(-om*(vxt_lw - 0.5*fn*(vr[VXt] - vl[VXt])) 
                                   - (1.0 - om)*vxt);
  emf->e[DIRt][DIRn](ti,tj,tk) =   -om*(vxb_lw - 0.5*fn*(vr[VXb] - vl[VXb])) 
                                   - (1.0 - om)*vxb;

  //emf->dL[DIRn](ti,tj,tk) = 0.5*(om*dtdx*vl[VXn]*vxn_lw + (1.0 - om)*lambda);
  //emf->dR[DIRn](ti,tj,tk) = 0.5*(om*dtdx*vr[VXn]*vxn_lw + (1.0 - om)*lambda);

  emf->dL[DIRn](ti,tj,tk) = 0.5*(om*fn*vl[VXn] + (1.0 - om)*lambda);
  emf->dR[DIRn](ti,tj,tk) = 0.5*(om*fn*vr[VXn] + (1.0 - om)*lambda);

}
#endif
#endif

#if CT_EMF_AVERAGE == CT_MAXWELL
// Function prototype
template<int dir>
inline void CT_Store_Maxwell(Data *, double *, double *, Ary2D &, 
                              int, int, int, int);

// Explicit template instantiations
/** \cond SuppressExplicitInstantiations */
template void CT_Store_Maxwell<IDIR>(Data *, double *, double *, Ary2D &, 
                                     int, int, int, int);
template void CT_Store_Maxwell<JDIR>(Data *, double *, double *, Ary2D &, 
                                     int, int, int, int);
template void CT_Store_Maxwell<KDIR>(Data *, double *, double *, Ary2D &, 
                                     int, int, int, int);
/** \endcond */

////////////////////////////////////////////////////////////////////////
#pragma acc routine seq
template<int dir>
inline void CT_Store_Maxwell(Data *data, double *fl, double *fr, Ary2D &flux, 
                             int i, int ti, int tj, int tk)
// 
// 
////////////////////////////////////////////////////////////////////////
{
  EMF  *emf = data->emf;
  
  GEN_INDEX_CYCLE (dir, BX1, BXn, BXt, BXb);
  GEN_INDEX_CYCLE (dir, IDIR, DIRn, DIRt, DIRb);
  #if (DIVE_CONTROL == CONSTRAINED_TRANSPORT)
  GEN_INDEX_CYCLE (dir, EX1, EXn, EXt, EXb);
  #endif

  emf->e[DIRb][DIRn](ti,tj,tk) = -0.5*(fl[BXt] + fr[BXt]);
  emf->e[DIRt][DIRn](ti,tj,tk) =  0.5*(fl[BXb] + fr[BXb]);
  #if (DIVE_CONTROL == CONSTRAINED_TRANSPORT)
  emf->B[DIRb][DIRn](ti,tj,tk) =  0.5*(fl[EXt] + fr[EXt]);
  emf->B[DIRt][DIRn](ti,tj,tk) = -0.5*(fl[EXb] + fr[EXb]);
  #endif

  #if (PHYSICS == ResRMHD) && (DIVE_CONTROL == CONSTRAINED_TRANSPORT)
  emf->Frho[DIRn](ti,tj,tk) = flux(i,RHO);
  #else
  #endif
}
#endif
