////////////////////////////////////////////////////////////////////////
/*!
  \file
  \brief   This function returns punctual values of
           both centered and staggered variables.

  Reckon variables' punctual values in their respective domains, if high order
  is enabled, this is done by high order quadrature techniques.
  High order quadrature uses Laplacian operators
  \f[
     \begin{array}{l}
      U_pc = \langle U_c \rangle - \frac{h^2}{24} \Delta U_c \\[0.5em]
      V_ps = \langle B_s \rangle - \frac{h^2}{24} \Delta^{\bot} B_s
     \end{array}
  \f]
  where \f$ U_c \f$ and \f$ U_s \f$ are 3D arrays containing, respectively,
  low order zone-centered and staggered (conservative) variables, \f$ \Delta \f$
  and \f$ \Delta^{\bot} \f$ are the Laplacian and transverse laplacian operators.
  This function performs the convolution of both centered and staggered variables,
  then converts conservative var. into primitive var.

- Set ho_lim (D)
- Exchange b.c. on ho_lim (D+B)
- Add buffer (D+B)


  \authors A. Mignone (andrea.mignone@unito.it)\n
           V. Berta   (vittoria.berta@unito.it)
           
  \date    Jan 26, 2026
*/
////////////////////////////////////////////////////////////////////////
#include "gpluto.hpp"

#if (HO_LAP_LIMITER == HO_JAMESON_LIM) || (HO_LAP_LIMITER == HO_RDER_LIM)
// int BitSet(int n, int k);

#if HO_LAP_LIMITER == HO_JAMESON_LIM
  #ifndef HO_LAP_LIMITER_EPS
    #define HO_LAP_LIMITER_EPS 0.07
  #endif
#else
  #ifndef HO_LAP_LIMITER_EPS
    #define HO_LAP_LIMITER_EPS 0.2
  #endif
#endif

#ifndef HO_LAP_LIMITER_BUF
  #define HO_LAP_LIMITER_BUF   0
#endif


////////////////////////////////////////////////////////////////////////
void DiscontinuityDetect(Data *data, Grid *grid, RBox *box)
/*!
  
  \param [in]      d        pointer to PLUTO Data structure
  \param [in]      grid     pointer to Grid structure
  \param [out]     eta      array of discontinuity detector switches
*/
////////////////////////////////////////////////////////////////////////
{
  int  i, j, k, nv;
  double umm, um, u, up, upp;

  #if HO_LAP_LIMITER_BUF > 0
  int nx1_tot = grid->ntot[IDIR];
  int nx2_tot = grid->ntot[JDIR];
  int nx3_tot = grid->ntot[KDIR];
  #endif
  int ibeg = box->ibeg, iend = box->iend;
  int jbeg = box->jbeg, jend = box->jend;
  int kbeg = box->kbeg, kend = box->kend;
  double tot = (iend-ibeg+1)*(jend-jbeg+1)*(kend-kbeg+1);

  // ------------------------------------------------------
  // 0. Allocate memory, define box size
  // ------------------------------------------------------

  #if (HO_LAP_LIMITER == NO) && (HO_ORDER_REDUCTION == NO)
  return;
  #elif HO_LAP_LIMITER == HO_JAMESON_LIM
  #elif HO_LAP_LIMITER == HO_RDER_LIM
  #else
    #error Invalid HO_LAP_LIMITER
  #endif

  RBox tot_box;
  
  RBoxDefine (0, grid->ntot[IDIR]-1, 
              0, grid->ntot[JDIR]-1, 
              0, grid->ntot[KDIR]-1, &tot_box);

  #pragma acc enter data copyin(tot_box)
  #pragma omp target enter data map(to: tot_box)

  // ------------------------------------------------------
  // 1. Evaluating discontonuity presence and flagging it
  //    by means of different discontinuity detectors
  // ------------------------------------------------------
  
  #pragma omp target teams distribute parallel for collapse(3)
  #pragma acc parallel loop collapse(3) present(data, tot_box)
  BOX_LOOP(tot_box,i,j,k) data->ho_lim(i,j,k) = 0;    // Default: no discontinuity 

  double xcount = 0.0;

  // ------------------------------------------------------
  // 2. Loop over domain + ghost zones (-1 or -2)
  // ------------------------------------------------------

  #pragma omp target teams distribute parallel for collapse(3)
  #pragma acc parallel loop collapse(3) present(data, box)
  BOX_LOOP(*box,i,j,k){

    double etax[NVAR], etay[NVAR], etaz[NVAR], eta[NVAR];
    double uref[NVAR], vref[NVAR];
    double eps[NVAR];

    NVAR_LOOP(nv) eps[nv] = HO_LAP_LIMITER_EPS;

    uref[RHO] = 1.0;
    uref[MX1] = uref[RHO]*0.3;
    uref[MX2] = uref[RHO]*0.3;
    uref[MX3] = uref[RHO]*0.3;
    uref[ENG] = 10.0;
    #if (PHYSICS == MHD) || (PHYSICS == RMHD) || (PHYSICS == ResRMHD)
    uref[BX1] = 1.0;
    uref[BX2] = 1.0;
    uref[BX3] = 1.0;
    #endif

    NVAR_LOOP(nv) {
      etax[nv] = etay[nv] = etaz[nv] = eta[nv] = 0.0;
    }

    NVAR_LOOP(nv) vref[nv] = fabs(data->Vc(i,j,k,nv));

    #if (PHYSICS == MHD) || (PHYSICS == RMHD) || (PHYSICS == ResRMHD)
    double pm  =   (    data->Vc(i,j,k,BX1)*data->Vc(i,j,k,BX1) 
                      + data->Vc(i,j,k,BX2)*data->Vc(i,j,k,BX2) 
                      + data->Vc(i,j,k,BX3)*data->Vc(i,j,k,BX3) );
    vref[BX1] = vref[BX2] = vref[BX3] = 0.5*sqrt(pm);
    #if PHYSICS == ResRMHD
    double pe  =   (    data->Vc(i,j,k,EX1)*data->Vc(i,j,k,EX1) 
                      + data->Vc(i,j,k,EX2)*data->Vc(i,j,k,EX2) 
                      + data->Vc(i,j,k,EX3)*data->Vc(i,j,k,EX3) );
    vref[EX1] = vref[EX2] = vref[EX3] = 0.5*sqrt(pe);
    #endif
    #endif

    #if HAVE_ENERGY 
    double cs2 = data->Vc(i,j,k,PRS)/data->Vc(i,j,k,RHO);
    double cs  = sqrt(cs2);
    #else
    double cs  = g_isoSoundSpeed;
    #endif

    vref[VX1] = vref[VX2] = vref[VX3] = 0.5*cs;
    #if PHYSICS == ResRMHD
    vref[EX1] = vref[EX2] = vref[EX3] = cs*sqrt(pm);
    #endif

    // ------------------------------------------------------
    // 2a. Jameson's switch function
    // ------------------------------------------------------

    #if HO_LAP_LIMITER == HO_JAMESON_LIM
    NVAR_LOOP(nv){
      double epsr = 1.e-2;
      double eps1 = 1.e-18;

      um   = data->Vc(i-1,j,k,nv);
      u    = data->Vc(i,j,k,nv);
      up   = data->Vc(i+1,j,k,nv);
      etax[nv] =   fabs(up - 2.0*u + um)
             /(fabs(up) + 2.0*fabs(u) + fabs(um) + epsr*vref[nv] + eps1);
      #if INCLUDE_JDIR
      um   = data->Vc(i,j-1,k,nv);
      u    = data->Vc(i,j,k,nv);
      up   = data->Vc(i,j+1,k,nv);
      etay[nv] =   fabs(up - 2.0*u + um)
            / (fabs(up) + 2.0*fabs(u) + fabs(um) + epsr*vref[nv] + eps1);
      #endif
      #if INCLUDE_KDIR
      um   = data->Vc(i,j,k-1,nv);
      u    = data->Vc(i,j,k,nv);
      up   = data->Vc(i,j,k+1,nv);
      etaz[nv] =   fabs(up - 2.0*u + um)
            / (fabs(up) + 2.0*fabs(u) + fabs(um) + epsr*vref[nv] + eps1);
      #endif

      eta[nv] = sqrt(etax[nv]*etax[nv] + etay[nv]*etay[nv] + etaz[nv]*etaz[nv]);

      #if PHYSICS == ResRMHD
      int selected_var =   (nv == RHO) || (nv == PRS) 
                        || (nv == VX1) || (nv == VX2) || (nv == VX3)
                        || (nv == BX1) || (nv == BX2) || (nv == BX3)
                        || (nv == EX1) || (nv == EX2) || (nv == EX3);
      #else 
      int selected_var = 1;
      #endif

      if (eta[nv] > eps[nv] && selected_var) data->ho_lim(i,j,k) = 1;
      
    } // end loop on nv 
    #endif

    // ------------------------------------------------------
    // 2b. Ratio of derivatives limiter
    // ------------------------------------------------------
 
    #if (HO_LAP_LIMITER == HO_RDER_LIM)
    double maxu;
    double d1u, d2u, d3u, d4u;
    double eta1, eta2;
    NVAR_LOOP(nv){

      #if DIVB_CONTROL == DIV_CLEANING  // Do not use the scalar 
      if (nv == PSI_GLM) continue;      // multiplier            
      #endif

      umm = data->Vc(i-2,j,k,nv);
      um  = data->Vc(i-1,j,k,nv);
      u   = data->Vc(i,j,k,nv);
      up  = data->Vc(i+1,j,k,nv);
      upp = data->Vc(i+2,j,k,nv);

      // maxL = MAX(fabs(umm), fabs(um));
      // maxR = MAX(fabs(upp), fabs(up));
      // maxu = MAX(fabs(maxL), fabs(maxR));
      // maxu = MAX(fabs(u), fabs(maxu));
      // maxu = (fabs(umm) + fabs(um) + fabs(u) + fabs(up) + fabs(upp))/5.0;
      maxu = vref[nv];

      d1u = 0.5*(up - um);
      d2u = (up - 2.0*u + um);
      d3u = (-0.5*umm + um - up + 0.5*upp);
      d4u = (umm - 4.0*um + 6.0*u - 4.0*up + upp);

      eta1 = fabs(d3u)/(fabs(d3u) + fabs(d1u) + maxu + 1.e-12);
      eta2 = fabs(d4u)/(fabs(d2u) + fabs(d4u) + maxu + 1.e-12);
      etax[nv] =  MAX(eta1, eta2);

      #if INCLUDE_JDIR
      umm = data->Vc(i,j-2,k,nv);
      um  = data->Vc(i,j-1,k,nv);
      u   = data->Vc(i,j,k,nv);
      up  = data->Vc(i,j+1,k,nv);
      upp = data->Vc(i,j+2,k,nv);

      // maxL = MAX(fabs(umm), fabs(um));
      // maxR = MAX(fabs(upp), fabs(up));
      // maxu = MAX(fabs(maxL), fabs(maxR));
      // maxu = MAX(fabs(u), fabs(maxu));
      // maxu = (fabs(umm) + fabs(um) + fabs(u) + fabs(up) + fabs(upp))/5.0;
      maxu = vref[nv];

      d1u = 0.5*(up - um);
      d2u = (up - 2.0*u + um);
      d3u = (-0.5*umm + um - up + 0.5*upp);
      d4u = (umm - 4.0*um + 6.0*u - 4.0*up + upp);

      eta1 = fabs(d3u)/(fabs(d3u) + fabs(d1u) + maxu + 1.e-12);
      eta2 = fabs(d4u)/(fabs(d2u) + fabs(d4u) + maxu + 1.e-12);

      etay[nv] = MAX(eta1, eta2);

      #endif

      #if INCLUDE_KDIR
      umm = data->Vc(i,j,k-2,nv);
      um  = data->Vc(i,j,k-1,nv);
      u   = data->Vc(i,j,k,nv);
      up  = data->Vc(i,j,k+1,nv);
      upp = data->Vc(i,j,k+2,nv);

      // maxL = MAX(fabs(umm), fabs(um));
      // maxR = MAX(fabs(upp), fabs(up));
      // maxu = MAX(fabs(maxL), fabs(maxR));
      // maxu = MAX(fabs(u), fabs(maxu));
      // maxu = (fabs(umm) + fabs(um) + fabs(u) + fabs(up) + fabs(upp))/5.0;
      maxu = vref[nv];

      d1u = 0.5*(up - um);
      d2u = (up - 2.0*u + um);
      d3u = (-0.5*umm + um - up + 0.5*upp);
      d4u = (umm - 4.0*um + 6.0*u - 4.0*up + upp);

      eta1 = fabs(d3u)/(fabs(d3u) + fabs(d1u) + maxu + 1.e-12);
      eta2 = fabs(d4u)/(fabs(d2u) + fabs(d4u) + maxu + 1.e-12);

      etaz[nv] = MAX(eta1, eta2);
      #endif

      eta[nv] = sqrt(etax[nv]*etax[nv] + etay[nv]*etay[nv] + etaz[nv]*etaz[nv]);

      #if PHYSICS == ResRMHD
      int selected_var =   (nv == RHO) || (nv == PRS) 
                        || (nv == VX1) || (nv == VX2) || (nv == VX3)
                        || (nv == BX1) || (nv == BX2) || (nv == BX3)
                        || (nv == EX1) || (nv == EX2) || (nv == EX3);
      #else
      int selected_var = 1;
      #endif
      
      if (eta[nv] > eps[nv] && selected_var) {
        data->ho_lim(i,j,k) = 1;
      }
    } // end loop on nv 
    #endif // HO_LAP_LIMITER == HO_RDER_LIM 

  } // end loop on i, j, k 

  FlagBoundary(data->ho_lim, FLAG_COPY, grid);

  // ------------------------------------------------------
  // 3. Now add buffer zones 
  // ------------------------------------------------------

  #pragma omp target teams distribute parallel for collapse(3)
  #pragma acc parallel loop collapse(3) present(data, tot_box)
  BOX_LOOP(tot_box,i,j,k){
    if (data->ho_lim(i,j,k) == 1){
      xcount += 1.0/tot;
      data->flag(i,j,k) |= FLAG_HO_LAP_LIMITER;
      #if HO_ORDER_REDUCTION == LINEAR
      data->flag(i,j,k) |= FLAG_MINMOD;
      #endif

      #if HO_LAP_LIMITER_BUF > 0
      int i1, j1, k1;
      ibeg = MAX(i - HO_LAP_LIMITER_BUF,0), iend = MIN(i + HO_LAP_LIMITER_BUF,nx1_tot-1);
      jbeg = MAX(j - HO_LAP_LIMITER_BUF,0), jend = MIN(j + HO_LAP_LIMITER_BUF,nx2_tot-1);
      kbeg = MAX(k - HO_LAP_LIMITER_BUF,0), kend = MIN(k + HO_LAP_LIMITER_BUF,nx3_tot-1);
      #pragma acc loop collapse(3) seq
      for (k1 = kbeg; k1 <= kend; k1++ ){
      for (j1 = jbeg; j1 <= jend; j1++ ){
      for (i1 = ibeg; i1 <= iend; i1++ ){
        data->flag(i1,j1,k1) |= FLAG_HO_LAP_LIMITER;
        #if HO_ORDER_REDUCTION == LINEAR
        data->flag(i1,j1,k1) |= FLAG_MINMOD;
        #endif
      }}}
      #endif
    }
  } // end BOX_LOOP() on i, j, k 

  // ------------------------------------------------------
  // Necessary when BUF > (nghost-3).
  // If the buffer is not exchanged between processors, then 
  // there's a potential inconsistency in situations like this:          
  //
  // Proc #1   ---|---|-B-|-F-|-B-|---|---|...|...|...|
  // Proc #2                  |...|...|...|---|---|---|---|---
  //
  // F = Flagged zone
  // B = Buffer zone 
  // ------------------------------------------------------

  //  FlagBoundary(data->flag, grid);
  #ifdef PARALLEL
  ParallelExchangeFlag<FLAG_OR> (data->flag, grid);
  #endif
  
  #pragma omp target exit data map(delete: tot_box)
  #pragma acc exit data delete(tot_box)

  return;

  // ------------------------------------------------------
  // 3. Clustering
  // ------------------------------------------------------
/*
  #pragma acc parallel loop collapse(3)
  for (k = kp_beg; k <= kp_end; k++){
  for (j = jp_beg; j <= jp_end; j++){
  for (i = ip_beg; i <= ip_end; i++){
    int nfx = data->ho_lim(i+1,j,k)*data->ho_lim(i-1,j,k);
    int nfy = data->ho_lim(i,j+1,k)*data->ho_lim(i,j-1,k);

    if (nfx*nfy > 0 && (data->ho_lim(i,j,k) == 0)){
      data->ho_lim(i,j,k) = 1;
    }

    if (data->ho_lim(i,j,k) == 1) {
      data->flag(i,j,k) |= FLAG_HO_LAP_LIMITER;
    }
  }}}
  */

}
#endif // (HO_LAP_LIMITER == HO_JAMESON_LIM) || (HO_LAP_LIMITER == HO_RDER_LIM) 
