////////////////////////////////////////////////////////////////////////
/*! 
  \file  
  \brief  Solves the linear transport step of the FARGO-MHD algorithm. 
 
  Implements a uniform shift of the conservative fluid variables in 
  the direction of orbital motion which is assumed to be 
  \f$ x_2=y,\phi\f$ (for Cartesian or polar coordinates) or 
  \f$ x_3 = \phi\f$  (for spherical coordinates).
  The shift is converted into an integer part "m" and a fractional
  remainder "eps".\n
 
  Parallelization is performed by adding extra data layers
  above and below the current data set so as to enlarge the required
  computational stencil in the orbital direction.
  The additional buffers are received from the processors lying,
  respectively, below and above and their size changes dynamically 
  from one step to the next.
 
  \b Reference
    - "A conservative orbital advection scheme for simulations 
       of magnetized shear flows with the PLUTO code"\n
       Mignone et al, A&A (2012) 545, A152  (--> Section 2.4)

  \author A. Mignone (andrea.mignone@unito.it)\n
          M. Rossazza (marco.rossazza@unito.it) \n
          G. Muscianisi (g.muscianisi@cineca.it)
  \date   Dec 06, 2025
*/
////////////////////////////////////////////////////////////////////////
#include "gpluto.hpp"

#ifdef PARALLEL
static void FARGO_ParallelExchange (Ary4D &U, Ary3D (&Us)[6],
                                    Ary4D_acc &recv_buf_lower, RBox *,
                                    Ary4D_acc &recv_buf_upper, RBox *, Grid *grid);
#endif

#define FARGO_MOD(s) ((s) < sbeg ? (s)+ns:((s) > send ? (s)-ns:(s)))

////////////////////////////////////////////////////////////////////////
/*!
  Compute the transport flux corresponding to a fractional 
  increment eps.
  
  \param [in]  q   1D array of a conserved quantity
  \param [out] flx the conservative flux
  \param [in]  eps the fractional increment
  
  The following MAPLE script may be used to check the 
  correctness of implementation
  \code
   restart;
   P0 := 3/2*q[i] - (q[p]+q[m])/4;
   P1 := q[p] - q[m];
   P2 := 3*(q[p] - 2*q[i] + q[m]);
 
   FS[p] := P0 + P1/2*(1 - epsilon) + P2*(1/4 - epsilon/2 + epsilon^2/3);
   FS[m] := P0 - P1/2*(1 - epsilon) + P2*(1/4 - epsilon/2 + epsilon^2/3);
 
   dq  := q[p] - q[m];
   d2q := q[p] - 2*q[i] + q[m];
   FF[p] := q[p] - epsilon/2*(dq + d2q*(3 - 2*epsilon));
   FF[m] := q[m] + epsilon/2*(dq - d2q*(3 - 2*epsilon));
  \endcode
*/
////////////////////////////////////////////////////////////////////////

inline void FARGO_Flux (Data *data, int stag[3], int nv_tot, Ary2DInt &m, Ary2D &eps, Grid *grid)
{
  int    s, i, nv, m_l;
  double dqc, d2q, dqp, dqm;
  int nx1_tot = grid->ntot[IDIR];
  int ibeg    = grid->lbeg[IDIR], iend = grid->lend[IDIR];
  int jbeg    = grid->lbeg[JDIR], jend = grid->lend[JDIR];
  int kbeg    = grid->lbeg[KDIR], kend = grid->lend[KDIR];
  #if GEOMETRY != SPHERICAL
  int k;
  int nx2_tot = grid->ntot[JDIR];
  int sbeg   = jbeg;  // Starting index in the orbital dir. 
  int send   = jend;  // Final index in the orbital dir. 
  int ns_tot = nx2_tot;
  #else
  int nx3_tot = grid->ntot[KDIR];
  int j;
  int sbeg   =  kbeg;
  int send   =  kend;
  int ns_tot =  nx3_tot;
  #endif
  int nproc_s = grid->nproc[SDIR];   // -- # of processors along shift dir -- 

  #pragma acc parallel loop collapse(3) present(data, eps, grid)
#if GEOMETRY == CARTESIAN || GEOMETRY == CYLINDRICAL
  for(nv=0; nv<nv_tot; nv++){
  for(k=kbeg; k<=kend+stag[2]; k++){
  for(i=ibeg; i<=iend+stag[0]; i++){
    long int offset = (ns_tot+2*FARGO_BUF)*(i + nx1_tot*k)*NVAR;
    double eps_l = eps(i,k);
    if (nproc_s > 1) m_l = m(i,k);
    else m_l = 0;
#else
  for(nv=0; nv<nv_tot; nv++){
  for(j=jbeg; j<=jend+stag[1]; j++){
  for(i=ibeg; i<=iend+stag[0]; i++){
    long int offset = (ns_tot+2*FARGO_BUF)*(i + nx1_tot*j)*NVAR;
    if (nproc_s > 1) m_l = m(i,j);
    else m_l = 0;
    double eps_l = eps(i,j);
#endif
    Ary2D q(&data->sweep->vC[offset],(ns_tot+2*FARGO_BUF),NVAR);
    Ary2D flux(&data->sweep->flux[offset],(ns_tot+2*FARGO_BUF),NVAR);
    Ary2D qm(&data->sweep->vL[offset],ns_tot,NVAR);
    Ary2D qp(&data->sweep->vR[offset],ns_tot,NVAR);
    Ary2D dq(&data->sweep->fL[offset],ns_tot,NVAR);
    Ary2D dql(&data->sweep->fR[offset],ns_tot,NVAR);
    
    #ifdef STAGGERED_MHD
    DIM_EXPAND(if (nv == BX1) continue;  ,
               if (nv == BX2) continue;  ,
               if (nv == BX3) continue;)
    #endif

      // ----------------------------------------------------------------------
      // 1. Compute the slopes dq, dql, dqp, dqm
      // Note: slopes are computed in s rather than in s1
      // for memory savings. 
      // ----------------------------------------------------------------------

    dq(0,nv) = q(1+FARGO_BUF-m_l,nv) - q(0+FARGO_BUF-m_l,nv);
    #pragma acc loop vector
    for (s = 1; s < ns_tot-1; s++){
      int s1 = s-m_l+FARGO_BUF;
      dq(s,nv) = q(s1+1,nv) - q(s1,nv);
    } 

    #pragma acc loop vector
    for (s = 1; s < ns_tot-1; s++){
      if (dq(s,nv)*dq(s-1,nv) > 0.0){
        dqc = 0.5*(dq(s,nv) + dq(s-1,nv));
        d2q = 2.0*(fabs(dq(s,nv)) < fabs(dq(s-1,nv)) ? dq(s,nv):dq(s-1,nv));
        dql(s,nv) = fabs(d2q) < fabs(dqc) ? d2q: dqc;
      }else dql(s,nv) = 0.0;
    } 

      // ----------------------------------------------------------------------
      // 2. Apply the Van Leer limiter or PPM
      // ----------------------------------------------------------------------

    #if FARGO_ORDER == 2  // MUSCL-Hancock, 2nd order 
    if (eps_l > 0.0){
      #pragma acc loop vector
      for (s=sbeg-1; s<=send; s++) {
        int s1 = s-m_l+FARGO_BUF;
        flux(s1+1,nv) = q(s1,nv) + 0.5*dql(s,nv)*(1.0 - eps_l);
      }
    }else { // -- remember: eps_l > 0 by definition -- 
      #pragma acc loop vector
      for (s=sbeg-1; s<=send; s++) {
        int s1 = s-m_l+FARGO_BUF;
        flux(s1+1,nv) = q(s1+1,nv) - 0.5*dql(s+1,nv)*(1.0 + eps_l);
      }
    }
    #elif FARGO_ORDER == 3  // PPM, 3rd order 
    #pragma acc loop vector
    for (s = sbeg-1; s <= send+1; s++){
      int s1 = s-m_l+FARGO_BUF;
      dqp =  0.5*dq(s,nv)   - (dql(s+1,nv) - dql(s,nv))/6.0;
      dqm = -0.5*dq(s-1,nv) + (dql(s-1,nv) - dql(s,nv))/6.0;
      if (dqp*dqm > 0.0) dqp = dqm = 0.0;
      else{
        if (fabs(dqp) >= 2.0*fabs(dqm)) dqp = -2.0*dqm;
        if (fabs(dqm) >= 2.0*fabs(dqp)) dqm = -2.0*dqp;
      }
      qp(s,nv) = q(s1,nv) + dqp;
      qm(s,nv) = q(s1,nv) + dqm;
    }
    if (eps_l > 0.0){
      #pragma acc loop vector
      for (s = sbeg-1; s <= send; s++){
        int s1 = s-m_l+FARGO_BUF;
        dqc = qp(s,nv) - qm(s,nv);
        d2q = qp(s,nv) - 2.0*q(s1,nv) + qm(s,nv);
        flux(s1+1,nv) = qp(s,nv) - 0.5*eps_l*(dqc + d2q*(3.0 - 2.0*eps_l));
      }
    }else{  
      #pragma acc loop vector
      for (s = sbeg-1; s <= send; s++){
        int s1 = s-m_l+FARGO_BUF;
        dqc = qp(s+1,nv) - qm(s+1,nv);
        d2q = qp(s+1,nv) - 2.0*q(s1+1,nv) + qm(s+1,nv);
        flux(s1+1,nv) = qm(s+1,nv) - 0.5*eps_l*(dqc - d2q*(3.0 + 2.0*eps_l));
      }
    }
    #else
    QuitPLUTO("! FARGO_ORDER should be either 2 or 3",1)
    #endif
  }}}

}

////////////////////////////////////////////////////////////////////////
/*!
  Shifts conserved variables in the orbital direction.
  
  \param [in,out]  U    a 3D array of conserved, zone-centered values
  \param [in,out]  Us   a 3D array of staggered magnetic fields
  \param [in]      grid  pointer to Grid structure;
  
  \return  This function has no return value.
*/
////////////////////////////////////////////////////////////////////////
void FARGO_ShiftSolution(Data *data, Ary4D &U, Ary3D (&Us)[6], Grid *grid, 
                         RunConfig &run_config, double dt)
#if GEOMETRY != SPHERICAL
 #define s j   // -- orbital direction index -- 
#else
 #define s k
#endif
{
mynvtxstart_("FARGO_ShiftSolution",BLUE);
  int    i, j, k, nv, ngh;
  int    sm;
  int    nproc_s = 1;
  double dphi, w, Lphi, dL;
  double **wA;
  double *x1  = grid->xC[IDIR];
  int nx1_tot = grid->ntot[IDIR];
  int nx2_tot = grid->ntot[JDIR];
  int nx3_tot = grid->ntot[KDIR];
  int ibeg    = grid->lbeg[IDIR], iend = grid->lend[IDIR];
  int jbeg    = grid->lbeg[JDIR], jend = grid->lend[JDIR];
  int kbeg    = grid->lbeg[KDIR], kend = grid->lend[KDIR];
  #if GEOMETRY != SPHERICAL
  int nx2     = grid->nint[JDIR];
  int sbeg   = jbeg;  // Starting index in the orbital dir. 
  int send   = jend;  // Final index in the orbital dir. 
  int ns     =  nx2;
  int ns_tot =  nx2_tot;
  #else
  double *x2 = grid->xC[JDIR];
  #ifdef STAGGERED_MHD
  double *x2p = grid->xR[JDIR], *x2m = grid->xL[JDIR];
  #endif
  int nx3     = grid->nint[KDIR];
  int sbeg   =  kbeg;
  int send   =  kend;
  int ns     =  nx3;
  int ns_tot =  nx3_tot;
  #endif
  #ifdef PARALLEL
  Ary4D_acc upper_buf;  // upper layer buffer  
  Ary4D_acc lower_buf;  // lower layer buffer
  RBox      lbox, ubox; // lower and upper box structures  
  #endif
  int stag[3] = {0,0,0}; // staggered directions for flux computation

  // ------------------------------------------------------
  // 0. In parallel mode, one-dimensional arrays q and flux
  //    must be large enough to contain data values coming 
  //    from neighbour processors.
  //    For this reason we augment them with an extra zone
  //    buffer containing FARGO_BUF cells on both sides.
  //  
  //    <..FARGO_BUF..>|-------| ... |-------|<..FARGO_BUF..>
   //                   sbeg-ng       send+ng
  // ------------------------------------------------------

  #if (GEOMETRY != SPHERICAL) && !INCLUDE_JDIR
  return;   // No shift needed in axisymmetric geometry 
  #endif

  // ------------------------------------------------------
  // 1. Allocate memory for static arrays.
  // ------------------------------------------------------

  #ifdef STAGGERED_MHD
  static Ary3D Ez, Ex;
  if (Ez.isAllocated() == 0){
    Ez.create(nx1_tot, nx2_tot, nx3_tot);
    Ex.create(nx1_tot, nx2_tot, nx3_tot);
    #pragma acc enter data copyin(Ez, Ex)
  }
  #endif

  // ------------------------------------------------------
  // 2. Get the pointer to the orbital speed array
  // ------------------------------------------------------

  wA = FARGO_Velocity();

  nproc_s = grid->nproc[SDIR];   // # of processors along shift dir 
  ngh     = grid->nghost[SDIR]; 
  Lphi    = run_config.domEnd[SDIR] - run_config.domBeg[SDIR];

  dphi = Lphi / grid->nint_glob[SDIR];

  // ------------------------------------------------------
  // 3. Fargo parallelization is performed by adding extra
  //    data layers above and below the current data set so
  //    as to enlarge the required computational stencil in
  //    the orbital direction:
  //
  //         ^
  //         |    |recv_lower|
  //         |    +----------+
  //         |    |          |
  //         |    |     U    |
  //         |    |          |
  //         |    +----------+
  //         |    |recv_upper|
  //
  //    The buffers recv_upper and recv_lower are received
  //    from the processors lying, respectively, below and
  //    above and their size changes dynamically from
  //    one step to the next.
  // ------------------------------------------------------
    
#ifdef PARALLEL

  int mmax, mmin;  // Used to determine buffer size from time to time 
  int nbuf_lower, nbuf_upper;

  // ------------------------------------------------------
  // 4. The value of mmax and mmin is kept from call to call
  //    because they are used to determine the size of the
  //    send/receive buffers each time step.
  //    At the very beginning, we initialize them to the
  //    largest possible value. 
  // ------------------------------------------------------

  mmax = MIN(ns, FARGO_BUF-1) - grid->nghost[SDIR] - 2;
  mmin = mmax;

  if (nproc_s > 1){ 
    nbuf_lower = abs(mmin) + ngh + 1;
    nbuf_upper =     mmax  + ngh + 1;

    // Check that buffer is not larger than processor size 

    if (nbuf_upper > ns || nbuf_lower > ns){
      QuitPLUTO("! FARGO_ShiftSolution(): buffer size exceeds processsor size",1);
    }

    // Check that buffer is less than FARGO_BUF
 
    if (nbuf_upper > FARGO_BUF || nbuf_lower > FARGO_BUF){
      QuitPLUTO ("! FARGO_ShiftSolution(): buffer size exceeds maximum length",1);
    }

    // Set grid indices of the lower & upper boxes

    #if GEOMETRY == CARTESIAN || GEOMETRY == CYLINDRICAL

    lbox.ibeg = ibeg; lbox.jbeg =            0; lbox.kbeg = kbeg; 
    lbox.iend = iend; lbox.jend = nbuf_lower-1; lbox.kend = kend;
    // lbox.di   = lbox.dj = lbox.dk = 1;

    ubox.ibeg = ibeg; ubox.jbeg =            0; ubox.kbeg = kbeg; 
    ubox.iend = iend; ubox.jend = nbuf_upper-1; ubox.kend = kend;
    // ubox.di = ubox.dj = ubox.dk = 1;

    #ifdef STAGGERED_MHD  // enlarge boxes to fit staggered fields 
    DIM_EXPAND(lbox.iend++; ubox.iend++;  ,
                                       ;  ,
               lbox.kend++; ubox.kend++;  )
    #endif 

    #elif GEOMETRY == SPHERICAL

    lbox.ibeg = ibeg; lbox.jbeg = jbeg; lbox.kbeg =            0; 
    lbox.iend = iend; lbox.jend = jend; lbox.kend = nbuf_lower-1; 
    // lbox.di = lbox.dj = lbox.dk = 1;

    ubox.ibeg = ibeg; ubox.jbeg = jbeg; ubox.kbeg =            0; 
    ubox.iend = iend; ubox.jend = jend; ubox.kend = nbuf_upper-1; 
    // ubox.di = ubox.dj = ubox.dk = 1;

    #ifdef STAGGERED_MHD  // -- enlarge boxes to fit staggered fields --
    DIM_EXPAND(lbox.iend++; ubox.iend++;  ,
               lbox.jend++; ubox.jend++;  ,
                                       ;  )
    #endif 

    #endif  // GEOMETRY 

    int uktot = ubox.kend - ubox.kbeg + 1;
    int ujtot = ubox.jend - ubox.jbeg + 1;
    int uitot = ubox.iend - ubox.ibeg + 1;
    int lktot = lbox.kend - lbox.kbeg + 1;    
    int ljtot = lbox.jend - lbox.jbeg + 1;
    int litot = lbox.iend - lbox.ibeg + 1;

  // -- allocate memory -- 

    upper_buf.create(uitot, ujtot, uktot, NVAR);
    lower_buf.create(litot, ljtot, lktot, NVAR);

  // -- exchange data values between neighbour processes -- 
 
    FARGO_ParallelExchange(U, Us, lower_buf, &lbox, upper_buf, &ubox, grid);
   }
#endif  // PARALLEL 

  // ------------------------------------------------------
  // 5. Shift cell-centered quantities
  // ------------------------------------------------------

  static Ary2DInt m;
  static Ary2D eps;
  if (m.isAllocated() == 0){
      #if GEOMETRY != SPHERICAL
      m.create(nx1_tot, nx3_tot);
      eps.create(nx1_tot, nx3_tot);
      #else
      m.create(nx1_tot, nx2_tot);
      eps.create(nx1_tot, nx2_tot);
      #endif
      #pragma acc enter data copyin(m, eps)
  }

  // 5a. Determine spatial shift, integer shift and fractional reminder

  #pragma acc parallel loop collapse(2) present(m, eps, wA, x1)
  #if GEOMETRY == CARTESIAN || GEOMETRY == CYLINDRICAL
  for(k=kbeg; k<=kend; k++){
  for(i=ibeg; i<=iend; i++){
    int *m_l = &m(i,k);
    double *eps_l = &eps(i,k); 
  #else
  for(j=jbeg; j<=jend; j++){
  for(i=ibeg; i<=iend; i++){
    int *m_l = &m(i,j);
    double *eps_l = &eps(i,j);
  #endif

    #if GEOMETRY == CARTESIAN
    w = wA[k][i];
    #elif GEOMETRY == CYLINDRICAL
    w = wA[k][i]/x1[i];
    #elif GEOMETRY == SPHERICAL
    w = wA[j][i]/(x1[i]*sin(x2[j]));
    #endif

    dL  = w*dt;                         // spatial shift 
    dL  = fmod(dL, Lphi);
    *m_l   = (int)floor(dL/dphi + 0.5); // nearest integer part 
    *eps_l = dL/dphi - *m_l;            // fractional remainder 
  }}

  // 5b. Copy to 1D buffers and use them to obtain Fargo Flux

  #pragma acc parallel loop collapse(3) present(U, data)
  #if GEOMETRY == CARTESIAN || GEOMETRY == CYLINDRICAL
  NVAR_LOOP(nv){
  for(k=kbeg; k<=kend; k++){
  for(i=ibeg; i<=iend; i++){
    long int offset = (ns_tot+2*FARGO_BUF)*(i + nx1_tot*k)*NVAR;
  #else
  NVAR_LOOP(nv){
  for(j=jbeg; j<=jend; j++){
  for(i=ibeg; i<=iend; i++){
    long int offset = (ns_tot+2*FARGO_BUF)*(i + nx1_tot*j)*NVAR;
  #endif

    #ifdef STAGGERED_MHD
    DIM_EXPAND(if (nv == BX1) continue;  ,
               if (nv == BX2) continue;  ,
               if (nv == BX3) continue;)
    #endif

    Ary2D q(&data->sweep->vC[offset],(ns_tot+2*FARGO_BUF),NVAR);
    Ary2D flux(&data->sweep->flux[offset],(ns_tot+2*FARGO_BUF),NVAR);
    
    for(s=sbeg; s<=send; s++) q(s+FARGO_BUF,nv) = U(i,j,k,nv);

    if (nproc_s > 1){ // -- copy values from lower and upper buffers -- 
      #ifdef PARALLEL
      for (s = sbeg - nbuf_upper; s < sbeg; s++){
        q(s+FARGO_BUF,nv) = FARGO_ARRAY_INDEX(upper_buf, nv, sbeg-1-s, k-kbeg, j-jbeg, i-ibeg);
      }
      for (s = send+1; s <= send + nbuf_lower; s++){
        q(s+FARGO_BUF,nv) = FARGO_ARRAY_INDEX(lower_buf, nv, s-send-1, k-kbeg, j-jbeg, i-ibeg);
      }
      #endif
    }else{               // -- impose periodic b.c. -- 
      #pragma acc loop vector
      for (s = 1; s <= ngh; s++) {
        q(sbeg - s, nv) = q(sbeg - s + ns, nv); 
        q(send + s, nv) = q(send + s - ns, nv); 
      }    
    }
  }}}

  FARGO_Flux (data, stag, NVAR, m, eps, grid);

  // 5c. Shift solution ahead in parallel mode

  if (nproc_s > 1){
    #ifdef PARALLEL
    #pragma acc parallel loop collapse(4) present(U, data, grid)
    #if GEOMETRY == CARTESIAN || GEOMETRY == CYLINDRICAL
    NVAR_LOOP(nv){
    for(k=kbeg; k<=kend; k++){
    for(i=ibeg; i<=iend; i++){
    for (s = sbeg; s <= send; s++){
      long int offset = (ns_tot+2*FARGO_BUF)*(i + nx1_tot*k)*NVAR;
      int *m_l = &m(i,k);
      double *eps_l = &eps(i,k);
    #else
    NVAR_LOOP(nv){
    for(j=jbeg; j<=jend; j++){
    for(i=ibeg; i<=iend; i++){
    for (s = sbeg; s <= send; s++){
      long int offset = (ns_tot+2*FARGO_BUF)*(i + nx1_tot*j)*NVAR;
      int *m_l = &m(i,j);
      double *eps_l = &eps(i,j);
    #endif

      #ifdef STAGGERED_MHD
      DIM_EXPAND(if (nv == BX1) continue;  ,
                 if (nv == BX2) continue;  ,
                 if (nv == BX3) continue;)
      #endif

      Ary2D q(&data->sweep->vC[offset],(ns_tot+2*FARGO_BUF),NVAR);
      Ary2D flux(&data->sweep->flux[offset],(ns_tot+2*FARGO_BUF),NVAR);

      sm = s - *m_l;
      U(i,j,k,nv) = q(sm+FARGO_BUF,nv) - *eps_l*(flux(sm+1+FARGO_BUF,nv) - flux(sm+FARGO_BUF,nv));
    }}}}
  #endif
  }else{
    #pragma acc parallel loop collapse(4) present(U, data, grid)
  #if GEOMETRY == CARTESIAN || GEOMETRY == CYLINDRICAL
    NVAR_LOOP(nv){
    for(k=kbeg; k<=kend; k++){
    for(i=ibeg; i<=iend; i++){
    for (s = sbeg; s <= send; s++){
      long int offset = (ns_tot+2*FARGO_BUF)*(i + nx1_tot*k)*NVAR;
      int *m_l = &m(i,k);
      double *eps_l = &eps(i,k);
  #else
    NVAR_LOOP(nv){
    for(j=jbeg; j<=jend; j++){
    for(i=ibeg; i<=iend; i++){
    for (s = sbeg; s <= send; s++){
      long int offset = (ns_tot+2*FARGO_BUF)*(i + nx1_tot*j)*NVAR;
      int *m_l = &m(i,j);
      double *eps_l = &eps(i,j);
  #endif

      #ifdef STAGGERED_MHD
      DIM_EXPAND(if (nv == BX1) continue;  ,
                 if (nv == BX2) continue;  ,
                 if (nv == BX3) continue;)
      #endif

      Ary2D q(&data->sweep->vC[offset],(ns_tot+2*FARGO_BUF),NVAR);
      Ary2D flux(&data->sweep->flux[offset],(ns_tot+2*FARGO_BUF),NVAR);

      sm = FARGO_MOD(s - *m_l);
      U(i,j,k,nv) = q(sm,nv) - *eps_l*(flux(sm+1,nv) - flux(sm,nv));
    }}}}
  }

  // ------------------------------------------------------
  // 6. Shift staggered magnetic fields
  // ------------------------------------------------------

#ifdef STAGGERED_MHD 
{
  double *dx1 = grid->dx[IDIR];
  double *dx2 = grid->dx[JDIR];
  double *dx3 = grid->dx[KDIR];
  double *x1p = grid->xR[IDIR];
  double *x1m = grid->xL[IDIR];
  int    ss;
  stag[BX1s] = 1;

  // ------------------------------------------------------
  // 6a. Compute:
  //      Ez   at x(i+1/2), y(j+1/2), z(k)  (Cartesian or Polar) 
  //     -Eth  at x(i+1/2), y(j), z(k+1/2)  (Spherical) 
  // ------------------------------------------------------

  #if GEOMETRY == CARTESIAN || GEOMETRY == CYLINDRICAL
  #pragma acc parallel loop collapse(2) present(wA, x1m, m, eps)
  for (k = kbeg; k   <= kend; k++){
  for (i = ibeg; i <= iend+1; i++){
    int *m_l = &m(i,k);
    double *eps_l = &eps(i,k);
  #else
  #pragma acc parallel loop collapse(2) present(wA, x1m, x2, m, eps)
  for (j = jbeg;   j <= jend; j++){
  for (i = ibeg; i <= iend+1; i++){
    int *m_l = &m(i,j);
    double *eps_l = &eps(i,j);
  #endif

    #if GEOMETRY == CARTESIAN
    w = 0.5*(wA[k][i-1] + wA[k][i]);
    #elif GEOMETRY == CYLINDRICAL
    w = 0.5*(wA[k][i-1] + wA[k][i])/x1m[i];
    #elif GEOMETRY == SPHERICAL
    w = 0.5*(wA[j][i-1] + wA[j][i])/(x1m[i]*sin(x2[j]));
    #endif

    // ----------------------------------------------------
    // 6b. Obtain shift in terms of an integer
    //     number of cells (m) plus a remainder eps. 
    // ----------------------------------------------------

    dL  = w*dt;                 // spatial shift 
    dL  = fmod(dL, Lphi);
    *m_l   = (int)floor(dL/dphi + 0.5); // nearest integer part 
    *eps_l = dL/dphi - *m_l;               // remainder in units of dphi 
  }}

  #if GEOMETRY == CARTESIAN || GEOMETRY == CYLINDRICAL
  #pragma acc parallel loop collapse(2) present(data, Us)
  for (k = kbeg; k   <= kend; k++){
  for (i = ibeg; i <= iend+1; i++){
    long int offset = (ns_tot+2*FARGO_BUF)*(i + nx1_tot*k)*NVAR;
  #else
  #pragma acc parallel loop collapse(2) present(data, Us)
  for (j = jbeg;   j <= jend; j++){
  for (i = ibeg; i <= iend+1; i++){
    long int offset = (ns_tot+2*FARGO_BUF)*(i + nx1_tot*j)*NVAR;
  #endif

    Ary2D q(&data->sweep->vC[offset],(ns_tot+2*FARGO_BUF),NVAR);
    Ary2D flux(&data->sweep->flux[offset],(ns_tot+2*FARGO_BUF),NVAR);

    // Copy 3D array into 1D buffer

    #pragma acc loop vector
    for(s=sbeg; s<=send; s++) q(s+FARGO_BUF,0) = Us[BX1s](i,j,k);

    if (nproc_s > 1){  // Copy values from lower and upper buffers
      #ifdef PARALLEL
      #pragma acc loop vector
      for (s = sbeg - nbuf_upper; s < sbeg; s++){
        q(s+FARGO_BUF,0) = FARGO_ARRAY_INDEX(upper_buf, BX1, sbeg-1-s, k-kbeg, j-jbeg, i-ibeg); 
      }
      #pragma acc loop vector
      for (s = send+1; s <= send + nbuf_lower; s++){
        q(s+FARGO_BUF,0) =  FARGO_ARRAY_INDEX(lower_buf, BX1, s-send-1, k-kbeg, j-jbeg, i-ibeg);
      }
      #endif
    }else{                  // Assign periodic b.c.
      #pragma acc loop vector
      for (s = 1; s <= ngh; s++) {
        q(sbeg - s, 0) = q(sbeg - s + ns, 0); 
        q(send + s, 0) = q(send + s - ns, 0); 
      }    
    }
  }}

  FARGO_Flux (data, stag, 1, m, eps, grid); 

  if (nproc_s > 1){ 
  #ifdef PARALLEL
  #if GEOMETRY == CARTESIAN || GEOMETRY == CYLINDRICAL
    #pragma acc parallel loop collapse(2) present(Ez, data, m, eps)
    for (k = kbeg; k   <= kend; k++){
    for (i = ibeg; i <= iend+1; i++){
      long int offset = (ns_tot+2*FARGO_BUF)*(i + nx1_tot*k)*NVAR;
      int m_l = m(i,k);
      double eps_l = eps(i,k);
  #else
    #pragma acc parallel loop collapse(2) present(Ez, data, m, eps)
    for (j = jbeg;   j <= jend; j++){
    for (i = ibeg; i <= iend+1; i++){
      long int offset = (ns_tot+2*FARGO_BUF)*(i + nx1_tot*j)*NVAR;
      int m_l = m(i,j);
      double eps_l = eps(i,j);
  #endif

      Ary2D q(&data->sweep->vC[offset],(ns_tot+2*FARGO_BUF),NVAR);
      Ary2D flux(&data->sweep->flux[offset],(ns_tot+2*FARGO_BUF),NVAR);

      // Update main solution staggered array 
      
      if (m_l >= 0){
        #pragma acc loop vector
        for (s = sbeg; s <= send+1; s++) {
          sm = s -m_l;
          Ez(i,j,k) = eps_l*flux(sm+FARGO_BUF,0);
          for (ss = s; ss > (s-m_l); ss--) Ez(i,j,k) += q(ss-1+FARGO_BUF,0);
          Ez(i,j,k) *= dphi;
        }
      }else{
        #pragma acc loop vector
        for (s = sbeg; s <= send+1; s++) {
          sm = s -m_l;
          Ez(i,j,k) = eps_l*flux(sm+FARGO_BUF,0);
          for (ss = s+1; ss < (s-m_l+1); ss++) Ez(i,j,k) -= q(ss-1+FARGO_BUF,0);
          Ez(i,j,k) *= dphi;
        }
      } 
    }}
  #endif
  } else {
  #if GEOMETRY == CARTESIAN || GEOMETRY == CYLINDRICAL
    #pragma acc parallel loop collapse(2) present(Ez, data, m, eps)
    for (k = kbeg; k   <= kend; k++){
    for (i = ibeg; i <= iend+1; i++){
      long int offset = (ns_tot+2*FARGO_BUF)*(i + nx1_tot*k)*NVAR;
      int m_l = m(i,k);
      double eps_l = eps(i,k);
  #else
    #pragma acc parallel loop collapse(2) present(Ez, data, m, eps)
    for (j = jbeg;   j <= jend; j++){
    for (i = ibeg; i <= iend+1; i++){
      long int offset = (ns_tot+2*FARGO_BUF)*(i + nx1_tot*j)*NVAR;
      int m_l = m(i,j);
      double eps_l = eps(i,j);
  #endif

      Ary2D q(&data->sweep->vC[offset],(ns_tot+2*FARGO_BUF),NVAR);
      Ary2D flux(&data->sweep->flux[offset],(ns_tot+2*FARGO_BUF),NVAR);

      // Update main solution staggered array

      if (m_l >= 0){
        #pragma acc loop vector
        for (s = sbeg; s <= send+1; s++) {
          sm = FARGO_MOD(s -m_l);
          Ez(i,j,k) = eps_l*flux(sm,0);
          for (ss = s; ss > (s-m_l); ss--) Ez(i,j,k) += q(FARGO_MOD(ss-1),0);
          Ez(i,j,k) *= dphi;
        }
      }else {
        #pragma acc loop vector
        for (s = sbeg; s <= send+1; s++) {
          sm = FARGO_MOD(s -m_l);
          Ez(i,j,k) = eps_l*flux(sm,0);
          for (ss = s+1; ss < (s-m_l+1); ss++) Ez(i,j,k) -= q(FARGO_MOD(ss-1),0);
          Ez(i,j,k) *= dphi;
        }
      }
    }}
  }

#if INCLUDE_KDIR
  
#if GEOMETRY == CARTESIAN || GEOMETRY == CYLINDRICAL
  stag[BX3s] = 1;
#else
  stag[BX2s] = 1;
#endif

  // ------------------------------------------------------
  // 6c. Compute:
  //      Ex  at x(i), y(j+1/2), z(k+1/2)  (Cartesian or Polar) 
  //     -Er  at x(i), y(j+1/2), z(k+1/2)  (Spherical) 
  // ------------------------------------------------------

#if GEOMETRY == CARTESIAN || GEOMETRY == CYLINDRICAL
  #pragma acc parallel loop collapse(2) present(wA, m, eps)
  for (k = kbeg; k <= kend+1; k++){
  for (i = ibeg; i <= iend; i++){
    int *m_l = &m(i,k);
    double *eps_l = &eps(i,k);
#else
  #pragma acc parallel loop collapse(2) present(wA, m, eps)
  for (j = jbeg; j <= jend+1; j++){
  for (i = ibeg; i <= iend; i++){
    int *m_l = &m(i,j);
    double *eps_l = &eps(i,j);
#endif

    #if GEOMETRY == CARTESIAN
    w = 0.5*(wA[k-1][i] + wA[k][i]);
    #elif GEOMETRY == CYLINDRICAL
    w = 0.5*(wA[k-1][i] + wA[k][i])/x1[i];
    #elif GEOMETRY == SPHERICAL
    w = 0.5*(wA[j-1][i] + wA[j][i])/(x1[i]*sin(x2p[j-1]));
    #endif

    // --------------------------------------------
    // 6d. Obtain shift in terms of an integer
    // number of cells (m) plus a remainder eps. 
    // --------------------------------------------

    dL  = w*dt;                 // spatial shift 
    dL  = fmod(dL, Lphi);
    *m_l   = (int)floor(dL/dphi + 0.5); // nearest integer part 
    *eps_l = dL/dphi -*m_l;               // remainder in units of dphi 
  }}

#if GEOMETRY == CARTESIAN || GEOMETRY == CYLINDRICAL
  #pragma acc parallel loop collapse(2) present(Ex, wA, x1, data, grid)
  for (k = kbeg; k <= kend+1; k++){
  for (i = ibeg; i <= iend; i++){
    long int offset = (ns_tot+2*FARGO_BUF)*(i + nx1_tot*k)*NVAR;
    int BS = BX3;
#else
  #pragma acc parallel loop collapse(2) present(Ex, wA, x1, data, grid)
  for (j = jbeg; j <= jend+1; j++){
  for (i = ibeg; i <= iend; i++){
    long int offset = (ns_tot+2*FARGO_BUF)*(i + nx1_tot*j)*NVAR;
    int BS = BX2;
#endif

    Ary2D q(&data->sweep->vC[offset],(ns_tot+2*FARGO_BUF),NVAR);
    Ary2D flux(&data->sweep->flux[offset],(ns_tot+2*FARGO_BUF),NVAR);

    // Copy 3D array into 1D buffer 

    #pragma acc loop vector
    #if GEOMETRY == CARTESIAN || GEOMETRY == CYLINDRICAL
    for(s=sbeg; s<=send; s++) q(s+FARGO_BUF,0) = -Us[BX3s](i,j,k); 
    #else
    for(s=sbeg; s<=send; s++) q(s+FARGO_BUF,0) = -Us[BX2s](i,j,k); 
    #endif

    if (nproc_s > 1){    // -- copy values from lower and upper buffers -- 
      #ifdef PARALLEL
      for (s = sbeg - nbuf_upper; s < sbeg; s++){
        q(s+FARGO_BUF,0) = -FARGO_ARRAY_INDEX(upper_buf, BS, sbeg-1-s, k-kbeg, j-jbeg, i-ibeg); 
      }
      for (s = send+1; s <= send + nbuf_lower; s++){
        q(s+FARGO_BUF,0) = -FARGO_ARRAY_INDEX(lower_buf, BS, s-send-1, k-kbeg, j-jbeg, i-ibeg);
      }
      #endif
    }else{               // -- assign periodic b.c. -- 
      #pragma acc loop vector
      for (s = 1; s <= ngh; s++) {                        
        q(sbeg - s,0) = q(sbeg - s + ns,0); 
        q(send + s,0) = q(send + s - ns,0); 
      }    
    }
  }}

  FARGO_Flux (data, stag, 1, m, eps, grid); 

  if (nproc_s > 1){ 
  #ifdef PARALLEL
  #if GEOMETRY == CARTESIAN || GEOMETRY == CYLINDRICAL
    #pragma acc parallel loop collapse(2) present(Ex, data, m, eps)
    for (k = kbeg; k <= kend+1; k++){
    for (i = ibeg; i <= iend; i++){
      long int offset = (ns_tot+2*FARGO_BUF)*(i + nx1_tot*k)*NVAR;
      int m_l = m(i,k);
      double eps_l = eps(i,k);
  #else
    #pragma acc parallel loop collapse(2) present(Ex, data, m, eps)
    for (j = jbeg; j <= jend+1; j++){
    for (i = ibeg; i <= iend; i++){
      long int offset = (ns_tot+2*FARGO_BUF)*(i + nx1_tot*j)*NVAR;
      int m_l = m(i,j);
      double eps_l = eps(i,j);
  #endif

      Ary2D q(&data->sweep->vC[offset],(ns_tot+2*FARGO_BUF),NVAR);
      Ary2D flux(&data->sweep->flux[offset],(ns_tot+2*FARGO_BUF),NVAR);
    
      if (m_l >= 0){
        #pragma acc loop vector
        for (s = sbeg; s <= send+1; s++) {
          sm = s -m_l;
          Ex(i,j,k) = eps_l*flux(sm+FARGO_BUF,0);
          for (ss = s; ss > (s-m_l); ss--) Ex(i,j,k) += q(ss-1+FARGO_BUF,0);
          Ex(i,j,k) *= dphi;
        }
      }else {
        #pragma acc loop vector
        for (s = sbeg; j <= send+1; s++) {
          sm = s -m_l;
          Ex(i,j,k) = eps_l*flux(sm+FARGO_BUF,0);
          for (ss = s+1; ss < (s-m_l+1); ss++) Ex(i,j,k) -= q(ss-1+FARGO_BUF,0);
          Ex(i,j,k) *= dphi;
        }
      }
    }}
  #endif
  } else {
  #if GEOMETRY == CARTESIAN || GEOMETRY == CYLINDRICAL
    #pragma acc parallel loop collapse(2) present(Ex, data, m, eps)
    for (k = kbeg; k <= kend+1; k++){
    for (i = ibeg; i <= iend; i++){
      long int offset = (ns_tot+2*FARGO_BUF)*(i + nx1_tot*k)*NVAR;
      int m_l = m(i,k);
      double eps_l = eps(i,k);
  #else
    #pragma acc parallel loop collapse(2) present(Ex, data, m, eps)
    for (j = jbeg; j <= jend+1; j++){
    for (i = ibeg; i <= iend; i++){
      long int offset = (ns_tot+2*FARGO_BUF)*(i + nx1_tot*j)*NVAR;
      int m_l = m(i,j);
      double eps_l = eps(i,j);
  #endif

      Ary2D q(&data->sweep->vC[offset],(ns_tot+2*FARGO_BUF),NVAR);
      Ary2D flux(&data->sweep->flux[offset],(ns_tot+2*FARGO_BUF),NVAR);

    // -- update main solution array -- 

      if (m_l >= 0){
        #pragma acc loop vector
        for (s = sbeg; s <= send+1; s++) {
          sm = FARGO_MOD(s -m_l);
          Ex(i,j,k) = eps_l*flux(sm,0);
          for (ss = s; ss > (s-m_l); ss--) Ex(i,j,k) += q(FARGO_MOD(ss-1),0);
          Ex(i,j,k) *= dphi;
        }
      }else {
        #pragma acc loop vector
        for (s = sbeg; j <= send+1; s++) {
          sm = FARGO_MOD(s -m_l);
          Ex(i,j,k) = eps_l*flux(sm,0);
          for (ss = s+1; ss < (s-m_l+1); ss++) Ex(i,j,k) -= q(FARGO_MOD(ss-1),0);
          Ex(i,j,k) *= dphi;
        }
      }
    }}
  }

#endif //  DIMENSIONS == 3 

  // ------------------------------------------------------
  // 6e. Update bx1 staggered magnetic field
  // ------------------------------------------------------

  #pragma acc parallel loop collapse(3) present(Us, Ez, grid)
  for(k=kbeg; k <= kend; k++){
  for(j=jbeg; j <= jend; j++){
  for (i = ibeg; i <= iend+1; i++){
    #if GEOMETRY == CARTESIAN || GEOMETRY == CYLINDRICAL
    Us[BX1s](i,j,k) -= (Ez(i,j+1,k) - Ez(i,j,k))/dphi;
    #elif GEOMETRY == SPHERICAL // in spherical coordinates Ez = -Eth 
//    Us[BX1s](i,j,k) -= (Ez(i,j,k) - Ez[k-1][j][i])/dphi;
    double dmu = cos(x2m[j]) - cos(x2p[j]);
    Us[BX1s](i,j,k) -= dx2[j]*sin(x2[j])/(dmu*dphi)*(Ez(i,j,k+1) - Ez(i,j,k));
    #endif
  }}}

  // ------------------------------------------------------
  // 6f. Update bx2 staggered magnetic field
  // ------------------------------------------------------

  #pragma acc parallel loop collapse(3) present(Us, Ez, Ex, grid)
  for(k=kbeg; k <= kend; k++){
  for (j = jbeg; j <= jend+1; j++){
  for(i=ibeg; i <= iend; i++){
    #if GEOMETRY == CARTESIAN
    Us[BX2s](i,j,k) += DIM_EXPAND(                                    ,
                           (Ez(i+1,j,k) - Ez(i,j,k))/dx1[i]  ,
                         - (Ex(i,j,k+1) - Ex(i,j,k))/dx3[k]);
    #elif GEOMETRY == CYLINDRICAL
    Us[BX2s](i,j,k) += DIM_EXPAND(                                             ,
                       (x1p[i]*Ez(i+1,j,k) -  x1m[i]*Ez(i,j,k))/dx1[i]  ,
                       - x1[i]*(Ex(i,j,k+1) - Ex(i,j,k))/dx3[k]); 
    #elif GEOMETRY == SPHERICAL // in spherical coordinates Ex = -Er 
    Us[BX2s](i,j,k) += (Ex(i,j,k+1) - Ex(i,j,k))/dphi;
    #endif
  }}}

  // ------------------------------------------------------
  // 6g. Update bx3 staggered magnetic field
  // ------------------------------------------------------

#if INCLUDE_KDIR
  #if GEOMETRY == SPHERICAL
  double A1p, A1m, A2p, A2m;
  #endif
  #pragma acc parallel loop collapse(3) present(Us, Ez, Ex, grid)
  for(k=kbeg; k <= kend+1; k++){
  for(j=jbeg; j <= jend; j++){
  for(i=ibeg; i <= iend; i++){
    #if GEOMETRY == CARTESIAN || GEOMETRY == CYLINDRICAL
    Us[BX3s](i,j,k) += (Ex(i,j+1,k) - Ex(i,j,k))/dphi;
    #elif GEOMETRY == SPHERICAL // Ez = -Eth, Ez = -Er 
    A1p = x1p[i]*x1p[i];
    A1m = x1m[i]*x1m[i];
    A2p = fabs(sin(x2p[j]));
    A2m = fabs(sin(x2m[j]));

    Us[BX3s](i,j,k) += sin(x2[j])*(   A1p*Ez(i+1,j,k) 
                                   -  A1m*Ez(i,j,k))/(x1[i]*dx1[i])
                      - (A2p*Ex(i,j+1,k) - A2m*Ex(i,j,k))/dx2[j];
    #endif
  }}}
#endif

  // ------------------------------------------------------
  // 7. Redefine zone-centered magnetic field
  // ------------------------------------------------------
  
  #pragma acc parallel loop collapse(3) present(Us, U, grid)
  for(k=kbeg; k <= kend; k++){
  for(j=jbeg; j <= jend; j++){
  for(i=ibeg; i <= iend; i++){
    double uc[NVAR];
    DIM_EXPAND(
      uc[BX1] = 0.5*(Us[BX1s](i+1,j,k) + Us[BX1s](i,j,k));  ,
      uc[BX2] = 0.5*(Us[BX2s](i,j+1,k) + Us[BX2s](i,j,k));  ,
      uc[BX3] = 0.5*(Us[BX3s](i,j,k+1) + Us[BX3s](i,j,k));)

    #if GEOMETRY != CARTESIAN
    Metric metric;
    MetricSet(metric, grid->rC[i], grid->sC[j]);
    OrthoToCons (uc, metric); 
    #endif
    
    DIM_EXPAND(
      U(i,j,k,BX1) = uc[BX1];  ,
      U(i,j,k,BX2) = uc[BX2];  ,
      U(i,j,k,BX3) = uc[BX3];)
   
  }}}
}
#endif // STAGGERED_MHD 

  #ifdef PARALLEL
  upper_buf.destroy();
  lower_buf.destroy();
  #endif

  mynvtxstop_();

}

#ifdef PARALLEL
////////////////////////////////////////////////////////////////////////
void FARGO_ParallelExchange (Ary4D &U, Ary3D (&Us)[6],
                             Ary4D_acc &recv_buf_lower, RBox *lbox,
                             Ary4D_acc &recv_buf_upper, RBox *ubox,
                             Grid *grid)
/*!
  Exchange data values between  adjacent processors lying in the
  the same direction (x2 for CARTESIAN/CYLINDRICAL or x3 for SPHERICAL)
  Values will be stored inside recv_buf_lower & recv_buf_upper.
  
  \param [in,out]  U    a 3D array of conserved, zone-centered values
  \param [in,out]  Us   a 3D array of staggered magnetic fields
  \param [out]     recv_buf_lower buffer array received from the 
  processor ahead
  \param [out]     recv_buf_upper buffer array received from the 
  processor behind
  \param [in]      ubox           pointer to RBox structure defining 
  the index range of the upper data layer
  \param [in]      lbox           pointer to RBox structure defining 
  the index range of the lower data layer
  \param [in]      grid  pointer to Grid structure;
  
  \return  This function has no return value.
*/
////////////////////////////////////////////////////////////////////////
{
  int   nv, i, j, k;
  int   ib, jb, kb;
  int   ie, je, ke;
  int   coords[3], dst, src;
  #if GEOMETRY == CARTESIAN || GEOMETRY == CYLINDRICAL
  int jbeg    = grid->lbeg[JDIR], jend = grid->lend[JDIR];
  #elif GEOMETRY == SPHERICAL
  int kbeg    = grid->lbeg[KDIR], kend = grid->lend[KDIR];
  #endif
  int sitot, sjtot, sktot;

  long int  ntot;
  Ary4D_acc lsend_buf,usend_buf;
  // double *usend_buf;
  double *usbuf, *urbuf;
  double *lsbuf, *lrbuf;
  MPI_Comm cartcomm;
  MPI_Request reqs[4];  
  int r = 0;

  // -------------------------------------------------------------------
  // get ranks of the processors lying above (upper) and below (lower)
  // ------------------------------------------------------------------

#if USE_AL == NO
cartcomm = new_SZ.cart_comm;
#else
  AL_Get_cart_comm(SZ, &cartcomm); 
#endif
  for (i = 0; i < DIMENSIONS; i++) coords[i] = grid->rank_coord[i];
  coords[SDIR] += 1;
  MPI_Cart_rank(cartcomm, coords, &dst);

  for (i = 0; i < DIMENSIONS; i++) coords[i] = grid->rank_coord[i];
  coords[SDIR] -= 1;
  MPI_Cart_rank(cartcomm, coords, &src);

  // -- determine buffer size & dimensions -- 

  ib = ubox->ibeg; jb = ubox->jbeg; kb = ubox->kbeg;
  ie = ubox->iend; je = ubox->jend; ke = ubox->kend;

  sitot = (ie - ib + 1);
  sjtot = (je - jb + 1);
  sktot = (ke - kb + 1);
  ntot  = (ie - ib + 1)*(je - jb + 1)*(ke - kb + 1)*NVAR;

  usend_buf.create(sitot, sjtot, sktot, NVAR);

  // -- determine buffer size & dimensions -- 

  ib = lbox->ibeg; jb = lbox->jbeg; kb = lbox->kbeg;
  ie = lbox->iend; je = lbox->jend; ke = lbox->kend;

  ntot  = (ie - ib + 1)*(je - jb + 1)*(ke - kb + 1)*NVAR;

  sitot = (ie - ib + 1);
  sjtot = (je - jb + 1);
  sktot = (ke - kb + 1);

  lsend_buf.create(sitot, sjtot, sktot, NVAR);

  usbuf = &(usend_buf(0,0,0,0));
  urbuf = &(recv_buf_upper(0,0,0,0));
  
  lsbuf = &(lsend_buf(0,0,0,0));
  lrbuf = &(recv_buf_lower(0,0,0,0));

  // -- MPI_Irecv -- 

  MPI_Irecv(urbuf, ntot, MPI_DOUBLE, src, 1, cartcomm, &reqs[r++]);
  MPI_Irecv(lrbuf, ntot, MPI_DOUBLE, dst, 100, cartcomm, &reqs[r++]);

  // -------------------------------------------------------------
  // Filling buffer the upper layer
  // -------------------------------------------------------------

  // -- copy data -- 

  #pragma acc parallel loop collapse(3) present(U, Us) async(1)
  BOX_LOOP(*ubox,i,j,k){
    #if GEOMETRY == CARTESIAN || GEOMETRY == CYLINDRICAL
// usend_buf[nv](i,j,k) = FARGO_ARRAY_INDEX(U, send-s, k,j,i); 
    for (nv = 0; nv < NVAR; nv++) usend_buf(i-ib,j-jb,k-kb,nv) = U(i,jend-j,k,nv); 

    #ifdef STAGGERED_MHD
    DIM_EXPAND(usend_buf(i-ib,j-jb,k-kb,BX1) = Us[BX1s](i,jend-j,k);  ,
                                                          ;  ,
                usend_buf(i-ib,j-jb,k-kb,BX3) = Us[BX3s](i,jend-j,k);)
    #endif

    #elif GEOMETRY == SPHERICAL
     for (nv = 0; nv < NVAR; nv++) usend_buf(i-ib,j-jb,k-kb,nv) = U(i,j,kend-k,nv);

     #ifdef STAGGERED_MHD
      DIM_EXPAND(usend_buf(i-ib,j-jb,k-kb,BX1) = Us[BX1s](i,j,kend-k);  ,
                 usend_buf(i-ib,j-jb,k-kb,BX2) = Us[BX2s](i,j,kend-k);  ,
                                                            ;)
     #endif

    #endif
  }

  // -------------------------------------------------------------
  // Filling buffer of the lower layer
  // -------------------------------------------------------------

  #pragma acc parallel loop collapse(3) present(U, Us) async(2)
  BOX_LOOP(*lbox,i,j,k) {
    #if GEOMETRY == CARTESIAN || GEOMETRY == CYLINDRICAL

     for (nv = 0; nv < NVAR; nv++) lsend_buf(i-ib,j-jb,k-kb,nv)= U(i,jbeg+j,k,nv); 

     #ifdef STAGGERED_MHD
      DIM_EXPAND(lsend_buf(i-ib,j-jb,k-kb,BX1) = Us[BX1s](i,jbeg+j,k);  ,
                                                                ;  ,
                 lsend_buf(i-ib,j-jb,k-kb,BX3) = Us[BX3s](i,jbeg+j,k);)
     #endif

    #elif GEOMETRY == SPHERICAL

     for (nv = 0; nv < NVAR; nv++) lsend_buf(i-ib,j-jb,k-kb,nv) = U(i,j,kbeg+k,nv); 

     #ifdef STAGGERED_MHD
      DIM_EXPAND(lsend_buf(i-ib,j-jb,k-kb,BX1) = Us[BX1s](i,j,kbeg+k);  ,
                 lsend_buf(i-ib,j-jb,k-kb,BX2) = Us[BX2s](i,j,kbeg+k);  ,
                                                                ;)
     #endif

    #endif
  }
  #pragma acc wait

  // -- MPI_Isend -- 

  // Post send to dst 
  MPI_Isend(usbuf, ntot, MPI_DOUBLE, dst, 1, cartcomm, &reqs[r++]);
  // Post send to src 
  MPI_Isend(lsbuf, ntot, MPI_DOUBLE, src, 100, cartcomm, &reqs[r++]);

  MPI_Waitall(r, reqs, MPI_STATUSES_IGNORE);

  lsend_buf.destroy();
  usend_buf.destroy();

  #ifndef STAGGERED_MHD
  #endif 
  
}
#endif // PARALLEL 

#if PARTICLES
////////////////////////////////////////////////////////////////////////
void FARGO_ShiftParticles(Data *data, Grid *grid, double dt)
/*!
  
  
*/
////////////////////////////////////////////////////////////////////////
{
//   int    i,  j,  k;
//   int    i1, j1, k1;
//   double w, v, **wA;
//   static double ***wp;
//   particleNode *curNode;
//   Particle *p;

//   #if !INCLUDE_JDIR
//   return;
//   #endif
  
// /* -----------------------------------------------------
//    0. Allocate memory
//    ----------------------------------------------------- */

//   if (wp == NULL){
//     wp = ARRAY_BOX (-1, 1, -1, 1, -1, 1, double);
//   }

// /* -----------------------------------------------------
//    1. Shift particles in the orbital direction
//    ----------------------------------------------------- */
  
//   wA = FARGO_Velocity();
//   PARTICLES_LOOP(curNode, d->PHead){
//     p = &(curNode->p);

//   /* -- Interpolate fargo shift velocity at particle position -- */

//     Particles_GetWeights(p, p->cell, wp, grid);  
//     i = p->cell[IDIR];
//     j = p->cell[JDIR];
//     k = p->cell[KDIR];

//     v = 0.0;
//     for (k1 = -INCLUDE_KDIR; k1 <= INCLUDE_KDIR; k1++){
//     for (j1 = -INCLUDE_JDIR; j1 <= INCLUDE_JDIR; j1++){
//     for (i1 = -INCLUDE_IDIR; i1 <= INCLUDE_IDIR; i1++){
//       v += wp[k1][j1][i1]*wA[k + k1][i + i1];
//     }}}
// // v = -SB_Q*SB_OMEGA*p->coord[IDIR];
//     p->coord[JDIR] += v*dt;
//   }

// /* -----------------------------------------------------
//    2. Shift particles in the orbital direction
//    ----------------------------------------------------- */

//   Particles_Boundary(data, grid);
//   Particles_BoundaryExchange(data, grid);
}
#endif
