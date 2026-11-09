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

  \authors V. Berta   (vittoria.berta@unito.it)\n
           A. Mignone (andrea.mignone@unito.it)

  \date   Jan 22, 2026
*/
////////////////////////////////////////////////////////////////////////
#include "gpluto.hpp"

////////////////////////////////////////////////////////////////////////
/*!
  \param [in,out]  d        pointer to PLUTO Data structure
  \param [in]      grid     pointer to Grid structure
  \param [in]      box      input box
*/
////////////////////////////////////////////////////////////////////////
void PointValue (Data *data, Grid *grid, RBox *box)
{
mynvtxstart_("POINT",YELLOW);

  int  i, j, k, nv;
  int ibeg = box->ibeg, iend = box->iend;
  int jbeg = box->jbeg, jend = box->jend;
  int kbeg = box->kbeg, kend = box->kend;

  uint32_t err_kernel = 0;    // Cumulative kernel error flag  

  double d2ux, d2uy, d2uz, du;
  #ifdef STAGGERED_MHD
  double d2bx, d2by, d2bz, db;
  #endif

  double *x1 = grid->xC[IDIR]; 
  double *x2 = grid->xC[JDIR]; 
  double *x3 = grid->xC[KDIR]; 
  #if GEOMETRY != CARTESIAN
  double *dx1 = grid->dx[IDIR];
  double *dx2 = grid->dx[JDIR];
  double *dx3 = grid->dx[KDIR];
  #endif
  
  RunConfig &run_config = *(data->run_config);

  // ------------------------------------------------------
  // 1. Revert to 2nd order when fallback ilsafe mode is on.
  //    (Point value and averages are the same)
  // ------------------------------------------------------

  // #if FALLBACK
  // if (data->integration_mode == FALLBACK){
  //   ConsToPrim3D (data, 0, grid, box);

  //   #pragma acc parallel loop collapse(4) present(data)
  //   NVAR_LOOP(nv) BOX_LOOP(*box,i,j,k){
  //     data->Upc(i,j,k,nv) = data->Uc(i,j,k,nv);
  //     data->Vpc(i,j,k,nv) = data->Vc(i,j,k,nv);
  //   }

  //   #ifdef STAGGERED_MHD
  //   #pragma acc parallel loop collapse(3) present(data)
  //   BOX_STAG_LOOP(*box,i,j,k){
  //     CT_EXPAND(nv, data->Vps[nv](i,j,k) = data->Vs[nv](i,j,k));
  //   }
  //   #endif

  //   return;
  // }
  // #endif

  // ------------------------------------------------------
  // 2. Recover point-value of staggered magnetic field. 
  //    Note: point-value evaluation extends to one extra 
  //          zone in the staggered direction in order to 
  //          later allow the cell-centered magentic field 
  //          to be interpolated with 4 values.
  // ------------------------------------------------------

  // 2a If STAGGERED_MHD -> Converting face-averaged staggered
  //    B field to high order point-value staggered B field

  #if (defined STAGGERED_MHD) && INCLUDE_IDIR // De-averaging along x 
  d2bx = d2by = d2bz = db = 0.0;
  #pragma omp target teams distribute parallel for collapse(3)
  #pragma acc parallel loop collapse(3) present(data)
  for (k = kbeg; k <= kend; k++){
  for (j = jbeg; j <= jend; j++){
  //#pragma acc loop vector
    for (i = ibeg - 1; i <= iend + 2; i++){

      DIM_EXPAND(                                                  ,
                 d2by = HO_DeAverage3D<JDIR>(data->Vs[BX1s],i,j,k);   ,
                 d2bz = HO_DeAverage3D<KDIR>(data->Vs[BX1s],i,j,k);)
      db = (d2by + d2bz)/24.;
      data->Vps[BX1s](i,j,k) = data->Vs[BX1s](i,j,k) - db;

      #if (PHYSICS == ResRMHD) && (DIVE_CONTROL == CONSTRAINED_TRANSPORT)
      DIM_EXPAND(                                     ,
                 d2ey = HO_LAPLY(data->Vs[EX1s],i,j,k);  ,
                 d2ez = HO_LAPLZ(data->Vs[EX1s],i,j,k);
      )
      de = (d2ey + d2ez)/24.;
      data->Vps[EX1s](i,j,k) = data->Vs[EX1s](i,j,k) - de;
      #endif
 
    }
  }} // end loop on j,k 
  #endif // INCLUDE_IDIR 

  #if (defined STAGGERED_MHD) && INCLUDE_JDIR // De-averaging along y 
  d2bx = d2by = d2bz = db = 0.0;
  #pragma omp target teams distribute parallel for collapse(3)
  #pragma acc parallel loop collapse(3) present(data)
  for (k = kbeg;     k <= kend;     k++){
  for (j = jbeg - 1; j <= jend + 2; j++){   
  //  #pragma acc loop vector
    for (i = ibeg;     i <= iend;     i++){

      DIM_EXPAND(d2bx = HO_DeAverage3D<IDIR>(data->Vs[BX2s],i,j,k);  ,
                                                                  ,
                 d2bz = HO_DeAverage3D<KDIR>(data->Vs[BX2s],i,j,k);)

      db = (d2bx + d2bz)/24.0;
      data->Vps[BX2s](i,j,k) = data->Vs[BX2s](i,j,k) - db;

      #if (PHYSICS == ResRMHD) && (DIVE_CONTROL == CONSTRAINED_TRANSPORT)
      DIM_EXPAND(d2ex = HO_LAPLX(data->Vs[EX2s],i,j,k);  ,
                                                      ,
                 d2ez = HO_LAPLZ(data->Vs[EX2s],i,j,k);
      )
      de = (d2ex + d2ez)/24.0;
      data->Vps[EX2s](i,j,k) = data->Vs[EX2s](i,j,k) - de;
      #endif
    }
  }} // end loop on i,k 
  #endif // INCLUDE_JDIR 

  #if (defined STAGGERED_MHD) &&  INCLUDE_KDIR // De-averaging along z 
  d2bx = d2by = d2bz = db = 0.0;
  #pragma omp target teams distribute parallel for collapse(3)
  #pragma acc parallel loop collapse(3) present(data)
  for (k = kbeg - 1; k <= kend + 2; k++){
  for (j = jbeg;     j <= jend;     j++){
   // #pragma acc loop vector
    for (i = ibeg;     i <= iend;     i++){

      DIM_EXPAND(d2bx = HO_DeAverage3D<IDIR>(data->Vs[BX3s],i,j,k); ,
                 d2by = HO_DeAverage3D<JDIR>(data->Vs[BX3s],i,j,k); ,
                )

      #if GEOMETRY == CYLINDRICAL
      d2bx += (data->Vs[BX3s](i+1,j,k) - data->Vs[BX3s](i-1,j,k))*dx1[i]/x1[i];
      #endif

      db = (d2bx + d2by)/24.0;
      data->Vps[BX3s](i,j,k) = data->Vs[BX3s](i,j,k) - db;

      #if (PHYSICS == ResRMHD) && (DIVE_CONTROL == CONSTRAINED_TRANSPORT)
      DIM_EXPAND(d2ex = HO_LAPLX(data->Vs[EX3s],i,j,k); ,
                 d2ey = HO_LAPLY(data->Vs[EX3s],i,j,k); ,
      )

      de = (d2ex + d2ey)/24.0;
      data->Vps[EX3s](i,j,k) = data->Vs[EX3s](i,j,k) - de;
      #endif
    }
  }} // end loop on i,j 
  #endif  // INCLUDE_KDIR 

  // ------------------------------------------------------
  // 3. Now obtain cell-centered point-values in conserved 
  //    variables data->Upc.
  // ------------------------------------------------------

  d2ux = d2uy = d2uz = du = 0.0;
  #pragma omp target teams distribute parallel for collapse(3)
  #pragma acc parallel loop collapse(3) present(data, run_config) copy(err_kernel)
  for (k = kbeg; k <= kend; k++){
  for (j = jbeg; j <= jend; j++){
  for (i = ibeg; i <= iend; i++){

    data->err(i,j,k) = 0; // Reset error

    // 3a. Check data integrity

    NVAR_LOOP(nv) {
      if (isnan(data->Uc(i,j,k,nv))) {
        ATOMIC_ERR_UPDATE(err_kernel, ERR_NAN, ERR_QUIT);
        KERNEL_ERR_UPDATE(data->err(i,j,k), ERR_NAN, ERR_QUIT);
      }
    }

    // 3b. Obtain point values of zone-centered var
  
    #if GEOMETRY != CARTESIAN
    double q1 = fabs(dx1[i]/x1[i]);
    #if GEOMETRY == SPHERICAL
    double q2 = dx2[j]*cos(x2[j])/sin(x2[j])*(1.0 + dx2[j]*dx2[j]/60.);
    #endif
    #endif

    #pragma acc loop
    NVAR_LOOP(nv) {

      DIM_EXPAND(d2ux = HO_DeAverage4D<IDIR>(data->Uc,i,j,k,nv);  ,
                 d2uy = HO_DeAverage4D<JDIR>(data->Uc,i,j,k,nv);  ,
                 d2uz = HO_DeAverage4D<KDIR>(data->Uc,i,j,k,nv);)

      #if GEOMETRY == CYLINDRICAL
      d2ux += (data->Uc(i+1,j,k,nv) - data->Uc(i-1,j,k,nv))*q1;
      #elif GEOMETRY == SPHERICAL
      DIM_EXPAND(d2ux *= 1.0 - 4.0/15.0*q1*q1;
                 d2ux += 2.0*(data->Uc(i+1,j,k,nv) - data->Uc(i-1,j,k,nv))*q1; ,
                 d2uy *= 1.0 + 2.0*dx2[j]*dx2[j]/15.0;
                 d2uy += (data->Uc(i,j+1,k,nv) - data->Uc(i,j-1,k,nv))*q2;     ,
      )
      #endif

      du = (d2ux + d2uy + d2uz)/24.0;
      #if HO_LAP_LIMITER != NO
      if (data->flag(i,j,k) & FLAG_HO_LAP_LIMITER) du = 0.0;
      #endif

      data->Upc(i,j,k,nv) = data->Uc(i,j,k,nv) - du;

    } // end NVAR_LOOP 

    // 3d. If STAGGERED_MHD -> High order reconstruction of
    //       cell-centered B field from primary staggered
    //       representation.

    #ifdef STAGGERED_MHD

    // Compute zone-centered point-value magnetic field energy 
    
    #if HO_FIELD_AVERAGE == 1
    
    DIM_EXPAND(
      data->Upc(i,j,k,BX1) = (-     (data->Vps[BX1s](i-1,j,k) + data->Vps[BX1s](i+2,j,k))
                              + 9.0*(data->Vps[BX1s](i,j,k)   + data->Vps[BX1s](i+1,j,k)))/16.0; ,
      data->Upc(i,j,k,BX2) = (-     (data->Vps[BX2s](i,j-1,k) + data->Vps[BX2s](i,j+2,k))
                              + 9.0*(data->Vps[BX2s](i,j,k)   + data->Vps[BX2s](i,j+1,k)))/16.0; ,
      data->Upc(i,j,k,BX3) = (-     (data->Vps[BX3s](i,j,k-1) + data->Vps[BX3s](i,j,k+2))
                              + 9.0*(data->Vps[BX3s](i,j,k)   + data->Vps[BX3s](i,j,k+1)))/16.0;
    )

    #if (PHYSICS == ResRMHD) && (DIVE_CONTROL == CONSTRAINED_TRANSPORT)
    DIM_EXPAND(
    data->Upc(i,j,k,EX1) = (-     (data->Vps[EX1s](i-1,j,k) + data->Vps[EX1s](i+2,j,k))
                            + 9.0*(data->Vps[EX1s](i,j,k)   + data->Vps[EX1s](i+1,j,k)))/16.0; ,
    data->Upc(i,j,k,EX2) = (-     (data->Vps[EX2s](i,j-1,k) + data->Vps[EX2s](i,j+2,k))
                            + 9.0*(data->Vps[EX2s](i,j,k)   + data->Vps[EX2s](i,j+1,k)))/16.0; ,
    data->Upc(i,j,k,EX3) = (-     (data->Vps[EX3s](i,j,k-1) + data->Vps[EX3s](i,j,k+2))
                            + 9.0*(data->Vps[EX3s](i,j,k)   + data->Vps[EX3s](i,j,k+1)))/16.0;
    )
    #endif

    #elif HO_FIELD_AVERAGE == 2

    double dBxR_dy, dBxR_dz, dBxL_dy, dBxL_dz;
    double dByR_dz, dByR_dx, dByL_dz, dByL_dx;
    double dBzR_dx, dBzR_dy, dBzL_dx, dBzL_dy;

    double d2BxR_dy2, d2BxR_dz2, d2BxL_dy2, d2BxL_dz2;
    double d2ByR_dz2, d2ByR_dx2, d2ByL_dz2, d2ByL_dx2;
    double d2BzR_dx2, d2BzR_dy2, d2BzL_dx2, d2BzL_dy2;

    // Compute 1st derivative correction 

    int ip = i+1, im = ip;
    int jp = j+1, jm = jp;
    int kp = k+1, km = kp;

    dBxR_dy = dBxL_dy = dByR_dx = dByL_dx = 0.0;
    dByR_dz = dByL_dz = dBzR_dy = dBzL_dy = 0.0;
    dBzR_dx = dBzL_dx = dBxR_dz = dBxL_dz = 0.0;

    d2BxR_dy2 = d2BxL_dy2 = d2ByR_dx2 = d2ByL_dx2 = 0.0;
    d2ByR_dz2 = d2ByL_dz2 = d2BzR_dy2 = d2BzL_dy2 = 0.0;
    d2BzR_dx2 = d2BzL_dx2 = d2BxR_dz2 = d2BxL_dz2 = 0.0;

    DIM_EXPAND(                                          
                                                              ,
      dByR_dx = 0.5*(data->Vs[BX2s](i+1,jp,k) - data->Vs[BX2s](i-1,jp,k));
      dByL_dx = 0.5*(data->Vs[BX2s](i+1,jm,k) - data->Vs[BX2s](i-1,jm,k));    

      dBxR_dy = 0.5*(data->Vs[BX1s](ip,j+1,k) - data->Vs[BX1s](ip,j-1,k));
      dBxL_dy = 0.5*(data->Vs[BX1s](im,j+1,k) - data->Vs[BX1s](im,j-1,k));  
                                                              ,
      dBzR_dx = 0.5*(data->Vs[BX3s](i+1,j,kp) - data->Vs[BX3s](i-1,j,kp));
      dBzL_dx = 0.5*(data->Vs[BX3s](i+1,j,km) - data->Vs[BX3s](i-1,j,km));  

      dBzR_dy = 0.5*(data->Vs[BX3s](i,j+1,kp) - data->Vs[BX3s](i,j-1,kp));
      dBzL_dy = 0.5*(data->Vs[BX3s](i,j+1,km) - data->Vs[BX3s](i,j-1,km));  

      dBxR_dz = 0.5*(data->Vs[BX1s](ip,j,k+1) - data->Vs[BX1s](ip,j,k-1));
      dBxL_dz = 0.5*(data->Vs[BX1s](im,j,k+1) - data->Vs[BX1s](im,j,k-1));  

      dByR_dz = 0.5*(data->Vs[BX2s](i,jp,k+1) - data->Vs[BX2s](i,jp,k-1));
      dByL_dz = 0.5*(data->Vs[BX2s](i,jm,k+1) - data->Vs[BX2s](i,jm,k-1));  
    )

    // Compute 2nd derivative correction

    DIM_EXPAND(                                          
                                                                            ,
      d2ByR_dx2 = data->Vs[BX2s](i+1,jp,k) - 2.0*data->Vs[BX2s](i,jp,k) + data->Vs[BX2s](i-1,jp,k);
      d2ByL_dx2 = data->Vs[BX2s](i+1,jm,k) - 2.0*data->Vs[BX2s](i,jm,k) + data->Vs[BX2s](i-1,jm,k);    

      d2BxR_dy2 = data->Vs[BX1s](ip,j+1,k) - 2.0*data->Vs[BX1s](ip,j,k) + data->Vs[BX1s](ip,j-1,k);
      d2BxL_dy2 = data->Vs[BX1s](im,j+1,k) - 2.0*data->Vs[BX1s](im,j,k) + data->Vs[BX1s](im,j-1,k);
                                                                            ,
      d2BzR_dx2 = data->Vs[BX3s](i+1,j,kp) - 2.0*data->Vs[BX3s](i,j,kp) + data->Vs[BX3s](i-1,j,kp);
      d2BzL_dx2 = data->Vs[BX3s](i+1,j,km) - 2.0*data->Vs[BX3s](i,j,km) + data->Vs[BX3s](i-1,j,km);

      d2BzR_dy2 = data->Vs[BX3s](i,j+1,kp) - 2.0*data->Vs[BX3s](i,j,kp) + data->Vs[BX3s](i,j-1,kp);
      d2BzL_dy2 = data->Vs[BX3s](i,j+1,km) - 2.0*data->Vs[BX3s](i,j,km) + data->Vs[BX3s](i,j-1,km);

      d2BxR_dz2 = data->Vs[BX1s](ip,j,k+1) - 2.0*data->Vs[BX1s](ip,j,k) + data->Vs[BX1s](ip,j,k-1);
      d2BxL_dz2 = data->Vs[BX1s](im,j,k+1) - 2.0*data->Vs[BX1s](im,j,k) + data->Vs[BX1s](im,j,k-1);

      d2ByR_dz2 = data->Vs[BX2s](i,jp,k+1) - 2.0*data->Vs[BX2s](i,jp,k) + data->Vs[BX2s](i,jp,k-1);
      d2ByL_dz2 = data->Vs[BX2s](i,jm,k+1) - 2.0*data->Vs[BX2s](i,jm,k) + data->Vs[BX2s](i,jm,k-1);
    )

    double dBt  = DIM_EXPAND(               0.0   ,
                            + dByR_dx - dByL_dx   , 
                            + dBzR_dx - dBzL_dx);
    double d2Bn = DIM_EXPAND(  0.0                    ,
                             + d2BxR_dy2 + d2BxL_dy2  , 
                             + d2BxR_dz2 + d2BxL_dz2);

    DIM_EXPAND(
      data->Upc(i,j,k,BX1) =   0.5 *(data->Vs[BX1s](ip,j,k) + data->Vs[BX1s](im,j,k)) 
                          + 0.125*dBt - d2Bn/48.0;  ,

      dBt  = DIM_EXPAND(                0.0  , 
                        + dBxR_dy - dBxL_dy  , 
                        + dBzR_dy - dBzL_dy);
      d2Bn = DIM_EXPAND(                  0.0  ,
                        + d2ByR_dx2 + d2ByL_dx2, 
                        + d2ByR_dz2 + d2ByL_dz2);
      data->Upc(i,j,k,BX2) =   0.5 *(data->Vs[BX2s](i,jp,k) + data->Vs[BX2s](i,jm,k)) 
                          + 0.125*dBt - d2Bn/48.0;      ,

      dBt  = dBxR_dz   + dByR_dz   - dBxL_dz   - dByL_dz;
      d2Bn = d2BzR_dx2 + d2BzL_dx2 + d2BzR_dy2 + d2BzL_dy2;
      data->Upc(i,j,k,BX3) =   0.5 *(data->Vs[BX3s](i,j,kp) + data->Vs[BX3s](i,j,km)) 
                          + 0.125*dBt - d2Bn/48.0;
    )

    #else 
      #error Invalid HO_FIELD_AVERAGE
    #endif

    #if (HO_LAP_LIMITER != NO) || (HO_DEAVERAGE == 0)
    if (data->flag(i,j,k) & FLAG_HO_LAP_LIMITER || (HO_DEAVERAGE == 0)){
      DIM_EXPAND(
        data->Upc(i,j,k,BX1) = 0.5*(data->Vs[BX1s](i,j,k) + data->Vs[BX1s](i+1,j,k)); ,
        data->Upc(i,j,k,BX2) = 0.5*(data->Vs[BX2s](i,j,k) + data->Vs[BX2s](i,j+1,k)); ,
        data->Upc(i,j,k,BX3) = 0.5*(data->Vs[BX3s](i,j,k) + data->Vs[BX3s](i,j,k+1));
      )
      #if (PHYSICS == ResRMHD) && (DIVE_CONTROL == CONSTRAINED_TRANSPORT)
      DIM_EXPAND(
        data->Upc(i,j,k,EX1) = 0.5*(data->Vs[EX1s](i,j,k) + data->Vs[EX1s](i+1,j,k)); ,
        data->Upc(i,j,k,EX2) = 0.5*(data->Vs[EX2s](i,j,k) + data->Vs[EX2s](i,j+1,k)); ,
        data->Upc(i,j,k,EX3) = 0.5*(data->Vs[EX3s](i,j,k) + data->Vs[EX3s](i,j,k+1));
      )
      #endif
    }
    #endif

    #endif // STAGGERED_MHD

    // 3e. Convert point values Upc --> Vpc

    {
      uint32_t err_loc = 0;
      double upc[NVAR], vpc[NVAR];
      State state;
      STATE_SET(state, vpc, upc, NULL, x1[i], x2[j], x3[k]);
      state.flag = &(data->flag(i,j,k));

      #pragma acc loop seq
      NVAR_LOOP(nv) upc[nv] = data->Upc(i,j,k,nv);

      err_loc = ConsToPrim (state, run_config);

      // data->err(i,j,k) |= err_loc;  // Do not flag error if a solution is 
                                       // is still possible

      // Conversion did not go through ? --> Use Uc instead of Upc

      if (err_loc){
        err_loc = 0;

        #pragma acc loop seq
        NVAR_LOOP(nv) upc[nv] = data->Uc(i,j,k,nv);

        err_loc = ConsToPrim (state, run_config);
        data->err(i,j,k) |= err_loc;
 
        // Update err & flag only using the latest value 

        if (err_loc){
          ATOMIC_ERR_UPDATE(err_kernel, err_loc, ERR_POINT_VALUE_FAIL);
          KERNEL_ERR_UPDATE(data->err(i,j,k), ERR_POINT_VALUE_FAIL);
          if (err_loc & ERR_CONS2PRIM_FAIL) {  // Need to quit ?
            ATOMIC_ERR_UPDATE(err_kernel, ERR_QUIT);
            KERNEL_ERR_UPDATE(data->err(i,j,k), ERR_QUIT);
          }
        }

        #pragma acc loop seq
        NVAR_LOOP(nv) data->Uc(i,j,k,nv) = upc[nv];
      }

      // Copy local values to main data array

      #pragma acc loop seq
      NVAR_LOOP(nv){
        data->Vpc(i,j,k,nv) = vpc[nv];  
        data->Upc(i,j,k,nv) = upc[nv];
      } 
    }

  }}} // end loop on i,j,k 

  // ------------------------------------------------------
  // 5. Error check
  // ------------------------------------------------------

  #if WARNING_MESSAGES > 0 
  Where ("PointValue()", data, err_kernel, box, grid);
  #endif

  #if WARNING_MESSAGES >= 0 
  int quit_condition = (err_kernel & ERR_QUIT);

  #if FALLBACK_MODE != NO
  quit_condition *= (data->integration_mode != DEFAULT);
  Fallback_Trigger(data->err, data->flag, box);
  #endif

  if (quit_condition) QuitPLUTO("from PointValue()",err_kernel);
  #endif

  mynvtxstop_();
}
