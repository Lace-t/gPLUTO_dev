////////////////////////////////////////////////////////////////////////
/*!
  \file
  \brief Average staggered magnetic field to zone center.

  \author A. Mignone (andrea.mignone@unito.it)
  \date   Nov 08, 2025
*/
////////////////////////////////////////////////////////////////////////
#include "gpluto.hpp"

// Similar to DOT_PRODUCT, but only for the included dimensions.
// The energy correction step requires only to replace cell-centered
// magnetic field with averaged staggered fields that exist
// only if the corresponding dimension is enabled.
#define DOT_PRODUCT_DIM(A,B) DIM_EXPAND(  (A)[IDIR]*(B)[IDIR], \
                                        + (A)[JDIR]*(B)[JDIR], \
                                        + (A)[KDIR]*(B)[KDIR]  )

////////////////////////////////////////////////////////////////////////
/*!
  Average staggered magnetic field components to form a zone-centered
  field, e.g., \f$ \av{B_i} = (B_{i+1/2} + B_{i-1/2})/2\f$.
  In cylindrical coordinates the volume-averaged radial component
  is obtained as
  \f[
  \av{B_R} = \frac{1}{\Delta V} \int B(R) R dR
  \qquad\mathrm{where}\qquad
  B(R) =   B_{R,i+\HALF} \frac{R - R_{i-\HALF}}{\Delta R}
  + B_{R,i-\HALF} \frac{R_{i+\HALF} - R}{\Delta R}
  \f]
  this yields
  \f[
  \av{B_R} =  \frac{R_{i+\HALF} + 2R_{i-\HALF}}{3(R_{i+\HALF}+R_{i-\HALF})}
  B_{R,i-\HALF}
  + \frac{2R_{i+\HALF} + R_{i-\HALF}}{3(R_{i+\HALF}+R_{i-\HALF})}
  B_{R,i-\HALF}
  \f]
  Simlar expressions hold in spherical coordinates.
  Check with the following MAPLE script:
  \code
   restart;
   J     := sin(xi); # J = xi, xi^2, sin(xi)
   Br    := (xi - a)/Delta*B[R] + (b - xi)/Delta*B[L];
   Bav   := int (Br*J, xi=a..b) / int (J, xi=a..b):
   Delta := b-a;
   Bav   := simplify(Bav);
   cp := coeff(Bav, B[R]);
   cm := coeff(Bav, B[L]);

   # Check symmetry propetries (sign must not change)

   eval([cm,cp],{a=0,b=1});
   eval([cm,cp],{a=-1,b=0});
  \endcode
  
  The averaging is done inside an arbitrary rectangular box.
  
  When the CT_ENERGY_CORRECTION flag is enabled, we also redefine the
  zone total energy using the newly formed cell-centered field, i.e.,
  
  \f[
  E_i \to E_i - \frac{\mathbf{B}_i^2}{2} + \frac{<\mathbf{B}_i^2>}{2}
  \f]
  
  \param [in]   d            pointer to PLUTO data strdata->Ucture
  \param [in]   ct_eng_corr  a flag to enable or disable the
                              CT energy correction step
  \param [in]   box          pointer to RBox structure
  \param [in]   grid         pointer to Grid structure
*/
////////////////////////////////////////////////////////////////////////
void CT_AverageStaggeredFields (Data *data, int ct_eng_corr, RBox *box, Grid *grid)
{
  int    i, j, k, nv;
  #if (GEOMETRY != CARTESIAN) || (PHYSICS != GRMHD)
  double *x1r = grid->xR[IDIR], *x1l = grid->xL[IDIR];
  #endif
  #if GEOMETRY == SPHERICAL
  double dp, dm, thp, thm;
  double *x2r = grid->xR[JDIR], *x2l = grid->xL[JDIR];
  #elif PHYSICS == GRMHD
  double *x1 = grid->xC[IDIR];
  double *x2 = grid->xC[JDIR];
  double *x3 = grid->xC[KDIR];
  double *x1r = grid->xR[IDIR], *x1l = grid->xL[IDIR];
  double *x2r = grid->xR[JDIR], *x2l = grid->xL[JDIR];
  double *x3r = grid->xR[KDIR], *x3l = grid->xL[KDIR];
  #endif
  double u[NVAR];
  #if PHYSICS == RMHD && (PCP_NR == YES || PCP_FOFC == YES || FALLBACK_MODE == FALLBACK_RECOMPUTE || CT_RMHD_CONSISTENT == YES)
  RunConfig &rc = *data->run_config;
  #endif
  static Ary4D Bc0, Ec0;

  // 0. Allocate extra memory for Bc (RMHD)

  #if (PHYSICS == RMHD) || (PHYSICS == GRMHD) || (PHYSICS == ResRMHD)
  int    nx1_tot = grid->ntot[IDIR];
  int    nx2_tot = grid->ntot[JDIR];
  int    nx3_tot = grid->ntot[KDIR];
  
  if (Bc0.isAllocated() == 0) {
    Bc0.create(nx1_tot,nx2_tot,nx3_tot,3);
    #pragma acc enter data copyin(Bc0)
    #pragma omp target enter data map(to: Bc0)
  }
  #endif

  #if (PHYSICS == ResRMHD) && (DIVE_CONTROL == CONSTRAINED_TRANSPORT)
  if (Ec0.isAllocated() == 0) {
    Ec0.create(nx1_tot,nx2_tot,nx3_tot,3);
    #pragma acc enter data copyin(Ec0)
    #pragma omp target enter data map(to: Ec0)
  }
  #endif

  // ------------------------------------------------------
  // 1. Loop over all zones in the given box.
  // ------------------------------------------------------

  mynvtxstart_("CT_AverageStaggeredFields", SILVER);

  #pragma omp target teams distribute parallel for collapse(3) private(u) \
              map(to:box[:1])
  #pragma acc parallel loop collapse(3) \
              present(grid, data) copyin(box,Bc0,Ec0) private(u[:NVAR])
  BOX_LOOP(*box, i,j,k){

    double Bf[3];
    #if (PHYSICS == ResRMHD) && (DIVE_CONTROL == CONSTRAINED_TRANSPORT)
    double Ef[3];   // Zone-averaged fields  
    #endif
    #if CT_ENERGY_CORRECTION == 1 && !(PHYSICS == RMHD && CT_RMHD_CONSISTENT)
    double Bc[3];
    #if PHYSICS == ResRMHD
    double Ec[3];   // Zone-centered fields
    #endif
    double B2c, B2f;
    #endif

    #if GEOMETRY != CARTESIAN
    Metric metric;
    MetricSet (metric, grid->rC[i], grid->sC[j]);
    #endif

    // 1a. Initialize staggered fields to zone-center values
    //     (e.g. Bf = Bc) before overwriting

    if (ct_eng_corr){
      #pragma acc loop seq
      NVAR_LOOP(nv) u[nv] = data->Uc(i,j,k,nv);

      #if GEOMETRY != CARTESIAN
      ConsToOrtho(u, metric);
      #endif
      Bf[IDIR] = u[BX1]; Bf[JDIR] = u[BX2]; Bf[KDIR] = u[BX3];
      #if CT_ENERGY_CORRECTION == 1 && !(PHYSICS == RMHD && CT_RMHD_CONSISTENT)
      Bc[IDIR] = u[BX1]; Bc[JDIR] = u[BX2]; Bc[KDIR] = u[BX3];
      #if PHYSICS == ResRMHD
      Ef[IDIR] = u[EX1]; Ef[JDIR] = u[EX2]; Ef[KDIR] = u[EX3]; 
      Ec[IDIR] = u[EX1]; Ec[JDIR] = u[EX2]; Ec[KDIR] = u[EX3];
      #endif
      #endif
    }

    // 1b. Average staggered fields to zone center [no GRMHD]
    //     For non-Cartesian geometries, the average process includes 
    //     some geometrical factors.

    #if ((GEOMETRY == CARTESIAN) && (PHYSICS != GRMHD)) || \
        ((PHYSICS == GRMHD) && (GR_METRIC == CARTESIAN))

    DIM_EXPAND(Bf[IDIR] = 0.5*(data->Vs[BX1s](i,j,k) + data->Vs[BX1s](i+1,j,k));  ,
               Bf[JDIR] = 0.5*(data->Vs[BX2s](i,j,k) + data->Vs[BX2s](i,j+1,k));  ,
               Bf[KDIR] = 0.5*(data->Vs[BX3s](i,j,k) + data->Vs[BX3s](i,j,k+1));)
    #if (PHYSICS == ResRMHD) && (DIVE_CONTROL == CONSTRAINED_TRANSPORT)
    DIM_EXPAND(Ef[IDIR] = 0.5*(data->Vs[EX1s](i,j,k) + data->Vs[EX1s](i+1,j,k));  ,
               Ef[JDIR] = 0.5*(data->Vs[EX2s](i,j,k) + data->Vs[EX2s](i,j+1,k));  ,
               Ef[KDIR] = 0.5*(data->Vs[EX3s](i,j,k) + data->Vs[EX3s](i,j,k+1));)
    #endif

    #elif GEOMETRY == CYLINDRICAL
    double cp, cm, rp, rm;

    rp = x1r[i];
    rm = x1l[i];
    cm = (rp + 2.0*rm)/(3.0*(rp + rm));
    cp = 1.0 - cm;
    DIM_EXPAND(Bf[IDIR] = cm*data->Vs[BX1s](i,j,k) + cp*data->Vs[BX1s](i+1,j,k);   ,
               Bf[JDIR] = 0.5*(data->Vs[BX2s](i,j,k) + data->Vs[BX2s](i,j+1,k));   ,
               Bf[KDIR] = 0.5*(data->Vs[BX3s](i,j,k) + data->Vs[BX3s](i,j,k+1));)

    #if (PHYSICS == ResRMHD) && (DIVE_CONTROL == CONSTRAINED_TRANSPORT)
    DIM_EXPAND(Ef[IDIR] = cm*data->Vs[EX1s](i,j,k) + cp*data->Vs[EX1s](i+1,j,k);   ,
               Ef[JDIR] = 0.5*(data->Vs[EX2s](i,j,k) + data->Vs[EX2s](i,j+1,k));   ,
               Ef[KDIR] = 0.5*(data->Vs[EX3s](i,j,k) + data->Vs[EX3s](i,j,k+1));)
    #endif
    

    #elif GEOMETRY == SPHERICAL
    double cp, cm, rp, rm;

    rp = x1r[i];
    rm = x1l[i];

    cm =  (3.0*rm*rm + 2.0*rm*rp + rp*rp)/(4.0*(rm*rm + rm*rp + rp*rp));
    cp = 1.0 - cm;

    thp = x2r[j];
    thm = x2l[j];

    dm  = (thm - thp)*cos(thm) - sin(thm) + sin(thp);
    dm /= (thm - thp)*(cos(thm) - cos(thp));
    dp  = 1.0 - dm;
    DIM_EXPAND(Bf[IDIR]  = cm*data->Vs[BX1s](i,j,k) + cp*data->Vs[BX1s](i+1,j,k);  ,
               Bf[JDIR]  = dm*data->Vs[BX2s](i,j,k) + dp*data->Vs[BX2s](i,j+1,k);  ,
               Bf[KDIR]  = 0.5*(data->Vs[BX3s](i,j,k) + data->Vs[BX3s](i,j,k+1));)

    #if (PHYSICS == ResRMHD) && (DIVE_CONTROL == CONSTRAINED_TRANSPORT)
    DIM_EXPAND(Ef[IDIR]  = cm*data->Vs[EX1s](i,j,k) + cp*data->Vs[EX1s](i+1,j,k);  ,
               Ef[JDIR]  = dm*data->Vs[EX2s](i,j,k) + dp*data->Vs[EX2s](i,j+1,k);  ,
               Ef[KDIR]  = 0.5*(data->Vs[EX3s](i,j,k) + data->Vs[EX3s](i,j,k+1));)    
    #endif

    #else // PHYSICS == GRMHD

    // To be optimized loading it to the Grid structrure
    // 1c. Average staggered fields to zone center [GRMHD]

    #if GR_FD
    // Not tested, probably missing some factors... 
    Metric Gl = DEFAULT_METRIC;
    Metric Gr = DEFAULT_METRIC;

    if (!SingularMetric(x1l[i],x2[j],x3[k])){
      ComputeLocalMetric(x1l[i], x2[j], x3[k], &Gl);
    }
    if (!SingularMetric(x1r[i],x2[j],x3[k])){
      ComputeLocalMetric(x1r[i], x2[j], x3[k], &Gr);
    }
    double cm = Gl.sqrtDet;
    double cp = Gr.sqrtDet;
    Bf[IDIR] = 0.5*(data->Vs[BX1s](i,j,k) + data->Vs[BX1s](i+1,j,k));

    #if INCLUDE_JDIR
    if (!SingularMetric(x1[i],x2l[j],x3[k])){
      ComputeLocalMetric(x1[i], x2l[j], x3[k], &Gl);
    }
    if (!SingularMetric(x1[i],x2r[j],x3[k])){
      ComputeLocalMetric(x1[i], x2r[j], x3[k], &Gr);
    }
    double cm = Gl.sqrtDet;
    double cp = Gr.sqrtDet;
    Bf[JDIR] = 0.5*(data->Vs[BX2s](i,j,k) + data->Vs[BX2s](i,j+1,k));
    #endif

    #if INCLUDE_KDIR
    if (!SingularMetric(x1[i],x2[j],x3l[k])){
      ComputeLocalMetric(x1[i], x2[j], x3l[k], &Gl);
    }
    if (!SingularMetric(x1[i],x2[j],x3r[k])){
      ComputeLocalMetric(x1[i], x2[j], x3r[k], &Gr);
    }
    double cm = Gl.sqrtDet;
    double cp = Gr.sqrtDet;
    Bf[KDIR] = 0.5*(data->Vs[BX3s](i,j,k) + data->Vs[BX3s](i,j,k+1));
    #endif

    #else

    // Compute extra metric factor 
    double x1b = grid->xgc[0](i,j,k);
    double dx1 = x1l[i] - x1r[i];
    double scrhm = (x1l[i] - x1b)/dx1;
    double scrhp = (x1b - x1r[i])/dx1;
    Bf[IDIR] = data->Vs[BX1s](i,j,k)*scrhm + data->Vs[BX1s](i+1,j,k)*scrhp;

    #if INCLUDE_JDIR
    double x2b = grid->xgc[1](i,j,k);
    double dx2 = x2l[j] - x2r[j];
    scrhm = (x2l[j] - x2b)/dx2;
    scrhp = (x2b - x2r[j])/dx2;
    Bf[JDIR] = data->Vs[BX2s](i,j,k)*scrhm + data->Vs[BX2s](i,j+1,k)*scrhp;
    #endif

    #if INCLUDE_KDIR
    double x3b = grid->xgc[2](i,j,k);
    double dx3 = x3l[k] - x3r[k];
    scrhm = (x3l[k] - x3b)/dx3;
    scrhp = (x3b - x3r[k])/dx3;
    Bf[KDIR] = data->Vs[BX3s](i,j,k)*scrhm + data->Vs[BX3s](i,j,k+1)*scrhp;
    #endif

    #endif
    #endif // (GEOMETRY == CARTESIAN) && (PHYSICS !=GRMHD) || ...

    #if PHYSICS == RMHD && CT_RMHD_CONSISTENT
    if (ct_eng_corr) {
      double prim[NVAR];
      if (CT_RMHD_CorrectCell(u, prim, Bf, rc)) {
        data->flag(i,j,k) |= FLAG_PCP_CT_FAIL;
      } else {
        #if GEOMETRY != CARTESIAN
        OrthoToCons(u,metric);
        #endif
        data->Uc(i,j,k,MX1) = u[MX1];
        data->Uc(i,j,k,MX2) = u[MX2];
        data->Uc(i,j,k,MX3) = u[MX3];
        data->Uc(i,j,k,ENG) = u[ENG];
        #if GEOMETRY != CARTESIAN
        ConsToOrtho(u,metric);
        #endif
      }
    }
    #endif

    // 1d. Apply energy correction if required. 
    //     Please note that, in ResRMHD, 
    //     Ef != Ec only when DIVE_CONTROL == CONSTRAINED_TRANSPORT

    #if CT_ENERGY_CORRECTION == 1 && !(PHYSICS == RMHD && CT_RMHD_CONSISTENT)
    if (ct_eng_corr){
    
      B2c = DOT_PRODUCT_DIM(Bc, Bc);
      B2f = DOT_PRODUCT_DIM(Bf, Bf);

      #if (PHYSICS == ResRMHD)
      {
        double Sc[3], Sf[3];

        #if DIVE_CONTROL == CONSTRAINED_TRANSPORT
        B2c += DOT_PRODUCT_DIM(Ec, Ec);
        B2f += DOT_PRODUCT_DIM(Ef, Ef);
        #endif

        CROSS_PRODUCT(Sc, Ec, Bc);
        CROSS_PRODUCT(Sf, Ef, Bf);

        data->Uc(i,j,k,MX1) += - Sc[IDIR] + Sf[IDIR];
        data->Uc(i,j,k,MX2) += - Sc[JDIR] + Sf[JDIR];
        data->Uc(i,j,k,MX3) += - Sc[KDIR] + Sf[KDIR];
      }
      #endif

      data->Uc(i,j,k,ENG) += 0.5*(B2f - B2c);

      // 1d1. For RMHD, save orthonormal field 
      //      before overwriting

      #if (PHYSICS == RMHD) || (PHYSICS == GRMHD) || (PHYSICS == ResRMHD)
      Bc0(i,j,k,0) = Bc[IDIR];
      Bc0(i,j,k,1) = Bc[JDIR];
      Bc0(i,j,k,2) = Bc[KDIR];
      #endif
      #if (PHYSICS == ResRMHD) && (DIVE_CONTROL == CONSTRAINED_TRANSPORT)
      Ec0(i,j,k,0) = Ec[IDIR];
      Ec0(i,j,k,1) = Ec[JDIR];
      Ec0(i,j,k,2) = Ec[KDIR];
      #endif
    }
    #endif // CT_ENERGY_CORRECTION == 1

    // 1e. Replace zone-centered field with staggered average

    DIM_EXPAND(u[BX1] = Bf[IDIR];  ,
               u[BX2] = Bf[JDIR];  ,
               u[BX3] = Bf[KDIR];)
    #if (PHYSICS == ResRMHD) && (DIVE_CONTROL == CONSTRAINED_TRANSPORT)
    DIM_EXPAND(u[EX1] = Ef[IDIR];  ,
               u[EX2] = Ef[JDIR];  ,
               u[EX3] = Ef[KDIR];)
    #endif

    #if GEOMETRY != CARTESIAN
    OrthoToCons(u, metric);
    #endif

    DIM_EXPAND(data->Uc(i,j,k,BX1) = u[BX1];  ,
               data->Uc(i,j,k,BX2) = u[BX2];  ,
               data->Uc(i,j,k,BX3) = u[BX3]; )

    #if (PHYSICS == ResRMHD) && (DIVE_CONTROL == CONSTRAINED_TRANSPORT)
    DIM_EXPAND(data->Uc(i,j,k,EX1) = u[EX1];  ,
               data->Uc(i,j,k,EX2) = u[EX2];  ,
               data->Uc(i,j,k,EX3) = u[EX3]; )
    #endif

  } // End BOX_LOOP()

mynvtxstop_();

  if (!ct_eng_corr) return;

  // 2. Implement the extra correction step as in Marti (2015),
  //    to complete the relativistic correction of momentum
  //    & energy.

  #if CT_ENERGY_CORRECTION == 1 && !(PHYSICS == RMHD && CT_RMHD_CONSISTENT)
  #if (PHYSICS == RMHD) || (PHYSICS == GRMHD)
  #if !(PHYSICS == RMHD && (PCP_NR == YES || PCP_FOFC == YES || FALLBACK_MODE == FALLBACK_RECOMPUTE || CT_RMHD_CONSISTENT == YES))
  ConsToPrim3D (data, 0, grid, box);
  #endif

  #pragma omp target teams distribute parallel for collapse(3) map(to:box[:1])
  #pragma acc parallel loop collapse(3) present(grid, data, Bc0) copyin(box)
  BOX_LOOP(*box, i,j,k){
    #if GEOMETRY != CARTESIAN
    Metric metric;
    MetricSet (metric, grid->rC[i], grid->sC[j]);
    #endif

    double vc[NVAR], uc[NVAR], Bc[3];

    Bc[IDIR] = Bc0(i,j,k,0); // Bc0 was saved as orthonormal
    Bc[JDIR] = Bc0(i,j,k,1);
    Bc[KDIR] = Bc0(i,j,k,2); 

    #pragma acc loop seq
    NVAR_LOOP(nv){
      vc[nv] = data->Vc(i,j,k,nv);
      uc[nv] = data->Uc(i,j,k,nv);
    }
    #if GEOMETRY != CARTESIAN
    ConsToOrtho(uc, metric);
    #endif

    #if PHYSICS == RMHD && (PCP_NR == YES || PCP_FOFC == YES || FALLBACK_MODE == FALLBACK_RECOMPUTE || CT_RMHD_CONSISTENT == YES)
    if (RMHD_PCP_Recover(uc, vc, rc, 0)) {
      data->flag(i,j,k) |= FLAG_PCP_CT_FAIL;
      continue;
    }
    #endif

    double Bc2 = DOT_PRODUCT(Bc, Bc);
    double Bf2 = DOT_PRODUCT(vc + BX1, vc + BX1);
    double vBc = DOT_PRODUCT(vc + VX1, Bc);
    double vBf = DOT_PRODUCT(vc + VX1, vc + BX1);
    double v2  = DOT_PRODUCT(vc + VX1, vc + VX1);

    uc[MX1] += -(Bc2 - Bf2)*vc[VX1] + Bc[0]*vBc - vc[BX1]*vBf;
    uc[MX2] += -(Bc2 - Bf2)*vc[VX2] + Bc[1]*vBc - vc[BX2]*vBf;
    uc[MX3] += -(Bc2 - Bf2)*vc[VX3] + Bc[2]*vBc - vc[BX3]*vBf;
    uc[ENG] += -0.5*v2*(Bc2 - Bf2) + 0.5*(vBc*vBc - vBf*vBf);

    #if GEOMETRY != CARTESIAN
    OrthoToCons(uc, metric);
    #endif
    #pragma acc loop seq
    NVAR_LOOP(nv) data->Uc(i,j,k,nv) = uc[nv];

  }
  #endif   // CT_ENERGY_CORRECTION == 1
  #endif   // (PHYSICS == RMHD) || (PHYSICS == GRMHD)

  #if (PHYSICS == MHD) && (GEOMETRY == CARTESIAN)
  #endif
}
