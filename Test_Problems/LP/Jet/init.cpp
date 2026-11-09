////////////////////////////////////////////////////////////////////////
/*!
  \file
  \brief 3D Jet propagation in a stratified halo atmosphere

  This problem models the injection of a relativistic/magnetized jet
  into a hydrostatic ambient medium (halo). The halo is constructed by
  solving the hydrostatic equilibrium equation in a fixed gravitational
  potential, assuming a prescribed temperature profile.

  The jet is injected at the X3_BEG boundary and smoothly blended with
  the ambient medium using radial profile functions and time-dependent
  ramping factors.

  ------------------------------------------------------------------------
  Halo model:

    - Density and pressure are obtained from hydrostatic equilibrium:
          dP/dr = -rho * dPhi/dr
    - Temperature profile: T(r) = Tc + (1 - 1/cosh(r/rc))*(Th - Tc)
    - Gravitational field: combination of bulge + dark matter halo
    - Normalization set by HALO_N0 (central number density)

  ------------------------------------------------------------------------
  Jet model:

    - Injection at X3_BEG boundary
    - Jet density:   rho_j = ETA * rho_halo(0)
    - Jet pressure:   JET_PRS (in cgs units)
    - Lorentz factor: JET_GAMMA
    - Magnetization:  JET_SIGMAB
    - Opening angle (optional): JET_ANG
    - Radial profile controlled by Profile(r, nv)

  ------------------------------------------------------------------------
  Time modulation:
    - Jet activation: JET_TON
    - Jet shutoff:    JET_TOFF
    - Ramp-up/down:   T_RAMP_UP, T_RAMP_DN

  ------------------------------------------------------------------------
  Units:
    - Length unit: UNIT_LENGTH
    - Density and pressure normalized using run_config units

  ------------------------------------------------------------------------
  Boundary conditions:
    - X1/X2/X3_END: outflow + halo refill if inflow detected
    - X3_BEG: jet injection + reflected outflow
    - Interior: density/pressure floors enforced


  \author A. Suriano \n
          A. Mignone (andrea.mignone@unito.it)

  \date   April 15, 2026
*/
////////////////////////////////////////////////////////////////////////
#include "gpluto.hpp"
#include "halo.hpp"
#include "jet.hpp"

// funtions for the jet
#pragma acc routine seq
double Profile(double, int, const RunConfig &);

////////////////////////////////////////////////////////////////////////
void Init (double *v, double x1, double x2, double x3, RunConfig &run_config)
/*
*/
////////////////////////////////////////////////////////////////////////
{
  run_config.unit_length = UNIT_LENGTH;
  run_config.small_rho = 5.e-6;
  run_config.small_prs = 5.e-9;

  v[VX1] = 0.0;
  v[VX2] = 0.0;
  v[VX3] = 0.0;

  v[BX1] = 0.0;
  v[BX2] = 0.0;
  v[BX3] = 0.0;

  v[TRC] = 0.0;
  v[TRC + 1] = 0.0;

}

////////////////////////////////////////////////////////////////////////
void InitDomain(Data *data, Grid *grid)
/*
  Assign initial condition by looping over the computational domain.
  Called after the usual Init() function to assign initial conditions
  on primitive variables.
  Value assigned here will overwrite those prescribed during Init().  
*/
////////////////////////////////////////////////////////////////////////
{
  const int nx1_tot = grid->ntot[IDIR];
  const int nx2_tot = grid->ntot[JDIR];
  const int nx3_tot = grid->ntot[KDIR];

  haloStruct halo_struct;
  RunConfig &run_config = *(data->run_config);

  HaloInit(&halo_struct, run_config);

  #pragma acc data copyin(halo_struct) 
  {
  #pragma acc parallel loop collapse(3) present(data, grid)
  for(int k = 0; k < nx3_tot; k++){
  for(int j = 0; j < nx2_tot; j++){
  for(int i = 0; i < nx1_tot; i++){
    const double x1 = grid->xC[IDIR][i], 
                 x2 = grid->xC[JDIR][j], 
                 x3 = grid->xC[KDIR][k];
    const double r = sqrt(x1*x1 + x2*x2 + x3*x3);
    double rho, prs;
    HaloInterp(&halo_struct, r, rho, prs);

    data->Vc(i,j,k,RHO) = rho;
    data->Vc(i,j,k,PRS) = prs;
    data->Vc(i,j,k,VX1) = 0.0;
    data->Vc(i,j,k,VX2) = 0.0;
    data->Vc(i,j,k,VX3) = 0.0;
    data->Vc(i,j,k,BX1) = 0.0;
    data->Vc(i,j,k,BX2) = 0.0;
    data->Vc(i,j,k,BX3) = 0.0;
    data->Vc(i,j,k,TRC) = 0.0;
    data->Vc(i,j,k,TRC + 1) = 0.0;
  }}}
  }
}

////////////////////////////////////////////////////////////////////////
void Analysis(Data *data, Grid *grid)
/*
  Perform runtime data analysis.
  
  \param [in] d the PLUTO Data structure
  \param [in] grid   pointer to array of Grid structures
*/
////////////////////////////////////////////////////////////////////////
{
}

////////////////////////////////////////////////////////////////////////
void UserDefBoundary(const Data *data, RBox *box, int side, Grid *grid)
/*
  Assign user-defined boundary conditions.
  
  \param [in,out] data  pointer to the PLUTO data structure
  \param [in]     box   pointer to a RBox structure containing boundary region
  \param [in]     side  boundary side where ghost zones need to be filled
  \param [in]     grid  pointer to an array of Grid structures.
*/
////////////////////////////////////////////////////////////////////////
{
  int i, j, k, nv;
  static int first_call = 1;
  const int ibeg = grid->lbeg[IDIR], iend = grid->lend[IDIR];
  const int jbeg = grid->lbeg[JDIR], jend = grid->lend[JDIR];
  const int kbeg = grid->lbeg[KDIR], kend = grid->lend[KDIR];

  double *x1 = grid->xC[IDIR], *x1L = grid->xL[IDIR];
  double *x2 = grid->xC[JDIR], *x2L = grid->xL[JDIR];
  double *x3 = grid->xC[KDIR];

  RunConfig &run_config = *(data->run_config);
  double inv_rj = 1.0 / (run_config.inputParam[JET_RAD]);

  static haloStruct halo_struct;
  static jetStruct  jet_struct;

  // ------------------------------------------------------
  // 0. Initialize halo & jet values
  // ------------------------------------------------------

  if (first_call) {
    double rhoh, prsh;
    HaloInit(&halo_struct, run_config);
    HaloInterp(&halo_struct, 0.0, rhoh, prsh);

    JetInit (&jet_struct,  run_config, rhoh,prsh);
    #pragma acc enter data copyin(halo_struct, jet_struct) 

    first_call = 0;
  }

  // ------------------------------------------------------
  // 1. Set boundary conditions: X1_BEG & X1_END 
  // ------------------------------------------------------

  if (side == X1_BEG || side == X1_END) { 
    int i1 = (side == X1_BEG ? ibeg:iend);
    int n  = (side == X1_BEG ? 1:-1);
    #pragma acc parallel loop collapse(3) present(box, grid, data, halo_struct) async(1)
    BOX_LOOP(*box, i,j,k) {
      const double vn = data->Vc(i1,j,k,VX1)*n;
      #pragma acc loop seq
      NVAR_LOOP(nv) data->Vc(i,j,k,nv) = data->Vc(i1,j,k,nv); // Outflow

      if (vn > 0.0) {
        const double r = sqrt(x1[i]*x1[i] + x2[j]*x2[j] + x3[k]*x3[k]);
        double rho, prs;
        HaloInterp(&halo_struct, r, rho, prs);

        data->Vc(i,j,k,RHO) = rho;
        data->Vc(i,j,k,PRS) = prs;
        data->Vc(i,j,k,VX1) = 0.0;
      } // if (vn > 0.0) 
    } // BOX_LOOP(*box, i,j,k) 

    #ifdef STAGGERED_MHD
      #if INCLUDE_JDIR
      #pragma acc parallel loop collapse(3) present(box, data) async(2)
      BOXJ_LOOP(*box, i,j,k) data->Vs[BX2s](i,j,k) = data->Vs[BX2s](i1,j,k);  
      #endif
      #if INCLUDE_KDIR
      #pragma acc parallel loop collapse(3) present(box, data) async(3)
      BOXK_LOOP(*box, i,j,k) data->Vs[BX3s](i,j, k) = data->Vs[BX3s](i1,j,k);
      #endif
    #endif
    #pragma acc wait
  } // if (side == X1_BEG || side == X1_END) 

  // ------------------------------------------------------
  // 2. Set boundary conditions: X2_BEG & X2_END 
  // ------------------------------------------------------

  if (side == X2_BEG || side == X2_END) { 
    int j1 = (side == X2_BEG ? jbeg:jend);
    int n  = (side == X2_BEG ? 1:-1);
    #pragma acc parallel loop collapse(3) present(box, grid, data, halo_struct) async(1)
    BOX_LOOP(*box, i,j,k) {
      const double vn = data->Vc(i,j1,k,VX2)*n;
      #pragma acc loop seq
      NVAR_LOOP(nv) data->Vc(i,j,k,nv) = data->Vc(i,j1,k,nv); // Outflow

      if (vn > 0.0) {
        const double r = sqrt(x1[i]*x1[i] + x2[j]*x2[j] + x3[k]*x3[k]);
        double rho, prs;
        HaloInterp(&halo_struct, r, rho, prs);

        data->Vc(i,j,k,RHO) = rho;
        data->Vc(i,j,k,PRS) = prs;
        data->Vc(i,j,k,VX2) = 0.0;
      } // if (vn > 0.0) 
    } // BOX_LOOP(*box, i,j,k) 

    #ifdef STAGGERED_MHD
      #if INCLUDE_IDIR
      #pragma acc parallel loop collapse(3) present(box, data) async(2)
      BOXI_LOOP(*box, i,j,k) data->Vs[BX1s](i,j,k) = data->Vs[BX1s](i,j1,k);  
      #endif
      #if INCLUDE_KDIR
      #pragma acc parallel loop collapse(3) present(box, data) async(3)
      BOXK_LOOP(*box, i,j,k) data->Vs[BX3s](i,j, k) = data->Vs[BX3s](i,j1,k);
      #endif
    #endif
    #pragma acc wait
  } // if (side == X2_BEG || side == X2_END) 

  // ------------------------------------------------------
  // 3a. Set boundary conditions: X3_BEG: Jet Launch
  // ------------------------------------------------------

  if (side == X3_BEG) { 
    const double t = run_config.time;
    double launch, dump;

    // 3a1. Set dump window

    const double t0    = run_config.inputParam[JET_TOFF];
    const double dt_dn = run_config.inputParam[T_RAMP_DN];
    if (t < t0) {
      dump = 1.0;
    } else if (t < t0 + dt_dn){
      dump = 1.0 / cosh(10.0 * (t0 - t) / (dt_dn));
    } else dump = 0.0;

    // 3a2. Set launch window

    const double t1    = run_config.inputParam[JET_TON];
    const double dt_up = run_config.inputParam[T_RAMP_UP];

    if (t < t1){
      launch = 0.0;
    } else if( t < t1 + dt_up) {
      launch = 1.0 / cosh(10.0 * (t - t1 - dt_up) / (dt_up));
    } else launch = 1.0;

    // 3a3. Assign b.c. on cell-centered variables

    #pragma acc parallel loop collapse(3) present(box, grid, data) async(1)
    BOX_LOOP(*box, i,j,k){

      double vout[NVAR], vjet[NVAR];
      #pragma acc cache(vout[:NVAR],vjet[:NVAR])  //[AM] What's this ?

      // Outflow portion
 
      #pragma acc loop seq
      NVAR_LOOP(nv) vout[nv] = data->Vc(i, j, 2*kbeg - k - 1, nv);

      // Reflect velocities and B
      vout[VX3] *= -1.0;
      vout[BX3] *= -1.0;

      // Compute profiles 
      const double r = sqrt(x1[i]*x1[i] + x2[j]*x2[j]);
      double prof[NVAR];
      #pragma acc loop seq
      NVAR_LOOP(nv) prof[nv] = Profile(r,nv,run_config)*dump*launch;

      // Assign jet values to inlet region, outflow outside
      JetValues(&jet_struct, vjet, x1[i], x2[j], x3[k], run_config);

      #pragma acc loop seq
      NVAR_LOOP(nv) {
        data->Vc(i,j,k,nv) = vout[nv] + (vjet[nv] - vout[nv])*prof[nv];
      } 

      #if PHYSICS == MHD || PHYSICS == RMHD
      data->Vc(i,j,k,BX1) = vout[BX1] + (vjet[BX1]*(r*inv_rj) - vout[BX1])*prof[BX1];
      data->Vc(i,j,k,BX2) = vout[BX2] + (vjet[BX2]*(r*inv_rj) - vout[BX2])*prof[BX2];
      data->Vc(i,j,k,BX3) = vout[BX3] + (vjet[BX3]*(r*inv_rj) - vout[BX3])*prof[BX3];
      #endif
    } // End BOX_LOOP()

    // 3a4. Assign b.c. on staggered fields

    #ifdef STAGGERED_MHD
    #if INCLUDE_IDIR
    #pragma acc parallel loop collapse(3) present(box, grid, data, run_config) async(2)
    BOXI_LOOP(*box,i,j,k) {
      double vjet[NVAR];
      const double r = sqrt(x1L[i]*x1L[i] + x2[j]*x2[j]);
      double prof = Profile(r, BX1, run_config)*dump*launch;
      const double vs = data->Vs[BX1s](i, j, 2*kbeg - k - 1);

      JetValues(&jet_struct, vjet, x1L[i], x2[j], x3[k], run_config);
      data->Vs[BX1s](i,j,k) = vs + (vjet[BX1]*(r*inv_rj) - vs)*prof;
    }     
    #endif // INCLUDE_IDIR 

    #if INCLUDE_JDIR
    #pragma acc parallel loop collapse(3) present(box, grid, data, run_config) async(3)
    BOXJ_LOOP(*box,i,j,k) {
      double vjet[NVAR];
      const double r = sqrt(x1[i]*x1[i] + x2L[j]*x2L[j]);
      double prof = Profile(r, BX2, run_config)*dump*launch;
      const double vs = data->Vs[BX2s](i, j, 2*kbeg - k - 1);

      JetValues(&jet_struct, vjet, x1[i], x2L[j], x3[k], run_config);
      data->Vs[BX2s](i,j,k) = vs + (vjet[BX2]*(r*inv_rj) - vs)*prof;
    }     
    #endif // INCLUDE_JDIR 
    #endif // #ifdef STAGGERED_MHD 
    #pragma acc wait
  }

  // ------------------------------------------------------
  // 3b. Set boundary conditions: X3_END
  // ------------------------------------------------------

  if (side == X3_END) { 
    #pragma acc parallel loop collapse(3) present(box, grid, data, halo_struct) async(1)
    BOX_LOOP(*box, i,j,k) {
      const double vn = -data->Vc(i,j,kend,VX3);
      #pragma acc loop seq
      NVAR_LOOP(nv) data->Vc(i,j,k,nv) = data->Vc(i,j,kend,nv); // Outflow

      if (vn > 0.0) {
        const double r = sqrt(x1[i]*x1[i] + x2[j]*x2[j] + x3[k]*x3[k]);
        double rho, prs;
        HaloInterp(&halo_struct, r, rho, prs);

        data->Vc(i,j,k,RHO) = rho;
        data->Vc(i,j,k,PRS) = prs;
        data->Vc(i,j,k,VX3) = 0.0;
      } // if (vn > 0.0) 
    } // BOX_LOOP(*box, i,j,k) 

    #ifdef STAGGERED_MHD
      #if INCLUDE_IDIR
      #pragma acc parallel loop collapse(3) present(box, data) async(2)
      BOXI_LOOP(*box, i,j,k) data->Vs[BX1s](i,j,k) = data->Vs[BX1s](i,j,kend);  
      #endif
      #if INCLUDE_JDIR
      #pragma acc parallel loop collapse(3) present(box, data) async(2)
      BOXJ_LOOP(*box, i,j,k) data->Vs[BX2s](i,j,k) = data->Vs[BX2s](i,j,kend);  
      #endif

    #endif
    #pragma acc wait
  } // if (side == X3_END) 

  // ------------------------------------------------------
  // 4. Set internal bounday conditions (floor values)
  // ------------------------------------------------------

  if (side == 0){
    #pragma acc parallel loop collapse(3) present(data, run_config)
    BOX_LOOP(*box,i,j,k){
      bool convert = false;
      double v[NVAR], u[NVAR];
      if(data->Vc(i,j,k, RHO) < run_config.small_rho)  {
        data->Vc(i,j,k, RHO) = run_config.small_rho;
        convert = true;
      }

      if(data->Vc(i,j,k, PRS) < run_config.small_prs) {
        data->Vc(i,j,k, PRS) = run_config.small_prs;
        convert = true;
      }
      // if(data->Vc(i,j,k, TRC) < SMALL_TRC)  {
      //   data->Vc(i,j,k, TRC) = 0.0;
      //   convert = true;
      // }

      if(convert){
        #pragma acc loop seq
        NVAR_LOOP(nv) v[nv] = data->Vc(i,j,k, nv);
        PrimToCons(v,u, run_config);
        #pragma acc loop seq
        NVAR_LOOP(nv) data->Uc(i,j,k, nv) = u[nv];
      }
    }
  }
}

////////////////////////////////////////////////////////////////////////
/*
*/
////////////////////////////////////////////////////////////////////////
#if CT_EMF_BOUNDARY == YES
void CT_EMF_Boundary(const Data *data, Grid *grid)
{
  int ibeg = grid->lbeg[IDIR], iend = grid->lend[IDIR] + INCLUDE_IDIR;
  int jbeg = grid->lbeg[JDIR], jend = grid->lend[JDIR] + INCLUDE_JDIR;
  int kbeg = grid->lbeg[KDIR], kend = grid->lend[KDIR] + INCLUDE_KDIR;
  EMF *emf = data->emf;

  if (grid->lbound[KDIR] != 0){
    const int k = kbeg;
    #pragma acc parallel loop collapse(2) present(emf) copyin(k)
    for (int j = jbeg; j <= jend; j++){
    for (int i = ibeg; i <= iend; i++){
      emf->Ex3e(i,j,k) = 0.0;
    }}
  }
}
#endif

#if ( BODY_FORCE == VECTOR) || (BODY_FORCE == (VECTOR + POTENTIAL))
////////////////////////////////////////////////////////////////////////
#pragma acc routine seq
void BodyForceVector(State const &state, RunConfig const &run_config, double *g)
/*
  Prescribe the acceleration vector as a function of the coordinates
  and the vector of primitive variables *v.
  
  \param [in] v  pointer to a cell-centered vector of primitive
  variables
  \param [out] g acceleration vector
  \param [in] x1  position in the 1st coordinate direction \f$x_1\f$
  \param [in] x2  position in the 2nd coordinate direction \f$x_2\f$
  \param [in] x3  position in the 3rd coordinate direction \f$x_3\f$
*/
////////////////////////////////////////////////////////////////////////
{
  const double x1 = state.x1;
  const double x2 = state.x2;
  const double x3 = state.x3;
  const double r = sqrt(x1*x1 + x2*x2 + x3*x3);
  const double unit_pressure =  run_config.unit_density
                               *run_config.unit_velocity
                               *run_config.unit_velocity;
  const double uphi = unit_pressure/(run_config.unit_density*run_config.unit_length);
  constexpr double kpc = 1.e3*CONST_pc;
  const double norm = 1.0/(r * uphi);

  g[IDIR] = -HaloGradPhi(r * kpc)*x1*norm;
  g[JDIR] = -HaloGradPhi(r * kpc)*x2*norm;
  g[KDIR] = -HaloGradPhi(r * kpc)*x3*norm;
}
#endif

#if (BODY_FORCE == POTENTIAL) || (BODY_FORCE == (VECTOR + POTENTIAL))
////////////////////////////////////////////////////////////////////////
double BodyForcePotential(double x1, double x2, double x3)
/*
  Return the gravitational potential as function of the coordinates.
  
  \param [in] x1  position in the 1st coordinate direction \f$x_1\f$
  \param [in] x2  position in the 2nd coordinate direction \f$x_2\f$
  \param [in] x3  position in the 3rd coordinate direction \f$x_3\f$
  
  \return The body force potential \f$ \Phi(x_1,x_2,x_3) \f$.
*/
////////////////////////////////////////////////////////////////////////
{
  return 0.0;
}
#endif

////////////////////////////////////////////////////////////////////////
/*
  Return the profile function for jet injection as a function of radius
*/
////////////////////////////////////////////////////////////////////////
double Profile(double r, int nv, const RunConfig &run_config) {
  const double rj = run_config.inputParam[JET_RAD]; // jet radius always given in code units

  if ((nv == BX1) || (nv == BX2) || (nv == BX3)) {
    if (r <= rj){
      return 1.0;
    } else {
      const double frac = r/rj;
      const double frac3 = frac*frac*frac;
      return 1.0/cosh(frac3*frac3);
    }
  }

  if (nv == VX3) {
    if (r <= (rj - 0.08)){
      return 1.0;
    } else {
      const double frac = 4.0 * r / (3.0 * rj);
      const double frac3 = frac*frac*frac;
      return 1.0/cosh(frac3*frac3);
    }
  }
  return (r <= rj) ? 1.0 : 0.0;
}
