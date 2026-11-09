////////////////////////////////////////////////////////////////////////
/*! 
  \file  
  \brief Relativistic jet propagation.

  This problem sets initial and boundary conditions for a jet
  propagating into a uniform medium with constant density and pressure.
  The computation can be carried using \c CYLINDRICAL or \c SPHERICAL
  coordinates.
  In the first case, the jet enters at the lower z-boundary with speed
  \f$ v_z = \beta\f$ while in the second case, a conical beam with a
  small aperture (\f$\theta = 5^\circ\f$) is injected with the same
  speed from the lower radial boundary.
  At the border of the nozzle, jet values are smoothly joined with
  ambient values using the Profile() function
  (in the current setting the profile is a sharp transition).
  The jet is pressure-matched so that the beam and ambient pressure
  coincide.
 
  The configuration is defined in terms of the following parameters:

  -# <tt>g_inputParam[BETA]</tt>:     the jet velocity;
  -# <tt>g_inputParam[RHO_IN]</tt>:   the jet density;
  -# <tt>g_inputParam[RHO_OUT]</tt>:  the ambient density.
  -# <tt>g_inputParam[PRESS_IN]</tt>: the jet pressure (also equal to ambient
                                     pressure)

  defined in \c pluto.ini.
  The \c TAUB equation of state is used.

  - Configurations #01 and #02 use \c CYLINDRICAL coordinates;
  - Configuration #03 employs \c SPHERICAL coordinates (see snapshot below)

  \author A. Mignone (andrea.mignone@unito.it)
  \date   Dec 08, 2025
*/
////////////////////////////////////////////////////////////////////////
#include "gpluto.hpp"

static double Profile(double r, int nv);
static void GetJetValues (double *, RunConfig &);

////////////////////////////////////////////////////////////////////////
/*!
*/
////////////////////////////////////////////////////////////////////////
void Init (double *v, double x1, double x2, double x3, RunConfig &run_config)
{
  double *inputParam = run_config.inputParam;

  #if EOS == IDEAL
  run_config.gamma = 5.0/3.0;
  #endif
  
  v[RHO] = inputParam[RHO_OUT];
  v[VX1] = 0.0;
  v[VX2] = 0.0;
  v[VX3] = 0.0;
  v[PRS] = inputParam[PRESS_IN];
  
  #if ENTROPY_SWITCH
  v[ENTR] = Entropy(v, run_config);
  #endif

}

////////////////////////////////////////////////////////////////////////
/*! 
  Assign initial condition by looping over the computational domain.
  Called after the usual Init() function to assign initial conditions
  on primitive variables.
  Value assigned here will overwrite those prescribed during Init().
  
*/
////////////////////////////////////////////////////////////////////////
void InitDomain (Data *data, Grid *grid)
{
}

////////////////////////////////////////////////////////////////////////
/*!
*/
////////////////////////////////////////////////////////////////////////
void Analysis (Data *data, Grid *grid)
{
}

////////////////////////////////////////////////////////////////////////
/*!
*/
////////////////////////////////////////////////////////////////////////
void UserDefBoundary (const Data *data, RBox *box, int side, Grid *grid) 
{
  int   nv, i, j, k;
  double  r, vjet[NVAR];
  RunConfig &run_config = *(data->run_config);

  #if GEOMETRY == SPHERICAL 
  int ibeg = grid->lbeg[IDIR];
  if (side == X1_BEG){
    GetJetValues(vjet, run_config);
    #pragma acc parallel loop collapse(3) present(data, box, grid) copyin(vjet[:NVAR])
    BOX_LOOP(*box,i,j,k){
      double vout[NVAR];
      NVAR_LOOP(nv) vout[nv] = data->Vc(2*ibeg-i-1,j,k,nv);
      vout[VX2] *= -1.0;
      //vout[VX3] *= -1.0;

      r = grid->xC[JDIR][j];
      NVAR_LOOP(nv)
        data->Vc(i,j,k,nv) = vout[nv] + (vjet[nv] - vout[nv])*Profile(r,nv);

      // Recompute entropy from primitive variables

      #if ENTROPY_SWITCH 
      NFLX_LOOP(nv) vout[nv] = data->Vc(i,j,k,nv);
      data->Vc(i,j,k,ENTR) = Entropy(vout, run_config);
      #endif
    }
  }
  #endif

  #if GEOMETRY == CYLINDRICAL || GEOMETRY == CARTESIAN
  int kbeg = grid->lbeg[KDIR];
  if (side == X3_BEG){
    GetJetValues(vjet, run_config);
    #pragma acc parallel loop collapse(3) present(data, box, grid) copyin(vjet[:NVAR])
    BOX_LOOP(*box,i,j,k){
      double vout[NVAR];
      NVAR_LOOP(nv) vout[nv] = data->Vc(i,j,2*kbeg-k-1,nv);
      vout[VX3] *= -1.0;

      r = grid->xC[IDIR][i];
      NVAR_LOOP(nv){
        data->Vc(i,j,k,nv) = vout[nv] + (vjet[nv] - vout[nv])*Profile(r,nv);
      }

      // Recompute entropy from primitive variables

      #if ENTROPY_SWITCH 
      NFLX_LOOP(nv) vout[nv] = data->Vc(i,j,k,nv);
      data->Vc(i,j,k,ENTR) = Entropy(vout, run_config);
      #endif

    }
  }
  #endif
}

////////////////////////////////////////////////////////////////////////
/*!
*/
////////////////////////////////////////////////////////////////////////
void GetJetValues (double *vjet, RunConfig &run_config)
{
  double *inputParam = run_config.inputParam;

  vjet[RHO] = inputParam[RHO_IN];
  vjet[VX1] = vjet[VX2] = vjet[VX3] = 0.0;
  #if GEOMETRY == CYLINDRICAL || GEOMETRY == CARTESIAN
   vjet[VX3] = inputParam[BETA];
  #elif GEOMETRY == SPHERICAL
   vjet[VX1] = inputParam[BETA];
  #endif
  vjet[PRS] = inputParam[PRESS_IN];

  #if ENTROPY_SWITCH
  vjet[ENTR] = Entropy(vjet, run_config);
  #endif
}
 
////////////////////////////////////////////////////////////////////////
/*!
*/
////////////////////////////////////////////////////////////////////////
double Profile(double r, int nv)
{
  int xn = 14;
  double r0 = 1.0;

  if (nv == RHO) r0 = 1.1;

  #if GEOMETRY == SPHERICAL
   r0 = 5.0/180.0*CONST_PI;
  #endif
  return 1.0/cosh(pow(r/r0,xn));
}

////////////////////////////////////////////////////////////////////////
/*!
*/
////////////////////////////////////////////////////////////////////////
#if MAPPED_GRID == YES
double MappedGrid (Grid *grid, int dir, double xi) {

  if (dir == IDIR) return MappedGrid_Exp(grid, dir, xi);
  else             return MappedGrid_Uniform(grid, dir, xi);
}
#endif
