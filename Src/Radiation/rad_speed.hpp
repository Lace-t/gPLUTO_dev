// /////////////////////////////////////////////////////////////////////
/*!
  \file 
  \brief   Compute the outermost wave speeds for the radiation block.

  Rad_Speed() computes an estimate to the leftmost and rightmost
  wave signal speeds of the radiation part included in the radiation
  module, using the information of the input states ::vR and ::vL.

  \authors J. D. Melon Fuksman (fuksman@mpia.de)
  \date    May 9, 2026

*/
// /////////////////////////////////////////////////////////////////////

template<int dir>
void MaxRadSpeed (State &, double &, double &);

#pragma omp declare target
#pragma acc routine seq
template<int dir>
////////////////////////////////////////////////////////////////////////
/*!
  Calculate the maximum and minimum wave speeds of the radiation part,
  and store it in cmin and cmax.
  Maximum and minimum eigenvalues are calculated explicitly as in, e.g.,
  M. A. Skinner & E. C. Ostriker, ApJS 206, 91 (2013).
  The eigenvalues depend only on the values of f=|F_{rad}|/E and
  cos(\theta)=Fn/|F_{rad}|.

  \tparam     dir    Sweep direction (IDIR, JDIR, or KDIR).
  \param[in]  state  Input radiation state.
  \param[out] cmin   Minimum wave speed.
  \param[out] cmax   Maximum wave speed.
*/
////////////////////////////////////////////////////////////////////////
inline void MaxRadSpeed (State &state, double &cmin, double &cmax)
{
  double *u = state.u;
  double fx, f, costheta;
  double sq1, sq1m1, sq2, qq, f2, fc;

  GEN_INDEX_CYCLE (dir, FR1, FRn, FRt, FRb);

  fx = u[FRn] ;
  f = sqrt( fx*fx + u[FRt]*u[FRt] + u[FRb]*u[FRb] ) ;

  if (f < 1e-150 ) {

    cmax = 1.0/sqrt(3) ;
    cmin = - cmax ;

  } else {

    costheta = fx/f ;

    f = f/u[ENR] ;
    if (f > 1.0) f = 1.0 ;
    f2 = f*f ;

    qq = fabs(4.0-3.0*f2) ;
    sq1 = sqrt(qq) ;
    sq1m1 = 1.0/sq1 ;

    sq2 = 0.666666667*(qq-sq1) + 2.0*costheta*costheta*(2.0-f2-sq1)  ;
    sq2 = sqrt(fabs(sq2)) ;

    fc = f*costheta ;

    cmax = sq1m1*(fc+sq2) ;
    cmin = sq1m1*(fc-sq2) ;

  }
}
#pragma omp end declare target

template<int dir>
double LimitRadWaveVelocity (State &, State &, Grid *, RunConfig &);

#pragma omp declare target
#pragma acc routine seq
template<int dir>
////////////////////////////////////////////////////////////////////////
/*!
  Takes an input the left and right states, the grid data, and the current
  position. Returns 4/(3*tau), where tau is the optical depth along the cell.

  \tparam    dir        Sweep direction (IDIR, JDIR, or KDIR).
  \param[in] stateL     Left state.
  \param[in] stateR     Right state.
  \param[in] grid       Grid structure.
  \param[in] run_config Run-time configuration.
*/
////////////////////////////////////////////////////////////////////////
inline double LimitRadWaveVelocity (State &stateL, State &stateR, Grid *grid,
                                    RunConfig &run_config)
{
  double *vL = stateL.v;
  double *vR = stateR.v;
  int i = stateL.i;
  int j = stateL.j;
  int k = stateL.k;

  double dx, q, lim_factor ;
  #if RADIATION_RELATIVISTIC
  double vel2 = 0.5*( vL[VX1]*vL[VX1] + vL[VX2]*vL[VX2] + vL[VX3]*vL[VX3]
                    + vR[VX1]*vR[VX1] + vR[VX2]*vR[VX2] + vR[VX3]*vR[VX3] ) ;
  double lgamma_1 = sqrt(1 - vel2) ;
  #endif

  // Define cell dx
  #if ( GEOMETRY == CARTESIAN ) || ( GEOMETRY == CYLINDRICAL )
    if constexpr (dir == IDIR) {
      dx = grid->dx[IDIR][i];
    } else if constexpr (dir == JDIR) {
      dx = grid->dx[JDIR][j];
    } else {
      dx = grid->dx[KDIR][k];
    }
  #elif GEOMETRY == SPHERICAL
    if constexpr (dir == IDIR) {
      dx = grid->dx[IDIR][i] ;
    } else if constexpr (dir == JDIR) {
      dx = grid->xC[IDIR][i] * grid->dx[JDIR][j] ;
    } else {
      dx = grid->xC[IDIR][i] * sin(grid->xC[JDIR][j]) * grid->dx[KDIR][k] ;
    }
  #endif

  // Compute optical depth of the cell
  #if RADIATION_VAR_OPACITIES == YES
    double abs_op, scat_op, tot_op ;
    UserDefOpacities (vL, &abs_op, &scat_op, run_config); 
    tot_op = 0.5*(abs_op + scat_op);
    UserDefOpacities (vR, &abs_op, &scat_op, run_config);
    tot_op += 0.5*(abs_op + scat_op);
  #else
    const double abs_op = run_config.rad_absorption_coeff ;
    const double tot_op = abs_op + run_config.rad_scattering_coeff ;
  #endif
  q = dx * tot_op * (vL[RHO] + vR[RHO]) ;

  // Small-q control. The signal speed is not limited for small opacities
  if ( q < 1e-100 ) return 2.0 ;

  #if RADIATION_NR
  lim_factor = 2.666666667 ;
  #else
  lim_factor = 2.666666667*lgamma_1 ;
  #endif

  return lim_factor/q ;
}
#pragma omp end declare target

template<int dir>
bool RadSpeed (State &, State &, double &, double &, Grid *, RunConfig&);

#pragma omp declare target
#pragma acc routine seq
template<int dir>
////////////////////////////////////////////////////////////////////////
/*!
  Compute leftmost (SL) and rightmost (SR) speed for the Riemann fan.
  Return true if diffusion limiter by Sadowski et al. 2013 is applied,
  return false otherwise.

  \tparam dir         Sweep direction (IDIR, JDIR, or KDIR).
  \param[in]  stateL     Left state.
  \param[in]  stateR     Right state.
  \param[out] aL         Leftmost signal speed.
  \param[out] aR         Rightmost signal speed.
  \param[in]  grid       Grid structure.
  \param[in]  run_config Runtime configuration.
*/
////////////////////////////////////////////////////////////////////////
inline bool RadSpeed (State &stateL, State &stateR, double &aL, double &aR,
                     Grid *grid, RunConfig &run_config)
{
  bool   diff_cond = false;
  double vLim;
  double sl_min, sl_max;
  double sr_min, sr_max;

  MaxRadSpeed<dir> (stateL, sl_min, sl_max);
  MaxRadSpeed<dir> (stateR, sr_min, sr_max);

  aL = MIN(sl_min, sr_min);
  aR = MAX(sl_max, sr_max);

  // Limit rad speed according to the cell optical depth
  #if RADIATION_DIFF_LIMITING
  vLim = LimitRadWaveVelocity<dir>(stateL, stateR, grid, run_config);
  diff_cond = aL < -vLim || aR > vLim;
  aL = MAX(aL, -vLim);
  aR = MIN(aR,  vLim);
  #endif

  return diff_cond ;
}
#pragma omp end declare target