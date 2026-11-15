#ifndef GPLUTO_RECONSTRUCTION_GEOMETRY_HPP
#define GPLUTO_RECONSTRUCTION_GEOMETRY_HPP

// Directional reconstruction in the physical coordinate, not cell index.
// No additional Grid storage: the same helpers run on CPU and accelerator.
#pragma omp declare target
struct LinearGeometry {
  double wp, wm, cp, cm, dp, dm;

  #pragma acc routine seq
  LinearGeometry(const Grid *grid, int dir, int i)
      : wp(1.0), wm(1.0), cp(2.0), cm(2.0), dp(0.5), dm(0.5) {
    // A mapped Cartesian grid must not take the uniform shortcut even if
    // UNIFORM_CARTESIAN_GRID was left at its Cartesian default.
#if MAPPED_GRID == YES || UNIFORM_CARTESIAN_GRID == NO
    const double h = grid->dx[dir][i];
    const double forward = grid->xC[dir][i+1] - grid->xC[dir][i];
    const double backward = grid->xC[dir][i] - grid->xC[dir][i-1];
    wp = h/forward;
    wm = h/backward;
    // xC is the geometrical midpoint in the current Grid representation.
    cp = 2.0/wp;
    cm = 2.0/wm;
#endif
  }

  // Bound extrapolations by adjacent cell values even at abrupt width
  // changes. Some smooth limiter formulas alone can exceed these bounds.
  #pragma acc routine seq
  double bound(double slope, double forward, double backward) const {
#if MAPPED_GRID == YES || UNIFORM_CARTESIAN_GRID == NO
    const double envelope = MINMOD_LIMITER(cp*forward,cm*backward);
    return MINMOD_LIMITER(slope,envelope);
#else
    return slope;
#endif
  }

  #pragma acc routine seq
  double minmod(double forward, double backward) const {
    const double p = forward*wp, m = backward*wm;
    return MINMOD_LIMITER(p,m);
  }
};

// Evaluate a polynomial of degree count-1 at a face from count cell
// averages. Use widths scaled by the adjacent cells to avoid cancellation
// at large coordinate offsets. count=4 gives the PPM face shared by cells
// face-1 and face; count=3 gives the biased auxiliary-array parabola.
#pragma acc routine seq
inline double ReconstructionFace(const Grid *grid, int dir, int first,
                                 int count, int face, const double *q) {
  const double scale = grid->dx[dir][face-1] + grid->dx[dir][face];
  double edge[5], node[4], derivative[5];
  const int at = face-first;
  edge[at] = 0.0;
  for (int k=at; k<count; ++k)
    edge[k+1] = edge[k] + grid->dx[dir][first+k]/scale;
  for (int k=at-1; k>=0; --k)
    edge[k] = edge[k+1] - grid->dx[dir][first+k]/scale;
#ifdef HIGH_ORDER
  // HIGH_ORDER reconstructs punctual values, rather than cell averages.
  double value = q[at-1];
  for (int k=0; k<count; ++k) node[k] = 0.5*(edge[k]+edge[k+1]);
  for (int k=0; k<count; ++k) {
    double weight = 1.0;
    for (int m=0; m<count; ++m)
      if (m!=k) weight *= -node[m]/(node[k]-node[m]);
    value += weight*(q[k]-q[at-1]);
  }
  return value;
#else
  // Derivative at this face of the Lagrange interpolant of the primitive
  // integral. Subtract a constant first so constant fields stay exact.
  double diagonal = 0.0;
  for (int k=0; k<=count; ++k) {
    if (k==at) continue;
    double numerator=1.0, denominator=1.0;
    for (int m=0; m<=count; ++m) {
      if (m!=k) denominator *= edge[k]-edge[m];
      if (m!=k && m!=at) numerator *= -edge[m];
    }
    derivative[k] = numerator/denominator;
    diagonal -= derivative[k];
  }
  derivative[at] = diagonal;
  double integral=0.0, value=q[at-1];
  for (int k=0; k<count; ++k) {
    integral += (edge[k+1]-edge[k])*(q[k]-q[at-1]);
    value += derivative[k+1]*integral;
  }
  return value;
#endif
}

#pragma acc routine seq
inline double PPMGeometryFace(const Grid *grid, int dir, int face,
                              double qm, double q0, double qp, double qpp) {
#if MAPPED_GRID == NO
  // Preserve the established uniform formula and operation order.
  return (-qm + 7.0*(q0+qp) - qpp)/12.0;
#else
  double q[4] = {qm,q0,qp,qpp};
  return ReconstructionFace(grid,dir,face-2,4,face,q);
#endif
}
#pragma omp end declare target
#endif
