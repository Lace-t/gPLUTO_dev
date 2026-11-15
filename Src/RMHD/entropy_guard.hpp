#pragma once

#include <cmath>

// Local orthonormal fluid entropy test shared by HLLC and HLLD.
// EW includes the full Powell correction. DC tests the fluid subproblem after
// GLM_Solve; CT tests the subproblem with a common staggered normal field.
// The excluded DC/CT flux components are preserved exactly.
//
// Switch to same-state HLL only when the candidate clearly violates the entropy
// inequality and HLL clearly satisfies it. Geometry factors and sources enter
// later in RHS: this test does not establish full reconstruction/CT/RK entropy
// stability. Requires ideal EOS and the reduced-energy convention.
#pragma omp declare target
namespace EntropyGuard {

// Components participating in the fluid entropy pair. In DC, the first eight
// entries are RMHD variables; the GLM scalar and normal B belong to cleaning.
#pragma acc routine seq
inline bool fluid_component(int n, int dir)
{
  #if DIVB_CONTROL == DIV_CLEANING
  return n < 8 && n != BX1 + dir;
  #elif DIVB_CONTROL == CONSTRAINED_TRANSPORT
  return n < NFLX && n != BX1 + dir;
  #else
  return n < NFLX;
  #endif
}

struct Result {
  double original = 0.,  // Candidate entropy residual; positive means violation.
         hll = 0.,       // Same-state HLL reference residual.
         after = 0.,     // Residual of the selected flux, when evaluated.
         lambda = 0.,    // 0: retain candidate; 1: switch completely to HLL.
         tolerance = 0.; // Roundoff allowance for the residual comparison.

  // 0: retained; 1: switched; 2: HLL cannot certify a correction;
  // 3: incompatible DC/CT states or nonfinite residuals.
  int status = 0;
};

// Mathematical entropy eta = -D*s/(gamma-1), where D is lab-frame density and
// s = log(p) - gamma*log(rho). Return w = d eta/d U, entropy potential
// psi = w dot F_physical - eta*v_n, and Powell coefficient phi = D*(v dot B)/p.
// scale accumulates absolute terms for the roundoff allowance. Input primitive
// and conservative states must already be physically admissible and consistent.
#pragma acc routine seq
inline void variables(const double *v, const double *u, const double *f, int dir,
                      const RunConfig &rc, double *w, double &psi, double &phi,
                      double &scale)
{
  const double p = v[PRS], D = u[RHO];
  double v2 = 0., vB = 0.;
  for (int n = 0; n < 3; n++) {
    v2 += v[VX1 + n]*v[VX1 + n];
    vB += v[VX1 + n]*v[BX1 + n];
  }
  const double invG = std::sqrt(1. - v2),
               h = 1. + rc.gamma/(rc.gamma - 1.)*p/v[RHO];
  const double s = std::log(p) - rc.gamma*std::log(v[RHO]);

  for (int n = 0; n < NVAR; n++)
    w[n] = 0.;

  // Reduced energy E_total-D changes the density entry of the entropy vector.
  w[RHO] = -D/p*(1. - h*invG) - s/(rc.gamma - 1.);
  for (int n = 0; n < 3; n++) {
    w[MX1 + n] = D/p*v[VX1 + n];
    w[BX1 + n] = D/p*(v[BX1 + n]*(1. - v2) + vB*v[VX1 + n]);
  }
  w[ENG] = -D/p;
  phi = D*vB/p;

  const double q = -D*s/(rc.gamma - 1.)*v[VX1 + dir];
  psi = -q;
  scale = std::fabs(q);
  for (int n = 0; n < NFLX; n++) {
    if (fluid_component(n, dir)) {
      psi += w[n]*f[n];
      scale += std::fabs(w[n]*f[n]);
    }
  }
}

// Compare the candidate flux with HLL using the same states and outer speeds.
// bn is the actual Powell interface field in EW; DC/CT do not use it in the
// residual. Only flux[0:NFLX] may change. The caller must synchronize pressure
// and compute tracer fluxes after a switch, before storing CT electric fields.
#pragma acc routine seq
inline Result apply(const double *vl, const double *vr,
                    const double *ul, const double *ur,
                    const double *fl, const double *fr,
                    double sl, double sr, double bn,
                    double *flux, int dir, const RunConfig &rc)
{
  Result result;

  #if DIVB_CONTROL == DIV_CLEANING
  // GLM_Solve must have assigned common Bn and psi. No entropy extension for
  // the full, unsplit cleaning system is assumed here.
  if (vl[BX1 + dir] != vr[BX1 + dir] || vl[PSI_GLM] != vr[PSI_GLM]) {
    result.status = 3;
    return result;
  }
  #endif
  #if DIVB_CONTROL == CONSTRAINED_TRANSPORT
  // Reconstruction must use the same staggered normal field on both sides.
  if (vl[BX1 + dir] != vr[BX1 + dir]) {
    result.status = 3;
    return result;
  }
  #endif

  double wl[NVAR], wr[NVAR], psil, psir, phil, phir, scale_l, scale_r;
  variables(vl, ul, fl, dir, rc, wl, psil, phil, scale_l);
  variables(vr, ur, fr, dir, rc, wr, psir, phir, scale_r);

  const double bl = MIN(sl, 0.), br = MAX(0., sr), inv = 1./(br - bl);
  double fh[NVAR], dw[NVAR];

  #if DIVB_CONTROL == EIGHT_WAVES
  const double source = phil*(bn - vl[BX1 + dir])
                      + phir*(vr[BX1 + dir] - bn);
  const double source_scale = std::fabs(phil*(bn - vl[BX1 + dir]))
                            + std::fabs(phir*(vr[BX1 + dir] - bn));
  #else
  const double source = 0., source_scale = 0.;
  #endif

  // R(F) = (wR-wL) dot F - (psiR-psiL) - Powell correction.
  // Fixed-Bn DC/CT have no Powell term. Both fluxes use the same offset.
  const double offset = -(psir - psil) - source;
  result.original = result.hll = offset;
  double scale = scale_l + scale_r + source_scale;
  for (int n = 0; n < NFLX; n++) {
    if (!fluid_component(n, dir)) {
      fh[n] = flux[n];
      dw[n] = 0.;
      continue;
    }
    fh[n] = bl*br*(ur[n] - ul[n]) + br*fl[n] - bl*fr[n];
    fh[n] *= inv;
    dw[n] = wr[n] - wl[n];
    result.original += dw[n]*flux[n];
    result.hll += dw[n]*fh[n];
    scale += std::fabs(dw[n]*flux[n]) + std::fabs(dw[n]*fh[n]);
  }

  // Use an absolute operation scale to avoid treating cancellation noise as a
  // violation. Require strict, opposite signs outside this band before switching.
  result.tolerance = 512.*2.220446049250313080847263336181640625e-16*scale + 1.e-25;
  result.after = result.original;
  if (!std::isfinite(result.original) || !std::isfinite(result.hll)) {
    result.status = 3;
    return result;
  }
  if (result.original <= result.tolerance)
    return result;
  if (result.hll >= -result.tolerance) {
    result.status = 2;
    return result;
  }

  result.lambda = 1.;
  for (int n = 0; n < NFLX; n++)
    flux[n] = fh[n];

  result.after = offset;
  for (int n = 0; n < NFLX; n++)
    result.after += dw[n]*flux[n];

  result.status = 1;
  return result;
}

} // namespace EntropyGuard
#pragma omp end declare target
