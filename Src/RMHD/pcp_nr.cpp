/* ///////////////////////////////////////////////////////////////////// */
/*
   Admissibility checks, conservative-state limiting, and the ideal-EOS
   physical-constraint-preserving primitive recovery for RMHD.

   The Wu admissibility and reconstruction-limiter routines are independent
   of the normal Mignone-McKinney recovery path and can also support stage
   limiting.  Only the Cai primitive-recovery solver is selected by PCP_NR.
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "gpluto.hpp"
#include <float.h>
#include "pcp_math.hpp"
#pragma omp declare target

#define PCP_MAX_ITER 500
#define PCP_LIMIT_ITER 48
#define PCP_TOL (8.0 * DBL_EPSILON)
#define PCP_ROOT_TOL 1.e-14
#define PCP_ABS_FLOOR (8.0 * DBL_EPSILON)
// Admission distance and cancellation detection serve different purposes.
// Keep PCP_TOL at rounding scale; reserve a larger energy distance for C2P.
#define PCP_CHECK_MARGIN (PCP_ADMISSIBILITY_EPS_FACTOR * DBL_EPSILON)
static_assert(PCP_ADMISSIBILITY_EPS_FACTOR >= 8.0 &&
              PCP_ADMISSIBILITY_EPS_FACTOR <= 1.0 / DBL_EPSILON,
              "PCP admission factor must be finite and at least 8");
static_assert(PCP_LIMIT_BACKOFF > 0.0 && PCP_LIMIT_BACKOFF < 1.0,
              "PCP limiter backoff must lie strictly between zero and one");

// NVIDIA devices implement double precision; CPU keeps PLUTO's extended
// precision boundary calculation. The fast checker path is always double.
#ifdef _OPENACC
using PCPWide = double;
#else
using PCPWide = long double;
#endif
#pragma acc routine seq
static PCPWide PCPWideSqrt(PCPWide x)
{
#ifdef _OPENACC
  return sqrt(x);
#else
  return sqrtl(x);
#endif
}
#pragma acc routine seq
static PCPWide PCPWideMax(PCPWide a, PCPWide b)
{
  return a > b ? a : b;
}

#pragma acc routine seq
static double PCPFullEnergy(const double *u)
{
#if RMHD_REDUCED_ENERGY == YES
  return u[ENG] + u[RHO];
#else
  return u[ENG];
#endif
}

#pragma acc routine seq
static int PCPFinite(double x)
{
  return isfinite(x);
}

/* Wu et al. G_0 admissibility test.  The small relative margins prevent a
   state on a floating-point boundary from being passed to the NR solver. */
#pragma acc routine seq
static int PCPCheck(const double *u)
{
  double D, E, m2, B2, tau, phi, rad, lhs, rhs, psi_d, psi_scale;
  PCPWide psi;
  double scale, eps_q, eps_E, a, Etest, q_d;

  if (u == NULL) {
    return 0;
  }
  D = u[RHO];
  E = PCPFullEnergy(u);
  m2 = u[MX1] * u[MX1] + u[MX2] * u[MX2] + u[MX3] * u[MX3];
  B2 = u[BX1] * u[BX1] + u[BX2] * u[BX2] + u[BX3] * u[BX3];
  tau = u[MX1] * u[BX1] + u[MX2] * u[BX2] + u[MX3] * u[BX3];

  if (!PCPFinite(D) || !PCPFinite(E) || !PCPFinite(m2) || !PCPFinite(B2) || !PCPFinite(tau) ||
      D <= PCP_ABS_FLOOR || E <= 0.0 || m2 < 0.0 || B2 < 0.0) {
    return 0;
  }

  a = D * D + m2;
  if (!PCPFinite(a) || a < 0.0) {
    return 0;
  }
  scale = fmax(1.0, fabs(E));
  /* Reserve an energy-scaled admission margin beyond rounding error. */
  eps_q = fmax(PCP_ABS_FLOOR, PCP_CHECK_MARGIN * scale);
  /* Most states are well inside q > 0.  Avoid an unconditional long-double
     square root in this hot checker; use extended precision only when the
     double result is close enough to the boundary for cancellation to
     matter. */
  q_d = E - sqrt(a);
  if (PCPFinite(q_d)) {
    if (q_d > eps_q) {
      /* Fast path: q has a clear positive margin. */
    } else if (q_d < -eps_q) {
      return 0;
    } else {
      PCPWide q = (PCPWide)E - PCPWideSqrt((PCPWide)a);
      if (!isfinite((double)q) || q <= (PCPWide)eps_q) {
        return 0;
      }
    }
  } else {
    return 0;
  }

  a = E * E - D * D - m2;
  if (!PCPFinite(a) || a < 0.0) {
    return 0;
  }
  /* Evaluate Psi on E - epsilon_E, as in the definition of G_epsilon.
     This avoids comparing Psi, which has a different physical scale, to an
     energy-sized tolerance. */
  /* Use the same admission distance for q and the shifted-energy Psi. */
  eps_E = fmax(PCP_ABS_FLOOR, PCP_CHECK_MARGIN * fmax(1.0, fabs(E)));
  Etest = E - eps_E;
  if (!PCPFinite(Etest) || Etest <= 0.0) {
    return 0;
  }
  a = B2 - Etest;
  phi = sqrt(a * a + 3.0 * (Etest * Etest - D * D - m2));
  rad = phi + a;
  if (!PCPFinite(phi) || !PCPFinite(rad) || rad < 0.0) {
    return 0;
  }
  /* Most states are well inside the admissible set.  Resolve those with a
     double-precision fast path; reserve PCPWide for the narrow region
     where the Psi subtraction is vulnerable to cancellation. */
  lhs = (phi - 2.0 * a) * sqrt(rad);
  rhs = sqrt(fmax(0.0, 13.5 * (D * D * B2 + tau * tau)));
  if (PCPFinite(lhs) && PCPFinite(rhs)) {
    psi_d = lhs - rhs;
    psi_scale = fmax(1.0, rhs);
    if (psi_d > PCP_TOL * psi_scale) {
      return 1;
    }
    if (psi_d < -PCP_TOL * psi_scale) {
      return 0;
    }
  }
  /* Psi is a difference of two nearly equal large terms in magnetically
     dominated states.  Evaluate only this checker expression in extended
     precision; the hot NR path remains double precision. */
  {
    PCPWide ldD = (PCPWide)D, ldE = (PCPWide)Etest;
    PCPWide ldB2 = (PCPWide)B2, ldm2 = (PCPWide)m2;
    PCPWide ldtau = (PCPWide)tau;
    PCPWide lda = ldB2 - ldE;
    PCPWide ldrad = lda * lda + 3.0 * (ldE * ldE - ldD * ldD - ldm2);
    PCPWide ldphi, ldradp;
    if (!(isfinite((double)ldrad) && ldrad >= 0.0)) {
      return 0;
    }
    ldphi = PCPWideSqrt(ldrad);
    ldradp = ldphi + lda;
    if (!(isfinite((double)ldphi) && isfinite((double)ldradp) && ldradp >= 0.0)) {
      return 0;
    }
    psi = (ldphi - 2.0 * lda) * PCPWideSqrt(ldradp) -
          PCPWideSqrt(PCPWideMax(0.0, 13.5 * (ldD * ldD * ldB2 + ldtau * ldtau)));
  }
  return isfinite((double)psi) && psi >= 0.0;
}

#pragma acc routine seq
int RMHD_PCP_Check(const double *u)
{
  return PCPCheck(u);
}

/* Scale only (D,m,E), keeping the face magnetic field untouched.  This is
   important for staggered CT: both states on a face retain the same B_n. */
#pragma acc routine seq
int RMHD_PCP_Limit(const double *u0, double *u)
{
  int k, n;
  double theta, lo, hi, mid, epsD;
  double cand[NVAR], anchor[NVAR], stage[NVAR];

  /* Return values: -1 = failed (input unchanged), 0 = unchanged, 1 = limited. */
  if (u == NULL) {
    return -1;
  }
  if (RMHD_PCP_Check(u)) {
    return 0;
  }
  if (u0 == NULL || !RMHD_PCP_Check(u0)) {
    return -1;
  }
  for (n = 0; n < NVAR; n++) {
    if (!PCPFinite(u[n])) {
      return -1;
    }
  }

  /* Keep the face magnetic field fixed for CT.  The actual bisection anchor
     is therefore the center fluid state paired with that face field, which
     must itself be admissible before convex limiting is meaningful. */
  for (n = 0; n < NVAR; n++) {
    anchor[n] = u0[n];
  }
  anchor[BX1] = u[BX1];
  anchor[BX2] = u[BX2];
  anchor[BX3] = u[BX3];
  if (!RMHD_PCP_Check(anchor)) {
    return -1;
  }

  /* Adapted to CT: first enforce the density floor at fixed full energy,
     then enforce q and Psi together along one fixed-B fluid-state segment. */
  for (n = 0; n < NVAR; n++) {
    cand[n] = u[n];
  }
  epsD = nextafter(PCP_ABS_FLOOR, INFINITY);
  if (u[RHO] < epsD) {
    theta = (anchor[RHO] - epsD) / (anchor[RHO] - u[RHO]);
    if (!PCPFinite(theta)) {
      return -1;
    }
    theta = fmax(0.0, fmin(1.0, theta));
    cand[RHO] = epsD;
  }
#if RMHD_REDUCED_ENERGY == YES
  /* ENG stores E-D: changing only D must leave the full energy unchanged. */
  cand[ENG] = u[ENG] + (u[RHO] - cand[RHO]);
#endif
  for (n = 0; n < NVAR; ++n) {
    stage[n] = cand[n];
  }

  /* The former q and Psi stages scaled the same fields toward the same
     anchor.  Combine them into one search using the full checker, which
     rejects failed q states before evaluating Psi.  Keep the admissible
     endpoint; neither B nor the other unscaled fields change. */
  if (!RMHD_PCP_Check(cand)) {
    lo = 0.0;
    hi = 1.0;
    for (k = 0; k < PCP_LIMIT_ITER; k++) {
      mid = 0.5 * (lo + hi);
      cand[RHO] = anchor[RHO] + mid * (stage[RHO] - anchor[RHO]);
      cand[MX1] = anchor[MX1] + mid * (stage[MX1] - anchor[MX1]);
      cand[MX2] = anchor[MX2] + mid * (stage[MX2] - anchor[MX2]);
      cand[MX3] = anchor[MX3] + mid * (stage[MX3] - anchor[MX3]);
      cand[ENG] = anchor[ENG] + mid * (stage[ENG] - anchor[ENG]);
      if (RMHD_PCP_Check(cand)) {
        lo = mid;
      } else {
        hi = mid;
      }
    }
    theta = lo;
    cand[RHO] = anchor[RHO] + theta * (stage[RHO] - anchor[RHO]);
    cand[MX1] = anchor[MX1] + theta * (stage[MX1] - anchor[MX1]);
    cand[MX2] = anchor[MX2] + theta * (stage[MX2] - anchor[MX2]);
    cand[MX3] = anchor[MX3] + theta * (stage[MX3] - anchor[MX3]);
    cand[ENG] = anchor[ENG] + theta * (stage[ENG] - anchor[ENG]);
  }
  /* Retreat from the admitted boundary to leave recovery room in double
     precision, including density-only repairs. Recheck the rounded output:
     interpolation itself can lose the margin for widely separated states.
     B and passive fields retain their input values throughout this guard. */
  for (int j = 0; j < NVAR; ++j) stage[j] = cand[j];
  theta = 1.0 - PCP_LIMIT_BACKOFF;
  for (int attempt = 0; attempt < 16; ++attempt) {
    const int index[5] = {RHO, MX1, MX2, MX3, ENG};
    for (int ix = 0; ix < 5; ++ix) {
      int j = index[ix];
      cand[j] = anchor[j] + theta * (stage[j] - anchor[j]);
    }
    if (RMHD_PCP_Check(cand)) break;
    theta *= 0.5;
  }
  if (!RMHD_PCP_Check(cand)) {
    cand[RHO] = anchor[RHO];
    cand[MX1] = anchor[MX1];
    cand[MX2] = anchor[MX2];
    cand[MX3] = anchor[MX3];
    cand[ENG] = anchor[ENG];
  }
  if (!RMHD_PCP_Check(cand)) {
    return -1;
  }
  for (n = 0; n < NVAR; ++n) {
    u[n] = cand[n];
  }
  return 1;
}

#if PCP_NR == YES
/* Produce a finite primitive state after an unsuccessful recovery.  The
   density uses the supplied lab-frame density when possible, so mass is not
   needlessly discarded; zero velocity makes it a valid rest-frame state. */
#pragma acc routine seq
void RMHD_PCP_SetFallback(const double *u, double *v, RunConfig &run_config)
{
  int n;
  double rho;

  for (n = 0; n < NVAR; n++) {
    v[n] = 0.0;
  }
  rho = (u != NULL && isfinite(u[RHO])) ? u[RHO] : run_config.small_rho;
  v[RHO] = fmax(run_config.small_rho, rho);
  v[PRS] = run_config.small_prs;
  if (u != NULL) {
    v[BX1] = isfinite(u[BX1]) ? u[BX1] : 0.0;
    v[BX2] = isfinite(u[BX2]) ? u[BX2] : 0.0;
    v[BX3] = isfinite(u[BX3]) ? u[BX3] : 0.0;
  }
}

#endif

// PCP_NR

typedef struct {
  double f, df, lor, rho, p, v2;
} PCPEval;

#pragma acc routine seq
static int PCPEvaluate(double xi, double D, double fluid_energy, double B2,
                       double beta1, double beta2, double g0, PCPEval *out)
{
  if (!(xi > 0.0) || !PCPFinite(xi)) {
    return 0;
  }
  const double eta = xi + B2;
  const double vp2 = (beta1 / xi) / xi;
  const double vt2 = (beta2 / eta) / eta;
  // Factor 1-v_parallel^2 before subtracting the transverse contribution.
  const double vp = sqrt(beta1) / xi;
  const double invw2 = PCPFma(1.0 - vp, 1.0 + vp, -vt2);
  if (!(invw2 > 0.0) || !PCPFinite(invw2)) {
    return 0;
  }
  const double invlor = sqrt(invw2);
  const double v2 = vp2 + vt2;
  const double phi_a = -((vp2 / xi) + (vt2 / eta));
  // Avoid xi/W^2 - D/W as two separately rounded products.
  const double vterm = invlor * PCPFma(xi, invlor, -D);
  out->lor = 1.0 / invlor;
  out->v2 = v2;
  out->rho = D * invlor;
  out->p = g0 * vterm;
  // E = xi-p+B^2/2+B^2*v_perp^2/2.  Parallel magnetic terms
  // cancel analytically, including in F', before floating-point evaluation.
  out->f = (xi - out->p) + PCPFma(0.5 * B2, vt2, -fluid_energy);
  out->df = 1.0 - (B2 / eta) * vt2 -
            g0 * (invw2 - 2.0 * xi * phi_a + D * phi_a / invlor);
  return PCPFinite(out->f) && PCPFinite(out->df) && PCPFinite(out->rho) &&
         PCPFinite(out->p) && PCPFinite(v2) && v2 >= 0.0;
}

#pragma acc routine seq
static double PCPCubicRoot(double E, double D, double B2, double tau)
{
  // Same positive root xi_c of x^3+(B^2-E)x^2-c=0 as Cai (2.11).
  // Scale the polynomial and use bracketed Newton to avoid cancellation
  // in Cardano/trigonometric formulas when the positive root is small.
  const double a = B2 - E;
  const double cmag = hypot(sqrt(B2) * D, tau) / sqrt(2.0);
  const double cr = pow(cmag, 2.0 / 3.0);
  const double scale = fmax(fabs(a), cr);
  if (!(scale > 0.0) || !PCPFinite(scale)) {
    return 0.0;
  }
  const double aa = a / scale;
  const double rr = cr / scale;
  const double cc = rr * rr * rr;
  double lo = fmax(0.0, -aa);
  double hi = lo + rr;
  double x = hi;
  for (int k = 0; k < 100; ++k) {
    const double f = PCPFma(x * x, x + aa, -cc);
    if (f > 0.0) {
      hi = x;
    } else {
      lo = x;
    }
    if (hi - lo <= 8.0 * DBL_EPSILON * hi) {
      break;
    }
    const double df = x * (3.0 * x + 2.0 * aa);
    const double trial = x - f / df;
    x = (PCPFinite(trial) && trial > lo && trial < hi)
            ? trial : lo + 0.5 * (hi - lo);
  }
  return scale * hi;
}

/* Cai et al. one-dimensional PCP Newton recovery for an ideal EOS. */
#pragma acc routine seq
int RMHD_PCP_NR(const double *u, double *v, RunConfig &run_config)
{
#if EOS != IDEAL
  (void)u;
  (void)v;
  return 1;
#else
  int k;
  int n, nosc = 0;
  double D, E, m2, B2, tau, bscale, bpar = 0.0, beta1;
  double g0, beta2, bhat[3] = {0.0, 0.0, 0.0};
  double fluid_energy;
  double xi, xnew, alpha1, phi, previous_f = 0.0, step;
  double v2, invw2, lor_check;
  PCPEval ev;

  /* Stage/limiter admission is separate from recovery. Allow NR to try
     states inside the checker's safety layer; validate the recovered state. */
  if (u == NULL || v == NULL) {
    return 1;
  }
  D = u[RHO];
  E = PCPFullEnergy(u);
  m2 = u[MX1] * u[MX1] + u[MX2] * u[MX2] + u[MX3] * u[MX3];
  g0 = (run_config.gamma - 1.0) / run_config.gamma;
  B2 = u[BX1] * u[BX1] + u[BX2] * u[BX2] + u[BX3] * u[BX3];
  tau = u[MX1] * u[BX1] + u[MX2] * u[BX2] + u[MX3] * u[BX3];
  if (!(D > 0.0) || !(E > 0.0) || !PCPFinite(D) || !PCPFinite(E) ||
      !PCPFinite(m2) || !PCPFinite(B2) || !PCPFinite(tau) ||
      !(g0 > 0.0) || !PCPFinite(g0)) {
    return 1;
  }
  fluid_energy = PCPFma(-0.5, B2, E);
  beta2 = m2;
  /* Normalize by the largest component BEFORE forming the norm.  Even
     hypot(B) has poor relative accuracy when B itself is subnormal.
     The scaled direction remains O(1) down to the smallest representable
     field, retaining the hydrodynamic limit without clipping/changing B. */
  bscale = fmax(fabs(u[BX1]), fmax(fabs(u[BX2]), fabs(u[BX3])));
  beta1 = 0.0;
  if (bscale > 0.0) {
    double b1 = u[BX1] / bscale, b2 = u[BX2] / bscale, b3 = u[BX3] / bscale;
    double bnorm = sqrt(b1 * b1 + b2 * b2 + b3 * b3);
    bhat[0] = b1 / bnorm;
    bhat[1] = b2 / bnorm;
    bhat[2] = b3 / bnorm;
    bpar = PCPFma(u[MX1], bhat[0], PCPFma(u[MX2], bhat[1], u[MX3]*bhat[2]));
    // |m x Bhat|^2 avoids m^2-(m.Bhat)^2 for nearly parallel flow.
    const double c1 = PCPFma(u[MX2], bhat[2], -u[MX3]*bhat[1]);
    const double c2 = PCPFma(u[MX3], bhat[0], -u[MX1]*bhat[2]);
    const double c3 = PCPFma(u[MX1], bhat[1], -u[MX2]*bhat[0]);
    beta2 = c1*c1 + c2*c2 + c3*c3;
    beta1 = bpar * bpar;
    if (!PCPFinite(beta1) || beta1 > m2 * (1.0 + 128.0 * DBL_EPSILON)) {
      return 1;
    }
    beta1 = fmin(m2, beta1);
  }
  /* Paper (2.9)-(2.10): use xi_d when F(xi_d) <= 0, otherwise xi_c.
     The proven safe initial value replaces the former bracket search. */
  alpha1 = B2 - E;
  phi = sqrt(alpha1 * alpha1 + 3.0 * (E * E - D * D - m2));
  if (alpha1 > 0.0) {
    // Rationalize phi-2a.  E^2-a^2 = 2 B^2 (E-B^2/2).
    const double numerator = PCPFma(2.0 * B2, fluid_energy, -D*D - m2);
    xi = numerator / (phi + 2.0 * alpha1);
  } else {
    xi = (phi - 2.0 * alpha1) / 3.0;
  }
  if (!PCPEvaluate(xi, D, fluid_energy, B2, beta1, beta2, g0, &ev)) {
    return 1;
  }
  if (ev.f > 0.0) {
    xi = PCPCubicRoot(E, D, B2, tau);
    if (!PCPEvaluate(xi, D, fluid_energy, B2, beta1, beta2, g0, &ev)) {
      return 1;
    }
  }
  for (k = 0; k < PCP_MAX_ITER; k++) {
    if (ev.p <= 0.0 || ev.rho <= 0.0 || ev.df <= 0.0) {
      return 1;
    }
    /* Algorithm A.1: unmodified Newton update, scale-relative step tolerance,
       and stop after more than three residual sign changes (round-off
       oscillation).  Compare signs without an underflow-prone product. */
    xnew = xi - ev.f / ev.df;
    if (!PCPFinite(xnew) || xnew <= 0.0) {
      return 1;
    }
    if ((previous_f < 0.0 && ev.f > 0.0) || (previous_f > 0.0 && ev.f < 0.0)) {
      nosc++;
    }
    previous_f = ev.f;
    step = fabs(xnew - xi);
    xi = xnew;
    /* Re-evaluate at the OUTPUT iterate: do not reuse rho,p at xi_n
       together with velocities calculated at xi_{n+1}. */
    if (!PCPEvaluate(xi, D, fluid_energy, B2, beta1, beta2, g0, &ev)) {
      return 1;
    }
    if (step <= PCP_ROOT_TOL * fabs(xi) || nosc > 3) {
      break;
    }
  }
  if (k == PCP_MAX_ITER) {
    return 1;
  }
  if (ev.p <= 0.0 || ev.rho <= 0.0 || ev.v2 >= 1.0) {
    return 1;
  }

  for (n = 0; n < NVAR; n++) {
    v[n] = 0.0;
  }
  v[RHO] = ev.rho;
  v[PRS] = ev.p;
  {
    double den = xi + B2;
    // Scaled Bhat avoids tau/xi overflow for large B.  Add the
    // parallel correction after division rather than multiplying by B.
    const double parallel = bscale > 0.0 ? (bpar / xi) * (B2 / den) : 0.0;
    v[VX1] = PCPFma(parallel, bhat[0], u[MX1] / den);
    v[VX2] = PCPFma(parallel, bhat[1], u[MX2] / den);
    v[VX3] = PCPFma(parallel, bhat[2], u[MX3] / den);
  }
  v[BX1] = u[BX1];
  v[BX2] = u[BX2];
  v[BX3] = u[BX3];
  v2 = v[VX1] * v[VX1] + v[VX2] * v[VX2] + v[VX3] * v[VX3];
  if (!PCPFinite(v[VX1]) || !PCPFinite(v[VX2]) || !PCPFinite(v[VX3]) || !PCPFinite(v2) ||
      v2 >= 1.0) {
    return 1;
  }
  invw2 = PCPFma(-v[VX1], v[VX1],
              PCPFma(-v[VX2], v[VX2], PCPFma(-v[VX3], v[VX3], 1.0)));
  if (!(invw2 > 0.0)) {
    return 1;
  }
  lor_check = 1.0 / sqrt(invw2);
  if (!PCPFinite(lor_check)) {
    return 1;
  }

  /* Paper (2.13): calculate the final thermodynamics from the recovered
     velocity and the final xi, so D = rho W and p = P(rho,h) use the
     same Lorentz factor as the returned primitive vector.  Do not require
     agreement with the internal Wcal to a fixed relative tolerance:
     reconstructing W from explicit velocities amplifies round-off by O(W^2).
     Retain finite, positive thermodynamics and subluminal velocity checks. */
  v[RHO] = D / lor_check;
  v[PRS] = g0 * PCPFma(xi, invw2, -v[RHO]);
  if (!PCPFinite(v[RHO]) || !PCPFinite(v[PRS]) || v[RHO] <= 0.0 || v[PRS] <= 0.0) {
    return 1;
  }

  // Disabled at user request: do not reject C2P using a P2C round trip.
#if 0
  double rel, ucheck[NVAR];
  /* Verify that the recovered primitive state maps back to the supplied
     conservative state.  This catches cancellation in extreme magnetized
     or ultra-relativistic states before they enter the flux calculation. */
  PrimToCons(v, ucheck, run_config);
  rel = fmax(1.0, fabs(E));
  if (fabs(PCPFullEnergy(ucheck) - E) > 1.e-8 * rel) {
    return 1;
  }
  rel = fmax(1.0, fabs(D));
  if (fabs(ucheck[RHO] - D) > 1.e-8 * rel) {
    return 1;
  }
  for (n = MX1; n <= MX3; n++) {
    rel = fmax(1.0, fabs(u[n]));
    if (fabs(ucheck[n] - u[n]) > 1.e-8 * rel) {
      return 1;
    }
  }
#endif
  return 0;
#endif
}

#undef PCP_MAX_ITER
#undef PCP_LIMIT_ITER
#undef PCP_TOL
#undef PCP_ROOT_TOL
#undef PCP_ABS_FLOOR

#pragma omp end declare target
