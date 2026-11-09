#include "gpluto.hpp"
#include <cassert>
#include <float.h>

int main()
{
  RunConfig rc;
  rc.gamma = 5.0 / 3.0;
  double cases[6][8] = {{1, 0, 0, 0, 0, 0, 0, 1},
                        {1, 0.7, 0.2, -0.1, 1, 2, -1, 0.1},
                        {0.01, 0.99, 0, 0, 0.1, 0.2, 0, 0.001},
                        {1, 0.1, -0.2, 0.3, 20, 10, 5, 0.01},
                        {1, 0, 0, 0, 1.e-200, 0, 0, 1},
                        {1, 0.9, 0, 0, 1.e-310, 0, 0, 0.1}};
  int failures = 0;
#pragma acc parallel loop copy(cases, failures) reduction(+ : failures)
  for (int c = 0; c < 6; ++c) {
    double p[NVAR] = {0}, u[NVAR], v[NVAR], u2[NVAR];
    for (int n = 0; n < 8; ++n) {
      p[n] = cases[c][n];
    }
#if NTRACER > 0
    p[TRC] = 0.37;
#endif
    PrimToCons(p, u, rc);
    if (!RMHD_PCP_Check(u) || RMHD_PCP_NR(u, v, rc)) {
      ++failures;
      continue;
    }
    PrimToCons(v, u2, rc);
    for (int n = 0; n < 8; ++n) {
      if (fabs(u[n] - u2[n]) > 1.e-8 * fmax(1.0, fabs(u[n]))) {
        ++failures;
      }
    }
    double bad[NVAR];
    for (int n = 0; n < NVAR; ++n) {
      bad[n] = u[n];
    }
    bad[ENG] = -1.0;
    if (RMHD_PCP_Check(bad) || RMHD_PCP_Limit(u, bad) != 1 || !RMHD_PCP_Check(bad)) {
      ++failures;
    }
    for (int n = 0; n < 3; ++n) {
      if (bad[BX1 + n] != u[BX1 + n]) {
        ++failures;
      }
    }
    double ct[NVAR];
    for (int n = 0; n < NVAR; ++n) {
      ct[n] = u[n];
    }
    double b[3] = {p[BX1] * 1.01, p[BX2] * 0.99, p[BX3]};
    if (CT_RMHD_CorrectCell(ct, v, b, rc) || ct[RHO] != u[RHO]) {
      ++failures;
    } else {
#if NTRACER > 0
      if (ct[TRC] != u[TRC]) {
        ++failures;
      }
#endif
      for (int n = 0; n < 3; ++n) {
        if (fabs(v[VX1 + n] - p[VX1 + n]) > 1.e-8) {
          ++failures;
        }
      }
      if (fabs(v[PRS] - p[PRS]) > 1.e-8 * fmax(1.0, p[PRS])) {
        ++failures;
      }
    }
    ct[RHO] = -1.0;
    double old = ct[ENG];
    if (!CT_RMHD_CorrectCell(ct, v, b, rc) || ct[ENG] != old) {
      ++failures;
    }
    bad[RHO] = NAN;
    if (RMHD_PCP_Check(bad)) {
      ++failures;
    }
  }
  // Parallel flow has exactly zero magnetic momentum.  The expanded
  // B^2*v - (v.B)*B formula loses this small fluid signal at large B.
  double parallel[NVAR] = {0}, mapped[NVAR];
  parallel[RHO] = 2.e-4;
  parallel[PRS] = 7.e-3;
  parallel[VX1] = 0.883193;
  parallel[BX1] = 5000.0;
  PrimToCons(parallel, mapped, rc);
  double expected_m = parallel[RHO] * Enthalpy(parallel, rc) *
                      parallel[VX1] / (1.0 - parallel[VX1]*parallel[VX1]);
  if (fabs(mapped[MX1] - expected_m) > 1.e-14 * fabs(expected_m)) {
    ++failures;
  }
  // First rejected FV state from the THETA=0, BMAG=5000 blast.
  double blast[NVAR] = {0}, recovered[NVAR];
  blast[RHO] = 0.00043371928267653024;
  blast[MX1] = 0.076011161476625455;
  blast[MX2] = -0.3056688594826531;
  blast[BX1] = 5000.0000524711213;
  blast[BX2] = 4.1126184909129552e-5;
  blast[ENG] = 12500000.340494687;
#if RMHD_REDUCED_ENERGY == YES
  blast[ENG] -= blast[RHO];
#endif
  if (!RMHD_PCP_Check(blast) || RMHD_PCP_NR(blast, recovered, rc)) {
    ++failures;
  }

  // Independent long-double forward states: near-parallel strong fields,
  // hydro, oblique fields, high Lorentz factors and changes of units.
  // Accuracy checks belong in tests, not the runtime recovery decision.
  int stress_count = 0;
  double max_velocity_error = 0.0, max_pressure_scaled_error = 0.0;
  const double field_values[] = {0.0, 1.e-200, 0.5, 2000.0, 5000.0};
  const double speeds[] = {0.0, 0.3, 0.9, 0.9999};
  const double angles[] = {0.0, 1.e-10, 0.7853981633974483};
  const double units[] = {1.e-20, 1.0, 1.e20};
  for (double field : field_values) {
    for (double speed : speeds) {
      for (double angle : angles) {
        for (double unit : units) {
          long double rho = 0.01L*unit, pressure = 0.001L*unit;
          long double vx = speed*cosl(angle), vy = speed*sinl(angle);
          long double bx = field*sqrtl(unit), by = 0.0L;
          long double vv = vx*vx + vy*vy, w2 = 1.0L/(1.0L-vv);
          long double enthalpy_density = rho + 2.5L*pressure;
          long double ez = vx*by-vy*bx;
          double input[NVAR] = {0}, output[NVAR] = {0};
          input[RHO] = rho*sqrtl(w2);
          input[MX1] = enthalpy_density*w2*vx + by*ez;
          input[MX2] = enthalpy_density*w2*vy - bx*ez;
          input[BX1] = bx;
          input[BX2] = by;
          input[ENG] = enthalpy_density*w2-pressure + 0.5L*(bx*bx+by*by+ez*ez);
#if RMHD_REDUCED_ENERGY == YES
          input[ENG] = enthalpy_density*w2-pressure-rho*sqrtl(w2) +
                       0.5L*(bx*bx+by*by+ez*ez);
#endif
          ++stress_count;
          if (RMHD_PCP_NR(input, output, rc)) {
            printf("NR stress recovery failed: B=%g v=%g angle=%g unit=%g\n",
                   field, speed, angle, unit);
            ++failures;
            continue;
          }
          double ve = fmax(fabs(output[VX1]-(double)vx), fabs(output[VX2]-(double)vy));
          // Thermal recovery is conditioned by stored energy precision.
          double tolerance = 1.e-8 + 256.0*DBL_EPSILON*
                             fabs(input[ENG]/(double)pressure);
          double pe = fabs(output[PRS]/(double)pressure - 1.0);
          max_velocity_error = fmax(max_velocity_error, ve);
          max_pressure_scaled_error = fmax(max_pressure_scaled_error, pe/tolerance);
          if (ve > tolerance || pe > tolerance || !(output[RHO] > 0.0) ||
              !(output[PRS] > 0.0)) {
            printf("NR stress accuracy failed: B=%g v=%g angle=%g unit=%g ve=%g pe=%g tol=%g\n",
                   field, speed, angle, unit, ve, pe, tolerance);
            ++failures;
          }
        }
      }
    }
  }
  printf("NR stress: %d states, max velocity error=%g, max pressure error/tolerance=%g\n",
         stress_count, max_velocity_error, max_pressure_scaled_error);

  // Strict FOFC recovery rejects this candidate without changing it.
  // The later, explicitly enabled pressure-repair phase can recover it.
  double bad_energy[NVAR] = {0};
  double repaired[NVAR] = {0};
  bad_energy[RHO] = 1.0;
  bad_energy[ENG] = -1.0;
  if (!RMHD_PCP_Recover(bad_energy, repaired, rc, 0) || bad_energy[ENG] != -1.0) {
    ++failures;
  }
  if (RMHD_PCP_Recover(bad_energy, repaired, rc, 1) || !RMHD_PCP_Check(bad_energy)) {
    ++failures;
  }

  printf("PCP kernel tests: %d failures\n", failures);
  return failures ? 1 : 0;
}
