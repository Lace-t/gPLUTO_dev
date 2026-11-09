#include "gpluto.hpp"

// Transactional CA2' correction: recover with FV B, rebuild with CT B.
// Mass and passive conservative scalars retain their exact original values.
#pragma acc routine seq
int CT_RMHD_CorrectCell(double *u, double *v, const double *b, RunConfig &run_config)
{
  double work[NVAR], primitive[NVAR], mapped[NVAR];
  for (int n = 0; n < NVAR; ++n) {
    work[n] = u[n];
    if (!isfinite(work[n])) {
      return 1;
    }
  }
  for (int n = 0; n < 3; ++n) {
    if (!isfinite(b[n])) {
      return 1;
    }
  }
  if (RMHD_PCP_Recover(work, primitive, run_config, 0)) {
    return 1;
  }
  double v2 = primitive[VX1] * primitive[VX1] + primitive[VX2] * primitive[VX2] +
              primitive[VX3] * primitive[VX3];
  if (!(v2 < 1.0) || primitive[RHO] <= 0.0 || primitive[PRS] <= 0.0) {
    return 1;
  }
  primitive[RHO] = u[RHO] * sqrt(1.0 - v2);
  for (int n = 0; n < 3; ++n) {
    primitive[BX1 + n] = b[n];
  }
  PrimToCons(primitive, mapped, run_config);
  mapped[RHO] = u[RHO];
  for (int n = 0; n < NVAR; ++n) {
    if (!isfinite(mapped[n]) || !isfinite(primitive[n])) {
      return 1;
    }
  }
  if (!RMHD_PCP_Check(mapped)) {
    return 1;
  }
  for (int n = 0; n < 3; ++n) {
    u[MX1 + n] = mapped[MX1 + n];
    u[BX1 + n] = b[n];
  }
  u[ENG] = mapped[ENG];
  for (int n = 0; n < NVAR; ++n) {
    v[n] = primitive[n];
  }
  return 0;
}
