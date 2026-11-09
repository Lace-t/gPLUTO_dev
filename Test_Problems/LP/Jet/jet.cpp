#include "gpluto.hpp"
#include "halo.hpp"
#include "jet.hpp"

////////////////////////////////////////////////////////////////////////
/*!
   Initialize Jet variables
*/
////////////////////////////////////////////////////////////////////////
void JetInit(jetStruct *jet,  const RunConfig &run_config,
             const double rhoh0, const double prsh0)
{
  double rad = run_config.inputParam[JET_RAD];
  double rhoj, prsj;
  rad *= run_config.unit_length;

  rhoj = jet->rhoj =  run_config.inputParam[JET_ETA]
                     *rhoh0
                     *run_config.unit_density;

  prsj = jet->prsj = run_config.inputParam[JET_PRS];

  const double area = CONST_PI*rad*rad;

  #if EOS == IDEAL
  jet->enth =   rhoj*CONST_c*CONST_c
              + run_config.gamma / (run_config.gamma - 1.0)*prsj;
  #elif EOS == TAUB
  jet->enth =   2.5*prsj
              + sqrt(2.25*prsj*prsj
              + (rhoj*CONST_c*CONST_c)*(rhoj*CONST_c*CONST_c));
  #else
  #error Wrong physics!
  #endif

  const double gam = run_config.inputParam[JET_GAMMA];
  jet->beta = sqrt(1.0 - 1.0 / (gam*gam));

  jet->unit_pressure =   run_config.unit_density
                        *run_config.unit_velocity
                        *run_config.unit_velocity;
  jet->unit_mag = sqrt(4.0*CONST_PI*jet->unit_pressure);

  const double v = jet->beta*CONST_c;
  const double pwr = fabs(v)*area*gam*(gam*jet->enth - rhoj*CONST_c*CONST_c);

  jet->B0 = sqrt((4.0*CONST_PI*pwr*run_config.inputParam[JET_SIGMAB]) /
                 (area*v));

  print("  --------------------------------------------------------\n");
  print("  o jet power : %e ergs-1\n", pwr);
  print("  o B0: %.12e\n", jet->B0);
  print("  o pj/pa = %lf \n", (prsj / jet->unit_pressure) / prsh0);
  print("  --------------------------------------------------------\n");

  #if (PERTURB_INLET == 1) && (DIMENSIONS == 3) && (GEOMETRY == CARTESIAN)
  int offset = 0;
  #ifdef PARLLEL
  offset = g_rank;
  #endif
  jet->seed[0] = 4*offset + 32768;
  jet->seed[1] = 16*offset + 32768;
  jet->seed[2] = 32*offset + 32768;
  #endif
}

////////////////////////////////////////////////////////////////////////
/*
   Compute jet values
*/
////////////////////////////////////////////////////////////////////////
void JetValues(const jetStruct *jet, double *vjet,
               const double x1, const double x2, const double x3,
               const RunConfig &run_config)
{
  const double rhoj = jet->rhoj;
  const double prsj = jet->prsj;

  const double v = jet->beta*CONST_c;

  double vx1 = 0.0;
  double vx2 = 0.0;
  double vx3 = v;

  const double R2  = x1*x1 + x2*x2;
  const double R   = sqrt(R2);

  double sphi = 0.0;
  double cphi = 1.0;
  if (R > 0.0) {
    sphi = x2 / R;
    cphi = x1 / R;
  }

  #if (DIV_JET == YES)
  {
    const double rad = run_config.inputParam[JET_RAD];
    const double h    = rad / tan(run_config.inputParam[JET_ANG]*CONST_PI / 180.0);
    const double diag = sqrt(h*h + R2);
    const double st   = R / diag;
    const double ct   = h / diag;

    vx1 = v*st*cphi;
    vx2 = v*st*sphi;
    vx3 = v*ct;
  }
  #endif

  double bx1 = -jet->B0*sphi;
  double bx2 =  jet->B0*cphi;
  double bx3 = 0.0;

  double csj;
  #if EOS == IDEAL
  csj = sqrt(run_config.gamma*prsj / jet->enth);
  #elif EOS == TAUB
  csj = prsj / (3.0*jet->enth)*(5.0*jet->enth - 8.0*prsj)
      / (jet->enth - prsj);
  csj = sqrt(csj);
  #else
  #error JetValues: EOS undefined
  #endif

  csj *= CONST_c / run_config.unit_velocity;

  #if (PERTURB_INLET == 1) && (DIMENSIONS == 3) && (GEOMETRY == CARTESIAN)
  {
    static constexpr double omega_coeff[8] = {
      0.50, 1.00, 2.00, 3.00, 0.03, 0.06, 0.12, 0.25
    };

    const double eps = 0.005;
    const double gam = run_config.inputParam[JET_GAMMA];
    const double vr  = sqrt(vx1*vx1 + vx2*vx2);

    double Anorm, th, omegal, bl, gam_new;
    int m, l;

    const double scrh2 = (1.0 + eps)*(1.0 + eps);
    Anorm = -vr + sqrt(vr*vr + (scrh2 - 1.0) / (gam*gam*scrh2));

    th = atan2(x2, x1);
    if (th < 0.0) th += 2.0*CONST_PI;

    const double cth = cos(th);
    const double sth = sin(th);

    #pragma acc loop seq
    for (m = 0; m <= 2; ++m) {
      #pragma acc loop seq
      for (l = 0; l < 8; ++l) {
        bl = -1.0 + 2.0*erand48(jet->seed);
        omegal = csj*omega_coeff[l]*g_time;

        const double phase = m*th + omegal + bl;
        const double amp   = Anorm / 24.0*cos(phase);

        vx1 += amp*cth*uvel;
        vx2 += amp*sth*uvel;
      }
    }

    gam_new = (vx1*vx1 + vx2*vx2 + vx3*vx3) / (CONST_c*CONST_c);
    gam_new = 1.0 / sqrt(1.0 - gam_new);

    if (gam_new > (1.0 + eps)*gam) {
      printLog("JetValues: Perturbed gmm[%lf] > (1+e)gmm[%lf]\n", gam_new, gam);
      QuitPLUTO("JetValues()",1);
    }
  }
  #endif

  const double inv_udens  = 1.0 / run_config.unit_density;
  const double inv_uvel   = 1.0 / run_config.unit_velocity;
  const double inv_upress = 1.0 / jet->unit_pressure;
  const double inv_umag   = 1.0 / jet->unit_mag;

  vjet[RHO]     = rhoj*inv_udens;
  vjet[PRS]     = prsj*inv_upress;
  vjet[VX1]     = vx1*inv_uvel;
  vjet[VX2]     = vx2*inv_uvel;
  vjet[VX3]     = vx3*inv_uvel;
  vjet[BX1]     = bx1*inv_umag;
  vjet[BX2]     = bx2*inv_umag;
  vjet[BX3]     = bx3*inv_umag;
  vjet[TRC]     = run_config.inputParam[JET_TRC];
  vjet[TRC + 1] = 1.0 - run_config.inputParam[JET_TRC];

}
