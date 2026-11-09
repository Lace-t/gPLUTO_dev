#include "gpluto.hpp"
#include "halo.hpp"

////////////////////////////////////////////////////////////////////////
/*
   Initialize Halo variables
*/
////////////////////////////////////////////////////////////////////////
void HaloInit(haloStruct *halo,  RunConfig &run_config)
{
  double nPointsGlob[3];

  nPointsGlob[0] =   0.5*static_cast<double>(run_config.Nx1)
                   + 2.0*static_cast<double>(GetNghost());
  nPointsGlob[1] =   0.5*static_cast<double>(run_config.Nx2)
                   + 2.0*static_cast<double>(GetNghost());
  nPointsGlob[2] =       static_cast<double>(run_config.Nx3)
                   + 2.0*static_cast<double>(GetNghost());

  halo->nPoints = sqrt(nPointsGlob[0]*nPointsGlob[0] +
                       nPointsGlob[1]*nPointsGlob[1] +
                       nPointsGlob[2]*nPointsGlob[2]);

  halo->res = 2.0;

  halo->Rm = sqrt((run_config.dx[0][0]*nPointsGlob[0])*
                  (run_config.dx[0][0]*nPointsGlob[0]) +
                  (run_config.dx[1][0]*nPointsGlob[1])*
                  (run_config.dx[1][0]*nPointsGlob[1]) +
                  (run_config.dx[2][0]*nPointsGlob[2])*
                  (run_config.dx[2][0]*nPointsGlob[2]));

  const int new_size = static_cast<int>(ceil(halo->res*halo->nPoints));
  halo->dr = halo->Rm/static_cast<double>(new_size);
  const double n0 = run_config.inputParam[HALO_N0];

  if (halo->rg.isAllocated() == 0 || halo->size != new_size) {
    if (halo->rg.isAllocated() != 0) {
      halo->rg.destroy();
      halo->prs.destroy();
      halo->rho.destroy();
    }

    halo->rg.create(new_size);
    halo->prs.create(new_size);
    halo->rho.create(new_size);

    halo->size = new_size;
  }

  double integ = 0.0;
  const double drLoc = halo->dr*run_config.unit_length;
  const double unit_pressure =  run_config.unit_density
                               *run_config.unit_velocity
                               *run_config.unit_velocity;
  double r = 0.0;

  for (int i = 0; i < halo->size; ++i) {
    const double integ2 = ((0.6*CONST_mp)/CONST_kB)*integ;
    const double p  = n0*CONST_kB*halo_Tc*exp(-integ2);
    const double de = (0.6*CONST_mp*p)/(CONST_kB*HaloTemperature(r));

    halo->rg(i)  = r/run_config.unit_length;
    halo->prs(i) = p/unit_pressure;
    halo->rho(i) = de/run_config.unit_density;

    if (i == 0) r += drLoc/100.0;
    integ += (HaloGradPhi(r)*drLoc)/HaloTemperature(r);
    if (i == 0) r -= drLoc/100.0;

    r += drLoc;
  }

  // Integrate with midpoint rule

  // for (int i = 0; i < halo->size; ++i) {
  //   const double integ2 = ((0.6*CONST_mp)/CONST_kB)*integ;
  //   const double p  = n0*CONST_kB*halo_Tc*exp(-integ2);
  //   const double de = (0.6*CONST_mp*p)/(CONST_kB*HaloTemperature(r));

  //   halo->rg(i)  = r/run_config.unit_length;
  //   halo->prs(i) = p/unit_pressure;
  //   halo->rho(i) = de/run_config.unit_density;

  //   double rm = r + 0.5*drLoc;  // Midpoint rule
  //   integ += (HaloGradPhi(rm)*drLoc)/HaloTemperature(rm);

  //   r += drLoc;
  // }

// if (g_rank == 0){
//   FILE *fp;
//   fp = fopen("halo.dat","w");
//   for (int i = 0; i < halo->size; ++i) {
//     fprintf (fp,"%12.6e  %12.6e  %12.6e\n",
//                 halo->rg(i), halo->rho(i), halo->prs(i));
//   }
//   fclose(fp);
// }

}

////////////////////////////////////////////////////////////////////////
/*! 
   Interpolate density and pressure at the required location
*/
////////////////////////////////////////////////////////////////////////
void HaloInterp(const haloStruct *halo, const double r, double &rho, double &prs)
{
  int i = static_cast<int>(floor(r / halo->dr));

  if (i < 0) i = 0;
  if (i >= halo->size - 1) i = halo->size - 2;

  const double r0 = halo->rg(i);
  const double r1 = halo->rg(i + 1);

  if (r1 == r0) {
    prs = halo->prs(i);
    rho = halo->rho(i);
    return;
  }else{
    const double w = (r - r0)/(r1 - r0);

    rho = halo->rho(i) + w*(halo->rho(i + 1) - halo->rho(i));
    prs = halo->prs(i) + w*(halo->prs(i + 1) - halo->prs(i));
  }
}

////////////////////////////////////////////////////////////////////////
/*
    Return grad(phi)
*/
////////////////////////////////////////////////////////////////////////
double HaloGradPhi(const double r)
{
  constexpr double kpc  = 1.e3*CONST_pc;
  constexpr double Mb   = 2.e11*CONST_Msun;
  constexpr double aH   = 2.0*kpc;
  constexpr double M200 = 2.1e47;
  constexpr double rs   = 1.e2*kpc;

  // Bulge contribution
  const double g_bulge = CONST_G*Mb / ((r + aH)*(r + aH));

  // Halo contribution
  const double inv_r2   = 1.0 / (r*r);
  const double log_term = log(1.0 + r / rs)*inv_r2;
  const double rational_term = 1.0 / (r*r + rs*r);

  const double g_halo = (CONST_G*M200 / 3.30698)*(log_term - rational_term);

  return g_bulge + g_halo;
}
