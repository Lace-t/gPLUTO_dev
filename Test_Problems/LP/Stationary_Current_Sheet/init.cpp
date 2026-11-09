////////////////////////////////////////////////////////////////////////
/*! 
  \file  
  \brief Particle acceleration in 2D reconnection layer.

  Setup a 2D reconnection layer as in sec. 4.1 of Suriano et al (2026). 
  The initial condition is kept constant and used only as a background 
  state for the particles. 
  It is not a stationary solution of the MHD equations.

  \authors A. Suriano (alessio.suriano@unito.it) \n
           A. Mignone (andrea.mignone@unito.it)

 \b References
    - "Particles acceleration with magnetic reconnection in large 
       scale RMHD simulations-II. Particle spectra"\n
      Suriano et al., MNRAS (2026) 546, 1-13

  \date   March 09, 2026
*/
////////////////////////////////////////////////////////////////////////
#include "gpluto.hpp"
//#include <gsl/gsl_sf_dawson.h>  /* special dawson function from gsl */

////////////////////////////////////////////////////////////////////////
void Init (double *v, double x, double y, double z, RunConfig &run_config)
/*! 
*/
////////////////////////////////////////////////////////////////////////
{
  const double v0    = run_config.inputParam[V0];
  const double a     = run_config.inputParam[SHEAR_WIDTH];
  const double eta   = run_config.inputParam[ETA];
  const double rho   = run_config.inputParam[RHO_0];
	const double sigma = run_config.inputParam[SIGMA];
  run_config.unit_velocity = CONST_c;
  constexpr double eps = 0.53;
  
  // Density constant 
  
  v[RHO] = rho;
  
  const double delta = sqrt(2.*eta*a/v0);
  const double E = (v0/a)*(delta/(2.*eps))*sqrt(v[RHO]*sigma);

  double theta = 0.0;
  double yp = -x*sin(theta) + y*cos(theta); // rotated system
  
  // Velocity with stagnation point at x,y=0 

  double vx = v0*x/a;
  double vy = v0*y/a;
  v[VX1] =   vx*cos(theta) - vy*sin(theta);     // Eq. (31)
  v[VX2] = -(vx*sin(theta) + vy*cos(theta));
  v[VX3] = 0.0;
  
  // Magnetic field 
  
  const double B    = 2.*E*a/(v0*delta);
  //  const double Bprof = B*gsl_sf_dawson(yp/delta) + run_config.inputParam[B0_PARAM];
  const double Bprof = B*Dawson(yp/delta) + run_config.inputParam[B0_PARAM];

  v[BX1] = (Bprof)*cos(theta);
  v[BX2] = (Bprof)*sin(theta);
  v[BX2]  = 0.0;
  v[BX3] = run_config.inputParam[B_GUIDE];
  
  // Pressure 
  
  const double c0 = v[RHO];

  v[PRS] = c0/100.0;

}

////////////////////////////////////////////////////////////////////////
void InitDomain (Data *data, Grid *grid)
{
}

////////////////////////////////////////////////////////////////////////
void Analysis (Data *data, Grid *grid)
/*! 
  Perform runtime data analysis.
  
  \param [in] d the PLUTO Data structure
  \param [in] grid   pointer to array of Grid structures  
*/
////////////////////////////////////////////////////////////////////////
{
  const uint64_t nptot = data->pc->Size();
  bool existing = false;
  double arr[12] = {};
  // #pragma acc parallel loop present(data,data->pc) copy(existing, arr[:12]) reduction(|:existing)

  #pragma acc update host(data->pc)

  for (uint64_t np = 0; np < nptot; np++){
    if(data->pc->id(np) == 32645){
      existing |= true;
      arr[0] = static_cast<double>(data->pc->id(np));
      arr[1] = data->pc->pos(np)[0];
      arr[2] = data->pc->pos(np)[1];
      arr[3] = data->pc->pos(np)[2];
      arr[4] = data->pc->color(np)[0];  //sig
      arr[5] = data->pc->color(np)[1];  //bet
      arr[6] = sqrt(data->pc->color(np)[2]);  //Brec
      arr[7] = sqrt(data->pc->color(np)[3]);  //Bguide
      arr[8] = data->pc->color(np)[4];  //p
      arr[9] = data->pc->eng(np)[PARTICLES_LP_NEBINS];
      break;
    }
  }
  
  FILE *fp = nullptr;

  static int first_call = 1; // Persistent across function calls

  char filename[64];
  snprintf(filename, sizeof(filename), "spectrum_check.csv");

  if(existing){
    fp = fopen(filename, (first_call ? "w" : "a"));
    if(first_call) fprintf(fp, "#t;id;posx;posy;posz;sigma;beta;Brec;Bguide;p;Emax\n");
  }
  first_call = 0;
  
  if (fp != nullptr && existing == true) {
    #pragma acc update host(data->pc)
    fprintf(fp, "%f", data->run_config->time);
    fprintf(fp, ";%llu", static_cast<uint64_t>(arr[0]));
    for(int i = 1; i < 10; i++) fprintf(fp, ";%.6e", arr[i]);
    fprintf(fp, "\n");
  }

  if (fp != nullptr) fclose(fp);

}

#if PHYSICS == MHD
////////////////////////////////////////////////////////////////////////
void BackgroundField (double x1, double x2, double x3, double *B0)
/*!
  Define the component of a static, curl-free background 
  magnetic field.
  
  \param [in] x1  position in the 1st coordinate direction \f$x_1\f$
  \param [in] x2  position in the 2nd coordinate direction \f$x_2\f$
  \param [in] x3  position in the 3rd coordinate direction \f$x_3\f$
  \param [out] B0 array containing the vector componens of the background
  magnetic field
*/
////////////////////////////////////////////////////////////////////////
{
   B0[0] = 0.0;
   B0[1] = 0.0;
   B0[2] = 0.0;

}

#endif

////////////////////////////////////////////////////////////////////////
void UserDefBoundary (const Data *data, RBox *box, int side, Grid *grid) 
/*! 
  Assign user-defined boundary conditions.
  
  \param [in,out] d  pointer to the PLUTO data structure containing
  cell-centered primitive quantities (data->Vc) and 
  staggered magnetic fields (data->Vs, when used) to 
  be filled.
  \param [in] box    pointer to a RBox structure containing the lower
  and upper indices of the ghost zone-centers/nodes
  or edges at which data values should be assigned.
  \param [in] side   specifies the boundary side where ghost zones need
  to be filled. It can assume the following 
  pre-definite values: X1_BEG, X1_END,
  X2_BEG, X2_END, 
  X3_BEG, X3_END.
  The special value side == 0 is used to control
  a region inside the computational domain.
  \param [in] grid  pointer to an array of Grid structures.
*/
////////////////////////////////////////////////////////////////////////
{
}

