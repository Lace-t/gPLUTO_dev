////////////////////////////////////////////////////////////////////////
/*!
  \file
  \brief  Header file for GLM Divergence Cleaning

  Contains function prototypes and global variable declaration
  for the GLM formulation to control the divergence-free condition
  of magnetic field.

  \authors A. Mignone (andrea.mignone@unito.it)
  \date    Nov 4, 2025

  \b References
     - "A Second-order unsplit Godunov scheme for cell-centered MHD:
        The CTU-GLM scheme"\n
        Mignone \& Tzeferacos, JCP (2010) 229, 2117

     - "High-order conservative finite difference GLM-MHD scheme for
        cell-centered MHD"\n
        Mignone, Tzeferacos \& Bodo, JCP (2010) 229, 5896
*/
////////////////////////////////////////////////////////////////////////

#ifndef GLM_ALPHA
  #define GLM_ALPHA    0.1  /**< Sets the damping rate of monopoles. */
#endif

#ifndef GLM_EXTENDED
 #define GLM_EXTENDED  NO  /**< The GLM_EXTENDED macro may be turned to YES to enable
                                the extended GLM formalism.
                                Although it breaks conservation of momentum and
                                energy, it has proven to be more robust in treating
                                low-beta plasma. */
#endif

#ifndef GLM_COMPUTE_DIVB
  #define GLM_COMPUTE_DIVB  NO
#endif

////////////////////////////////////////////////////////////////////////
/* Function prototyping                                               */
////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////
/*!
  Add source terms to the right hand side of the conservative equations, 
  momentum and energy equations only. 
  This yields the extended GLM equations given by Eq. (24a)--(24c) in 
  
  "Hyperbolic Divergence cleaning for the MHD Equations"
  Dedner et al. (2002), JcP, 175, 645
*/
////////////////////////////////////////////////////////////////////////
#pragma acc routine seq
inline void GLM_ExtendedSource (const double divB, const double grad_psi, 
                                const double* v, double* src, const int dir)
{
  double bx, by, bz;

// #if BACKGROUND_FIELD == YES
//   GetBackgroundField (stateC, beg - 1, end, CELL_CENTER, grid);
// #endif

  bx = v[BX1];
  by = v[BX2];
  bz = v[BX3];

  // #if BACKGROUND_FIELD == YES
  // bx += stateC->Bbck[i][IDIR];
  // by += stateC->Bbck[i][JDIR];
  // bz += stateC->Bbck[i][KDIR];
  // #endif
  
  src[MX1] -= divB*bx;
  src[MX2] -= divB*by;
  src[MX3] -= divB*bz;

  #if HAVE_ENERGY
  int BXn = BX1 + dir;
  src[ENG] -= v[BXn]*grad_psi;
  #else
  #endif
}

////////////////////////////////////////////////////////////////////////
#pragma acc routine seq
inline void GLM_Solve (double *vL, double *vR, RunConfig &run_config, const int dir)
/*!
*/
////////////////////////////////////////////////////////////////////////
{
  int    BXn = BX1 + dir;
  
  // -------------------------------------------------
  // Solve the Riemann problem for the 2x2 linear
  // GLM system. Use the solution to assign a single
  // value to both left and right (Bn,psi).
  // -------------------------------------------------

  double glm_ch = run_config.glm_ch;
  double dB     = vR[BXn]     - vL[BXn];
  double dpsi   = vR[PSI_GLM] - vL[PSI_GLM];
  double Bm     = 0.5*(vL[BXn]     + vR[BXn])     - 0.5*dpsi/glm_ch;
  double psim   = 0.5*(vL[PSI_GLM] + vR[PSI_GLM]) - 0.5*glm_ch*dB;

  vL[BXn]     = vR[BXn]     = Bm;
  vL[PSI_GLM] = vR[PSI_GLM] = psim;

  #if (PHYSICS == ResRMHD) && (DIVE_CONTROL == DIV_CLEANING)
  int    EXn = EX1 + dir;
  double dE   = vR[EXn] - vL[EXn];
  double dphi = vR[PHI_GLM] - vL[PHI_GLM];
  double Em   = 0.5*(vL[EXn] + vR[EXn]) - 0.5*dphi/glm_ch;
  double phim = 0.5*(vL[PHI_GLM] + vR[PHI_GLM]) - 0.5*glm_ch*dE;

  vL[EXn]     = vR[EXn]     = Em;
  vL[PHI_GLM] = vR[PHI_GLM] = phim;
  #endif
}

void  GLM_Init  (const Data *, Grid *);
void  GLM_Source (Data *, double, Grid *);
