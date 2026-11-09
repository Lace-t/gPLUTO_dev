////////////////////////////////////////////////////////////////////////
/*!
  \file
  \brief FARGO-MHD module header file.

  Contains contains basic definitions and declarations used by the
  FARGO-MHD module.

  \authors A. Mignone    (andrea.mignone@unito.it)\n
           G. Muscianisi (g.muscianisi@cineca.it)\n
           M. Rossazza   (marco.rossazza@unito.it)
  \date    Apr 28, 2025
*/
////////////////////////////////////////////////////////////////////////

#ifdef PARALLEL
  #define FARGO_BUF 16
 #else
  #define FARGO_BUF 0
 #endif

/**<  Set the order of interpolation during the linear transport step.
      Either 2 or 3. Default is 3. */
#ifndef FARGO_ORDER
 #define FARGO_ORDER    3  
#endif
        
/*! Set how often (in number of steps) the total azimuthal 
    velocity should be averaged.                           */
#ifndef FARGO_NSTEP_AVERAGE
  #ifdef SHEARINGBOX
    #define FARGO_NSTEP_AVERAGE  -1   // Never average with ShearinBox model 
  #else 
    #define FARGO_NSTEP_AVERAGE  10   // Default is 10 
  #endif  
#endif

/*! Used to write total/residual velocity */
#ifndef FARGO_OUTPUT_VTOT
  #define FARGO_OUTPUT_VTOT    NO
#endif

// 1,2,3 already taken by vec. potential 
#define  FARGO_W  (NVAR+4)

  // ----------------------------------------------------------------
  // Note: when both FARGO and SHEARINGBOX modules are used, *ALL* 
  // source terms are accounted for by the body_force function 
  // ----------------------------------------------------------------

// -- SBEG and SEND are the initial and final 
//    indices in the direction of the orbital velocity -- 

#if GEOMETRY != SPHERICAL
 #define SDIR   JDIR  // Orbital direction. 
//  #define SBEG   JBEG  /* Starting index in the orbital dir. */
//  #define SEND   JEND  /* Final index in the orbital dir. */
//  #define NS     NX2
//  #define NS_TOT NX2_TOT
#else
 #define SDIR   KDIR
//  #define SBEG   KBEG
//  #define SEND   KEND
//  #define NS     NX3
//  #define NS_TOT NX3_TOT
#endif

  // ------------------------------------------------------------------
  // Handy macros for dimensional loops
  // ------------------------------------------------------------------

#define SDOM_LOOP(s) for ((s) = SBEG; (s) <= SEND; (s)++)
  
#if GEOMETRY == SPHERICAL
  #define FARGO_ARRAY_INDEX(A,nv,s,k,j,i)  A(i,j,s,nv)
#else
  #define FARGO_ARRAY_INDEX(A,nv,s,k,j,i)  A(i,s,k,nv)
#endif 

#if (GEOMETRY == CARTESIAN) || (GEOMETRY == CYLINDRICAL)
  #define FARGO_SPEED(w, i,j,k)  w[k][i]
#elif GEOMETRY == SPHERICAL
  #define FARGO_SPEED(w, i,j,k)  w[j][i]
#endif

  // ----------------------------------------------------------------- 
  // Prevent compilation when Fargo and AMR are given simultaneously
  // -----------------------------------------------------------------

void     FARGO_AverageVelocity(Data *, Grid *);
void     FARGO_ComputeResidualVelocity(const Data *, double ***, Grid *);
void     FARGO_ComputeTotalVelocity   (const Data *, double ***, Grid *);
void     FARGO_Initialize(Grid *);
void     FARGO_RestartOld(const Data *, Grid *);
void     FARGO_Restart(const Data *, char *, int, int, Grid *);
void     FARGO_ShiftSolution(Data *, Ary4D &, Ary3D (&Us)[6], Grid *, RunConfig &, double);
#if PARTICLES
void     FARGO_ShiftParticles(Data *, Grid *, double);
#endif
void     FARGO_Source(Ary4D &, Ary3D (&Us)[6], double, Grid *);
double **FARGO_Velocity(void);
void     FARGO_Write(const Data *, char *, int, Grid *);


