#include "grid_spec.hpp"
////////////////////////////////////////////////////////////////////////
/*!
  \file
  \brief PLUTO header file for structure declarations.

  \author A. Mignone (andrea.mignone@unito.it)
  \date   Dec 1, 2025
*/
////////////////////////////////////////////////////////////////////////

struct DomainDecomp {

  // -- Communicators --
  ParallelComm comm;             // base communicator (ParallelComm_WORLD)
  ParallelComm cart_comm;        // Cartesian topology communicator
  ParallelComm oned_comm[3];  // 1D sub-communicators along each direction

  // -- Datatype --
  FileDataType type;         // MPI_DOUBLE or MPI_FLOAT

  // -- Topology info --
  int ndim;                  // number of dimensions
  int rank;                  // this process's rank in cart_comm
  int size;                  // total number of processes

  int coords[3];          // this process's coordinates in Cartesian grid
  int cart_dims[3];       // number of processes along each direction
  int periods[3];         // periodic flags per direction
  int isparallel[3];      // which dimensions are decomposed

  // -- Global domain --
  int gsize[3];           // global number of cells per direction

  // -- Local domain (no ghosts) --
  int lsize[3];           // local number of cells per direction
  int lbeg[3];            // local start index (no ghosts)
  int lend[3];            // local end index (no ghosts)

  // -- Local domain (with ghosts) --
  int beg[3];             // local start index (with ghosts)
  int end[3];             // local end index (with ghosts)
  int ghosts[3];          // ghost zone size per direction

  // -- Global position of this subdomain --
  int gbeg[3];            // global start index of local domain
  int gend[3];            // global end index of local domain

  // -- Boundary flags --
  int is_gbeg[3];         // 1 if this rank is at the global start boundary
  int is_gend[3];         // 1 if this rank is at the global end boundary

  // -- Neighbor ranks (for ghost exchange) --
  int left[3];            // left neighbor rank per direction
  int right[3];           // right neighbor rank per direction

  // -- Staggered flag (for SZ_stagx/y/z) --
  int isstaggered[3];     // 1 if staggered in this direction

  // -- MPI derived types for ghost exchange --
  FileDataType type_lr[3]; // strided type for Left->Right send
  FileDataType type_rl[3]; // strided type for Right->Left send

  // -- MPI-IO subarrays --
  FileDataType gsubarr;       // global subarray for MPI-IO
  FileDataType lsubarr;       // local subarray for MPI-IO
};

////////////////////////////////////////////////////////////////////////
/*! The cmdLine structure includes information related to command
    line arguments given to PLUTO                                     */
////////////////////////////////////////////////////////////////////////
#pragma omp declare target
typedef struct cmdLine_{
  char restart;           /**< Enable restart from double precision binary files */
  char h5restart;         /**< Enable restart from hdf5 files   */
  char prestart;          /**< Enable / disable partice restart */
  char makegrid;          /**< Write grid only                  */
  char write;             /**< Set it 1 or 0 to enable or disable writing. */
  char parallel_dim[3];   /**< Enable/disable domain decomp. in a given direction */
  char input_file[32];    /**< Input file name */
  int  nrestart;          /**< The file number used for restart */
  int  maxsteps;          /**< The maximum number of steps (unless negative) */
  int  jet;               /**< Follow jet evolution in a given direction */
  int  nproc[3];          /**< User supplied number of processors */
  int  xres;              /**< Change the resolution via command line */
  char fill[26];            /**< useless, it makes the struct a power of 2 */ 
} cmdLine;
#pragma omp end declare target

////////////////////////////////////////////////////////////////////////
/*! The EMF structure is used to pull together all the information
    necessary to build / use the electromotive force used to update
    the staggered components of magnetic field.                       */
////////////////////////////////////////////////////////////////////////

#pragma omp declare target
typedef struct ElectroMotiveForce{

/*! \name Face-centered electric field components.
    Three-dimensional arrays storing the emf components computed
    at cell faces during the dimensional sweeps. */

/**@{ */
  Ary3D e[3][3]; /**< Store either flux or velocity components. 
                      The 1st index refers to the vector component while
                      the 2nd index indicates the sweep direction.
                      For instance, e[IDIR][JDIR] is the y-component of 
                      electric field stored at x1 faces. (see also ct_store.hpp). */
  Ary3D e_dff[3][3];
/**@} */

#if PHYSICS == ResRMHD
/*! \name Face-centered magnetic field components.
    Three-dimensional arrays storing the emf components computed
    at cell faces during the dimensional sweeps. */

/**@{ */
  Ary3D B[3][3];  
  Ary3D Frho[3];
/**@} */
#endif

  Ary3Dschar svx, svy, svz;
  Ary3Dschar sv[3];

/*! \name Signal velocities and emf coefficients */
 
/**@{ */
  Ary3D SxL, SxR;
  Ary3D SyL, SyR;
  Ary3D SzL, SzR;
  Ary3D dL[3], dR[3], aL[3], aR[3];
/**@} */

Ary3D Ee[3];  // Edge-averaged electric fields  ("e" = edge)  */
#if (PHYSICS == ResRMHD) || (PHYSICS == GRMHD)
Ary3D Be[3];  // Edge-averaged magnetic fields  ("e" = edge)  */
#endif

#if CT_ENERGY_CORRECTION > 1
  Ary3D Bx1ey, Bx1ez;
  Ary3D Bx2ex, Bx2ez;
  Ary3D Bx3ex, Bx3ey;
  Ary3D Ex1e_phi, Ex2e_phi, Ex3e_phi;
  Ary3D Ex3e_phix, Ex3e_phiy;
  Ary3D Pfx, Pfy, Pfz;
#endif

} EMF;
#pragma omp end declare target

// #pragma omp declare target
// typedef struct FLOAT_VECT{
//   float v1, v2, v3;
// } Float_Vect;
// #pragma omp end declare target

////////////////////////////////////////////////////////////////////////
/*! The PLUTO Grid structure contains information pertaining to the
    computational mesh in a specific 1D coordinate direction.
    Since PLUTO assumes a logically rectangular system of coordinates,
    the whole computational domain is obtained as the cartesian product
    of 2 or 3 grid structures.\n

    In parallel, each processor owns a different portion of the domain
    and the grid structures will be different.
    For this reason, in the following member description, we use the
    word "global" or "local" to refer the the whole computational domain
    or to the sub-domain owned by a single processor.

    Similarly, variables ending with a "glob" suffix are intended to be
    global, i.e., they refer to the whole computational stencil and
    not to the local processor sub-domain.                            */
////////////////////////////////////////////////////////////////////////
#pragma omp declare target
typedef struct Grid_{
  double xbeg[3], xbeg_glob[3]; /**< Leftmost acive domain interface (local, global)  */
  double xend[3], xend_glob[3]; /**< Rightmost acive domain interface (local, global) */
  double *xC[3], *xC_glob[3];   /**< Cell geometrical central points. */
  double *xL[3], *xL_glob[3];   /**< Left interface  */
  double *xR[3], *xR_glob[3];   /**< Right interface */
  double *dx[3], *dx_glob[3];   /**< Cell width      */

  #if MAPPED_GRID == YES
  double dxi[3];
  double *dxi_dx[3], *dxi_dx_glob[3];
  #endif

  double dl_min[3];   /**< Minimum cell length */

  Ary3D xgc[3];       /**< Cell volumetric centroid 
                           (!= x when geometry != CARTESIAN). */

  Ary3D dV;           /**< Cell volume.  */
  Ary3D A[3];         /**< Left interface area, A[i] = \f$A_{i+\HALF}\f$. */
  Ary3D L[3];         /**< Cell edge length, L[i] = \f$L_{i,j+\HALF,k+\HALF}\f$. */
  double  **dx_dl[3];      /**< dx/dl (dl = length), used for gradient-like operators */

  double *rt;         /**< In spherical coordinates, gives \tilde{r} */

  double *rC;         /**< Radius in cylindrical / spherical coords, |r|  */
  double *rL;         /**< Same as rC, but at left interface interfaces  */
  double *rC_inv;     /**< Inverse radius in cylindrical / spherical coords, 1/|r|  */
  double *rL_inv;     /**< Same as rC_inv, but at left interface interfaces  */

  double *sC;         /**< In spherical coordinates, gives |sin(th)|
                           at the cell center */
  double *sL;         /**< In spherical coordinates, gives |sin(th)|
                            at a j-1/2 interface */
  double *dmu;        /** < In spherical coordinates, gives the \theta
                            volume = fabs(cos(th_m) - cos(th_p)) */
  double *dsn;        /** < In spherical coordinates, gives 
                            volume = fabs(sin(th_p) - sin(th_m)) */
  double *cot;        /** < In spherical coordinates, gives the volume
                            integrated cotangent: d\sin /d\mu */
  int end_ring;       /**< In RING_AVERAGE mode, it gives the
                           index of the last active ring */
  double *inv_dx[3];  /**<      */
  double *inv_dxi[3]; /**< inverse of the distance between the center of 
                           two cells, inv_dxi = \f$\DS \frac{2}{\Delta x_i +
                           \Delta x_{i+1}}\f$.     */
  int lbound[3];      /**< When different from zero, it specifies the boundary
                           condition to be applied at leftmost grid side where
                           the physical boundary is located.
                           Otherwise, it equals zero if the current
                           processor does not touch the leftmost physical boundary.
                           This evantuality (lbound = 0) is possible only
                           in PARALLEL mode.  */
  int rbound[3];      /**< Same as lbound, but for the right edge of the grid. */
  int gbeg[3];        /**< Global start index for the global array */
  int gend[3];        /**< Global end   index for the global array */
  int beg[3];         /**< Global start index for the local  array */
  int end[3];         /**< Global end   index for the local  array */
  int lbeg[3];        /**< Local start  index for the local  array */
  int lend[3];        /**< Local end    index for the local  array */

  int ntot_glob[3];    /**< Total number of points in the global domain
                             (boundaries included)  */
  int nint_glob[3];    /**< Total number of points in the global domain
                             (boundaries excluded)  */
  int ntot[3];         /**< Total number of points in the local domain
                             (boundaries included)  */
  int ntot_stag[3][3]; /**< Total number of points in the local 
                            for staggered arrays (boundaries included). */
  int nint[3];         /**< Total number of points in the local domain
                            (boundaries excluded)  */

  int nghost[3];       /**< Number of ghost zones. */
  int nproc[3];        /**< Number of processors per direction            */
  int mg_type[3];      /**< Mapped grid type                              */
  int nmax;            /**< Max number of zones in the three directions   */
  int rank_coord[3];   /**< Parallel coordinate in a Cartesian topology   */
  int *ring_av_csize;  /**< The chunk size when RING_AVERAGE is turned on */
  char padding[24];    /**< Align to multiple of 128 */

} Grid;
#pragma omp end declare target

////////////////////////////////////////////////////////////////////////
/*! The List defines a collection of integer values typically used
    as argument to the FOR_EACH() macro.
    Variables are stored sequentially in such a way that structure
    elements can be conveniently initialized upon declaration:

    intList list = {4, MX1, MX2, MX3, ENG};

    list will have 4 elements so that, using the FOR_EACH macro
    a loop will be done on MX1, MX2, MX3, ENG.                        */
////////////////////////////////////////////////////////////////////////
#pragma omp declare target
typedef struct intList_{
  int nvar;       /**< Number of variables. */
  int indx[2046]; /**< Array of integers containg variables indices. */
  int i;          /**< Internal counter. */
} intList;
#pragma omp end declare target

#if PHYSICS != GRMHD
////////////////////////////////////////////////////////////////////////
/*! Metric structure: compute the scale factors                       */
////////////////////////////////////////////////////////////////////////
typedef struct Metric_{
    double h1 = 1.0;
    double h2 = 1.0;
    double h3 = 1.0;
    double h1_inv = 1.0;
    double h2_inv = 1.0;
    double h3_inv = 1.0;
} Metric;
#endif

#if PHYSICS == GRMHD
////////////////////////////////////////////////////////////////////////
/*! Structures to store the local (i.e. for a given point in space)
    metric coeffients and their spatial derivatives.                  */
////////////////////////////////////////////////////////////////////////
typedef struct Metric_{
    double lapse;
    double shift[3];
    double gCov[3][3];
    double gCon[3][3];
    double sqrtDet;
    double sqrtDet_1;
} Metric;

typedef struct MetricDer_{
    double lapseDer[3];
    double shiftDer[3][3];
    double gCovDer[3][3][3];
} MetricDer;
#endif

////////////////////////////////////////////////////////////////////////
/*! The Output structure contains essential information for I/O.         */
////////////////////////////////////////////////////////////////////////

typedef struct Output_{
  Ary4D *Vary4D;       /**< Pointer to the Ary4D class that includes the 
                            the first NVAR variables. */    
  Ary3D *Vary3D[MAX_OUTPUT_VARS]; /**< Array of pointers to different Ary3D */
  double ***Vc_buf_AL;    /**< Storage buffer for cell-centered data */ 
  double ***Vs_buf_AL[3]; /**< Storage buffer for staggered field data */
  double *Vc_buf;    /**< Storage buffer for cell-centered data */ 
  double *Vs_buf[3]; /**< Storage buffer for staggered field data */
  double dt;           /**< Time increment between outputs. */
  double dclock;       /**< Time increment in clock hours. */

  int    type;         /**< Output data format (DBL, FLT, VTK, ...). */
  int    nvar;         /**< (Fluid only) Total # of vars that can potentially 
                            be written. This is the same for all kind of outputs  */
  int    cgs;          /**< (Fluid only) When set to 1, save data in c.g.s units  */
  int    nfile;        /**< Current number being saved. */
  int    nfileParticles; /**< Current number being saved. */
  int    dn;           /**< Step increment between outputs. */
  int    stag_var[MAX_OUTPUT_VARS];  /**< (Fluid only). Centered or staggered
                                           variable - same for all outputs. */
  int    dump_var[MAX_OUTPUT_VARS];  /**< (Fluid only) Include/exclude variables
                                           being written.  */
  int    field_dim[MAX_OUTPUT_VARS]; /**< (Particle only) The dimensionality
                                          of the field being written
                                          (scalar = 1, array > 1) */
  int    user_var;     /**< Number of user-defined variables */
  char   mode[32];     /**< (Fluid only) Single or multiple files. */
  char   **var_name;   /**< (Fluid only) Variable names. Same for all output types.  */
  char   ext[8];       /**< File extension (.flt, .dbl, etc...)           */
  char   dir[256];     /**< Output directory name                        */

  char   fill[140];    /**< Useless, just to make the structure size a power of 2 */
} Output;

////////////////////////////////////////////////////////////////////////
/*! The RBox (= Rectangular Box) defines a rectangular portion of the
    domain in terms of the grid indices <tt>[ibeg,jbeg,kbeg]</tt> corresponding
    to the lower corner and <tt>[iend,jend,kend]</tt> corresponding to the
    upper corner.                                                     */
////////////////////////////////////////////////////////////////////////
#pragma omp declare target
typedef struct RBox_{
  int ibeg; /**< Lower corner index in the x1 direction */
  int iend; /**< Upper corner index in the x1 direction */
  int jbeg; /**< Lower corner index in the x2 direction */
  int jend; /**< Upper corner index in the x2 direction */
  int kbeg; /**< Lower corner index in the x3 direction */
  int kend; /**< Upper corner index in the x3 direction */

  int nx1;  /**< Size in the x1 direction */
  int nx2;  /**< Size in the x2 direction */
  int nx3;  /**< Size in the x3 direction */
} RBox;
#pragma omp end declare target

////////////////////////////////////////////////////////////////////////
/*! The restart structure contains restart information that must be
    read from disk to restart PLUTO.

    Important: The restart structure should be aligned to a power of 2 to
    prevent (for some compilers) from changing the alignment of the
    structure and therefore troubleshooting when restarting
    from files written on different architectures.                    */
////////////////////////////////////////////////////////////////////////
typedef struct Restart_{
  double t;
  double dt;
  double param[MAX_RESTART_PARAMS];
  int    nstep;
  int    nfile[MAX_OUTPUT_TYPES];
  int    nfileParticles[MAX_OUTPUT_TYPES];
  char   fill[10];  // Align the structure to power of 2 
} Restart;

////////////////////////////////////////////////////////////////////////
/*! The runConfig structure contains relevant quantities for the   
    specific problem.                                                 */
////////////////////////////////////////////////////////////////////////
#pragma omp declare target
typedef struct RunConfig_{ 
  double inputParam[32];  /**< Array of input parameter values read from 
                               pluto.ini */

  double gamma;               /**< The adiabtic index for the ideal EoS  */
  double isoSoundSpeed = 1.0; /**< The isothermal sound speed (for isothermal EOS) */
  double time;                /**< Simulation time */
  double domBeg[3];           /**< Lower limits of the computational domain. */
  double domEnd[3];           /**< Upper limits of the computational domain. */
  double OmegaZ;              /**< Rotating frame angular frequency */
  double small_prs = 1.e-12;
  double small_rho = 1.e-12;
  double unit_density  = CONST_mp;  /**< Unit density (in gr/cm^3) */
  double unit_velocity = CONST_c;   /**< Unit velocity (in cm/sec) */
  double unit_length   = CONST_au;  /**< Unit length   (in cm)     */
  double unit_time;                 /**< Unit time (in s) computed in startup */
  double unit_Bmag;                 /**< Unit magnetic field (in cgs)
                                         (compute in startup)        */
  double *x[3];
  double *dx[3];

  double glm_ch = -1.0;

  #if COOLING != NO
    /** < Values for COOLING functions */
    double minCoolingTemp;        /** Minimum Cooling Temperature */
    double maxCoolingRate = 0.1;  /** Maximum Cooling Rate (default) */
    double smallPressure;         /** minimum pressure value (used as control value) */
    double smallDensity;          /** minimum density value (used as control value)*/

    double *L_tab, *T_tab;
    
    int     ntab; /** "static" number of rows read in radiat.cpp */
  #endif

  #if RADIATION
  double rad_constant = 1.0;         /**< Radiation constant in code units */
  double rad_ideal_gas_constant = 1.0;   /**< Pressure-to-temperature conversion for radiation */
  double rad_absorption_coeff = 0.0; /**< Constant absorption opacity for radiation */
  double rad_scattering_coeff = 0.0; /**< Constant scattering opacity for radiation */
  double rad_c = 1.0;                /**< Speed of light in code units */
  double rad_reduced_c = 1.0;        /**< Reduced speed of light (c_r <= c) */
  int    rad_max_iter;         /**< Max iterations for the radiation implicit step (\c rad_max_iter) */
  double rad_tolerance;        /**< Convergence tolerance for the radiation implicit step (\c rad_tolerance) */
  double small_enr;         /**< Minimum radiation energy density if enr<0 (\c small_enr) */
  #endif

  int     Nx1;
  int     Nx2;
  int     Nx3;
} RunConfig;
#pragma omp end declare target

////////////////////////////////////////////////////////////////////////
/*! The Runtime structure contains runtime initialization parameters
    read from pluto.ini (or equivalent).                              */
////////////////////////////////////////////////////////////////////////
#pragma omp declare target
typedef struct Runtime_{
  GridSpec grid_spec[3];    /**< Optional PLUTO-style grid patches. */
  int    npoint[3];           /**< Global number of zones in the interior domain */
  int    left_bound[3];       /**< Array of left boundary types */
  int    right_bound[3];      /**< Array of right boundary types */
  int    grid_is_uniform[3];  /**< = 1 when grid is uniform, 0 otherwise */
  int    log_freq;            /**< The log frequency (\c log) */
  int    user_var;            /**< The number of additional user-variables being
                                 held in memory and written to disk */
  int    anl_dn;              /**< Number of step increment for ANALYSIS */
  int    Panl_dn;             /**< Number of step increment for PARTICLE ANALYSIS */
  char   solv_type[MAX_STRING_SIZE];     /**< The Riemann solver (\c Solver) */
  char   rad_solv_type[MAX_STRING_SIZE]; /**< The radiation Riemann solver (\c Solver) */
  char   user_var_name[MAX_OUTPUT_VARS][MAX_STRING_SIZE]; /**< Names of user-defined
                                                               variables */
  char   output_dir[MAX_STRING_SIZE];  /**< The name of the output directory.
                                             Default is current directory.
                                            (\c output_dir for static PLUTO,
                                             \c Output_dir for PLUTO-Chombo)  */
  char   log_dir[MAX_STRING_SIZE]; /**< The name of the output directory
                                        where log files will be written to.
                                        Default is output_dir. */
  Output output[MAX_OUTPUT_TYPES];
  double xbeg[3];            /**< self-expl. */
  double xend[3];            /**< self-expl. */
  double  cfl;             /**< Hyperbolic cfl number (\c CFL) */
  double  cfl_max_var;     /**< Maximum increment between consecutive time
                                  steps (\c CFL_max_var). */
  double  cfl_par;         /**< (STS) parabolic  cfl number */
  double  cfl_grav;        /**< gravity CFL > for Self Gravity module **/
  double  rmax_par;        /**< (STS) max ratio between current time
                                  step and parabolic time step */
  double particles_tstart; /**< A copy of runtime->particles_tstart */
  
  double  tstop;           /**< The final integration time (\c tstop) */
  double  tfreeze;         /**< The fluid freezing time  (\c tfreeze) */
  double  first_dt;        /**< The initial time step (\c first_dt) */
  double  rad_first_dt;    /**< The initial radiation sub-step (\c rad_first_dt) */
  double  rad_cfl;         /**< Radiation CFL number (\c rad_CFL) */
  double  rad_cfl_max_var; /**< Maximum increment between consecutive radiation
                                sub-steps (\c rad_CFL_max_var). */
  int     rad_max_iter;    /**< Max iterations for the radiation implicit step (\c rad_max_iter) */
  double  rad_tolerance;   /**< Convergence tolerance for the radiation implicit step (\c rad_tolerance) */
  double  small_enr;    /**< Minimum radiation energy density (\c small_enr) */
  double  anl_dt;          /**< Time step increment for Analysis()
                                ( <tt> analysis (double) </tt> ) */
  double  Panl_dt;         /**< Time step increment for ParticleAnalysis()
                                ( <tt> analysis (double) </tt> ) */

  int     Nparticles_per_block; /**< Total number of particles per cell */
  int     Ncells_per_block;     /**< Initial grid particle size 
                                     (multiple of base grid) */
  double  aux[32];         /**< We keep aux inside this structure,
                               since in parallel execution it has
                               to be comunicated to all processors */ 
} Runtime;
#pragma omp end declare target

////////////////////////////////////////////////////////////////////////
/*! The State_ struct defines a fluid state at a given point in space.
    Useful when passing several arguments between functions.          */
////////////////////////////////////////////////////////////////////////
#pragma omp declare target
typedef struct State_{
  double *v;       /**< Array of primitive variables */
  double *u;       /**< Array of cons. variable */
  double *f;       /**< Array of fluxes */
  double  p;       /**< Pressure */
  double *Bgf;     /**< Background field         */
  double x1;       /**< Pointer to x1 coordinate */
  double x2;       /**< Pointer to x2 coordinate */
  double x3;       /**< Pointer to x3 coordinate */
  double a2;       /**< Sound speed squared      */
  double h;        /**< Enthalpy                 */
  double pm;       /**< Magnetic energy (testing)   */ 
  double w = 0.0;  /**< Orbital Velocity (FARGO + ROTATIONS)   */ 
  double eta;      /**< Magnetic resistivity */
  double divB;     /**< Local divergence of B */
  uint32_t *flag;  /**< Pointer to local flag */
  int    i = -1;   /**< Zone x1-grid index  */
  int    j = -1;   /**< Zone x2-grid index  */
  int    k = -1;   /**< Zone x3-grid index  */
} State; 
#pragma omp end declare target

////////////////////////////////////////////////////////////////////////
/*! This structure contains one-dimensional vectors of conserved
    variables, primitive variables, fluxes and so on, used during
    the reconstruct-Solve-Average strategy.
    It is a frequently passed to the Riemann solver routines, source and
    flux functions, etc.                                              */
////////////////////////////////////////////////////////////////////////
#pragma omp declare target
typedef struct Sweep_{
  double *flux;     /**< Upwind flux computed with the Riemann solver */
  double *vC;       /**< Array of primitive variables at cell center   */
  double *uC;       /**< Array of primitive variables at cell center   */

  double *vL;       /**< Array of primitive variables at left of interf.   */
  double *vR;       /**< Array of primitive variables at right of interf.   */
  double *fL;
  double *fR;

  double *cs2;      /**< Sound speed squared */
  double *rhs;      /**< Conservative right hand side */
  double *press;    /**< Upwind pressure term computed with the Riemann solver */
  double *pL;
  double *pR;

  double *Bn;       /**< Face-centered magentic field, e.g., Bn = Bx(i+1/2) */
  double *En;       /**< Face-centered electric field, e.g., En = Ex(i+1/2) */
  double *PSIn;     /**< */
  double *SL;       /**< Leftmost  velocity in the Riemann fan at i+1/2 */
  double *SR;       /**< Rightmost velocity in the Riemann fan at i+1/2 */

  #if RADIATION
  double *SrL;      /**< Leftmost interface velocity for radiation fields */
  double *SrR;      /**< Rightmost interface velocity for radiation fields */
  #endif

  #if (SELF_GRAV == YES)
  double* g_c;     /**<gravitational accleration */
  #endif

  double *SaL, *SaR, *Sc; /**< MHD alfven waves, contact wave */
  double *dL, *dR;        /**< Diffusion coefficient for EMF  */
  double *aL, *aR;        /**< Flux averaging coefficients    */
  double **u;
  double **uXB;
  double **lambda;        /**< Conservative right hand side */

  uint32_t *flag;
  char fill[29];
} Sweep;
#pragma omp end declare target

////////////////////////////////////////////////////////////////////////
/*! The timeStep structure contains essential information for
    determining the time step.                                        */
////////////////////////////////////////////////////////////////////////
#pragma omp declare target
typedef struct timeStep_{
  double *cmax;       /**< Maximum signal velocity for hyperbolic eqns */
  double *Dmax;       /**< Maximum diffusion coefficient for parabolic terms */
  double invDt_hyp;   /**< Inverse of hyperbolic time step,
                           \f$ \lambda/\Delta l\f$.*/
  double invDt_par;   /**< Inverse of parabolic (diffusion)  time step
                           \f$ \eta/\Delta l^2\f$. */
  #if COOLING != NO
  double dt_cool;     /**< Cooling time step. */
  #endif
  
  double cfl;         /**< Courant number for advection. */
  double cfl_par;     /**< Courant number for diffusion (STS only). */
  double rmax_par;
  double clock_particles;
  double clock_particles_bound;
  double cfl_grav;    /**< Courant number for gravitational stability */
  double clock_hyp;
  double clock_par;
  double clock_cooling;
  double clock_tot;

  double time;        /**< simulation time */
  double dt;          /**< simulation time increment */
  double dt_physics = 0.0; /**< Latest NextTimeStep bound before growth/event limits; 0 until computed. */
  int    stepNumber;  /**< simulation time step */
  int    intStage;    /**< simulation integration stage */
  int    Nsub_particles; /**< Number of sub-cycles in particles */
  int    Nsts;        /**< Maximum number of substeps used in STS. */
  int    Nrkc;        /**< Maximum number of substeps used in RKC. */
  int    Nrkl;        /**< Maximum number of substeps used in RKL. */

  double particles_tstart; /**< A copy of runtime->particles_tstart */
  double invDt_particles;  /**< Max inverse dt for particles */
  double omega_particles;  /**< Max Larmor frequency for particles */

} timeStep;
#pragma omp end declare target

////////////////////////////////////////////////////////////////////////
/*! The Data structure contains the main solution 3D arrays used by
    the code.                                                         */
////////////////////////////////////////////////////////////////////////
#pragma omp declare target
struct Data{
  Ary4D Vc;      /**< The main four-index data array used for cell-centered
                      primitive variables.  */
  Ary4D Uc;      /**< The main four-index data array used for cell-centered
                      conservative variables. */
  
  Ary4D Uc0;    /**< Copy of Uc at the beginning of the step */
  Ary4D Usafe;   /**< Safe solution                           */

  Ary4D Upc;     /**< The main four-index data array used for cell-centered
                      punctual conservative variables.  */
  Ary4D Vpc;     /**< The main four-index data array used for cell-centered
                      punctual primitive variables.     */
  Ary3D Vs[6];   /**< Main arrays for staggered data (B and E for ResRMHD)  */
  Ary3D Vps[6];  /**< The main four-index data array used for face-centered
                         punctual staggered variables.*/
  Ary3Duint32_t ho_lim;   /**< Discontinuity detector switch */

  Ary3D A[3];     /**< Vector potential */
  Ary3D J[3];     /**< Electric current */
  Ary3D divB;     /**< Divergence of the magnetic field.*/
  Ary3D q;        /**< Electric charge density (only for ResRMHD)    */
  Ary3D C_dt;
  Ary3D C_dtpar[3];
  Ary3D phi_m;
  Ary3Duint32_t flag; /**< 3D Array of integration flags   */

  // EMF 
  Ary3D Ex1;      /**< cell-centered emf used in CT averaging or CR particles */
  Ary3D Ex2;      /**< cell-centered emf used in CT averaging or CR particles */
  Ary3D Ex3;      /**< cell-centered emf used in CT averaging or CR particles */

  Ary3D eta[3];   /**< resistivity */
  Ary3D Ex1Prv;   /**< cell-centered emf used in CT averaging or CR particles */
  Ary3D Ex2Prv;   /**< cell-centered emf used in CT averaging or CR particles */
  Ary3D Ex3Prv;   /**< cell-centered emf used in CT averaging or CR particles */

  // Error diagnostic members
  Ary3Duint32_t err;      /**< 3D Array used for cell-by-cell error detection & 
                               diagnostic in each kernel.  */

  #if (POISSON_SOLVER == YES)
    Ary3D poisson_sol; /**< Array storing the solution of the Poisson solver **/
    #if (SELF_GRAV == YES)
    Ary4D acc_grav;    /**< Array storing the gravitational acceleration */
    #endif
    #if (POISSON_DEBUG == YES) 
    Ary3D err_res;
    #endif 
  #endif

  #if RING_AVERAGE != NO
  #if RING_AVERAGE_FOURIER == YES
  Ary2D A0[3], A1[3], B1[3]; 
  #endif
  Ary3D Phi[3], PhiM[3];
  #endif

  uint64_t err_count = 0; /**< Count total number of errors during the step  */

  // -- Particles-related quantities -- 

  class ParticlesContainer *pc;

  double ****Fdust; /**< Drag force (dust particles only)   */

  struct ElectroMotiveForce *emf;

  struct timeStep_  *Dts;
  struct RunConfig_ *run_config;
  struct Output_    *output;

  void (*fluidRiemannSolver) (struct Data *, timeStep *, Grid *, RBox *, int);
  #if RADIATION
  void (*radiationRiemannSolver)  (struct Data *, timeStep *, Grid *, RBox *, int);
  #endif
  struct Sweep_ *sweep;

  double max_mach;

  double restartParam[MAX_RESTART_PARAMS];
  signed char integration_mode;  /**< Sets the current integration mode (DEFAULT/FALLBACK) */
  signed char use_flat = 0;      /**< Reduce the reconstruction at runtime 
                                      by switching to flat reconstruction */
  char padding[175];  // Align to multiple of 128 
};
#pragma omp end declare target
