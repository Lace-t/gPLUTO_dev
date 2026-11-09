////////////////////////////////////////////////////////////////////////
/*!
  \file
  \brief PLUTO header file for function-like macros.

  \author A. Mignone (andrea.mignone@unito.it)\n
          M. Rossazza
  \date   Dec 25, 2025
*/
////////////////////////////////////////////////////////////////////////

#define ABS(a)  ((a) > 0.0 ? (a):(-1.0*a))
/*! Return the number with the smaller absolute value. */
#define ABS_MIN(a,b)  (fabs(a) < fabs(b) ? (a):(b))

/*!
    ARRAY_STRIDES(grid, dir, nX, nY)

    Set stride parameters used to linearize multidimensional grid data
    into contiguous 1D sweep arrays.

    The macro defines:
      nX = total number of points along the sweep direction `dir`
           (length of a single 1D sweep)
      nY = total number of points along the first transverse direction

    These parameters are used by sweep routines to compute the linear
    offset of a given sweep line as:

      offset = nX * ( t1 + nY * t2 )

    where (t1, t2) are the two indices transverse to the sweep direction.

    Current ordering used in sweep buffers:

      IDIR (sweep along i):
        offset = NX * ( j + NY * k )

      JDIR (sweep along j):
        offset = NY * ( i + NX * k )

      KDIR (sweep along k):
        offset = NZ * ( i + NX * j )

    This layout stores all 1D sweeps contiguously in memory, enabling
    direction-independent algorithms and efficient execution on CPUs,
    GPUs, and MPI buffers.

    NOTE:
    The macro does not perform the full 3D-to-1D mapping by itself;
    it only provides stride information required by the calling code. */

#define ARRAY_STRIDES(grid, dir, nX, nY)   \
    nX = grid->ntot[dir];                  \
    nY = (dir == IDIR ? grid->ntot[JDIR] : grid->ntot[IDIR])

/*!
  Convert 3D grid coordinates (i,j,k) into a 1D linear array index
  suitable for flattened contiguous storage.
  The mapping assumes a structured box with row-major ordering:
  i varies fastest, followed by j, then k.
 
  Equivalent mapping <tt> 
  index = (i - ibeg) + nx1 * ( (j - jbeg) + nx2 * (k - kbeg) )
  </tt>
 */
#define BOX_INDEX(box,i,j,k) ( (i - (box).ibeg) + (box).nx1*  \
                             ( (j - (box).jbeg) + (box).nx2*(k - (box).kbeg)))

/*! The BOX_LOOP() macro implements a loop over (i,j,k) in a rectangular
    portion of the domain with indices defined by the (pointer to)
    RBox structure B.
    The BOX_STAG_LOOP() loops over all the edges of the domain and it can be
    also used to loop over staggered fields with a minimal extra cost.
    The BOX{I,J,K}_LOOP() allows to travel across the box by including one extra
    point in the direction of staggering. */

#define BOX_LOOP(box,i,j,k) for (k = (box).kbeg; k <= (box).kend; k++) \
                             for (j = (box).jbeg; j <= (box).jend; j++) \
                             for (i = (box).ibeg; i <= (box).iend; i++)

#define BOX_STAG_LOOP(box,i,j,k)  for (k = (box).kbeg; k <= (box).kend + INCLUDE_KDIR; k++) \
                                   for (j = (box).jbeg; j <= (box).jend + INCLUDE_JDIR; j++) \
                                   for (i = (box).ibeg; i <= (box).iend + INCLUDE_IDIR; i++)

#define BOXI_LOOP(box,i,j,k)  for (k = (box).kbeg; k <= (box).kend;   k++) \
                               for (j = (box).jbeg; j <= (box).jend;   j++) \
                               for (i = (box).ibeg; i <= (box).iend+1; i++)

#define BOXJ_LOOP(box,i,j,k)  for (k = (box).kbeg; k <= (box).kend;   k++) \
                               for (j = (box).jbeg; j <= (box).jend+1; j++) \
                               for (i = (box).ibeg; i <= (box).iend;   i++)

#define BOXK_LOOP(box,i,j,k)  for (k = (box).kbeg; k <= (box).kend+1; k++) \
                               for (j = (box).jbeg; j <= (box).jend;   j++) \
                               for (i = (box).ibeg; i <= (box).iend;   i++)

/*! The FOR_EACH(p, intList) macro implements a loop over the elements
    of the array \c intList->indx (in a similar way to Python lists).

    Example:
    \code
      intList list;
      list.nvar = 3;
      list.indx[0] = 2;
      list.indx[1] = 5;
      list.indx[2] = 17;
      FOR_EACH(nv, &list) printf ("value is = %d\n",nv);
    \endcode
*/
#define FOR_EACH(nv, list)  \
  for ((list)->i = 0, nv = (list)->indx[0]; \
       (list)->i < (list)->nvar; \
       nv = (list)->indx[++((list)->i)])

/*! 
  @brief A one-line conditional macro.

  The IF(id, a) macro expands to 'a' if id == 1, and to nothing if id == 0.
  This is useful for enabling or disabling code fragments at compile time 
  based on preprocessor definitions.
 */
#define IF_1(a) a
#define IF_0(a) 
#define IF(id, a)  IF_I(id, a)
#define IF_I(id, a) IF_##id(a)

/*!
  @brief Cycles three indices based on a direction value.
 
  This macro reassigns the values of three variables (`e1`, `e2`, and `e3`) 
  based on the value of `dir`, using a modulo-3 increment.
  
  @param dir The direction value used to offset the indices.
  @param e1 The first index value (input and output).
  @param e2 The second index value (output).
  @param e3 The third index value (output).
 
  @endcode
 */
#define INDEX_CYCLE(dir, e1, e2, e3) \
        e3 = e1 + (dir+2)%3;  \
        e2 = e1 + (dir+1)%3;  \
        e1 = e1 + (dir+0)%3;

#define GEN_INDEX_CYCLE(dir, base, e1, e2, e3) \
        constexpr int e3 = base + (dir+2)%3;  \
        constexpr int e2 = base + (dir+1)%3;  \
        constexpr int e1 = base + (dir+0)%3;

// -------------------------------------------------------------
// 1. Macros that can be placed at any point in the code
// -------------------------------------------------------------
   
// Math Macros fot Self Gravity Module 

#define POW2(x)   ((x)*(x))
#define POW3(x)   ((x)*(x)*(x))
#define POW4(x)   ((x)*(x)*(x)*(x))
#define SECH(x) ((x) > 30.0 ? 0.0:1.0/cosh(x))

////////////////////////////////////////////////////////////////////////
/*! Variable looping macros                                           */
////////////////////////////////////////////////////////////////////////

/**@{ */
#define NFLX_LOOP(n)  for ((n) = 0;        (n) < NFLX;               (n)++)
#define NIONS_LOOP(n) for ((n) = IONS_BEG; (n) < (IONS_BEG + NIONS); (n)++)
#define NSCL_LOOP(n)  for ((n) = NFLX;     (n) < (NFLX + NSCL);      (n)++)
#define NTRC_LOOP(n)  for ((n) = TRC;      (n) < (TRC  + NTRACER);   (n)++)
#define NVAR_LOOP(n)  for ((n) = 0;        (n) < NVAR; (n)++)

// NFLUID_LOOP loops over radiation fields only in the relativistic case
#if RADIATION
#define NRAD_LOOP(n)   for ((n) = ENR;      (n) < (ENR + 4);       (n)++)
#if RADIATION_RELATIVISTIC
#define NFLUID_LOOP(n) NVAR_LOOP(n)
#else
#define NFLUID_LOOP(n) for ((n) = 0;        (n) < (NFLX + NSCL);      (n)++)
#endif
#else
#define NFLUID_LOOP(n) NVAR_LOOP(n)
#endif
/**@} */

/*! Faster implementation than stdlib floor() function.
    It returns the largest integer value less than or equal to z. */
#define INT_FLOOR(z)   ((int)((z) + 32768.) - 32768)

/*! Return the maximum between two numbers. */
#define MAX(a,b)  ( (a) >= (b) ? (a) : (b) )

/*! Return the minimum between two numbers. */
#define MIN(a,b)  ( (a) <= (b) ? (a) : (b) )

/*! Return the sign of x. */
#define DSIGN(x)      ( (x) >= 0.0 ? (1.0) : (-1.0))

/*! Quick limiting macros (more general ones are found in States/plm_coeffs.hpp) */
#define MINMOD_LIMITER(a,b)  ((a)*(b) > 0.0 ? (fabs(a) < fabs(b) ? (a):(b)):0.0)
#define VANLEER_LIMITER(a,b) ((a)*(b) > 0.0 ? 2.0*(a)*(b)/((a)+(b)):0.0)
#define MC_LIMITER(a,b)      (MINMOD_LIMITER(0.5*((a)+(b)), 2.0*MINMOD_LIMITER((a),(b))))
#define SWAP_VAR(x) SwapEndian(&x, sizeof(x));

// -------------------------------------------------------------
// 2. Macros that must be placed only *AFTER* definitions.hpp has been
//    included
// -------------------------------------------------------------

////////////////////////////////////////////////////////////////////////
/* Expand dimension- or component-dependent expressions               */
////////////////////////////////////////////////////////////////////////

/*! \def DIM_EXPAND(a,b,c)
    Allows to write dimension-ndependent expressions involving vectors
    by evaluating as many arguments as the number of DIMENSIONS.
    The result is that only the first argument will be compiled in 1D,
    the first two arguments in 2D and all of them in 3D.
    As an example:
    \code
    DIM_EXPAND(v[VX1] =  1.0;  ,
               v[VX2] =  5.0;  ,
               v[VX3] = -4.0; )
    \endcode
    becomes
    \code
     v[VX1] = 1.0;
    \endcode
    when \c DIMENSIONS is equal to 1 or
    \code
     v[VX1] = 1.0;
     v[VX2] = 5.0;
    \endcode
    when \c DIMENSIONS is equal to 2 or
    \code
     v[VX1] =  1.0;
     v[VX2] =  5.0;
     v[VX3] = -4.0;
    \endcode
    when \c DIMENSIONS is equal to 3. */

/*! \def DIM_EXPAND(a,b,c), DIM_SELECT(a,b,c) \& INCLUDE_XDIR
    Expand only the 1st, 2nd or 3rd argument based on the value of
    DIMENSIONS. */

/*! \def INCLUDE_IDIR, INCLUDE_JDIR \& INCLUDE_KDIR macros */

#ifndef INCLUDE_IDIR
  #define INCLUDE_IDIR  (DIMENSIONS >= 1)
#endif

#ifndef INCLUDE_JDIR
  #define INCLUDE_JDIR  (DIMENSIONS >= 2)
#endif

#ifndef INCLUDE_KDIR
  #define INCLUDE_KDIR  (DIMENSIONS == 3)
#endif

#if INCLUDE_IDIR
  #define DIM1(a)     a
  #define NOT_DIM1(a) 
#else
  #define DIM1(a)
  #define NOT_DIM1(a) a
#endif

#if INCLUDE_JDIR
  #define DIM2(a)     a
  #define NOT_DIM2(a) 
#else
  #define DIM2(a)
  #define NOT_DIM2(a) a
#endif

#if INCLUDE_KDIR
  #define DIM3(a)     a
  #define NOT_DIM3(a) 
#else
  #define DIM3(a)
  #define NOT_DIM3(a) a
#endif

/*! \def DIM_EXPAND(a,b,c) \& DIM_SELECT(a,b,c)
    Expand only the 1st, 2nd or 3rd argument based on the value of
    DIMENSIONS. */

#define DIM_EXPAND(a, b, c)    DIM1(a) DIM2(b) DIM3(c)
#define DIM_EXPAND_NOT(a,b,c)  NOT_DIM1(a) NOT_DIM2(b) NOT_DIM3(c)

#if DIMENSIONS == 1
  #define DIM_SELECT(a, b, c) a
#elif DIMENSIONS == 2
  #define DIM_SELECT(a, b, c) b
#elif DIMENSIONS == 3
  #define DIM_SELECT(a, b, c) c
#else
  #error "Unsupported DIMENSIONS value"
#endif

// #if DIMENSIONS == 3
//   #if !INCLUDE_JDIR
//   #define DIM_LOOP(d) for ((d)=0; (d)<DIMENSIONS; (d) += 2)
//   #else
//   #define DIM_LOOP(d) for ((d)=0; (d)<DIMENSIONS; (d) += 1)
// #endif

#define DIM_LOOP(d) for ((d)=0; (d)<DIMENSIONS; (d)+=(DIMENSIONS==3 && INCLUDE_JDIR==0 ? 2 : 1))

#define DOT_PRODUCT(a,b)  ((a)[0]*(b)[0] + (a)[1]*(b)[1] + (a)[2]*(b)[2])

#define CROSS_PRODUCT(S, A, B) \
        S[IDIR] = ((A)[JDIR] * (B)[KDIR] - (A)[KDIR] * (B)[JDIR]); \
        S[JDIR] = ((A)[KDIR] * (B)[IDIR] - (A)[IDIR] * (B)[KDIR]); \
        S[KDIR] = ((A)[IDIR] * (B)[JDIR] - (A)[JDIR] * (B)[IDIR]);

/*! \name Forward and central finite differences macros.
    The following set of macros provide a compact way to perform
    two-point, undivided finite difference operations in a specified
    direction.
    Differences can be either \c forward or \c central.
    For instance, \c FDIFF_X2(Q,i,j,k) will compute a forward
    difference of \c Q in the \c y direction, \c (Q[j+1]-Q[j]), while
    \c CDIFF_X3(Q,i,j,k) will compute a central different approximation
    of \c Q in the \c z direction: \c (Q[k+1] - Q[k-1])/2. */

/**@{ */
#define OoT 0.08333333333
#if INCLUDE_IDIR
  #define FDIFF_X1(Q,i,j,k,var)      (Q(i+1,j,k,var) - Q(i,j,k,var))
  #define CDIFF_X1(Q,i,j,k,var)  (0.5*(Q(i+1,j,k,var) - Q(i-1,j,k,var)))
  #define CDIFF4_X1(Q,i,j,k,var)  (OoT*(Q(i-2,j,k,var) + 8.0*(Q(i+1,j,k,var) - Q(i-1,j,k,var)) - Q(i+2,j,k,var)))
#else
  #define FDIFF_X1(Q,i,j,k,var)  (0.0)
  #define CDIFF_X1(Q,i,j,k,var)  (0.0)
  #define CDIFF4_X1(Q,i,j,k,var)  (0.0)
#endif

#if INCLUDE_JDIR
  #define FDIFF_X2(Q,i,j,k,var)       (Q(i,j+1,k,var) - Q(i,j,k,var))
  #define CDIFF_X2(Q,i,j,k,var)  (0.5*(Q(i,j+1,k,var) - Q(i,j-1,k,var)))
  #define CDIFF4_X2(Q,i,j,k,var)  (OoT*(Q(i,j-2,k,var) + 8.0*(Q(i,j+1,k,var) - Q(i,j-1,k,var)) - Q(i,j+2,k,var)))
#else
  #define FDIFF_X2(Q,i,j,k,var)  (0.0)
  #define CDIFF_X2(Q,i,j,k,var)  (0.0)
  #define CDIFF4_X2(Q,i,j,k,var)  (0.0)
#endif

#if INCLUDE_KDIR
  #define FDIFF_X3(Q,i,j,k,var)       (Q(i,j,k+1,var) - Q(i,j,k,var))
  #define CDIFF_X3(Q,i,j,k,var)  (0.5*(Q(i,j,k+1,var) - Q(i,j,k-1,var)))
  #define CDIFF4_X3(Q,i,j,k,var)  (OoT*(Q(i,j,k-2,var) + 8.0*(Q(i,j,k+1,var) - Q(i,j,k-1,var)) - Q(i,j,k+2,var)))
#else
  #define FDIFF_X3(Q,i,j,k,var)  (0.0)
  #define CDIFF_X3(Q,i,j,k,var)  (0.0)
  #define CDIFF4_X3(Q,i,j,k,var)  (0.0)
#endif
/**@} */

/*! \name Spatial averages macros.
    The following set of macros provide a compact way to perform multi-D
    averages from cell centered values to interfaces.
    For instance, \C AVERAGE_X(q,k,j,i) will simply take the
    arithmetic average betwen q(i) and q(i+1) at the i+1/2 interface.
    Likewise, AVERAGE_YZ(q,k,j,i) will produce an average at the
    j+1/2 and k+1/2 edge. */ 

/**@{ */
#if INCLUDE_IDIR
  #define AVERAGE_X(q,i,j,k)   (0.5*(q(i,j,k) + q(i+1,j,k)))
#else
  #define AVERAGE_X(q,i,j,k)   (q(i,j,k))
#endif

#if INCLUDE_JDIR
  #define AVERAGE_Y(q,i,j,k)   (0.5*(q(i,j,k) + q(i,j+1,k)))
  #define AVERAGE_XY(q,i,j,k)  (0.5*(AVERAGE_X(q,i,j,k) + AVERAGE_X(q,i,j+1,k)))
#else
  #define AVERAGE_Y(q,i,j,k)    (q(i,0,k))
  #define AVERAGE_XY(q,i,j,k)   AVERAGE_X(q,i,0,k)
#endif

#if INCLUDE_KDIR
  #define AVERAGE_Z(q,i,j,k)    (0.5*(q(i,j,k) + q(i,j,k+1)))
  #define AVERAGE_XZ(q,i,j,k)   (0.5*(AVERAGE_X(q,i,j,k) + AVERAGE_X(q,i,j,k+1)))
  #define AVERAGE_YZ(q,i,j,k)   (0.5*(AVERAGE_Y(q,i,j,k) + AVERAGE_Y(q,i,j,k+1)))
  #define AVERAGE_XYZ(q,i,j,k)  (0.5*(AVERAGE_XY(q,i,j,k) + AVERAGE_XY(q,i,j,k+1)))
#else
  #define AVERAGE_Z(q,i,j,k)    (q(i,j,0))
  #define AVERAGE_XZ(q,i,j,k)   AVERAGE_X(q,i,j,0)
  #define AVERAGE_YZ(q,i,j,k)   AVERAGE_Y(q,i,j,0)
  #define AVERAGE_XYZ(q,i,j,k)  AVERAGE_XY(q,i,j,0)
#endif
/**@} */

/*! Set members of the state structure. */
#define STATE_SET(state, vv, uu, bgf, xx, yy, zz)          \
        do {                                               \
           state.v  = vv; state.u  = uu; state.Bgf = bgf;  \
           state.x1 = xx; state.x2 = yy; state.x3 = zz;    \
        } while (0)

/*! \name Warning & error detection macros.
    \C The macros \C ATOMIC_ERR_UPDATE() and \C KERNEL_ERR_UPDATE()
    provide one-line compact construct to perform atomic or non-atomic
    flagging of errors. */

/**@{ */
//
// WARNING_MESSAGES    QUIT Only   Print cumulative err     Print local err
// -----------------+------------+------------------------+-----------------
//      -1          |    NO      |      NO                |      NO     
//       0          |   YES      |      NO                |      NO 
//       1          |   YES      |     YES                |      NO
//       2          |   YES      |     YES                |     YES
// 

#if WARNING_MESSAGES >= 0
  // Define the main macro to handle multiple expressions
  #define ATOMIC_ERR_UPDATE(err, ...) \
      do { \
        uint32_t _arr[] = {__VA_ARGS__}; \
        for (int _i = 0; _i < sizeof(_arr)/sizeof(_arr[0]); ++_i) { \
          _Pragma("acc atomic update") \
          err |= _arr[_i];             \
        }  \
      } while (0)
#else

  // Empty macro
  #define ATOMIC_ERR_UPDATE(flag, ...)   {}

#endif

#if WARNING_MESSAGES >= 1
  // Define the main macro to handle multiple expressions
  #define KERNEL_ERR_UPDATE(flag, ...) \
      do { \
        uint32_t _arr[] = {__VA_ARGS__}; \
        for (int _i = 0; _i < sizeof(_arr)/sizeof(_arr[0]); ++_i) flag |= _arr[_i];  \
      } while (0)

#else

  // Empty macro
  #define KERNEL_ERR_UPDATE(flag, ...)  {}

#endif

/*! Return 1 if an err "err_code" is turned on inside err  */
#define ERR_CHECK(err, err_code) (((err) & (err_code)) ? 1 : 0)
/**@} */

#ifdef _OPENACC
 #define MALLOC_DOUBLE(size)     (double *)acc_malloc(size * sizeof(double))
#else
 #define MALLOC_DOUBLE(size)     (double *)malloc(size * sizeof(double))
#endif

#ifdef _OPENACC
 #define FREE(ptr)   acc_free(ptr)
#else
 #define FREE(ptr)   free(ptr)
#endif

#define PROC_EXIST(bound, proc) ((bound == 0 || bound == PERIODIC || bound == SHEARING) && proc > 1)
