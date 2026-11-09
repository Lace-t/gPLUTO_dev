////////////////////////////////////////////////////////////////////////
/*!
  \file
  \brief Collection of functions for debugging purposes.

  \author A. Mignone (andrea.mignone@unito.it)\n
          S. Truzzi  
  \date   March 28, 2026
*/
////////////////////////////////////////////////////////////////////////
#include "gpluto.hpp"

////////////////////////////////////////////////////////////////////////
/*
 * Check whether the array \c U, \c V or \c Vs contain
 * Not-a-Number (NaN), negative density or pressures.
 * NOTE: this function works on the CPU and should be never called for 
 *       performance issues.
 *
  \param [in]  U       an array of conserved variables of the type
 *                      U(i,j,k,nv); Ignored when U == NULL
  \param [in]  V       an array of primitive variables of the type
 *                      V(i,j,k,nv); Ignored when V == NULL
  \param [in]  Vs      an array of staggered fields of the type
 *                      Vs[nv](i,j,k); Ignored when Vs == NULL
  \param [in]  positivity_check   an integer (0/1) that enables 
 *                                 check of density / energy / pressure
  \param [in]  box     the box to be checked
  \param [in]  str     a reference string to be printed.
 *
 * \return Return 0 if no problem has been found. Otherwise return
 *         1: for nan, 2: for negative density and 3: for negative
 *         pressure or energy.
*/
////////////////////////////////////////////////////////////////////////
int CheckData (Ary4D &U, Ary4D &V, Ary3D (&Vs)[6], int positivity_check,
               RBox *box, const char *str, Grid *grid)
{
  int nv, i,j,k;
  uint32_t err = 0;
  static Ary3Duint32_t flag;
  
  mynvtxstart_("DBG_TOOLS", BLUE);
  // 0. Static storage for local 3D flag

  if (flag.isAllocated() == 0){
    flag.create(grid->ntot[IDIR],grid->ntot[JDIR],grid->ntot[KDIR]);
    #pragma acc enter data copyin(flag)
  }

  // 1. Check nan in conserved variables

  if (U.isAllocated() != 0) {
    #pragma acc parallel loop collapse(3) present(U, flag, box)
    BOX_LOOP(*box, i,j,k){
      flag(i,j,k) = 0;

      NVAR_LOOP(nv) {
        if (isnan(U(i,j,k,nv))) {
          #pragma acc atomic update
          err         |= ERR_NAN;
          flag(i,j,k) |= ERR_NAN;
        }
      }

      if (positivity_check && U(i,j,k,RHO) < 0.0) {
        #pragma acc atomic update
        err         |= ERR_NEGATIVE_RHO;
        flag(i,j,k) |= ERR_NEGATIVE_RHO;
      }

      #ifdef ENG
      if (positivity_check && U(i,j,k,ENG) < 0.0) {
        #pragma acc atomic update
        err         |= ERR_NEGATIVE_ENG;
        flag(i,j,k) |= ERR_NEGATIVE_ENG;
      }
      #endif
    }
  }  

  // 2. Check nan in primitive variables

  if (V.isAllocated() != 0) {
    #pragma acc parallel loop collapse(3) present(V, box)
    BOX_LOOP(*box, i,j,k){
    
      NVAR_LOOP(nv) {
        if (isnan(V(i,j,k,nv))) {
          #pragma acc atomic update
          err         |= ERR_NAN;
          flag(i,j,k) |= ERR_NAN;
        }
      }

      if (positivity_check && V(i,j,k,RHO) < 0.0) {
        #pragma acc atomic update
        err         |= ERR_NEGATIVE_RHO;
        flag(i,j,k) |= ERR_NEGATIVE_RHO;
      }

      #ifdef ENG
      if (positivity_check && V(i,j,k,ENG) < 0.0) {
        #pragma acc atomic update
        err         |= ERR_NEGATIVE_ENG;
        flag(i,j,k) |= ERR_NEGATIVE_ENG;
      }
      #endif
    }
  }

  // 3. Check nan in staggered fields

#ifdef STAGGERED_MHD
  if (Vs[0].isAllocated() != 0) {

    // 3a. Check nan in BX1s, EX1s 

    #pragma acc parallel loop collapse(3) present(Vs[:6], box)
    BOXI_LOOP(*box, i,j,k){
      if (isnan(Vs[BX1s](i,j,k))) {
        #pragma acc atomic update
        err         |= ERR_NAN;
        flag(i,j,k) |= ERR_NAN;
      }
      #if (PHYSICS == ResRMHD) && (DIVE_CONTROL == CONSTRAINED_TRANSPORT)
      if (isnan(Vs[EX1s](i,j,k))) {
        #pragma acc atomic update
        err         |= ERR_NAN;
        flag(i,j,k) |= ERR_NAN;
      }
      #endif
    }

    // 3b. Check nan in BX2s, EX2s 

    #if INCLUDE_JDIR
    #pragma acc parallel loop collapse(3) present(Vs[:6], box)
    BOXJ_LOOP(*box, i,j,k){
      if (isnan(Vs[BX2s](i,j,k))) {
        #pragma acc atomic update
        err         |= ERR_NAN;
        flag(i,j,k) |= ERR_NAN;
      }
      #if (PHYSICS == ResRMHD) && (DIVE_CONTROL == CONSTRAINED_TRANSPORT)
      if (isnan(Vs[EX2s](i,j,k))){
        #pragma acc atomic update
        err         |= ERR_NAN;
        flag(i,j,k) |= ERR_NAN;
      }
      #endif

    }
    #endif

    // 3c. Check nan in BX3s, EX3s 
    
    #if INCLUDE_KDIR
    #pragma acc parallel loop collapse(3) present(Vs[:6], box)
    BOXK_LOOP(*box, i,j,k){
      if (isnan(Vs[BX3s](i,j,k))) {
        #pragma acc atomic update
        err         |= ERR_NAN;
        flag(i,j,k) |= ERR_NAN;
      }
      #if (PHYSICS == ResRMHD) && (DIVE_CONTROL == CONSTRAINED_TRANSPORT)
      if (isnan(Vs[EX3s](i,j,k))) {
        #pragma acc atomic update
        err         |= ERR_NAN;
        flag(i,j,k) |= ERR_NAN;
      }
      #endif

    }
    #endif   // INCLUDE_KDIR 
  }       // Vs != NULL 
#endif // STAGGERED_MHD 

  // 4. Error Check
  // #if WARNING_MESSAGES > 0 
  //  Where ("CheckData()", flag, err, box, grid);
  // #endif

  // 5. Need to Quit ? 
  // if (err & ERR_NAN) {
  //   QuitPLUTO("from CheckData(): NaN found in U",err);
  // }

  mynvtxstop_();


  return 0;
}

////////////////////////////////////////////////////////////////////////
/*!
  Show which bits of the error flag are set.
  The variable filter can be used to select only certain bits.
  
  \param [in] err        an error variable
*/
////////////////////////////////////////////////////////////////////////
string ErrorShow(uint32_t err)
{
  int nbits = 19;
  int bit_position;
  // int isBitSet;
  int show_err = 0;
  string fstring = "";

  // 1. Check bits to see which of them are on. 
  //    Enable show_flag only if they correspond to those
  //    specified by filter.

  for (bit_position = 0; bit_position < nbits; bit_position++){
    // isBitSet  = err & (1 << bit_position);
    show_err = 1;
  }

  if (show_err){
    if (err & ERR_NEGATIVE_RHO)   fstring += "/ERR_NEGATIVE_RHO"; 
    if (err & ERR_NEGATIVE_ENG)   fstring += "/ERR_NEGATIVE_ENG"; 
    if (err & ERR_NEGATIVE_ENR)   fstring += "/ERR_NEGATIVE_ENR";
    if (err & ERR_NEGATIVE_PRS)   fstring += "/ERR_NEGATIVE_PRS"; 
    if (err & ERR_NEGATIVE_ENTR)  fstring += "/ERR_NEGATIVE_ENTR"; 
    if (err & ERR_CONS2PRIM_FAIL) fstring += "/ERR_CONS2PRIM_FAIL"; 
    if (err & ERR_PRESSURE_FIX)   fstring += "/ERR_PRESSURE_FIX"; 
    if (err & ERR_ENG_SOLVE_INACCURATE)  fstring += "/ERR_ENG_SOLVE_INACCURATE"; 
    if (err & ERR_ENG_SOLVE_NO_SOLUTION) fstring += "/ERR_ENG_SOLVE_NO_SOLUTION"; 
    if (err & ERR_SUPERLUMINAL_SPEED)    fstring += "/ERR_SUPERLUMINAL_SPEED"; 
    if (err & ERR_TOO_MANY_ITERATIONS)   fstring += "/ERR_TOO_MANY_ITERATIONS"; 
    if (err & ERR_RIEMANN_SOLVE_FAIL)    fstring += "/ERR_RIEMANN_SOLVE_FAIL"; 
    if (err & ERR_NAN )                  fstring += "/ERR_NAN"; 
    if (err & ERR_POINT_VALUE_FAIL)      fstring += "/ERR_POINT_VALUE_FAIL"; 
    if (err & ERR_IMEX_GUESS_FAIL)       fstring += "/ERR_IMEX_GUESS_FAIL"; 
    if (err & ERR_QUIT)                  fstring += "/ERR_QUIT"; 
    if (err & ERR_SUPERLUMINAL_RADIATION_FLUX) fstring += "/ERR_SUPERLUMINAL_RADIATION_FLUX";
    if (err & ERR_GAUSSIAN_SOLVE)        fstring += "/ERR_GAUSSIAN_SOLVE";
  }
  return fstring;
}

////////////////////////////////////////////////////////////////////////
/*!
  Show which bits of contain active flags.
  The variable filter can be used to select only certain bits.
  
  \param [in] flag        a flag variable
*/
////////////////////////////////////////////////////////////////////////
string FlagShow(uint32_t flag)
{
  int nbits = 16;
  int bit_position;
  // int isBitSet;
  int show_flag = 0;
  string fstring = "";

  // 1. Check bits to see which of them are on. 
  //    Enable show_flag only if they correspond to those
  //    specified by filter.

  for (bit_position = 0; bit_position < nbits; bit_position++){
    // isBitSet  = flag & (1 << bit_position);
    show_flag = 1;
  }

  if (show_flag){
    if (flag & FLAG_FLAT)              fstring += "/FLAG_FLAT"; 
    if (flag & FLAG_MINMOD)            fstring += "/FLAG_MINMOD"; 
    if (flag & FLAG_HLL)               fstring += "/FLAG_HLL"; 
    if (flag & FLAG_ENTROPY)           fstring += "/FLAG_ENTROPY"; 
    if (flag & FLAG_ENERGY)            fstring += "/FLAG_ENERGY"; 
    if (flag & FLAG_INTERNAL_BOUNDARY) fstring += "/FLAG_INTERNAL_BOUDARY"; 
    if (flag & FLAG_HO_LAP_LIMITER)    fstring += "/FLAG_HO_LAP_LIMITER"; 
    if (flag & FLAG_FALLBACK)          fstring += "/FLAG_FALLBACK"; 
    if (flag & FLAG_RAD_MAX_ITER)      fstring += "/FLAG_RAD_MAX_ITER";
    if (flag & FLAG_GAUSSIAN_SOLVE)    fstring += "/FLAG_GAUSSIAN_SOLVE";
    if (flag & FLAG_SHOCK)             fstring += "/FLAG_SHOCK"; 
  }
  return fstring;
}

////////////////////////////////////////////////////////////////////////
/*
 * Print explicit values of U, V or Vs in the box.
 *
  \param [in]  U       an array of conserved variables of the type
 *                      U[k][j][i][nv]; Ignored when U == NULL
  \param [in]  V       an array of primitive variables of the type
 *                      V[nv][k][j][i]; Ignored when V == NULL
  \param [in]  Vs      an array of staggered fields of the type
 *                      Vs[nv][k][j][i]; Ignored when Vs == NULL
  \param [in]  box     the box to be checked
  \param [in]  str     a reference string to be printed.
 * 
 * Example usage:
  {  ////////////////////////////////////////////////////////////
    RBox tbox;
    Ary4D Qnull;
    Ary3D Vsnull[6];
    RBoxDefine (grid->lbeg, grid->lend, &tbox);
    ShowData (Qnull, Qnull, Bs0, &tbox, "RKStep(stage=1)");
    exit(0);
  }  ///////////////////////////////////////////////////////////////
*/
////////////////////////////////////////////////////////////////////////
void ShowData (Ary4D &U, Ary4D &V, Ary3D (&Vs)[6], RBox *box, const char *str)
{
  int i,j,k,nv;
  double q[256];

  // 1. Print value of conserved variables

  if (U.isAllocated()){
    BOX_LOOP(*box,i,j,k){
      NVAR_LOOP(nv) q[nv] = U(i,j,k, nv);
      printLog ("ShowData() [%s]: U(%d, %d, %d) = \n", str,i,j,k);
      ShowState(q, 2);
    }
  }

  // 2. Print value of primitive variables

  if (V.isAllocated()){
    BOX_LOOP(*box,i,j,k){
      NVAR_LOOP(nv) q[nv] = V(i,j,k, nv);
      printLog ("ShowData() [%s]: V(%d, %d, %d) = \n", str,i,j,k);
      ShowState(q, 1);
    }
  }

  // 3. Print value of staggered fields

#ifdef STAGGERED_MHD
  if (Vs[BX1s].isAllocated()){

    BOXI_LOOP(*box,i,j,k){
      CT_ASSIGN( printLog ("ShowData() [%s]: Vs[BX1s](%d, %d, %d) = %12.6e\n",
                            str, i,j,k,Vs[BX1s](i,j,k));   ,
                 printLog ("                 Vs[EX1s](%d, %d, %d) = %12.6e\n",
                            i,j,k,Vs[EX1s](i,j,k));) 
    } 

    #if INCLUDE_JDIR
    BOXJ_LOOP(*box,i,j,k){
      CT_ASSIGN( printLog ("ShowData() [%s]: Vs[BX2s](%d, %d, %d) = %12.6e\n",
                            str, i,j,k,Vs[BX2s](i,j,k));   ,
                 printLog ("                 Vs[EX2s](%d, %d, %d) = %12.6e\n",
                             i,j,k,Vs[EX2s](i,j,k));)
    } 
    #endif

    #if INCLUDE_KDIR
    BOXK_LOOP(*box,i,j,k){
      CT_ASSIGN( printLog ("ShowData() [%s]: Vs[BX3s](%d, %d, %d) = %12.6e\n",
                            str, i,j,k,Vs[BX3s](i,j,k));   ,
                 printLog ("                 Vs[EX3s](%d, %d, %d) = %12.6e\n",
                            i,j,k,Vs[EX3s](i,j,k));)
    } 
    #endif

  }
#endif // STAGGERED_MHD 
}

////////////////////////////////////////////////////////////////////////
/*!
  Make a nice printLoging of a 2D square  matrix \c A[0..n-1][0..n-1]
  Entries with values below eps will display "0.0"
*/
////////////////////////////////////////////////////////////////////////
void ShowMatrix(double **A, int n, double eps)
{
  int k1,k2;

  printLog ("----------------------------------------------------------------\n");
  for (k1 = 0; k1 < n; k1++){
    for (k2 = 0; k2 < n; k2++){
      printLog ("%12.3e   ", fabs(A[k1][k2]) > eps ? A[k1][k2]:0.0);
    }
    printLog ("\n");
  }
  printLog ("----------------------------------------------------------------\n");
}

////////////////////////////////////////////////////////////////////////
/*! 
  Print the components of the array v.
  Here mode = 1 (-1) if the array contains primitive variables;
  Here mode = 2 (-2) if the array contains conservative variables;
  Here mode = 3      if the array contains fluxes;
  
  A positive value of mode will produce code-ready vertical layout 
  while a negative value will produce horizontal layout.
*/
////////////////////////////////////////////////////////////////////////
void ShowState (double *vc, int mode)
{
  if (mode == 1){  // Primitive variables, vertical layout 
    printLog ("  vc[RHO] = %18.12e;\n", vc[RHO]);
    printLog ("  vc[VX1] = %18.12e;\n", vc[VX1]);
    printLog ("  vc[VX2] = %18.12e;\n", vc[VX2]);
    printLog ("  vc[VX3] = %18.12e;\n", vc[VX3]);
    #if HAVE_ENERGY
    printLog ("  vc[PRS] = %18.12e;\n", vc[PRS]);
    #endif
    #if ENTROPY_SWITCH != NO
    printLog ("  vc[ENTR] = %18.12e;\n", vc[ENTR]);
    #endif
  
    #if PHYSICS == MHD || PHYSICS == RMHD || PHYSICS == ResRMHD || PHYSICS == GRMHD
    printLog ("  vc[BX1] = %18.12e;\n", vc[BX1]);
    printLog ("  vc[BX2] = %18.12e;\n", vc[BX2]);
    printLog ("  vc[BX3] = %18.12e;\n", vc[BX3]);
    #endif
    #if PHYSICS == ResRMHD
    printLog ("  vc[EX1] = %18.12e;\n", vc[EX1]);
    printLog ("  vc[EX2] = %18.12e;\n", vc[EX2]);
    printLog ("  vc[EX3] = %18.12e;\n", vc[EX3]);
    #ifdef CRG
    printLog ("  vc[CRG] = %18.12e;\n", vc[CRG]);
    #endif
    #endif
   
    #ifdef PHI_GLM
    printLog ("  vc[PHI_GLM] = %18.12e;\n", vc[PHI_GLM]);
    #endif
    #ifdef PSI_GLM
    printLog ("  vc[PSI_GLM] = %18.12e;\n", vc[PSI_GLM]);
    #endif
    printLog ("\n");

  } else if (mode == -1) { // Primitive variables, horizontal layout 
    printLog ("rho = %+8.3e; vx1 = %+8.3e; vx2 = %+8.3e; vx3 = %+8.3e",
               vc[RHO], vc[VX1], vc[VX2], vc[VX3]);
    #if HAVE_ENERGY
    printLog ("; prs = %+8.3e", vc[PRS]);
    #endif
    #if ENTROPY_SWITCH != NO
    printLog ("; entr = %18.12e", vc[ENTR]);
    #endif
    #if PHYSICS == MHD || PHYSICS == RMHD || PHYSICS == ResRMHD || PHYSICS == GRMHD
    printLog ("; Bx1 = %+8.3e; Bx2 = %+8.3e; Bx3 = %+8.3e",
                 vc[BX1], vc[BX2], vc[BX3]);  
    #endif             
    printLog ("\n");

  }else if (mode == 2){
    printLog ("  uc[RHO] = %18.12e;\n", vc[RHO]);
    printLog ("  uc[MX1] = %18.12e;\n", vc[MX1]); 
    printLog ("  uc[MX2] = %18.12e;\n", vc[MX2]); 
    printLog ("  uc[MX3] = %18.12e;\n", vc[MX3]);
    #if HAVE_ENERGY
    printLog ("  uc[ENG] = %18.12e;\n", vc[ENG]);
    #endif
    #if ENTROPY_SWITCH != NO
    printLog ("  uc[ENTR] = %18.12e;\n", vc[ENTR]);
    #endif

    #if PHYSICS == MHD || PHYSICS == RMHD || PHYSICS == ResRMHD || PHYSICS == GRMHD
    printLog ("  uc[BX1] = %18.12e;\n", vc[BX1]);
    printLog ("  uc[BX2] = %18.12e;\n", vc[BX2]);
    printLog ("  uc[BX3] = %18.12e;\n", vc[BX3]);
    #endif
    #if PHYSICS == ResRMHD
    printLog ("  uc[EX1] = %18.12e;\n", vc[EX1]);
    printLog ("  uc[EX2] = %18.12e;\n", vc[EX2]);
    printLog ("  uc[EX3] = %18.12e;\n", vc[EX3]);
    #ifdef CRG
    printLog ("  uc[CRG] = %18.12e;\n", vc[CRG]);
    #endif
    #endif
   
    #ifdef PHI_GLM
    printLog ("  uc[PHI_GLM] = %18.12e;\n", vc[PHI_GLM]);
    #endif
    #ifdef PSI_GLM
    printLog ("  uc[PSI_GLM] = %18.12e;\n", vc[PSI_GLM]);
    #endif
    printLog ("\n");
    
  }else if (mode == -2){

    printLog ("  rho = %+8.3e; mx1 = %+8.3e; mx2 = %+8.3e; mx3 = %+8.3e",
               vc[RHO], vc[MX1], vc[MX2], vc[MX3]);
    #if HAVE_ENERGY
    printLog ("  eng = %+8.3e", vc[ENG]);
    #endif
    #if ENTROPY_SWITCH != NO
    printLog ("  entr = %+8.3e", vc[ENTR]);
    #endif
    #if PHYSICS == MHD || PHYSICS == RMHD || PHYSICS == ResRMHD || PHYSICS == GRMHD
    printLog ("  Bx1 = %+8.3e; Bx2 = %+8.3e; Bx3 = %+8.3e;",
                 vc[BX1], vc[BX1], vc[BX3]);  
    #endif
    printLog("\n");

  }else if (mode == 3){   // -- Flux mode -- 

    printLog ("  F[RHO] = %18.12e;\n", vc[RHO]);
    printLog ("  F[MX1] = %18.12e;\n", vc[MX1]); 
    printLog ("  F[MX2] = %18.12e;\n", vc[MX2]); 
    printLog ("  F[MX3] = %18.12e;\n", vc[MX3]);
    #if HAVE_ENERGY
    printLog ("  F[ENG] = %18.12e;\n", vc[ENG]);
    #endif
    #if ENTROPY_SWITCH != NO
    printLog ("  F[ENTR] = %18.12e;\n", vc[ENTR]);
    #endif

    #if PHYSICS == MHD || PHYSICS == RMHD || PHYSICS == ResRMHD || PHYSICS == GRMHD
    printLog ("  F[BX1] = %18.12e;\n", vc[BX1]);
    printLog ("  F[BX2] = %18.12e;\n", vc[BX2]);
    printLog ("  F[BX3] = %18.12e;\n", vc[BX3]);
    #endif
    #if PHYSICS == ResRMHD
    printLog ("  F[EX1] = %18.12e;\n", vc[EX1]);
    printLog ("  F[EX2] = %18.12e;\n", vc[EX2]);
    printLog ("  F[EX3] = %18.12e;\n", vc[EX3]);
    #ifdef CRG
    printLog ("  F[CRG] = %18.12e;\n", vc[CRG]);
    #endif
    #endif
   
    #ifdef PHI_GLM
    printLog ("  F[PHI_GLM] = %18.12e;\n", vc[PHI_GLM]);
    #endif
    #ifdef PSI_GLM
    printLog ("  F[PSI_GLM] = %18.12e;\n", vc[PSI_GLM]);
    #endif
    printLog ("\n");
    
  }else{
    printLog ("! Invalid mode specified in ShowState()\n");
  }
}

////////////////////////////////////////////////////////////////////////
/*!
  Print the first n components of the vector v[]
*/
////////////////////////////////////////////////////////////////////////
void ShowVector (double *v, int n)
{
  int k;

  for (k = 0; k < n; k++)  printLog ("%+12.6e  ", v[k]);
  printLog ("\n");
}

////////////////////////////////////////////////////////////////////////
/*!
  When WARNING_MESSAGES is set to 1, only the existence of the error 
  and its type are recorded (but not zone index).
  When WARNING_MESSAGES is set to 2, more diagnostic is printed 
  indicating also the zone.
  
  \param [in] func_name   the caller function name
  \param [in] err         the flag() array, used when WARNING_MESSAGES == 2
  \param [in] err_tot     the cumulative error (used when WARNING_MESSAGES == 1)
  \param [in] box         the box on which this function operates
  \param [in] grid        pointer to Grid structure
*/
////////////////////////////////////////////////////////////////////////
void Where (const string func_name, Data *data, uint32_t err_tot, 
            RBox *box, Grid *grid)
{
#if WARNING_MESSAGES == 0 
#else
  int i,j,k;

  if (err_tot == 0) return;

  double time = data->Dts->time;

#if WARNING_MESSAGES == 1 
#endif

#if WARNING_MESSAGES == 1 || WARNING_MESSAGES == 2

  // 1. Count how many zones have been flagged with errors
  //    [Since there's no print involved, this loop can be 
  //     carried out on the GPU]
  
  uint64_t tmp_count = 0;

  #pragma acc parallel loop collapse(3) present(data, box) reduction(+:tmp_count)
  BOX_LOOP(*box, i,j,k){
    if (data->err(i,j,k) != 0) tmp_count++; 
  } 

  data->err_count += tmp_count;

  #if WARNING_MESSAGES == 1
  if (data->err_count > 0){
    printLog ("! %llu cells reported errors ",data->err_count);
    printLog ("[step: %d; t = %10.4e; g_intStage = %d]\n", 
                g_stepNumber, time , g_intStage);
  }
  #endif

#endif // WARNING_MESSAGES == 1 || WARNING_MESSAGES == 2

#if WARNING_MESSAGES == 2

  // 2. Print a cumulative error (not zone-by-zone) using 
  //    just the err_tot int.

  if (err_tot){ 
    printLog ("! %s: step: %d; t = %10.4e; g_intStage = %d", 
            func_name.c_str(), g_stepNumber, time, g_intStage);
    #if FALLBACK_MODE != NO
    printLog ("; integration mode = %s", 
              data->integration_mode == DEFAULT ? "DEFAULT":"FALLBACK");
    #endif
    printLog("\n");

  // Align left 

    for (int l = 0; l < func_name.length(); l++) printLog (" ");  
    printLog ("    %llu cells reported errors: ", data->err_count);

    printLog ("%s\n", ErrorShow(err_tot).c_str());
  }
  printLog ("\n");  

#elif WARNING_MESSAGES >= 3

  // 3. Perform a more accurate zone by zone search and 
  //    print local errors.

  BOX_LOOP(*box, i,j,k){

    // 3a. Print errors only if flag is > 0

    if (data->err(i,j,k) > 0) {

      double x1 = grid->xC[IDIR][i];
      double x2 = grid->xC[JDIR][j];
      double x3 = grid->xC[KDIR][k];

      data->err_count++;    

      // 3b. Print error string names when enabled. 

      uint32_t local_err  = data->err(i,j,k);  // Copy errors 
      uint32_t local_flag = data->flag(i,j,k); // Copy flags  

      if (local_err) {   // Error found ? 
        printLog ("! %s: zone: ", func_name.c_str());
        DIM_SELECT(
          printLog ("i = %d, x1 = %f\n", i, x1);  ,
          printLog ("(i,j) = (%d, %d), (x1,x2) = (%f, %f)\n", i, j, x1, x2);  ,
          printLog ("(i,j,k) = (%d, %d, %d), (x1,x2,x3) = (%f, %f, %f)\n", 
                      i, j, k, x1, x2, x3);
        )

        for (int l = 0; l < func_name.length(); l++) printLog (" "); // Align left
        printLog ("    step: %d; t = %10.4e; g_intStage = %d", 
                  g_stepNumber, time, g_intStage);
        #if FALLBACK_MODE != NO
        printLog ("; integration mode = %s", 
                    data->integration_mode == DEFAULT ? "DEFAULT":"FALLBACK");
        #endif
        printLog("\n");

        for (int l = 0; l < func_name.length(); l++) printLog (" "); // Align left
        printLog ("    error(s) found: %s\n", ErrorShow(local_err).c_str());
        for (int l = 0; l < func_name.length(); l++) printLog (" "); // Align left  
        printLog ("    active flag(s): %s\n", FlagShow (local_flag).c_str());
        for (int l = 0; l < func_name.length(); l++) printLog (" "); // Align left  
        printLog ("    err_count = %llu\n", data->err_count);
      }
    }  // End if flag() 
  } // End BOX_LOOP()  
#endif  // WARNING_MESSAGES >= 3 
#endif  // WARNING_MESSAGES == 0 
}
