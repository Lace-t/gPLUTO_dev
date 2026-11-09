////////////////////////////////////////////////////////////////////////
/*!
  \file
  \brief Solve the stiff-implicit part of the RK-Imex method.

  Main driver for including the stiff source terms of the resistive
  RMHD equation within the IMEX formalism.

  - "General relativistic magnetohydrodynamic dynamo in thick accretion discs: 
     fully non-linear simulations"\n
     Tomei \& al. MNRAS (2020) 491 2346-2359

  - "A constrained transport method for the solution of the
     resistive relativistic MHD equations"\n
     Mignone \& al. MNRAS (2019) 486, 4252

  \authors A. Mignone  (andrea.mignone@unito.it)\n
           M. Rossazza (marco.rossazza@unito.it)\n
           G. Mattia\n          
           
  \date    Jan 28, 2026
*/
////////////////////////////////////////////////////////////////////////
#include "gpluto.hpp"

int IdealRMHD_Guess(double *, double *, RunConfig *);

#define  ENERGY_SOLVE   1
#define  ENTROPY_SOLVE  2

#define MAX_ITER  20
#define IMEX_DEBUG  0  // == NO (no debug             
                       // == 1  (only at crash + verify solution
                       // == 2  (every step) 

////////////////////////////////////////////////////////////////////////
void IMEXImplicitUpdate (Data *data, double dt1, Ary4D &S, Ary4D &Jf, 
                         Grid *grid, RBox *box)
/*!
  Perform an implicit update of the conservative variables by adding
  the stiff source term to the electric field.
  Only the electric field is actually updated al.hppough this requires
  primitive variables to change.
  
  \f[
  U^{(1)} = U^{(0)} + \frac{\Gamma}{\eta}\Delta tS(U^{(1)}
  \qquad{\rm where}\qquad \Gamma = \left(1 - \frac{1}{\sqrt{2}}\right)
  \f]
  
  \param [in,out]  data   PLUTO data array.
  On input, Vc is an array of primitive quantities
  at the previous time level.
  On output Vc is converted from Uc at the
  next time level.
  \param [in]      dt1    The time step to be taken
  \param [out]     S      A 4D array containing the three component of the
  source term vector
  \param [out]     Jf     ...
  \param [in]      grid   Pointer to grid structure.
  \param [out]     box    The box on which the imex step is taken
*/
////////////////////////////////////////////////////////////////////////
{
  int    i, j, k, nv;
  #if IMEX_DEBUG != NO
  int debug_condition = 1;
  #endif
  uint32_t err_kernel  = 0;   // Cumulative kernel error flag 
  RunConfig &run_config = *(data->run_config);

  #if EOS == IDEAL
  double gmm  = run_config.gamma;
  double gmm1 = gmm/(gmm - 1.0);
  #endif   
  double tol_f = 1.e-14, tol_g = 1.e-14; // Tolerance used in the iteration cycle
  double eta;

  #ifdef HIGH_ORDER
  int iarr = 1;
  const Ary4D& Vc = data->Vpc; 
  const Ary4D& Uc = data->Upc; 
  #else
  int iarr = 0;
  const Ary4D& Vc = data->Vc; 
  const Ary4D& Uc = data->Uc; 
  #endif

  // 0. Allocate static arrays

  // 1. Set boundary conditions when CT on electric field
  //    is used. In this case, ghost zone values are required
  //    to be advanced as well.

  #if DIVE_CONTROL == CONSTRAINED_TRANSPORT
  Boundary (data, 0, grid);
  PrimToCons3D(data, iarr, grid, box);
  #endif

  #ifdef HIGH_ORDER
  PrimToCons3D(data, 1, grid, box);
  #endif

  mynvtxstart_("Imex_Source",YELLOW);

  // ------------------------------------------------------
  // 2, 3. Loop over computational zones (2) 
  //       Apply Newton-Broyden method to each cell (3) 
  // ------------------------------------------------------

  #pragma acc parallel loop collapse(3) present(data, data->run_config, box, grid)\
                                        copy(err_kernel)
  #pragma omp target teams distribute parallel for collapse (3) map (to: data->run_config, S)  //map (alloc:err_kernel)
  BOX_LOOP(*box, i,j,k) {
    int    done;
    int    method = ENERGY_SOLVE;  // Default method for now. 
                                   // No major improvements is obeserved when
                                   // we switch to ENTROPY_SOLVE (parts below
                                   // are commented)

    double vc[NVAR], uc[NVAR], vc0[NVAR], uc0[NVAR];
    double func[3], J[3][3];
    double u[3], uxB[3], du[3];
    double err_f, err_g; 
    double glor, glor_1;
    double Emag2, Bmag2, uE;
    double rho, rhoh, h;
    double Pem, D, Dh, theta;
    double scrh1, scrh2;
    double A1, A5, A0, etap;
    double A0d, A1d;
    double dE[3][3];
    double dE_2[3], dh[3], x[3]; 

    State state;

    data->err(i,j,k) = 0;

    NVAR_LOOP(nv) uc[nv] = uc0[nv] = Uc(i,j,k,nv);

    #if GEOMETRY != CARTESIAN
    Metric metric;
    MetricSet(metric, grid->rC[i], grid->sC[j]);
    ConsToOrtho(uc, metric);
    #endif

    // 2b. Check for Nan in input state

    for (nv   = 0; nv < NVAR; nv++){
      vc[nv] = vc0[nv] = Vc(i,j,k,nv);

      if (isnan(vc[nv])){
        ATOMIC_ERR_UPDATE(err_kernel, ERR_NAN);
        KERNEL_ERR_UPDATE(data->err(i,j,k), ERR_NAN, ERR_QUIT);
      }
    }

    STATE_SET (state, vc, uc, NULL, grid->xC[IDIR][i], 
                                    grid->xC[JDIR][j], 
                                    grid->xC[KDIR][k]);
    eta = Resistive_eta(state, *(data->run_config)); 

    #if IMEX_DEBUG == 2
    debug_condition = (g_stepNumber == 84) && (i == 203) && (j == 9);
    if (debug_condition) printf ("   Beg loop: i,j,k = %d, %d, %d\n",i,j,k);
    #endif
  
    // 2c. Provide a suitable guess to the 4-velocity
    //     to start the iteration cycle.
     
    if (dt1/eta > 10.0*(1.0 - 1.0/sqrt(2.0))) {
      uint32_t status;
      status = IdealRMHD_Guess(uc, vc, data->run_config);
      if (status != 0) {
        ATOMIC_ERR_UPDATE(err_kernel, ERR_IMEX_GUESS_FAIL);
        KERNEL_ERR_UPDATE(data->err(i,j,k), ERR_IMEX_GUESS_FAIL);
        NVAR_LOOP(nv) vc[nv] = vc0[nv];  // Guess is input state
      }
    }

    glor = 1.0/sqrt(1.0 - vc[VX1]*vc[VX1] - vc[VX2]*vc[VX2] - vc[VX3]*vc[VX3]);

    u[0] = vc[VX1]*glor;
    u[1] = vc[VX2]*glor;
    u[2] = vc[VX3]*glor;
    D    = uc[RHO];

    // ----------------------------------------------------
    // 3. Start iterating
    // ----------------------------------------------------

    err_g  = err_f = 1.0;
    done   = 0;

    #if IMEX_DEBUG == 2
    if (debug_condition) {
      printLog ("     Beg NB iteration cycle\n"); 
      printLog ("     Input state (vc):\n"); ShowState(vc,1);
      printLog ("     Input state (uc):\n"); ShowState(uc,1);
    }
    #endif

   for(int kcycle = 0; kcycle < MAX_ITER; kcycle++){
        
      #if IMEX_DEBUG == 2
      if (debug_condition) printf ("      ----> kcycle = %d\n",kcycle);
      #endif

      // 3a. Using the current four-velocity estimate compute primitive 
      //     variables and electric field.
      //     Important: while vc[i][EX#] changes through the iteration cycle, 
      //     uc[i+nx1_tot*EX#] must be kept constant to the input value.

      glor   = sqrt(1.0 + u[0]*u[0] + u[1]*u[1] + u[2]*u[2]);
      glor_1 = 1.0/glor;
      rho    = D*glor_1;
      
      uxB[0] = u[1]*uc[BX3] - u[2]*uc[BX2];
      uxB[1] = u[2]*uc[BX1] - u[0]*uc[BX3];
      uxB[2] = u[0]*uc[BX2] - u[1]*uc[BX1];
      uE     = vc0[EX1]*u[0] + vc0[EX2]*u[1] + vc0[EX3]*u[2];  

      etap = dt1/(eta + dt1*glor);    
      
      A1   = eta/(dt1 + eta*glor);    // Eq. [A7] of Tomei et al.
      A5   = -1.0;                    // Eq. [A11] of Tomei et al.
      A0 = etap;                    // Inverse of Eq. [A6] of Tomei et al.

      // 3b. Now we find the new E as a function of the old E (Eq. [52] of 
      //     Mignone et al 2019).
      //     Perhaps it would be worth trying inverting the
      //     linear system.
     
      vc[EX1] = eta/(eta + dt1*glor)*uc[EX1] + A0*(A1*uE*u[0] + A5*uxB[0]);
      vc[EX2] = eta/(eta + dt1*glor)*uc[EX2] + A0*(A1*uE*u[1] + A5*uxB[1]);
      vc[EX3] = eta/(eta + dt1*glor)*uc[EX3] + A0*(A1*uE*u[2] + A5*uxB[2]);
  
      Emag2 = vc[EX1]*vc[EX1] + vc[EX2]*vc[EX2] + vc[EX3]*vc[EX3];
      Bmag2 = uc[BX1]*uc[BX1] + uc[BX2]*uc[BX2] + uc[BX3]*uc[BX3];
      Pem   = 0.5*(Emag2 + Bmag2);
  
      if (method == ENERGY_SOLVE){
        double Eg = uc[ENG] - Pem; // Hydro energy 
        #if EOS == IDEAL
        double denp  = gmm1*glor*glor - 1.0;
        #if RMHD_REDUCED_ENERGY == YES
        vc[PRS] = (Eg + D - D*glor)/denp;

        #error "RMHD_REDUCED_ENERGY"  not compatible with the ResRMHD / RHD module

        #else
        vc[PRS] = (Eg - D*glor)/denp;  // Better off using primitive vars 
                                       // to avoid catastrophic cancellation ? 
        #endif
        #endif

        #if EOS == TAUB
        {
          double glor2 = glor*glor;
          double denp  = 1.0 + glor2*(-5.0 + 4.0*glor2);
          #if RMHD_REDUCED_ENERGY == YES
            #error in Imex_Source
          #else
          scrh1   =  4.0*glor2*D*D*denp + 9.0*glor2*glor2*Eg*Eg;
          vc[PRS] =  2.0*(Eg*Eg - glor2*D*D)
                    /(Eg*(5.0*glor2 - 2.0) + sqrt(scrh1));
          #endif
        }
        #endif
      }
      
      // 3c. If a negative pressure is found, reset it and correct energy.

      if (vc[PRS] < 0.0){
        ATOMIC_ERR_UPDATE(err_kernel,       ERR_NEGATIVE_PRS, ERR_QUIT);
        KERNEL_ERR_UPDATE(data->err(i,j,k), ERR_NEGATIVE_PRS, ERR_QUIT);   

        vc[PRS] = run_config.small_prs;

        #if EOS == IDEAL
        rhoh    = rho + gmm1*vc[PRS];
        #elif EOS == TAUB
        rhoh = 2.5*vc[PRS] + sqrt(2.25*vc[PRS]*vc[PRS] + rho*rho);
        #endif
        #if RMHD_REDUCED_ENERGY == YES
        uc[ENG] = rhoh*glor*glor - vc[PRS] - D + Pem;
        #else
        uc[ENG] = rhoh*glor*glor - vc[PRS] + Pem;
        #endif
      }

      if (done) break;

      theta = vc[PRS]/rho;
      #if EOS == IDEAL
      h  = 1.0 + gmm1*theta;
      #elif EOS == TAUB
      h  = 2.5*theta + sqrt(2.25*theta*theta + 1.0);
      #endif
      Dh = D*h;
      
      // 3d. func(u) = m - (w*v + ExB) is the actual function
      //     we need to invert to find u (4-velocity)
      
      func[0] = uc[MX1] - Dh*u[0] - vc[EX2]*uc[BX3] + vc[EX3]*uc[BX2];
      func[1] = uc[MX2] - Dh*u[1] - vc[EX3]*uc[BX1] + vc[EX1]*uc[BX3];
      func[2] = uc[MX3] - Dh*u[2] - vc[EX1]*uc[BX2] + vc[EX2]*uc[BX1];

      // 3e. Compute Jacobian J = dfunc/du.
      //     These are found in the Appendix of Tomei et al (2020) (A12 - A17),
      //     which are more compact.
      
      A0d = 1.0;
      //      A1d = -eta*eta/(dt1 + eta*glor)/(dt1 + eta*glor); //[AM] Isn't this A1^2 ?
      A1d = -A1*A1;

      dE[0][0]  = uE*(A1d*u[0]*u[0]    + glor*A1);        
      dE[1][0]  = uE*(A1d*u[1]*u[0])   - glor*A5*uc[BX3];
      dE[2][0]  = uE*(A1d*u[2]*u[0])   + glor*A5*uc[BX2];
      dE[0][1]  = uE*(A1d*u[0]*u[1])   + glor*A5*uc[BX3];
      dE[1][1]  = uE*(A1d*u[1]*u[1]    + glor*A1);
      dE[2][1]  = uE*(A1d*u[2]*u[1])   - glor*A5*uc[BX1];
      dE[0][2]  = uE*(A1d*u[0]*u[2])   - glor*A5*uc[BX2];
      dE[1][2]  = uE*(A1d*u[1]*u[2])   + glor*A5*uc[BX1];
      dE[2][2]  = uE*(A1d*u[2]*u[2]    + glor*A1);
      
      dE[0][0] += glor*A1*uc[EX1]*u[0] - A0d*vc[EX1]*u[0];
      dE[1][0] += glor*A1*uc[EX1]*u[1] - A0d*vc[EX2]*u[0];
      dE[2][0] += glor*A1*uc[EX1]*u[2] - A0d*vc[EX3]*u[0];
      dE[0][1] += glor*A1*uc[EX2]*u[0] - A0d*vc[EX1]*u[1];
      dE[1][1] += glor*A1*uc[EX2]*u[1] - A0d*vc[EX2]*u[1];
      dE[2][1] += glor*A1*uc[EX2]*u[2] - A0d*vc[EX3]*u[1];
      dE[0][2] += glor*A1*uc[EX3]*u[0] - A0d*vc[EX1]*u[2];
      dE[1][2] += glor*A1*uc[EX3]*u[1] - A0d*vc[EX2]*u[2];
      dE[2][2] += glor*A1*uc[EX3]*u[2] - A0d*vc[EX3]*u[2];
      
      scrh1 = A0/glor;
      dE[0][0] *= scrh1;   // d(E_x)/d(u_x) 
      dE[1][0] *= scrh1;   // d(E_y)/d(u_x)
      dE[2][0] *= scrh1;
      dE[0][1] *= scrh1;
      dE[1][1] *= scrh1;
      dE[2][1] *= scrh1;
      dE[0][2] *= scrh1;
      dE[1][2] *= scrh1;
      dE[2][2] *= scrh1;   

      dE_2[0] = vc[EX1]*dE[0][0] + vc[EX2]*dE[1][0] + vc[EX3]*dE[2][0];   // WARNING!!!
      dE_2[1] = vc[EX1]*dE[0][1] + vc[EX2]*dE[1][1] + vc[EX3]*dE[2][1];   // This is 0.5dE^2/du_i !
      dE_2[2] = vc[EX1]*dE[0][2] + vc[EX2]*dE[1][2] + vc[EX3]*dE[2][2];

      #if EOS == IDEAL
      double hdot = gmm1;
      #elif EOS == TAUB
      double hdot = 2.5 + 2.25*theta/(h - 2.5*theta);
      #endif
      scrh1 = (Dh + vc[PRS]/glor);
      scrh2 = 1.0/(uc[RHO]*(1.0/hdot - glor*glor));
      dh[0] = (u[0]*scrh1 + glor*dE_2[0])*scrh2;
      dh[1] = (u[1]*scrh1 + glor*dE_2[1])*scrh2;
      dh[2] = (u[2]*scrh1 + glor*dE_2[2])*scrh2;

      J[0][0] = -D*u[0]*dh[0] - dE[1][0]*uc[BX3] + dE[2][0]*uc[BX2] - Dh;
      J[1][0] = -D*u[1]*dh[0] - dE[2][0]*uc[BX1] + dE[0][0]*uc[BX3];
      J[2][0] = -D*u[2]*dh[0] - dE[0][0]*uc[BX2] + dE[1][0]*uc[BX1];
      J[0][1] = -D*u[0]*dh[1] - dE[1][1]*uc[BX3] + dE[2][1]*uc[BX2];
      J[1][1] = -D*u[1]*dh[1] - dE[2][1]*uc[BX1] + dE[0][1]*uc[BX3] - Dh;
      J[2][1] = -D*u[2]*dh[1] - dE[0][1]*uc[BX2] + dE[1][1]*uc[BX1];
      J[0][2] = -D*u[0]*dh[2] - dE[1][2]*uc[BX3] + dE[2][2]*uc[BX2];
      J[1][2] = -D*u[1]*dh[2] - dE[2][2]*uc[BX1] + dE[0][2]*uc[BX3];
      J[2][2] = -D*u[2]*dh[2] - dE[0][2]*uc[BX2] + dE[1][2]*uc[BX1] - Dh;

      // 3f. Solve  J.du = -f using Gaussian elimination

      du[0] = -func[0];
      du[1] = -func[1];
      du[2] = -func[2];

      GaussianElimination(J, x, du, 3);

      u[0] += x[0];
      u[1] += x[1];
      u[2] += x[2];

      #if IMEX_DEBUG == 2
      if (debug_condition) {
        printf("      u[k=%d] = {%12.6e, %12.6e, %12.6e}\n",kcycle, u[0], u[1], u[2]);
      }
      #endif

      // 3g. Compute error. Stop computation if kcycle
      //     exceeds a critical value.

      err_g  = sqrt(du[0]*du[0] + du[1]*du[1] + du[2]*du[2])/glor;
      err_f  = sqrt(func[0]*func[0] + func[1]*func[1] + func[2]*func[2]);
      err_f /= sqrt(uc[MX1]*uc[MX1] + uc[MX2]*uc[MX2] + uc[MX3]*uc[MX3]) + 1.e-6;

      if (err_g < tol_g || err_f < tol_f) done = 1;

      #if IMEX_DEBUG == 2
      if (debug_condition) {
        printf ("      err_g = %12.6e, err_f = %12.6e\n", err_g, err_f);
      }
      #endif 

      if (kcycle == MAX_ITER-1) {
        ATOMIC_ERR_UPDATE(err_kernel,       ERR_TOO_MANY_ITERATIONS, ERR_QUIT);
        KERNEL_ERR_UPDATE(data->err(i,j,k), ERR_TOO_MANY_ITERATIONS, ERR_QUIT);
      } 
    } // end for ( MAX_ITER ) - Newton-Broyden iteration cycle 
    
    // 4. Update main array of primitive variables &
    //    electric field in conservative array.

    vc[VX1] = u[0]*glor_1;
    vc[VX2] = u[1]*glor_1;
    vc[VX3] = u[2]*glor_1;
    vc[RHO] = rho;   // D*glor_1;
    #if ENTROPY_SWITCH      // Recompute primitive entropy
    vc[ENTR] = Entropy(vc,run_config);
    uc[ENTR] = uc[RHO]*vc[ENTR];
    #endif

    NVAR_LOOP(nv) Vc(i,j,k,nv) = vc[nv];
   
    uc[EX1] = vc[EX1];
    uc[EX2] = vc[EX2];
    uc[EX3] = vc[EX3];
    
    #if GEOMETRY != CARTESIAN
    OrthoToCons (uc, metric); 
    #endif

    Uc(i,j,k,EX1) = uc[EX1];
    Uc(i,j,k,EX2) = uc[EX2];
    Uc(i,j,k,EX3) = uc[EX3];

    // 5. Compute the source term.

    S(i,j,k,EX1) = (uc[EX1] - uc0[EX1])/dt1;
    S(i,j,k,EX2) = (uc[EX2] - uc0[EX2])/dt1;
    S(i,j,k,EX3) = (uc[EX3] - uc0[EX3])/dt1;


/*
    S(i,j,k,EX1) = -(glor*vc[EX1] + uxB[IDIR] - Ev*u[IDIR])/eta;
    S(i,j,k,EX2) = -(glor*vc[EX2] + uxB[JDIR] - Ev*u[JDIR])/eta;
    S(i,j,k,EX3) = -(glor*vc[EX3] + uxB[KDIR] - Ev*u[KDIR])/eta;
*/
    #if ENTROPY_SWITCH

    double Ev = vc[EX1]*vc[VX1] + vc[EX2]*vc[VX2] + vc[EX3]*vc[VX3];

    Uc(i,j,k,ENTR) = uc[ENTR]; // Recompute conservative entropy
    
    double q[3], q2, chi_1;
    q[IDIR] = (vc[EX1] + vc[VX2]*vc[BX3] - vc[VX3]*vc[BX2]);
    q[JDIR] = (vc[EX2] + vc[VX3]*vc[BX1] - vc[VX1]*vc[BX3]);
    q[KDIR] = (vc[EX3] + vc[VX1]*vc[BX2] - vc[VX2]*vc[BX1]);
    q2 = DOT_PRODUCT(q,q);

    #if EOS == IDEAL
    chi_1 = gmm - 1.0;
    #elif EOS == TAUB
    chi_1 = 2.0/3.0;
    #endif

    S(i,j,k,ENTR) = glor*glor*chi_1*(q2 - Ev*Ev)*vc[ENTR]*vc[RHO]/(vc[PRS]*eta);

    #endif // ENTROPY_SWITCH 

  } // End BOX_LOOP()  

  // ------------------------------------------------------
  // 5. Update staggered electric field
  // ------------------------------------------------------

  #if DIVE_CONTROL == CONSTRAINED_TRANSPORT
  CT_IMEXImplicitUpdate(data, Jf, dt1, grid);
  ConsToPrim3D(data, iarr, grid, box);
  #else
  #endif

  // 6. Check errors and quit condition(s)

  #if WARNING_MESSAGES > 0 
  Where ("ImexImplicitUpdate()", data, err_kernel, box, grid);
  #endif

  #if WARNING_MESSAGES >= 0 
  if (err_kernel & ERR_QUIT) QuitPLUTO("from ImexImplicitUpdate()", err_kernel);
  #endif

  mynvtxstop_();
}

////////////////////////////////////////////////////////////////////////
void IMEXExplicitSource (Data *data, Ary4D &S, Grid *grid)
/*!
  Compute the explicit source term in Ampere's law.
*/
////////////////////////////////////////////////////////////////////////
{
  int i,j,k, nv;
  int ibeg = grid->lbeg[IDIR], iend = grid->lend[IDIR];
  int jbeg = grid->lbeg[JDIR], jend = grid->lend[JDIR];
  int kbeg = grid->lbeg[KDIR], kend = grid->lend[KDIR];

  double Ev, glor, eta;
  double *x1 = grid->xC[IDIR];
  double *x2 = grid->xC[JDIR];
  double *x3 = grid->xC[KDIR];

  double vc[NVAR], uc[NVAR], u[3], uXB[3];
  RBox   box;

  RBoxDefine (ibeg-1, iend+1, jbeg-1, jend+1, kbeg-1, kend+1, &box);

  #pragma omp target teams distribute parallel for collapse(3) \
                     firstprivate(vc, uc, u, uXB) map(to:S, data->run_config) //map (to:data[:1],S[:1])
  #pragma acc parallel loop collapse(3) private(vc[:NVAR], uc[:NVAR], u[:3], uXB[:3]) \
                       present(data, data->run_config, S) copyin(box)
  BOX_LOOP(box, i,j,k){

    State state;

  // -- Copy to local array -- 
    #pragma acc loop
    NVAR_LOOP(nv){
      #ifdef HIGH_ORDER             // With high order we         
      vc[nv] = data->Vpc(i,j,k,nv); // operate on the point value 
      uc[nv] = data->Upc(i,j,k,nv);
      #else
      vc[nv] = data->Vc(i,j,k,nv);
      uc[nv] = data->Uc(i,j,k,nv);
      #endif
    }

    STATE_SET(state, vc, uc, NULL, x1[i], x2[j], x3[k]);

    eta = Resistive_eta(state, *(data->run_config));  

    Ev = vc[EX1]*vc[VX1] + vc[EX2]*vc[VX2] + vc[EX3]*vc[VX3];
    glor = 1.0/sqrt(1.0 - vc[VX1]*vc[VX1] - vc[VX2]*vc[VX2] - vc[VX3]*vc[VX3]);

    u[0]  = vc[VX1]*glor;  // Four velocity 
    u[1]  = vc[VX2]*glor;
    u[2]  = vc[VX3]*glor;

    uXB[0] = u[1]*vc[BX3] - u[2]*vc[BX2];
    uXB[1] = u[2]*vc[BX1] - u[0]*vc[BX3];
    uXB[2] = u[0]*vc[BX2] - u[1]*vc[BX1];
    
    #pragma acc loop
    NVAR_LOOP(nv) S(i,j,k,nv) =  0.0;

    S(i,j,k,EX1) = -(glor*vc[EX1] + uXB[IDIR] - Ev*u[IDIR])/eta;
    S(i,j,k,EX2) = -(glor*vc[EX2] + uXB[JDIR] - Ev*u[JDIR])/eta;
    S(i,j,k,EX3) = -(glor*vc[EX3] + uXB[KDIR] - Ev*u[KDIR])/eta;
  }
}

////////////////////////////////////////////////////////////////////////
#pragma omp declare target
int IdealRMHD_Guess(double *uc0, double *vg, RunConfig *run_config)
/*!
  \return Error codes: see RMHD_EneergySolve()
*/
////////////////////////////////////////////////////////////////////////
{
  uint32_t err;

  err = RMHD_EnergySolve (uc0, vg, *run_config);

  vg[BX1] = uc0[BX1];
  vg[BX2] = uc0[BX2];
  vg[BX3] = uc0[BX3];

  vg[EX1] = -(vg[VX2]*vg[BX3] - vg[VX3]*vg[BX2]);
  vg[EX2] = -(vg[VX3]*vg[BX1] - vg[VX1]*vg[BX3]);
  vg[EX3] = -(vg[VX1]*vg[BX2] - vg[VX2]*vg[BX1]);

  return err;
}
#pragma omp end declare target
#undef IMEX_DEBUG
#undef  ENERGY_SOLVE
#undef  ENTROPY_SOLVE
