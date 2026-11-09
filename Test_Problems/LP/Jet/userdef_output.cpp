#include "gpluto.hpp"

////////////////////////////////////////////////////////////////////////
void UserDefOutput (const Data * data, Grid * grid)
/*!
*/
////////////////////////////////////////////////////////////////////////
{

  int i,j,k;
  const int nelem = 1;
  static Ary4D Q;
  Ary2D deposit_arr;
  int nx1_tot = grid->ntot[IDIR];
  int nx2_tot = grid->ntot[JDIR];
  int nx3_tot = grid->ntot[KDIR];
  RBox box;
  RBoxDefine (grid->lbeg, grid->lend, &box);
  #pragma acc enter data copyin(box)
  // 0a. Allocate arrays 
  uint64_t nparticles = data->pc->Size();

  static Ary3D J250MHz, J1GHz, J5GHz, J10GHz, J20GHz;
  if (J250MHz.isAllocated() == NO){
    J250MHz.create (nx1_tot,nx2_tot,nx3_tot);
    J1GHz.create (nx1_tot,nx2_tot,nx3_tot);
    J5GHz.create (nx1_tot,nx2_tot,nx3_tot);
    J10GHz.create (nx1_tot,nx2_tot,nx3_tot);
    J20GHz.create (nx1_tot,nx2_tot,nx3_tot);
    Q.create (nx1_tot,nx2_tot,nx3_tot,nelem);
    #pragma acc enter data copyin(J250MHz, J1GHz, J5GHz, J10GHz, J20GHz, Q)
  }

  deposit_arr.create(nparticles, nelem); 
  #pragma acc enter data copyin(deposit_arr)
  // 0b. Link arrays to a particular output

  LinkUserVar(&J250MHz, "J250MHz", FLT_OUTPUT);
  LinkUserVar(&J1GHz, "J1GHz", FLT_OUTPUT);
  LinkUserVar(&J5GHz, "J5GHz", FLT_OUTPUT);
  LinkUserVar(&J10GHz, "J10GHz", FLT_OUTPUT);
  LinkUserVar(&J20GHz, "J20GHz", FLT_OUTPUT);

  // 1. Define quantities to be deposited  

  //#pragma acc parallel loop present(data, data->pc, data->run_config, deposit_arr)
  #pragma acc update host(data->pc)

  // [AS]: Note that this loop has to be done on CPU because of how the bessel funcions are treated in Particles_LP_Sync_Emissivity.
  // This means that we need to call #pragma acc update host(data->pc) and then #pragma acc update device(deposit_arr)
  // To run on GPU we need to allocate the bessel functions before calling Particles_LP_Sync_Emissivity
  for (uint64_t np = 0; np < nparticles; np++){

    double nu    = 0.25; // Frequency in Ghz
    double theta =  90.0; // View angle 
    double JnuSync, StkQ, StkU;
    Particles_LP_Sync_Emissivity(*data->pc, *data->run_config, nu, theta,
                                 &JnuSync, &StkQ, &StkU, np);

    // Copy relevant quantities on deposit_arr()

    deposit_arr(np, 0) = JnuSync;
    // deposit_arr(np, 1) = StkQ;  // Need this ? //onnly if you savealso stokes parameters Q and U
    // deposit_arr(np, 2) = StkU;
  }

  #pragma acc update device(deposit_arr)
  Particles_Deposit (*data->pc, Q, deposit_arr, nelem, grid);

  // 2. Copy deposited quantities on 3D Arrays before writing
  
  #pragma acc parallel loop collapse(3) present(Q, J250MHz, box)
  BOX_LOOP(box,i,j,k){
    J250MHz(i,j,k) = Q(i,j,k,0);
  }

  // 1. Define quantities to be deposited  

  // #pragma acc parallel loop present(data, data->pc, data->run_config, deposit_arr)
  for (uint64_t np = 0; np < nparticles; np++ ){

    double nu    = 1.; // Frequency in Ghz
    double theta =  90.0; // View angle 
    double JnuSync, StkQ, StkU;
    Particles_LP_Sync_Emissivity(*data->pc, *data->run_config, nu, theta,
                                 &JnuSync, &StkQ, &StkU, np);

    // Copy relevant quantities on deposit_arr()

    deposit_arr(np, 0) = JnuSync;
    // deposit_arr(np, 1) = StkQ;  // Need this ? 
    // deposit_arr(np, 2) = StkU;
  }

  #pragma acc update device(deposit_arr)
  Particles_Deposit (*data->pc, Q, deposit_arr, nelem, grid);

  // 2. Copy deposited quantities on 3D Arrays before writing

  #pragma acc parallel loop collapse(3) present(Q, J1GHz, box)
  BOX_LOOP(box,i,j,k){
    J1GHz(i,j,k) = Q(i,j,k,0);
  }
  
  // 1. Define quantities to be deposited  
  
  // #pragma acc parallel loop present(data, data->pc, data->run_config, deposit_arr)
  for (uint64_t np = 0; np < nparticles; np++){

    double nu    = 5.; // Frequency in Ghz
    double theta =  90.0; // View angle 
    double JnuSync, StkQ, StkU;
    Particles_LP_Sync_Emissivity(*data->pc, *data->run_config, nu, theta,
                                 &JnuSync, &StkQ, &StkU, np);

    // Copy relevant quantities on deposit_arr()

    deposit_arr(np, 0) = JnuSync;
    // deposit_arr(np, 1) = StkQ;  // Need this ? 
    // deposit_arr(np, 2) = StkU;
  }

  #pragma acc update device(deposit_arr)
  Particles_Deposit (*data->pc, Q, deposit_arr, nelem, grid);

  // 2. Copy deposited quantities on 3D Arrays before writing

  #pragma acc parallel loop collapse(3) present(Q, J5GHz, box)
  BOX_LOOP(box,i,j,k){
    J5GHz(i,j,k) = Q(i,j,k,0);
  }
  
  // 1. Define quantities to be deposited  

  // #pragma acc parallel loop present(data, data->pc, data->run_config, deposit_arr)
  for (uint64_t np = 0; np < nparticles; np++ ){

    double nu    = 10.; // Frequency in Ghz
    double theta =  90.0; // View angle 
    double JnuSync, StkQ, StkU;
    Particles_LP_Sync_Emissivity(*data->pc, *data->run_config, nu, theta,
                                 &JnuSync, &StkQ, &StkU, np);

    // Copy relevant quantities on deposit_arr()

    deposit_arr(np, 0) = JnuSync;
    // deposit_arr(np, 1) = StkQ;  // Need this ? 
    // deposit_arr(np, 2) = StkU;
  }

  #pragma acc update device(deposit_arr)

  Particles_Deposit (*data->pc, Q, deposit_arr, nelem, grid);

  // 2. Copy deposited quantities on 3D Arrays before writing

  #pragma acc parallel loop collapse(3) present(Q, J10GHz, box)
  BOX_LOOP(box,i,j,k){
    J10GHz(i,j,k) = Q(i,j,k,0);
  }  
  
   // 1. Define quantities to be deposited  
  // #pragma acc parallel loop present(data, data->pc, data->run_config, deposit_arr)
  for (uint64_t np = 0; np < nparticles; np++ ){

    double nu    = 20.; // Frequency in Ghz
    double theta =  90.0; // View angle 
    double JnuSync, StkQ, StkU;
    Particles_LP_Sync_Emissivity(*data->pc, *data->run_config, nu, theta,
                                 &JnuSync, &StkQ, &StkU, np);

    // Copy relevant quantities on deposit_arr()

    deposit_arr(np, 0) = JnuSync;
    // deposit_arr(np, 1) = StkQ;  // Need this ? 
    // deposit_arr(np, 2) = StkU;
  }

  #pragma acc update device(deposit_arr)

  Particles_Deposit (*data->pc, Q, deposit_arr, nelem, grid);

  // 2. Copy deposited quantities on 3D Arrays before writing

  #pragma acc parallel loop collapse(3) present(Q, J20GHz, box)
  BOX_LOOP(box,i,j,k){
    J20GHz(i,j,k) = Q(i,j,k,0);
  } 
 
  #pragma acc update host(J250MHz, J1GHz, J5GHz, J10GHz, J20GHz)
  #pragma acc exit data delete(deposit_arr, box)
  
}


