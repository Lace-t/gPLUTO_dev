////////////////////////////////////////////////////////////////////////
/*!
  \file
  \brief Conservative CT correction

  Compute the Poynting vector at zone faces using the upwinded emf. 
  This function is called only when:

  - CT_ENERGY_CORRECTION == 2 (average of Poynting flux): 
    in this case the Poynting flux is computed by averaging the Poynting 
    vector at edges:

    F(i-1/2) = < EB > = [ (EB)^NW + EB^(SW) ]/2

  - CT_ENERGY_CORRECTION == 3 (Poynting flux of average):
    in this case the Poynting flux is computed  
    by averaging the electric and magnetic fields and the computing the
    Poynting flux:

    F(i-1/2) = < E > < B > = [ (E)^NW + EB^(SW) ][ (E)^NW + EB^(SW) ]/4

  \authors A. Mignone    (andrea.mignone@unito.it)
  \date    Apr 1, 2026

  \b Reference
*/
////////////////////////////////////////////////////////////////////////
#include "gpluto.hpp"

#define QUAD_RULE(Elo, Ehi, Blo, Bhi) (0.25*(Elo + Ehi)*(Blo + Bhi))
//#define QUAD_RULE(Elo, Ehi, Blo, Bhi) ( (  (2.0*Elo + Ehi)*Blo \
//                                         + (Elo + 2.0*Ehi)*Bhi)/6.0)

////////////////////////////////////////////////////////////////////////
/*! 
  \param [in]     data  pointer to data structure
  \param [in,out] dUc   the conservative variable increment
  \param [in]     dt    the current time step 
  \param [in]     grid  pointer to grid structure
  \param [in]     rbox  Box on the operation is performed
*/
////////////////////////////////////////////////////////////////////////
void CT_EnergyCorrection (Data *data, Ary4D &dUc, double dt, Grid *grid, RBox *box)
{
#if CT_ENERGY_CORRECTION > 1
  int i,j,k;
  EMF &emf = *(data->emf);
  const Ary3D &Ex1e = emf.Ee[IDIR];
  const Ary3D &Ex2e = emf.Ee[JDIR];
  const Ary3D &Ex3e = emf.Ee[KDIR];

  // ----------------------------------------------
  // 1. Compute Poynting vector from edges to center
  // 
  // By at i+1/2,j+1/2,k
  // Bx at i+1/2,j,k+1/2
  // ----------------------------------------------

  BOX_STAG_LOOP(*box, i,j,k){
    double ExBy = 0.0;
    double ExBz = 0.0;
    double EyBx = 0.0;
    double EyBz = 0.0;
    double EzBy = 0.0;
    double EzBx = 0.0;

    #if CT_ENERGY_CORRECTION == 2
    DIM_EXPAND(                                                     ,
               EzBy = 0.5*(  Ex3e(i,j,k)  *emf.Bx2ez(i,j,k)
                           + Ex3e(i,j+1,k)*emf.Bx2ez(i,j+1,k)); ,
               EyBz = 0.5*(  Ex2e(i,j,k)  *emf.Bx3ey(i,j,k)
                           + Ex2e(i,j,k+1)*emf.Bx3ey(i,j,k+1)) ); 

    emf.Pfx(i,j,k) = EyBz - EzBy;
                       
    DIM_EXPAND(                                                     ,
               EzBx = 0.5*(  Ex3e(i,j,k)  *emf.Bx1ez(i,j,k)
                           + Ex3e(i+1,j,k)*emf.Bx1ez(i+1,j,k)); ,
               ExBz = 0.5*(  Ex1e(i,j,k)  *emf.Bx3ex(i,j,k)
                           + Ex1e(i,j,k+1)*emf.Bx3ex(i,j,k+1)); ) 

    emf.Pfy(i,j,k) = EzBx - ExBz;

    #if INCLUDE_KDIR
    DIM_EXPAND(                                                     ,
               ExBy = 0.5*(  Ex1e(i,j,k)  *emf.Bx2ex(i,j,k)
                           + Ex1e(i,j+1,k)*emf.Bx2ex(i,j+1,k)); ,
               EyBx = 0.5*(  Ex2e(i,j,k)  *emf.Bx1ey(i,j,k)
                           + Ex2e(i+1,j,k)*emf.Bx1ey(i+1,j,k)); ) 
    emf.Pfz(i,j,k) = ExBy - EyBx;
    #endif

    #endif // CT_ENERGY_CORRECTION == 2

    // --------------------------------------------
    // 2. Second version
    // --------------------------------------------

    #if CT_ENERGY_CORRECTION == 3
    double Elo, Ehi, Blo, Bhi;

  // -- 2a. Compute Sx -- 

#if 0
    DIM_EXPAND(                                                     ,
               EzBy =   0.5*(emf.Ex3e(i,j,k) + emf.Ex3e(i,j+1,k))
                       *0.5*(emf.Bx2ez(i,j,k)+ emf.Bx2ez(i,j+1,k)); ,
               EyBz =   0.5*(emf.Ex2e(i,j,k) + emf.Ex2e(i,j,k+1))
                       *0.5*(emf.Bx3ey(i,j,k)+ emf.Bx3ey(i,j,k+1)); )
    emf.Pfx(i,j,k) = EyBz - EzBy;
#else
    DIM_EXPAND(                                                       ,
               Elo  = Ex3e (i,j,k); Ehi = Ex3e (i,j+1,k); 
               Blo  = emf.Bx2ez(i,j,k); Bhi = emf.Bx2ez(i,j+1,k);  
               EzBy = QUAD_RULE(Elo, Ehi, Blo, Bhi); 
                                                                      ,
               Elo  = Ex2e (i,j,k); Ehi = Ex2e (i,j,k+1);
               Blo  = emf.Bx3ey(i,j,k); Bhi = emf.Bx3ey(i,j,k+1);
               EyBz = QUAD_RULE(Elo, Ehi, Blo, Bhi);)

//    EzBy = (2.0*Elo + Ehi)*Blo/6.0 + Bhi*(Elo + 2.0*Ehi)/6.0;
    // DIM_EXPAND(                                             ,
    //            EzBy = ((2.0*Ez_lo + Ez_hi)*By_lo + (Ez_lo + 2.0*Ez_hi)*By_hi)/6.0; ,
    //            EyBz = ((2.0*Ey_lo + Ey_hi)*Bz_lo + (Ez_lo + 2.0*Ez_hi)*By_hi)/6.0;)

    emf.Pfx(i,j,k) = EyBz - EzBy;
#endif

  // -- 2b. Compute Sy -- 

#if 0
    DIM_EXPAND(                                                     ,
               EzBx =  0.5*(Ex3e(i,j,k)  + Ex3e(i+1,j,k))
                      *0.5*(emf.Bx1ez(i,j,k) + emf.Bx1ez(i+1,j,k)); ,
               ExBz =  0.5*(Ex1e(i,j,k)  + Ex1e(i,j,k+1))
                      *0.5*(emf.Bx3ex(i,j,k) + emf.Bx3ex(i,j,k+1)); ) 
    emf.Pfy(i,j,k) = EzBx - ExBz;
#else
    DIM_EXPAND(                                                       ,
               Elo  = Ex3e (i,j,k); Ehi = Ex3e (i+1,j,k); 
               Blo  = emf.Bx1ez(i,j,k); Bhi = emf.Bx1ez(i+1,j,k);  
               EzBx = QUAD_RULE(Elo, Ehi, Blo, Bhi);                  ,
               Elo  = Ex1e (i,j,k); Ehi = Ex1e (i,j,k+1);
               Blo  = emf.Bx3ex(i,j,k); Bhi = emf.Bx3ex(i,j,k+1);
               ExBz = QUAD_RULE(Elo, Ehi, Blo, Bhi);  )

    emf.Pfy(i,j,k) = EzBx - ExBz;
#endif

    #if INCLUDE_KDIR
#if 0
    DIM_EXPAND(                                                     ,
               ExBy =  0.5*(Ex1e(i,j,k)  + Ex1e(i,j+1,k))
                      *0.5*(emf.Bx2ex(i,j,k) + emf.Bx2ex(i,j+1,k)); ,
               EyBx =  0.5*(Ex2e(i,j,k)  + Ex2e(i+1,j,k))
                      *0.5*(emf.Bx1ey(i,j,k) + emf.Bx1ey(i+1,j,k)); )

    emf.Pfz(i,j,k) = ExBy - EyBx;
#else
    DIM_EXPAND(                                                       ,
               Elo  = Ex1e (i,j,k); Ehi = Ex1e (i,j+1,k); 
               Blo  = emf.Bx2ex(i,j,k); Bhi = emf.Bx2ex(i,j+1,k);  
               ExBy = QUAD_RULE(Elo, Ehi, Blo, Bhi);                   ,
               Elo  = Ex2e (i,j,k); Ehi = Ex2e (i+1,j,k);
               Blo  = emf.Bx1ey(i,j,k); Bhi = emf.Bx1ey(i+1,j,k);
               EyBx = QUAD_RULE(Elo, Ehi, Blo, Bhi);)

    emf.Pfz(i,j,k) = ExBy - EyBx;
#endif

    #endif

    #endif // CT_ENERGY_CORRECTION == 3

    // --------------------------------------------
    // 3. Third version
    // --------------------------------------------

    #if CT_ENERGY_CORRECTION == 4
    DIM_EXPAND(                                                         ,
               EzBy = 0.5*(  emf.Ex3e_phi(i,j,k)  *emf.Bx2ez(i,j,k)
                           + emf.Ex3e_phi(i,j+1,k)*emf.Bx2ez(i,j+1,k)); ,
               EyBz = 0.5*(  emf.Ex2e_phi(i,j,k)  *emf.Bx3ey(i,j,k)
                           + emf.Ex2e_phi(i,j,k+1)*emf.Bx3ey(i,j,k+1)) ); 
    emf.Pfx(i,j,k) = EyBz - EzBy;
                       
    DIM_EXPAND(                                                     ,
               EzBx = 0.5*(  emf.Ex3e_phi(i,j,k)  *emf.Bx1ez(i,j,k)
                           + emf.Ex3e_phi(i+1,j,k)*emf.Bx1ez(i+1,j,k)); ,
               ExBz = 0.5*(  emf.Ex1e_phi(i,j,k)  *emf.Bx3ex(i,j,k)
                           + emf.Ex1e_phi(i,j,k+1)*emf.Bx3ex(i,j,k+1)) ); 
    emf.Pfy(i,j,k) = EzBx - ExBz;

    #if INCLUDE_KDIR
    DIM_EXPAND(                                                     ,
               ExBy = 0.5*(  emf.Ex1e_phi(i,j,k)  *emf.Bx2ex(i,j,k)
                           + emf.Ex1e_phi(i,j+1,k)*emf.Bx2ex(i,j+1,k)); ,
               EyBx = 0.5*(  emf.Ex2e_phi(i,j,k)  *emf.Bx1ey(i,j,k)
                           + emf.Ex2e_phi(i+1,j,k)*emf.Bx1ey(i+1,j,k)) ); 
    emf.Pfz(i,j,k) = ExBy - EyBx;
    #endif

    #endif // CT_ENERGY_CORRECTION == 4

  }

  // --------------------------------------------------------
  // 2. Update energy
  // --------------------------------------------------------

  BOX_LOOP(*box, i,j,k){
    double dtdx = dt/grid->dx[IDIR][i];
    double dtdy = dt/grid->dx[JDIR][j];
    double dtdz = dt/grid->dx[KDIR][k];

    dUc(i,j,k,ENG) -= DIM_EXPAND(  dtdx*(emf.Pfx(i+1,j,k) - emf.Pfx(i,j,k)), 
                                 + dtdy*(emf.Pfy(i,j+1,k) - emf.Pfy(i,j,k)), 
                                 + dtdz*(emf.Pfz(i,j,k+1) - emf.Pfz(i,j,k))); 


  }
#endif
}


#define VEL_SQUARE(v) (v[VX1]*v[VX1] + v[VX2]*v[VX2] + v[VX3]*v[VX3])
#define SQUARE(x)  ((x)*(x))
////////////////////////////////////////////////////////////////////////
/*!
  Compute the hydro energy flux and the hydrodynamic energy.
*/
////////////////////////////////////////////////////////////////////////
double HD_EnergyFlux(double *v, double *Ehd, int dir, RunConfig &run_config)
{
#if CT_ENERGY_CORRECTION > 1
#if HAVE_ENERGY
  #if PHYSICS == MHD
  double *u = v + VX1;
  double *B = v + BX1;

  double Ek  = 0.5*v[RHO]*VEL_SQUARE(v);  
  (*Ehd) = Ek + v[PRS]/(run_config.gamma - 1.0);

/*
In 2.5D we also have to include

div(F) = div(Fg) + (Ey*Bz - Ez*By)_x + (Ez*Bx - Ex*Bz)_y 
                   [Cent]   [CT]        [CT]    [Cent]

  The two terms [Cent] must be included here in the standard cell-centered
  way:

   Ey*Bz = -(vz*Bx - vx*Bz)*Bz = vx*Bz^2 - (vz*Bz)*Bx  
  -Ex*Bz =  (vy*Bz - vz*By)*Bz = vy*Bz^2 - (vz*Bz)*By

  The extra terms include out of the plane components.
*/
double pt = v[PRS], Fhd;
#if (INCLUDE_IDIR+INCLUDE_JDIR+INCLUDE_KDIR) < 3
double pm = 0.5*DIM_EXPAND_NOT(v[BX1]*v[BX1], + v[BX2]*v[BX2], + v[BX3]*v[BX3]);
(*Ehd) += pm;
pt     += pm;
#endif

Fhd     = ((*Ehd) + pt)*v[VX1 + dir];
#if (INCLUDE_IDIR+INCLUDE_JDIR+INCLUDE_KDIR) < 3
double vB  = DIM_EXPAND_NOT(v[VX1]*v[BX1], + v[VX2]*v[BX2], + v[VX3]*v[BX3]);
Fhd   -= v[BX1 + dir]*vB;
#endif
  #if (CT_ENERGY_CORRECTION == 2) || (CT_ENERGY_CORRECTION == 3)
  return Fhd;
  #endif

  // #if (CT_ENERGY_CORRECTION == 2) || (CT_ENERGY_CORRECTION == 3)
  // return ((*Ehd) + v[PRS])*v[VX1 + dir];
  // #endif

  #elif PHYSICS == ResMHD

  double  h = Enthalpy(v, run_config);

  double vel2 = v[VX1]*v[VX1] + v[VX2]*v[VX2] + v[VX3]*v[VX3];
  double g2   = 1.0/(1.0 - vel2);
  double g    = sqrt(g2);
  double rhoh = v[RHO]*h;

  (*Ehd) = rhoh*g2 - v[PRS];
  return rhoh*g2*v[VX1+dir];

  #endif

#endif
#endif

  return 0;
}


////////////////////////////////////////////////////////////////////////
/*!
  Compute magnetic energy depending on the value of CT_MAGNETIC_ENERGY:
  
  1:    Pm = ((Bxc)^2 + (Byc)^2)/2, where Bxc = (BxE + BxW)/2
  Byc = (ByN + ByS)/2   
  2:    Pm = (BxE^2 + BxW^2 + ByN^2 + ByS^2)/4 
  3:    Pm =   (Bxp^2 + Bxm^2 + Bxp*Bxm)/6
  + (Byp^2 + Bym^2 + Byp*Bym)/6
  
*/
////////////////////////////////////////////////////////////////////////
double CT_MagneticEnergy(const Data *data, int i, int j, int k)
{
#if CT_MAGNETIC_ENERGY
  double Bxp = 0.0, Bxm = 0.0;
  double Byp = 0.0, Bym = 0.0;
  double Bzp = 0.0, Bzm = 0.0;

  DIM_EXPAND(Bxp = data->Vs[BX1s](i+1,j,k);
             Bxm = data->Vs[BX1s](i,j,k);    ,
             Byp = data->Vs[BX2s](i,j+1,k);
             Bym = data->Vs[BX2s](i,j,k);    ,
             Bzp = data->Vs[BX3s](i,j,k+1);
             Bzm = data->Vs[BX3s](i,j,k);)

  #if CT_MAGNETIC_ENERGY == 1
  double Pm1  =   0.5*(DIM_EXPAND(  SQUARE(0.5*(Bxp + Bxm))  ,
                                  + SQUARE(0.5*(Byp + Bym))  , 
                                  + SQUARE(0.5*(Bzp + Bzm)) ));
  return Pm1;
  #elif CT_MAGNETIC_ENERGY == 2
  double Pm2  = 0.25*(DIM_EXPAND(  Bxp*Bxp + Bxm*Bxm,
                                 + Byp*Byp + Bym*Bym,
                                 + Bzp*Bzp + Bzm*Bzm));
  return Pm2;
  #elif CT_MAGNETIC_ENERGY == 3
  double Pm3  = (DIM_EXPAND(  Bxp*Bxp + Bxm*Bxm + Bxp*Bxm,
                            + Byp*Byp + Bym*Bym + Byp*Bym, 
                            + Bzp*Bzp + Bzm*Bzm + Bzp*Bzm)) /6.0;
  return Pm3;
  #else
  QuitPLUTO ("! Invalid CT_MAGNETIC_ENERGY value. Stop in CT_MagneticEergy",1);
  #endif
#endif

  return -1;
}