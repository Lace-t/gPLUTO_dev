////////////////////////////////////////////////////////////////////////
/*!
  \file
  \brief Shared header file used for velocity limiting / conversion


  \author A. Mignone (andrea.mignone@unito.it) 
  \date    Feb 17, 2025
*/
////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////
#pragma acc routine seq
inline void VelocityLimiter (Ary2D &v, Ary2D &vL, Ary2D &vR, const int i)
/*!
  Check whether the total reconstructed velocity is > 1
  If a superluminal value occurs, flatten distribution.
*/
////////////////////////////////////////////////////////////////////////
{
#if RECONSTRUCT_4VEL == NO
  int    nv;
  double v2m, v2p;

  v2m = vR(i,VX1)*vR(i,VX1)  + vR(i,VX2)*vR(i,VX2) + vR(i,VX3)*vR(i,VX3);
  v2p = vL(i+1,VX1)*vL(i+1,VX1) + vL(i+1,VX2)*vL(i+1,VX2) + vL(i+1,VX3)*vL(i+1,VX3);
  if (v2m >= 1.0 || v2p >= 1.0){
    #pragma acc loop seq
    for (nv = 0; nv < NVAR; nv++) vR(i,nv) = vL(i+1,nv) = v(i,nv);
  }
#endif

  #if !RECONSTRUCT_4VEL
  (void)v; (void)vL; (void)vR; (void)i;
  #endif
}

////////////////////////////////////////////////////////////////////////
#pragma acc routine vector
inline void ConvertTo4vel (Ary2D &v, int i)
/*!
*/
////////////////////////////////////////////////////////////////////////
{
  double lor;

  lor = v(i,VX1)*v(i,VX1) + v(i,VX2)*v(i,VX2) + v(i,VX3)*v(i,VX3);
  lor = 1.0/sqrt(1.0 - lor);

  v(i,VX1) *= lor;
  v(i,VX2) *= lor;
  v(i,VX3) *= lor;
}

////////////////////////////////////////////////////////////////////////
#pragma acc routine vector
inline void   ConvertTo3vel (Ary2D &v, int i)
/*!
*/
////////////////////////////////////////////////////////////////////////
{
  double lor;

  lor = v(i,VX1)*v(i,VX1) + v(i,VX2)*v(i,VX2) + v(i,VX3)*v(i,VX3);
  lor = sqrt(1.0 + lor);
  
  v(i,VX1) /= lor;
  v(i,VX2) /= lor;
  v(i,VX3) /= lor;
}
