////////////////////////////////////////////////////////////////////////
/*!
  \file
  \brief  Miscellaneous math functions.

  \author A. Mignone (andrea.mignone@unito.it)
  \date   March 09, 2026
*/
////////////////////////////////////////////////////////////////////////
#include "gpluto.hpp"

////////////////////////////////////////////////////////////////////////
/*!
  Sort an integer array x in ascending order.
  Taken from http://www.cquestions.com/2008/01/c-program-for-quick-sort.html
  
  \param [in,out] *x     the array to be sorted.
  \param [in]     first  the starting element of the array
  \param [in]     last   the final element of the array
*/
////////////////////////////////////////////////////////////////////////
void QuickSort(int *x,int first,int last)
{
  int pivot,j,temp,i;

  if (first < last){
    pivot = first;
    i     = first;
    j     = last;

    while(i<j){
      while(x[i]<=x[pivot]&&i<last)
        i++;
      while(x[j]>x[pivot])
        j--;
      if(i<j){
        temp=x[i];
        x[i]=x[j];
        x[j]=temp;
      }
    }

    temp=x[pivot];
    x[pivot]=x[j];
    x[j]=temp;
    QuickSort(x,first,j-1);
    QuickSort(x,j+1,last);
  }
}

////////////////////////////////////////////////////////////////////////
/*!
  Transform the vector v = (v_{x1}, v_{x2}, v_{x3}) from
  the input coordinate system to the actual (= GEOMETRY) system.
  
  - Input:  v[0..2] = components in the "input_geom" basis
  - Output: v[IDIR..KDIR] = components in the "GEOMETRY" basis
  
  \param [in,out]  *v     The original vector (overwritten by transformed one).
  \param [in] x1,x2,x3    Grid coordinates in the target geometry system 
  given by GEOMETRY
  \param [in] input_geom  The geometry of the input vector
*/
////////////////////////////////////////////////////////////////////////
void VectorTransform(double *v, double x1, double x2, double x3, int input_geom)
{
  double vx1 = v[0];
  double vx2 = v[1];
  double vx3 = v[2];

#if GEOMETRY == CARTESIAN
  // Target: Cartesian
  //         {x1, x2, x3} = {x, y, z}

  if (input_geom == CARTESIAN) {
    // Already Cartesian
    v[IDIR] = vx1;
    v[JDIR] = vx2;
    v[KDIR] = vx3;
  }

  if (input_geom == CYLINDRICAL) {
    // (x, y, z) -> (r, phi, z)
    double sphi = x1/sqrt(x1*x1 + x2*x2);
    double cphi = x2/sqrt(x1*x1 + x2*x2);

    v[IDIR] = vx1*cphi - vx2*sphi;
    v[JDIR] = vx1*sphi + vx2*cphi;
    v[KDIR] = vx3;
  }

  if (input_geom == SPHERICAL) {
    double R    = sqrt(x1*x1 + x2*x2);          // Cylindrical radius
    double r    = sqrt(x1*x1 + x2*x2 + x3*x3);  // Spherical radius
    double cth  = x3/r;     // cos(theta)
    double sth  = R/r;      // sin(theta)
    double tphi = x2/x1;    // tan(phi)
    double cphi = 1.0/sqrt(1.0 + tphi*tphi);
    double sphi = tphi*cphi;

    v[IDIR] = vx1*sth*cphi + vx2*cth*cphi - vx3*sphi;
    v[JDIR] = vx1*sth*sphi + vx2*cth*sphi + vx3*cphi;
    v[KDIR] = vx1*cth - vx2*sth;
  }

#elif GEOMETRY == CYLINDRICAL
  // Target: Cylindrical 
  //         {x1, x2, x3} = {R, phi, z}

  if (input_geom == CARTESIAN) {
    v[IDIR] =  vx1*cos(x2) + vx2*sin(x2);   // v_r
    v[JDIR] = -vx1*sin(x2) + vx2*cos(x2);   // v_phi
    v[KDIR] =  vx3;                         // v_z
  }

  if (input_geom == CYLINDRICAL) {
    // Already cylindrical
    v[IDIR] = vx1;
    v[JDIR] = vx2;
    v[KDIR] = vx3;
  }

  if (input_geom == SPHERICAL) {
    double r    = sqrt(x1*x1 + x3*x3);  // Spherical radius
    double R    = fabs(x1);             // Cylindrical radius
    double cth  = R/r;
    double sth  = x3/r;
    double cphi = cos(x2);
    double sphi = sin(x2);

    // First: spherical -> Cartesian
    double vx = vx1*sth*cphi + vx2*cth*cphi - vx3*sphi;
    double vy = vx1*sth*sphi + vx2*cth*sphi + vx3*cphi;
    double vz = vx1*cth - vx2*sth;

    // Then Cartesian -> cylindrical

    v[IDIR] =  vx*cphi + vy*sphi;   // v_r
    v[JDIR] = -vx*sphi + vy*cphi;   // v_phi
    v[KDIR] = vz;                   // v_z
  }

#elif GEOMETRY == SPHERICAL
  // Target: Spherical (r, th, phi)
  //         {x1, x2, x3} = {r, theta, phi}

  if (input_geom == CARTESIAN) {
    double th  = x2;
    double phi = x3;

    v[IDIR] =  vx1*sin(th)*cos(phi) + vx2*sin(th)*sin(phi) + vx3*cos(th);  // v_r
    v[JDIR] =  vx1*cos(th)*cos(phi) + vx2*cos(th)*sin(phi) - vx3*sin(th);  // v_theta
    v[KDIR] = -vx1*sin(phi) + vx2*cos(phi);                                // v_phi
  }

  if (input_geom == CYLINDRICAL) {
    double th  = x2;
    double phi = x3;

    // cylindrical -> Cartesian
    double vx = vx1*cos(phi) - vx2*sin(phi);
    double vy = vx1*sin(phi) + vx2*cos(phi);
    double vz = vx3;

    // Cartesian -> spherical
    v[IDIR] =  vx*sin(th)*cos(phi) + vy*sin(th)*sin(phi) + vz*cos(th);  // v_r
    v[JDIR] =  vx*cos(th)*cos(phi) + vy*cos(th)*sin(phi) - vz*sin(th);  // v_theta
    v[KDIR] = -vx*sin(phi) + vy*cos(phi);                               // v_phi
  }

  if (input_geom == SPHERICAL) {
    // Already spherical
    v[IDIR] = vx1;
    v[JDIR] = vx2;
    v[KDIR] = vx3;
  }
  
#endif
}

