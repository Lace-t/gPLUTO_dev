#pragma once

#include <cmath>

#if defined(__GNUC__) && !defined(_OPENACC) && !defined(_OPENMP)
#define ENTROPY_INLINE inline __attribute__((always_inline))
#else
#define ENTROPY_INLINE inline
#endif

// Local orthonormal fluid entropy test shared by HLLC and HLLD.
// EW includes the full Powell correction. DC tests the fluid subproblem after
// GLM_Solve; CT tests the subproblem with a common staggered normal field.
// The excluded DC/CT flux components are preserved exactly.
//
// Switch to same-state HLL only when the candidate clearly violates the entropy
// inequality and HLL clearly satisfies it. Geometry factors and sources enter
// later in RHS: this test does not establish full reconstruction/CT/RK entropy
// stability. Requires ideal EOS and the reduced-energy convention.
#pragma omp declare target
namespace EntropyGuard {

// Components participating in the fluid entropy pair. In DC, the first eight
// entries are RMHD variables; the GLM scalar and normal B belong to cleaning.
#pragma acc routine seq
ENTROPY_INLINE bool fluid_component(int n, int dir)
{
  #if DIVB_CONTROL == DIV_CLEANING
  return n < 8 && n != BX1 + dir;
  #elif DIVB_CONTROL == CONSTRAINED_TRANSPORT
  return n < NFLX && n != BX1 + dir;
  #else
  return n < NFLX;
  #endif
}

// Finite zero-field states have zero magnetic entropy variables. Keep all
// excluded flux components unchanged; diagnostics always use the general path.
#pragma acc routine seq
template<bool hydro>
ENTROPY_INLINE bool active_component(int n, int dir)
{
  return fluid_component(n, dir) && (!hydro || n < BX1 || n > BX3);
}

struct CompactPair {
 double jump[NFLX]{}, offset=0., base_scale=0.;
};
#pragma acc routine seq
template<int dir>
ENTROPY_INLINE void compact_pair(const double *vl, const double *vr,
 const double *ul, const double *ur, double bn, double pl, double pr,
 const RunConfig &rc, CompactPair &c)
{
 #if DIVB_CONTROL != EIGHT_WAVES
 bn=vl[BX1+dir];
 #endif
 const double al=ul[RHO]/vl[PRS],ar=ur[RHO]/vr[PRS];
 const double betal=vl[RHO]/vl[PRS],betar=vr[RHO]/vr[PRS];
 double v2l=0.,v2r=0.,vbl=0.,vbr=0.;
 for(int j=0;j<3;++j){v2l+=vl[VX1+j]*vl[VX1+j];v2r+=vr[VX1+j]*vr[VX1+j];vbl+=vl[VX1+j]*vl[BX1+j];vbr+=vr[VX1+j]*vr[BX1+j];}
 const double rp=vr[PRS]/vl[PRS],rr=vr[RHO]/vl[RHO];
 const double lp=std::fabs(rp-1.)<.125?std::log1p((vr[PRS]-vl[PRS])/vl[PRS]):std::log(rp);
 const double lr=std::fabs(rr-1.)<.125?std::log1p((vr[RHO]-vl[RHO])/vl[RHO]):std::log(rr);
 c.jump[RHO]=betar-betal-(ar-al)-(lp-rc.gamma*lr)/(rc.gamma-1.);
 c.jump[ENG]=-(ar-al);
 const double kl=al*pl*vl[VX1+dir],kr=ar*pr*vr[VX1+dir];
 const double phil=al*vbl,phir=ar*vbr;
 c.offset=-(kr-kl)+bn*(phir-phil);
 c.base_scale=std::fabs(kl)+std::fabs(kr)+std::fabs(bn)*(std::fabs(phil)+std::fabs(phir));
 for(int j=0;j<3;++j){
  c.jump[MX1+j]=ar*vr[VX1+j]-al*vl[VX1+j];
  if(fluid_component(BX1+j,dir)) c.jump[BX1+j]=ar*((1.-v2r)*vr[BX1+j]+vbr*vr[VX1+j])-al*((1.-v2l)*vl[BX1+j]+vbl*vl[VX1+j]);
 }

}
#pragma acc routine seq
ENTROPY_INLINE double compact_residual(const CompactPair &c,const double *flux,double &scale){
 double r=c.offset;scale=c.base_scale;
 for(int n=0;n<NFLX;++n){const double t=c.jump[n]*flux[n];r+=t;scale+=std::fabs(t);}
 return r;
}


#pragma acc routine seq
template<int dir>
ENTROPY_INLINE bool apply_flux(const double *vl,const double *vr,
 const double *ul,const double *ur,const double *fl,const double *fr,
 double sl,double sr,double bn,double *flux,const RunConfig &rc,double pl,double pr)
{
 // DC/CT require a common normal field; incompatible input has no fluid entropy pair.
 #if DIVB_CONTROL != EIGHT_WAVES
 if(vl[BX1+dir]!=vr[BX1+dir])return false;
 #if DIVB_CONTROL == DIV_CLEANING
 if(vl[PSI_GLM]!=vr[PSI_GLM])return false;
 #endif
 #endif
 CompactPair c;
 compact_pair<dir>(vl,vr,ul,ur,bn,pl,pr,rc,c);
 double scale;const double r=compact_residual(c,flux,scale);
 constexpr double eps=512.*2.220446049250313080847263336181640625e-16;
 // The same compact residual is used throughout, including near cancellation.
 if(!std::isfinite(r) || !(r>eps*scale+1.e-25))return false;
 const double bl=MIN(sl,0.),br=MAX(sr,0.),inv=1./(br-bl);
 double hll[NFLX];
 for(int n=0;n<NFLX;++n)hll[n]=fluid_component(n,dir)?
  (bl*br*(ur[n]-ul[n])+br*fl[n]-bl*fr[n])*inv:flux[n];
 double hs;const double rh=compact_residual(c,hll,hs);
 const double tol=eps*(scale+hs-c.base_scale)+1.e-25;
 if(!std::isfinite(rh) || !(r>tol && rh < -tol))return false;
 for(int n=0;n<NFLX;++n)flux[n]=hll[n];
 return true;
}
ENTROPY_INLINE double total_pressure(const double *v) {
 double b2=0.,v2=0.,vb=0.;for(int j=0;j<3;++j){b2+=v[BX1+j]*v[BX1+j];v2+=v[VX1+j]*v[VX1+j];vb+=v[BX1+j]*v[VX1+j];}
 return v[PRS]+.5*(b2*(1.-v2)+vb*vb);
}
#pragma acc routine seq
ENTROPY_INLINE bool apply_flux(const double *vl,const double *vr,
 const double *ul,const double *ur,const double *fl,const double *fr,
 double sl,double sr,double bn,double *flux,int dir,const RunConfig &rc)
{
 const double pl=total_pressure(vl),pr=total_pressure(vr);
 if(dir==0)return apply_flux<0>(vl,vr,ul,ur,fl,fr,sl,sr,bn,flux,rc,pl,pr);
 if(dir==1)return apply_flux<1>(vl,vr,ul,ur,fl,fr,sl,sr,bn,flux,rc,pl,pr);
 return apply_flux<2>(vl,vr,ul,ur,fl,fr,sl,sr,bn,flux,rc,pl,pr);
}
} // namespace EntropyGuard
#pragma omp end declare target
#undef ENTROPY_INLINE
