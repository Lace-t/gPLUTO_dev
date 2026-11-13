#include "gpluto.hpp"
#include <cstring>
#if PCP_FOFC == YES
#include <cstdlib>
#include <cstdio>
#include <cassert>
#if !defined(_OPENACC) && !defined(_OPENMP)
#include <vector>
#include <algorithm>
#endif

namespace fofc_cache {
bool active = false;
int pass = 0;
Ary4D *rhs_target = nullptr;

Stage::Stage()
{
  assert(!active);
#ifdef _OPENACC
  const char *disable = "GPU_FOFC_CACHE_DISABLE";
#else
  const char *disable = "FOFC_CACHE_DISABLE";
#endif
  active = reusable() && std::getenv(disable) == nullptr;
  pass = 0;
}
Stage::~Stage() { active = false; }

bool stream_enabled()
{
#if defined(_OPENACC) || defined(_OPENMP)
  return false;
#else
  return active && std::getenv("FOFC_STREAM_DISABLE") == nullptr;
#endif
}

bool rk_enabled()
{
#ifdef STAGGERED_MHD
  // CT must materialize U before averaging/correcting the staggered field.
  return false;
#else
  return std::getenv("FOFC_RK_FUSE_DISABLE") == nullptr;
#endif
}

RHSTarget::RHSTarget(Ary4D &u) : previous(rhs_target)
{
#ifdef _OPENACC
  const bool enabled = reusable() && std::getenv("GPU_RHS_FUSE_DISABLE") == nullptr;
#else
  const bool enabled = active && std::getenv("FOFC_FUSE_DISABLE") == nullptr;
#endif
  rhs_target = enabled ? &u : nullptr;
}
RHSTarget::~RHSTarget() { rhs_target = previous; }

#ifdef _OPENACC
void transfer(Data *data, timeStep *Dts, Grid *grid, RBox *box, int dir, bool save)
{
  static double *storage[3] = {};
  const long total = (long)grid->ntot[0] * grid->ntot[1] * grid->ntot[2];
  const long count = total * (NVAR + 7);
  if (!storage[dir]) {
    storage[dir] = new double[count];
    std::fprintf(stderr, "GPUCache rank=%d dir=%d bytes=%ld\n",
                 g_rank, dir, count * (long)sizeof(double));
  }
  double *cache = storage[dir];
  int ni, nj, ib, ie, jb, je, kb, ke;
  ARRAY_STRIDES(grid, dir, ni, nj);
  if      (dir == IDIR) RBoxPermuteIndices<IDIR>(box, ib, ie, jb, je, kb, ke);
  else if (dir == JDIR) RBoxPermuteIndices<JDIR>(box, ib, ie, jb, je, kb, ke);
  else                 RBoxPermuteIndices<KDIR>(box, ib, ie, jb, je, kb, ke);

  // Preserve rounded center velocities over full reconstructed sweep rows.
#pragma acc parallel loop collapse(3) present(data, Dts, cache[0:count])
  for (int k = kb; k <= ke; ++k) {
    for (int j = jb; j <= je; ++j) {
      for (int i = 0; i < ni; ++i) {
        const long row = (long)ni * (j + (long)nj * k);
        const long q = row + i;
#pragma acc loop seq
        for (int n = 0; n < 3; ++n) {
          const long at = NVAR * row + (long)ni * (VX1 + n) + i;
          if (save) cache[(NVAR + 4 + n) * total + q] = data->sweep->vC[at];
          else data->sweep->vC[at] = cache[(NVAR + 4 + n) * total + q];
        }
      }
    }
  }

  // Save before LLF corrections; restore before applying the current mask.
#pragma acc parallel loop collapse(3) present(data, Dts, cache[0:count])
  for (int k = kb; k <= ke; ++k) {
    for (int j = jb; j <= je; ++j) {
      for (int i = ib; i <= ie + 1; ++i) {
        const long row = (long)ni * (j + (long)nj * k);
        const long q = row + i;
#pragma acc loop seq
        for (int n = 0; n < NVAR; ++n) {
          const long at = NVAR * row + (long)ni * n + i;
          if (save) cache[at] = data->sweep->flux[at];
          else data->sweep->flux[at] = cache[at];
        }
        if (save) {
          cache[NVAR * total + q]       = data->sweep->press[q];
          cache[(NVAR + 1) * total + q] = data->sweep->Bn[q];
          cache[(NVAR + 2) * total + q] = Dts->cmax[q];
          cache[(NVAR + 3) * total + q] = data->sweep->PSIn[q];
        } else {
          data->sweep->press[q] = cache[NVAR * total + q];
          data->sweep->Bn[q]    = cache[(NVAR + 1) * total + q];
          Dts->cmax[q]          = cache[(NVAR + 2) * total + q];
          data->sweep->PSIn[q]  = cache[(NVAR + 3) * total + q];
        }
      }
    }
  }
}
#endif

#if !defined(_OPENACC) && !defined(_OPENMP)
// PSIn is produced only by GLM solvers. Do not read uninitialized storage
// in EW/CT builds; Bn is needed by both Powell sources and CT's LLF faces.
#ifdef GLM_MHD
constexpr int auxiliary_fields=4;
#else
constexpr int auxiliary_fields=3;
#endif

#ifdef STAGGERED_MHD
// Snapshot the directional outputs of CT_Store_* before FOFC changes them.
// Restore before applying the current pass's LLF mask. Edge EMFs are rebuilt
// by CT_ComputeEMF, so caching them would retain the previous pass's mask.
static void transfer_ct(Data *d,Grid *g,RBox *box,int dir,bool save){
 static std::vector<double> cache[3];
#if CT_EMF_AVERAGE == CT_CONTACT
 static std::vector<signed char> sign[3];
#endif
 EMF *emf=d->emf;
 const int transverse=(dir+1)%3, binormal=(dir+2)%3;
 Ary3D *fields[]={&emf->e[binormal][dir],&emf->e[transverse][dir],
#if CT_EMF_AVERAGE == CT_FLUX
                 &emf->e_dff[binormal][dir],&emf->e_dff[transverse][dir],
#elif CT_EMF_AVERAGE == UCT_HLL
                 &emf->aL[dir],&emf->aR[dir],&emf->dL[dir],&emf->dR[dir],
#endif
 };
 constexpr size_t count=sizeof(fields)/sizeof(fields[0]);
 int ni,nj,ib,ie,jb,je,kb,ke;
 ARRAY_STRIDES(g,dir,ni,nj);
 if(dir==IDIR)RBoxPermuteIndices<IDIR>(box,ib,ie,jb,je,kb,ke);
 else if(dir==JDIR)RBoxPermuteIndices<JDIR>(box,ib,ie,jb,je,kb,ke);
 else RBoxPermuteIndices<KDIR>(box,ib,ie,jb,je,kb,ke);
 const size_t faces=(size_t)(ie-ib+2)*(je-jb+1)*(ke-kb+1);
 cache[dir].resize(count*faces);
#if CT_EMF_AVERAGE == CT_CONTACT
 sign[dir].resize(faces);
#endif
 size_t face=0;
 for(int k=kb;k<=ke;++k)for(int j=jb;j<=je;++j)for(int i=ib;i<=ie+1;++i,++face){
  const int x=dir==IDIR?i:j;
  const int y=dir==JDIR?i:(dir==IDIR?j:k);
  const int z=dir==KDIR?i:k;
  for(size_t n=0;n<count;++n){
   double &value=(*fields[n])(x,y,z);
   if(save)cache[dir][n*faces+face]=value;
   else value=cache[dir][n*faces+face];
  }
#if CT_EMF_AVERAGE == CT_CONTACT
  if(save)sign[dir][face]=emf->sv[dir](x,y,z);
  else emf->sv[dir](x,y,z)=sign[dir][face];
#endif
 }
}
#endif

static void transfer_legacy(Data *d,timeStep *dt,Grid *g,RBox *box,int dir,bool save){
 static std::vector<double> cache[3], center_velocity[3];
 const size_t total=(size_t)g->ntot[0]*g->ntot[1]*g->ntot[2];
 const size_t fields=NVAR+auxiliary_fields;
 if(cache[dir].size()!=total*fields)cache[dir].resize(total*fields);
 if(center_velocity[dir].size()!=total*3)center_velocity[dir].resize(total*3);
 int ni,nj,ib,ie,jb,je,kb,ke;
 ARRAY_STRIDES(g,dir,ni,nj);
 if(dir==IDIR)RBoxPermuteIndices<IDIR>(box,ib,ie,jb,je,kb,ke);
 else if(dir==JDIR)RBoxPermuteIndices<JDIR>(box,ib,ie,jb,je,kb,ke);
 else RBoxPermuteIndices<KDIR>(box,ib,ie,jb,je,kb,ke);
 // Reconstruction converts vC to four-velocity and back in place. Preserve
 // its rounded three-velocity as well: FOFC and RHS consume this buffer.
 for(int k=kb;k<=ke;++k)for(int j=jb;j<=je;++j)for(int i=0;i<ni;++i){
  size_t q=i+(size_t)ni*(j+(size_t)nj*k);
  for(int n=0;n<3;++n){
   double &v=d->sweep->vC[(size_t)NVAR*ni*(j+(size_t)nj*k)+i+(size_t)ni*(VX1+n)];
   if(save)center_velocity[dir][3*q+n]=v;else v=center_velocity[dir][3*q+n];
  }
 }
 for(int k=kb;k<=ke;++k)for(int j=jb;j<=je;++j)for(int i=ib;i<=ie+1;++i){
  size_t q=i+(size_t)ni*(j+(size_t)nj*k);
  double *c=cache[dir].data()+q*fields;
  for(int n=0;n<NVAR;++n){
   double &f=d->sweep->flux[(size_t)NVAR*ni*(j+(size_t)nj*k)+i+(size_t)ni*n];
   if(save)c[n]=f;else f=c[n];
  }
  double *values[]={d->sweep->press+q,d->sweep->Bn+q,dt->cmax+q,
#ifdef GLM_MHD
                    d->sweep->PSIn+q,
#endif
  };
  for(int n=0;n<auxiliary_fields;++n){if(save)c[NVAR+n]=*values[n];else *values[n]=c[NVAR+n];}
 }
}

// Cache matches each sweep's variable-major row layout. Only initialized
// faces and reconstructed central velocities are saved.
void transfer(Data *d,timeStep *dt,Grid *g,RBox *box,int dir,bool save){
#ifdef STAGGERED_MHD
 transfer_ct(d,g,box,dir,save);
#endif
 if(!stream_enabled()){transfer_legacy(d,dt,g,box,dir,save);return;}
 static std::vector<double> flux[3],aux[3],vel[3];
 size_t total=(size_t)g->ntot[0]*g->ntot[1]*g->ntot[2];
 flux[dir].resize(total*NVAR);aux[dir].resize(total*auxiliary_fields);vel[dir].resize(total*3);
 int ni,nj,ib,ie,jb,je,kb,ke;
 ARRAY_STRIDES(g,dir,ni,nj);
 if(dir==IDIR)RBoxPermuteIndices<IDIR>(box,ib,ie,jb,je,kb,ke);
 else if(dir==JDIR)RBoxPermuteIndices<JDIR>(box,ib,ie,jb,je,kb,ke);
 else RBoxPermuteIndices<KDIR>(box,ib,ie,jb,je,kb,ke);
 for(int k=kb;k<=ke;++k)for(int j=jb;j<=je;++j){
  size_t row=(size_t)ni*(j+(size_t)nj*k);
  for(int n=0;n<NVAR;++n){
   if(save){
    if(n>=VX1 && n<=VX3)
     std::copy_n(d->sweep->vC+NVAR*row+ni*n,ni,vel[dir].data()+3*row+ni*(n-VX1));
   }else{
    for(int i=0;i<ni;++i){
     int x=dir==IDIR?i:j;
     int y=dir==JDIR?i:(dir==IDIR?j:k);
     int z=dir==KDIR?i:k;
     d->sweep->vC[NVAR*row+ni*n+i]=(n>=VX1 && n<=VX3)
       ?vel[dir][3*row+ni*(n-VX1)+i]:d->Vc(x,y,z,n);
    }
   }
   double *a=d->sweep->flux+NVAR*row+ni*n+ib;
   double *b=flux[dir].data()+NVAR*row+ni*n+ib;
   if(save)std::copy_n(a,ie-ib+2,b);else std::copy_n(b,ie-ib+2,a);
  }
  if(!save)for(int i=0;i<ni;++i){
   int x=dir==IDIR?i:j;
   int y=dir==JDIR?i:(dir==IDIR?j:k);
   int z=dir==KDIR?i:k;
   d->sweep->flag[row+i]=d->flag(x,y,z);
  }
  double *a[]={d->sweep->press,d->sweep->Bn,dt->cmax,
#ifdef GLM_MHD
               d->sweep->PSIn,
#endif
  };
  for(int n=0;n<auxiliary_fields;++n){
   double *b=aux[dir].data()+n*total+row+ib;
   if(save)std::copy_n(a[n]+row+ib,ie-ib+2,b);
   else std::copy_n(b,ie-ib+2,a[n]+row+ib);
  }
 }
}

#elif defined(_OPENMP)
void transfer(Data *, timeStep *, Grid *, RBox *, int, bool) {}
#endif
} // namespace fofc_cache
#endif // PCP_FOFC

#pragma acc routine seq
int RMHD_PCP_Recover(double *u, double *v, RunConfig &run_config, int allow_fix)
{
  double work[NVAR], primitive[NVAR];
  for (int n = 0; n < NVAR; ++n) {
    work[n] = u[n];
    primitive[n] = 0.0;
    if (!isfinite(work[n])) {
      return 1;
    }
  }
  // EnergySolve recovers fluid variables; magnetic fields map directly.
  primitive[BX1] = work[BX1];
  primitive[BX2] = work[BX2];
  primitive[BX3] = work[BX3];

  int admissible = RMHD_PCP_Check(work);
  uint32_t err = 0;
  if (!admissible) {
    if (!allow_fix || work[RHO] <= 0.0) {
      return 1;
    }
    primitive[BX1] = work[BX1];
    primitive[BX2] = work[BX2];
    primitive[BX3] = work[BX3];
    err = RMHD_PressureFix(work, primitive, run_config);
  } else {
    State state;
    STATE_SET(state, primitive, work, NULL, 0.0, 0.0, 0.0);
    uint32_t flag = 0;
    state.flag = &flag;
    if (allow_fix) {
      err = ConsToPrim(state, run_config);
    } else {
      // FOFC and limiter trials must not commit a pressure-floor repair.
      err = RMHD_EnergySolve(work, primitive, run_config);
#if PCP_NR == YES
      if (err && !RMHD_PCP_NR(work, primitive, run_config)) {
        err = 0;
      }
#endif
      if (err) {
        return 1;
      }
    }
  }
  if (err & (ERR_CONS2PRIM_FAIL | ERR_PRESSURE_FIX)) {
    return 1;
  }
  if (!allow_fix) {
    // The strict path cannot alter U: identical input retains its admission.
    // Reuse the entrance criterion instead of evaluating it a second time.
    for (int n = 0; n < NVAR; ++n) {
      if (work[n] != u[n]) {
        return 1;
      }
    }
  } else if (!RMHD_PCP_Check(work)) {
    // PressureFix can change energy, so repaired states need a fresh check.
    return 1;
  }
#if NSCL > 0
  for (int n = NFLX; n < NVAR; ++n) {
    primitive[n] = work[n] / work[RHO];
  }
#endif
#ifdef GLM_MHD
  primitive[PSI_GLM] = work[PSI_GLM];
#endif
  double v2 = primitive[VX1] * primitive[VX1] + primitive[VX2] * primitive[VX2] +
              primitive[VX3] * primitive[VX3];
  if (!(v2 < 1.0) || primitive[RHO] <= 0.0 || primitive[PRS] <= 0.0) {
    return 1;
  }
  for (int n = 0; n < NVAR; ++n) {
    if (!isfinite(primitive[n])) {
      return 1;
    }
  }
  for (int n = 0; n < NVAR; ++n) {
    u[n] = work[n];
    v[n] = primitive[n];
  }
  return 0;
}

// The face at sweep index i lies between cells i-1 and i.
void RMHD_PCP_LimitFaces(Data *data, Grid *grid, RBox *box, int dir)
{
  int ib, ie, jb, je, kb, ke, ni, nj;
  ARRAY_STRIDES(grid, dir, ni, nj);
  if (dir == IDIR) {
    RBoxPermuteIndices<IDIR>(box, ib, ie, jb, je, kb, ke);
  } else if (dir == JDIR) {
    RBoxPermuteIndices<JDIR>(box, ib, ie, jb, je, kb, ke);
  } else {
    RBoxPermuteIndices<KDIR>(box, ib, ie, jb, je, kb, ke);
  }
  RunConfig &run_config = *data->run_config;
#pragma acc parallel loop collapse(3) present(data, grid, run_config)
  for (int k = kb; k <= ke; ++k) {
    for (int j = jb; j <= je; ++j) {
      for (int i = ib; i <= ie + 1; ++i) {
        long off = ni * (j + nj * k);
        Ary2D vc(&data->sweep->vC[NVAR * off], ni, NVAR);
        Ary2D vL(&data->sweep->vL[NVAR * off], ni, NVAR);
        Ary2D vR(&data->sweep->vR[NVAR * off], ni, NVAR);
        double trial_v[2][NVAR];
        bool changed[2] = {false, false};
        bool failed = false;
        for (int side = 0; side < 2; ++side) {
          double face_u[NVAR], center_u[NVAR], center_v[NVAR];
          for (int n = 0; n < NVAR; ++n) {
            trial_v[side][n] = side == 0 ? vL(i, n) : vR(i, n);
          }
          PrimToCons(trial_v[side], face_u, run_config);
          if (RMHD_PCP_Check(face_u)) {
            continue;
          }
          for (int n = 0; n < NVAR; ++n) {
            center_v[n] = vc(i - 1 + side, n);
          }
          PrimToCons(center_v, center_u, run_config);
          int limited = RMHD_PCP_Limit(center_u, face_u);
          if (limited < 0 || RMHD_PCP_Recover(face_u, trial_v[side], run_config, 0)) {
            failed = true;
            break;
          }
          changed[side] = (limited > 0);
        }
        // Keep interface updates atomic: discard both trials if either fails.
        if (!failed) {
          if (changed[0]) {
            for (int n = 0; n < NVAR; ++n) {
              vL(i, n) = trial_v[0][n];
            }
          }
          if (changed[1]) {
            for (int n = 0; n < NVAR; ++n) {
              vR(i, n) = trial_v[1][n];
            }
          }
        }
      }
    }
  }
}

#if (PCP_NR == YES || PCP_FOFC == YES || FALLBACK_MODE == FALLBACK_RECOMPUTE || CT_RMHD_CONSISTENT == YES)
#if PCP_FOFC == YES
// Device snapshots include the upper face in each component's own direction.
struct PCPSnapshot {
  Ary4D u, v;
  Ary3D b[6];
  void create(Grid *g)
  {
    int nx = g->ntot[0], ny = g->ntot[1], nz = g->ntot[2];
    u.create(nx, ny, nz, NVAR);
    v.create(nx, ny, nz, NVAR);
#ifdef STAGGERED_MHD
    for (int n = 0; n < DIMENSIONS; ++n) {
      b[n].create(nx + (n == 0), ny + (n == 1), nz + (n == 2));
    }
#endif
#pragma acc enter data copyin(this[ : 1])
  }
  void transfer(Data *d, Grid *g, int save)
  {
    int nx = g->ntot[0], ny = g->ntot[1], nz = g->ntot[2];
#pragma acc parallel loop collapse(3) present(d, u, v)
    for (int k = 0; k < nz; ++k) {
      for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
          for (int n = 0; n < NVAR; ++n) {
            if (save) {
              u(i, j, k, n) = d->Uc(i, j, k, n);
              v(i, j, k, n) = d->Vc(i, j, k, n);
            } else {
              d->Uc(i, j, k, n) = u(i, j, k, n);
              d->Vc(i, j, k, n) = v(i, j, k, n);
            }
          }
        }
      }
    }
#ifdef STAGGERED_MHD
    for (int n = 0; n < DIMENSIONS; ++n) {
      int ex = nx + (n == 0), ey = ny + (n == 1), ez = nz + (n == 2);
#pragma acc parallel loop collapse(3) present(d, b)
      for (int k = 0; k < ez; ++k) {
        for (int j = 0; j < ey; ++j) {
          for (int i = 0; i < ex; ++i) {
            if (save) {
              b[n](i, j, k) = d->Vs[n](i, j, k);
            } else {
              d->Vs[n](i, j, k) = b[n](i, j, k);
            }
          }
        }
      }
    }
#endif
  }
};
#endif

// Only the CPU backend owns candidate storage. Recovery callers share this
// interface; device and non-FOFC builds recover directly without a cache.
namespace candidate_cache {
#if PCP_FOFC == YES && !defined(_OPENACC) && !defined(_OPENMP)
static std::vector<double> input, primitive;
static std::vector<unsigned char> failed, ct_failed;

inline bool enabled()
{
  return fofc_cache::active && std::getenv("FOFC_CANDIDATE_DISABLE") == nullptr;
}

inline void prepare(Grid *grid)
{
  const size_t cells = (size_t)grid->ntot[0] * grid->ntot[1] * grid->ntot[2];
  input.resize(cells * NVAR);
  primitive.resize(cells * NVAR);
  failed.resize(cells);
  ct_failed.resize(cells);
}

inline int recover(double *u, double *v, RunConfig &config, int allow_fix,
                   size_t cell, bool ct_failure, bool use_cache)
{
  bool same = use_cache && fofc_cache::pass > 0 && ct_failed[cell] == ct_failure;
  if (same) {
    for (int n = 0; n < NVAR; ++n) {
      // Preserve signed zero and every other bit of the orthonormal input.
      if (std::memcmp(&u[n], &input[cell * NVAR + n], sizeof(double))) {
        same = false;
        break;
      }
    }
  }
  if (same) {
    if (!failed[cell]) {
      for (int n = 0; n < NVAR; ++n) v[n] = primitive[cell * NVAR + n];
    }
    return failed[cell];
  }
  if (use_cache) {
    for (int n = 0; n < NVAR; ++n) input[cell * NVAR + n] = u[n];
  }
  const int err = ct_failure ? 1 : RMHD_PCP_Recover(u, v, config, allow_fix);
  if (use_cache) {
    failed[cell] = (err != 0);
    ct_failed[cell] = ct_failure;
    if (!err) {
      for (int n = 0; n < NVAR; ++n) primitive[cell * NVAR + n] = v[n];
    }
  }
  return err;
}
#else
inline bool enabled() { return false; }
inline void prepare(Grid *) {}

#pragma acc routine seq
inline int recover(double *u, double *v, RunConfig &config, int allow_fix,
                   size_t, bool ct_failure, bool)
{
  return ct_failure ? 1 : RMHD_PCP_Recover(u, v, config, allow_fix);
}
#endif
} // namespace candidate_cache

// Both entry points use the same metric conversion, recovery and commit path.
// Only the fused variant constructs U in the loop instead of reading it.
template<bool fused>
static int PCPRecoverStageImpl(Data *data, Grid *grid, RBox &box,
                               RunConfig &run_config, int allow_pressure_fix,
                               const Ary4D &Ubase, const Ary4D &stage_u,
                               const Ary4D &rhs, double initial_weight,
                               double stage_weight)
{
  int i, j, k;
  int bad = 0;
  const bool use_candidate = candidate_cache::enabled() && !allow_pressure_fix;

#pragma acc parallel loop collapse(3) present(data, grid, run_config, box, Ubase, stage_u, rhs) reduction(+ : bad)
  BOX_LOOP (box, i, j, k) {
    double u[NVAR], v[NVAR];
    for (int n = 0; n < NVAR; ++n) {
      if constexpr (fused) {
        u[n] = initial_weight * Ubase(i, j, k, n) +
               stage_weight * (stage_u(i, j, k, n) + rhs(i, j, k, n));
      } else {
        u[n] = data->Uc(i, j, k, n);
      }
    }
#if GEOMETRY != CARTESIAN
    Metric metric;
    MetricSet(metric, grid->rC[i], grid->sC[j]);
    ConsToOrtho(u, metric);
#endif
    const size_t cell = i + (size_t)grid->ntot[0] * (j + (size_t)grid->ntot[1] * k);
    const bool ct_failure = (data->flag(i, j, k) & FLAG_PCP_CT_FAIL) != 0;
    const int err = candidate_cache::recover(u, v, run_config, allow_pressure_fix,
                                             cell, ct_failure, use_candidate);
    if (err) {
      if constexpr (fused) {
        // Keep the original coordinate representation on failure. Neither a
        // metric round trip nor a failed recovery may alter the candidate U.
        for (int n = 0; n < NVAR; ++n) {
          data->Uc(i, j, k, n) = initial_weight * Ubase(i, j, k, n) +
              stage_weight * (stage_u(i, j, k, n) + rhs(i, j, k, n));
        }
      }
      data->flag(i, j, k) |= FLAG_PCP_FOFC | FLAG_FALLBACK;
      ++bad;
    } else {
      data->flag(i, j, k) &= ~FLAG_FALLBACK;
#if GEOMETRY != CARTESIAN
      OrthoToCons(u, metric);
#endif
      for (int n = 0; n < NVAR; ++n) {
        data->Uc(i, j, k, n) = u[n];
        data->Vc(i, j, k, n) = v[n];
      }
    }
  }
#ifdef PARALLEL
  MPI_Allreduce(MPI_IN_PLACE, &bad, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
#endif
  return bad;
}

static int PCPRecoverStage(Data *data, Grid *grid, RBox &box,
                           RunConfig &config, int allow_fix)
{
  // The unused inputs disappear from the non-fused template instantiation.
  return PCPRecoverStageImpl<false>(data, grid, box, config, allow_fix,
                                    data->Uc, data->Uc, data->Uc, 0.0, 0.0);
}

#if PCP_FOFC == YES
static int PCPRecoverStageFused(Data *data, Grid *grid, RBox &box,
                                RunConfig &config, const Ary4D &Ubase,
                                const Ary4D &stage_u, const Ary4D &rhs,
                                double initial_weight, double stage_weight)
{
  return PCPRecoverStageImpl<true>(data, grid, box, config, 0, Ubase, stage_u,
                                   rhs, initial_weight, stage_weight);
}
#endif

void RMHD_PCP_SetStageTime(Data *data, timeStep *Dts, int stage_number)
{
  RunConfig &run_config = *data->run_config;
  double offset = stage_number == 2 ? Dts->dt
                  : stage_number == 3 ? 0.5 * Dts->dt
                                      : 0.0;
  run_config.time = Dts->time + offset;
#pragma acc update device(run_config.time)
}

uint32_t RMHD_PCP_FinishStage(Data *data, Grid *grid, RBox *box)
{
  int bad = PCPRecoverStage(data, grid, *box, *data->run_config, 0);
#if PCP_PRESSURE_FIX == YES
  if (bad) {
    bad = PCPRecoverStage(data, grid, *box, *data->run_config, 1);
    //printLog("> PCP PressureFix: stage=%d rejected=%d\n", g_intStage, bad);
  }
#endif
  if (bad) {
#if FALLBACK_MODE == FALLBACK_RECOMPUTE
    return ERR_CONS2PRIM_FAIL;
#else
    QuitPLUTO("ERR_CONS2PRIM_FAIL", 1);
#endif
  }
  return 0;
}

#if PCP_FOFC == YES
// The RK driver owns boundaries, particles, radiation and source scheduling.
// Ubase/Bbase are its initial RK state; rhs/rs are its reusable work arrays.
uint32_t RMHD_PCP_RunStage(Data *data, timeStep *Dts, double dt,
                           const Ary4D &Ubase, Ary3D (&Bbase)[6],
                           Ary4D &rhs, Ary3D (&rs)[6], double initial_weight,
                           double stage_weight, RBox *box, Grid *grid)
{
  fofc_cache::Stage cache_stage;
  if (candidate_cache::enabled()) candidate_cache::prepare(grid);
  static PCPSnapshot stage;
  static bool allocated = false;
  if (!allocated) {
    stage.create(grid);
    allocated = true;
  }
  stage.transfer(data, grid, 1);
  RunConfig &run_config = *data->run_config;
  int nx = grid->ntot[0], ny = grid->ntot[1], nz = grid->ntot[2];
  const int step_retry = (data->integration_mode == FALLBACK_RECOMPUTE);
  int i, j, k;
  if (g_intStage == 1) {
#pragma acc parallel loop collapse(3) present(data, Ubase, box)
    BOX_LOOP ((*box), i, j, k) {
      for (int n = 0; n < NVAR; ++n) {
        Ubase(i, j, k, n) = data->Uc(i, j, k, n);
      }
    }
#ifdef STAGGERED_MHD
    for (int n = 0; n < DIMENSIONS; ++n) {
      int ex = box->iend + (n == 0);
      int ey = box->jend + (n == 1);
      int ez = box->kend + (n == 2);
#pragma acc parallel loop collapse(3) present(data, Bbase, box)
      for (int k = box->kbeg; k <= ez; ++k) {
        for (int j = box->jbeg; j <= ey; ++j) {
          for (int i = box->ibeg; i <= ex; ++i) {
            Bbase[n](i, j, k) = data->Vs[n](i, j, k);
          }
        }
      }
    }
#endif
  }
#pragma acc parallel loop collapse(3) present(data)
  for (int k = 0; k < nz; ++k) {
    for (int j = 0; j < ny; ++j) {
      for (int i = 0; i < nx; ++i) {
        data->flag(i, j, k) &= ~(FLAG_PCP_FOFC | FLAG_PCP_CT_FAIL);
        // Whole-step fallback selects only troubled roots, not every face.
        if (step_retry > 0 && (data->flag(i, j, k) & FLAG_HLL)) {
          data->flag(i, j, k) |= FLAG_PCP_FOFC;
        }
      }
    }
  }
  const int local_passes = PCP_STAGE_LOCAL_PASSES;
  int bad = 0;
  for (int pass = 0; pass < local_passes; ++pass) {
    fofc_cache::pass = pass;
    if (pass > 0) {
      stage.transfer(data, grid, 0);
    }
    // Keep COPY on later RK stages for other flag consumers; retain the
    // established first-stage optimization after the initial FlagShock.
#if SHOCK_FLATTENING == MULTID && ENTROPY_SWITCH == NO
    if (pass > 0 || g_intStage != 1) {
      FlagBoundary(data->flag, FLAG_COPY, grid);
    }
#else
    FlagBoundary(data->flag, FLAG_COPY, grid);
#endif
    const bool pcp_flux_possible = (step_retry > 0 || pass > 0);
    UpdateStage(data, rhs, rs, grid, Dts, box, dt, pcp_flux_possible);
    const bool fuse_rk = fofc_cache::rk_enabled();
    if(!fuse_rk){
#pragma acc parallel loop collapse(3) present(data, rhs, box, Ubase, stage.u)
    BOX_LOOP ((*box), i, j, k) {
      for (int n = 0; n < NVAR; ++n) {
        data->Uc(i, j, k, n) = initial_weight * Ubase(i, j, k, n) +
                               stage_weight * (stage.u(i, j, k, n) + rhs(i, j, k, n));
      }
    }
    }
#ifdef STAGGERED_MHD
    for (int n = 0; n < DIMENSIONS; ++n) {
      int ex = box->iend + (n == 0);
      int ey = box->jend + (n == 1);
      int ez = box->kend + (n == 2);
#pragma acc parallel loop collapse(3) present(data, box, rs, Bbase, stage.b)
      for (int k = box->kbeg; k <= ez; ++k) {
        for (int j = box->jbeg; j <= ey; ++j) {
          for (int i = box->ibeg; i <= ex; ++i) {
            data->Vs[n](i, j, k) = initial_weight * Bbase[n](i, j, k) +
                                   stage_weight * (stage.b[n](i, j, k) + rs[n](i, j, k));
          }
        }
      }
    }
#pragma acc parallel loop collapse(3) present(data, box)
    BOX_LOOP ((*box), i, j, k) {
      data->flag(i, j, k) &= ~FLAG_PCP_CT_FAIL;
    }
    CT_AverageStaggeredFields(data, 1, box, grid);
#endif
    if (fuse_rk) {
      bad = PCPRecoverStageFused(data, grid, *box, run_config, Ubase, stage.u, rhs,
                                 initial_weight, stage_weight);
    } else {
      bad = PCPRecoverStage(data, grid, *box, run_config, 0);
    }
    if (!bad) {
      return 0;
    }
    // printLog("> PCP: stage=%d pass=%d retry=%d rejected=%d dt=%.16e\n",
    //          g_intStage, pass, step_retry, bad, dt);
  }
#if PCP_PRESSURE_FIX == YES
  // Exhaust FOFC before pressure repair; a remaining failure goes to main.
  bad = PCPRecoverStage(data, grid, *box, run_config, 1);
  // printLog("> PCP PressureFix: stage=%d retry=%d rejected=%d\n",
  //          g_intStage, step_retry, bad);
#endif
  if (bad) {
#if FALLBACK_MODE == FALLBACK_RECOMPUTE
    return ERR_CONS2PRIM_FAIL;
#else
    QuitPLUTO("ERR_CONS2PRIM_FAIL", 1);
#endif
  }
  return 0;
}
#endif // PCP_FOFC
#endif // stage checks / recovery
