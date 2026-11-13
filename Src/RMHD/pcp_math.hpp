#ifndef GPLUTO_PCP_MATH_HPP
#define GPLUTO_PCP_MATH_HPP

#include <cmath>
// The device-intrinsic workaround is needed from NVHPC 26.3 onward.
// Older SDKs use std::fma and avoid cudadevice.h conflicts with CUDA half types.
#if defined(__NVCOMPILER) && defined(_OPENACC) && \
    (__NVCOMPILER_MAJOR__ * 100 + __NVCOMPILER_MINOR__ >= 2603)
#if __has_include(<cudadevice.h>) && __has_include(<nv/target>)
#define GPLUTO_PCP_NV_FMA 1
#include <cudadevice.h>
#include <nv/target>
#endif
#endif

#pragma omp declare target
#pragma acc routine seq
#ifdef GPLUTO_PCP_NV_FMA
// NVHPC 26.3 rejects the DFMA lowered from std::fma in OpenACC routines.
// Keep host/device selection until device compilation: early host inlining
// would reintroduce DFMA before the compiler discards the host branch.
__attribute__((noinline))
#endif
static double PCPFma(double a, double b, double c)
{
#ifdef GPLUTO_PCP_NV_FMA
  if target (nv::target::is_device) {
    return __fma_rn(a, b, c);
  } else {
    return std::fma(a, b, c);
  }
#else
  return std::fma(a, b, c);
#endif
}
#pragma omp end declare target

#ifdef GPLUTO_PCP_NV_FMA
#undef GPLUTO_PCP_NV_FMA
#endif

#endif
