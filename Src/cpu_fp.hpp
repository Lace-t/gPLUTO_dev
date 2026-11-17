#pragma once

// CPU startup policy. FTZ flushes underflowing results; DAZ treats subnormal
// inputs as zero. This deliberately trades gradual underflow for throughput.
// GPU arithmetic is unaffected. Set GPLUTO_FTZ_DAZ=0 for an IEEE control run.
#if defined(__x86_64__) && !defined(_OPENACC) && !defined(_OPENMP)
#include <immintrin.h>
#include <cstdlib>
namespace CpuFp {
inline unsigned configure() {
  constexpr unsigned mask = 0x8040u; // MXCSR: FTZ (bit 15), DAZ (bit 6).
  const char *mode = std::getenv("GPLUTO_FTZ_DAZ");
  const bool enabled = !(mode && mode[0] == '0' && mode[1] == '\0');
  const unsigned csr = _mm_getcsr();
  _mm_setcsr(enabled ? csr | mask : csr & ~mask);
  return _mm_getcsr();
}
}
#else
namespace CpuFp { inline unsigned configure() { return 0; } }
#endif
