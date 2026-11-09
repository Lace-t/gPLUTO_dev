// Standalone regression: compile with g++ or nvc++ -acc, without -cuda.
#include "../../../Src/RMHD/pcp_math.hpp"
#include <cstdio>
#include <cstring>
#include <cstdint>
#include <limits>
#ifdef _OPENACC
#include <openacc.h>
#endif

int main()
{
  const double huge = std::numeric_limits<double>::max();
  double inputs[][3] = {
    {1.0 + 0x1p-27, 1.0 - 0x1p-27, -1.0},
    {1.0 + 0x1p-52, 1.0 - 0x1p-52, -1.0},
    {huge, 2.0, -huge}, // Intermediate overflow must not overflow FMA.
    {0x1p-600, 0x1p-600, 0x1p-1000}, // Underflowing product, normal sum.
    {0.0, 1.0, -0.0}, {-0.0, 1.0, -0.0},
    {0x1p-500, 0x1p-500, -0x1p-1000},
    {0x1p500, 0x1p500, -0x1p1000}
  };
  constexpr int count = sizeof(inputs) / sizeof(inputs[0]);
  double expected[count], host[count], device[count];
  for (int n = 0; n < count; ++n) {
    expected[n] = std::fma(inputs[n][0], inputs[n][1], inputs[n][2]);
    host[n] = PCPFma(inputs[n][0], inputs[n][1], inputs[n][2]);
  }
  int on_gpu = 0;
#pragma acc serial copyin(inputs) copyout(device, on_gpu)
  {
#ifdef _OPENACC
    on_gpu = acc_on_device(acc_device_nvidia);
#endif
    for (int n = 0; n < count; ++n) {
      device[n] = PCPFma(inputs[n][0], inputs[n][1], inputs[n][2]);
    }
  }
  int failures = 0;
  for (int n = 0; n < count; ++n) {
    if (std::memcmp(&host[n], &expected[n], sizeof(double)) ||
        std::memcmp(&device[n], &expected[n], sizeof(double))) {
      std::printf("FMA case %d failed: host=%a device=%a expected=%a\n",
                  n, host[n], device[n], expected[n]);
      ++failures;
    }
  }
#ifdef _OPENACC
  // Catch NVHPC's silent host fallback when a GPU region rejects DFMA.
  if (acc_get_device_type() == acc_device_nvidia && !on_gpu) {
    std::printf("FMA kernel fell back to host execution\n");
    ++failures;
  }
#endif
  std::printf("PCP FMA: %d cases, %d failures, NVIDIA GPU=%d\n",
              count, failures, on_gpu);
  return failures ? 1 : 0;
}
