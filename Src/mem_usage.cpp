#include "gpluto.hpp"
#include <stdint.h>
#include <stdio.h>
#if _OPENACC
 #include <openacc.h>
#endif

// # define ACC_MEASURE_MEM
#ifdef ACC_MEASURE_MEM
int64_t max_mem_usage = 0;

void measure_acc_mem_usage() {
    auto dev_ty = acc_get_device_type();
    auto dev_mem = acc_get_property(0, dev_ty, acc_property_memory);
    auto dev_free_mem = acc_get_property(0, dev_ty, acc_property_free_memory);
    auto mem = dev_mem - dev_free_mem;
    if (mem > max_mem_usage)
        max_mem_usage = mem;
}
void print_mem_usage() {
    printf("mem_usage %f GB\n", (double)max_mem_usage * 1e-9);
}
#else
// #else
void measure_acc_mem_usage() {}
void print_mem_usage() {}
#endif
// #endif
