# PCP regression tests

`pcp_math.hpp` provides fused, double-precision multiply-add without changing
solver equations. NVHPC OpenACC uses `__fma_rn` only on the device; host calls,
CPU builds and other offload compilers use `std::fma`. The NVHPC workaround is
selected by compiler and header availability, without SDK-version or GPU-model
checks. Older SDKs without the target-selection headers retain the standard
FMA implementation. No `-cuda` flag or case-local source override is required.

The NVHPC helper must retain its `noinline` attribute: otherwise early host
inlining can reintroduce an unsupported DFMA before device selection, causing
missing device routines or silent host fallback. Both implementations retain
fused rounding; replacing FMA with `a*b+c` loses cancellation accuracy.

Run the standalone regression from the gPLUTO root:

```sh
g++ -std=c++17 -O3 Test_Problems/RMHD/PCP/test_fma.cpp -o /tmp/pcp-fma-cpu
/tmp/pcp-fma-cpu
nvc++ -std=c++17 -acc -O3 -Minline=1000 Test_Problems/RMHD/PCP/test_fma.cpp -o /tmp/pcp-fma-gpu
/tmp/pcp-fma-gpu
```

The test compares bit patterns with standard FMA, exercises cancellation,
intermediate overflow/underflow and signed zero, and rejects silent host
fallback when NVIDIA GPU execution is selected. Subnormal host arithmetic can
be affected by the compiler/runtime's flush-to-zero mode; the cases do not
assume a different floating-point environment from the application's default.

To test a different installed SDK, use its `nvc++` on PATH. To check other GPU
architectures, add that SDK's supported `-gpu=ccXX,...` targets at both compile
and link time. Cross-compilation verifies code generation and linking, while
execution requires matching hardware or supported PTX JIT.

Run solver, blast-wave, stage-rejection and rollback regressions in a fresh
output directory (the runner also executes `test_fma`):

```sh
python3 Test_Problems/RMHD/PCP/run_tests.py /tmp/pcp-regression-new serial nvc
```

Optional modes are `mpi` and `mpi-nvc`; load the corresponding MPI/compiler
environment first. Existing output directories are not overwritten.
