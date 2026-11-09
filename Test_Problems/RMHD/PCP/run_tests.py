"""Build and exercise PCP kernels, stage rejection and rollback.

Usage: python run_tests.py /absolute/new/output [serial mpi nvc mpi-nvc]
Load the NVHPC environment (nvhpc on this workstation) before GPU builds.
Uses only Python's standard library. Existing results are never deleted.
"""
from pathlib import Path
import json
import os
import shutil
import struct
import subprocess
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
OUT = Path(sys.argv[1]).resolve()
MODES = sys.argv[2:] or ['serial', 'mpi', 'nvc', 'mpi-nvc']
ARCH = {'serial':'linux.g++.defs', 'mpi':'linux.mpicxx.defs',
        'nvc':'linux.nvc++.acc.defs', 'mpi-nvc':'linux.mpicxx.acc.defs'}
results = {}

def command(args, cwd, logfile):
    with (cwd/logfile).open('w') as stream:
        subprocess.run(args, cwd=cwd, stdout=stream, stderr=subprocess.STDOUT,
                       check=True, timeout=600)

def run(cwd, mode, label):
    args = ['./pluto']
    if mode in ('mpi', 'mpi-nvc'):
        args = [os.environ.get('CPU_MPIRUN', '/usr/bin/mpirun') if mode == 'mpi' else 'mpirun', '-np', '2'] + args
    command(args, cwd, label+'.stdout.log')
    text = (cwd/(label+'.stdout.log')).read_text()
    if mode in ('mpi', 'mpi-nvc'):
        shutil.copy2(cwd/'pluto.0.log', cwd/(label+'.rank0.log'))
        text += (cwd/(label+'.rank0.log')).read_text()
    if '> Done.' not in text:
        raise RuntimeError(mode+': missing completion')
    return text

for mode in MODES:
    if mode not in ARCH:
        raise ValueError(mode)
    work = OUT/mode
    work.mkdir(parents=True, exist_ok=False)
    defs = (ROOT/'Test_Problems/RMHD/Blast/definitions_01.hpp').read_text()
    defs += '\n#define PCP_NR YES\n#define PCP_FOFC YES\n#define FALLBACK_MODE FALLBACK_RECOMPUTE\n#define CT_RMHD_CONSISTENT YES\n'
    (work/'definitions.hpp').write_text(defs)
    shutil.copy2(ROOT/'Test_Problems/RMHD/Blast/init.cpp', work/'init.cpp')
    ini = (ROOT/'Test_Problems/RMHD/Blast/pluto_01.ini').read_text()
    ini = ini.replace('256', '32').replace('tstop            0.4', 'tstop            0.002')
    ini = ini.replace('first_dt         1.e-6', 'first_dt         1.e-4')
    (work/'pluto.ini').write_text(ini)
    make = (ROOT/'Src/Templates/makefile').read_text()
    make = '\n'.join('ARCH = '+ARCH[mode] if line.startswith('ARCH ') else
                     'gPLUTO_DIR = '+str(ROOT) if line.startswith('gPLUTO_DIR ') else line
                     for line in make.splitlines())
    make = make.replace('#   Additional_object_files_here',
        'include $(SRC)/RMHD/makefile\ninclude $(SRC)/MHD/CT/makefile\n'
        'OBJ += plm_states.o rk_step.o update_stage.o\n#   Additional_object_files_here')
    make += '\nunit: $(filter-out main.o,$(OBJ))\n'
    make += '\t$(CC) $(CFLAGS) $(INCLUDE_DIRS) '+str(HERE/'test_fma.cpp')+' -o test_fma.o\n'
    make += '\t$(CC) test_fma.o $(LDFLAGS) -o test_fma\n'
    make += '\t$(CC) $(CFLAGS) $(INCLUDE_DIRS) '+str(HERE/'test_pcp.cpp')+' -o test_pcp.o\n'
    make += '\t$(CC) $(filter-out main.o,$(OBJ)) test_pcp.o $(LDFLAGS) -o test_pcp\n'
    (work/'makefile').write_text(make)
    make_args = ['make', '-j6'] + (['CC='+os.environ.get('CPU_MPICXX', '/usr/bin/mpicxx')] if mode == 'mpi' else [])
    stage = (ROOT/'Src/RMHD/pcp_stage.cpp').read_text()
    stage = stage.replace('// printLog(', 'printLog(').replace('//          g_intStage', '           g_intStage')
    (work/'pcp_stage.cpp').write_text(stage)
    command(make_args + ['pluto', 'unit'], work, 'build.log')
    command(['./test_fma'], work, 'fma.log')
    command(['./test_pcp'], work, 'unit.log')
    run(work, mode, 'normal')
    shutil.copy2(work/'data.0001.dbl', work/'normal.dbl')
    for fault in ['fofc', 'failsafe']:
        command([sys.executable, str(HERE/'inject_stage.py'),
                 str(ROOT/'Src/update_stage.cpp'), fault], work, fault+'.prepare.log')
        command(make_args, work, fault+'.build.log')
        text = run(work, mode, fault)
        if 'PCP: stage=1 pass=0 retry=0 rejected=1' not in text:
            raise RuntimeError(mode+': rejection not exercised')
        if fault == 'failsafe' and ('pass=3 retry=0' not in text or 'PCP PressureFix: stage=1 retry=0' not in text or '[Fallback] retry=1' not in text or 't = 5.000000e-05' not in text):
            raise RuntimeError(mode+': rollback/halving not exercised')
        if fault == 'failsafe' and text.index('pass=3 retry=0') > text.index('PCP PressureFix: stage=1 retry=0'):
            raise RuntimeError(mode+': PressureFix ran before FOFC exhausted its passes')
        if fault == 'fofc' and 'PCP PressureFix:' in text:
            raise RuntimeError(mode+': successful FOFC unexpectedly reached PressureFix')
    results[mode] = {'kernel':'passed', 'normal':'passed', 'fofc':'passed', 'failsafe':'passed'}
    print(mode, results[mode], flush=True)

# Compare physical arrays from the normal runs, before fault injection.
if 'serial' in MODES:
    raw = (OUT/'serial/normal.dbl').read_bytes()
    reference = struct.unpack('='+str(len(raw)//8)+'d', raw)
    for mode in MODES:
        raw = (OUT/mode/'normal.dbl').read_bytes()
        values = struct.unpack('='+str(len(raw)//8)+'d', raw)
        if len(values) != len(reference):
            raise RuntimeError('output shape differs: '+mode)
        err = max(abs(a-b)/max(1.0,abs(a)) for a,b in zip(reference,values))
        results[mode]['max_scaled_difference_from_serial'] = err
        if not err < 1.e-8:
            raise RuntimeError('output differs: '+mode+' '+str(err))
(OUT/'results.json').write_text(json.dumps(results, indent=2)+'\n')
