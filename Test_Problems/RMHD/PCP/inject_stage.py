"""Create a test-only UpdateStage override; never alters production sources."""
from pathlib import Path
import sys
src=Path(sys.argv[1]).read_text()
mode=sys.argv[2]
condition='injected++ == 0' if mode=='fofc' else 'Dts->dt > 7.5e-5 && g_stepNumber == 0'
code='''
  // Deterministic rejection on one side of a rank boundary.
  static int injected=0;
  if(g_rank==0 && g_intStage==1 && (CONDITION)){
    #pragma acc serial present(data,dUc,box)
    {
      dUc(box->iend,box->jbeg,box->kbeg,RHO)=-2.0*data->Uc(box->iend,box->jbeg,box->kbeg,RHO);
    }
  }
'''.replace('CONDITION',condition)
needle='  mynvtxstop_();\n}\n\n////////////////////////////////////////////////////////////////////////\n/*! \n  Note:'
assert needle in src
src=src.replace(needle,code+needle,1)
Path('update_stage.cpp').write_text(src)
