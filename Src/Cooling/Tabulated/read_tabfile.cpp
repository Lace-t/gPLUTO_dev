#include "gpluto.hpp"

void ReadTabFile(RunConfig &run_conf){
  static double *L_tab, *T_tab;
  static int ntab;
  
  if (T_tab == NULL){
    FILE *fcool;
    printLog (" > Reading table from disk...\n");
    fcool = fopen("cooltable.dat","r");
    if (fcool == NULL){
      QuitPLUTO ("! from ReadTableFile(): cooltable.dat could not be found",1);
    }
    L_tab = ARRAY_1D(20000, double);
    T_tab = ARRAY_1D(20000, double);

    ntab = 0;
    while (fscanf(fcool, "%lf  %lf\n", T_tab + ntab, 
                                      L_tab + ntab)!=EOF) {
      ntab++;
    }
    
    fclose(fcool);
  }
  run_conf.ntab = ntab;
  run_conf.L_tab = L_tab;
  run_conf.T_tab = T_tab;
  
}