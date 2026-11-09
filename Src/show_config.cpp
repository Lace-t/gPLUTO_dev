////////////////////////////////////////////////////////////////////////
/*! 
  \file  
  \brief Print useful information about the current computations.

  \author  A. Mignone (andrea.mignone@unito.it)
  \date    Apr 19, 2026
*/
////////////////////////////////////////////////////////////////////////
#include "gpluto.hpp"

////////////////////////////////////////////////////////////////////////
/*!
  Write a summary of the selected options
*/
////////////////////////////////////////////////////////////////////////
void ShowConfig (char *ini_file, RunConfig &run_config)
{
  FILE *fp;
  time_t time_now;
  char  str1[128], str2[128], str3[128], sline[512];
  double *dbl_pnt;  // For printLoging size of pointer to double 
  int    *int_pnt;  // For printLoging size of pointer to int 

  print ("\n");
  print ("          ____  __    __  ____________ \n");
  print ("   ____  / __ \\/ /   / / / /_  __/ __ \\ \n");
  print ("  / __ `/ /_/ / /   / / / / / / / / / / \n");
  print (" / /_/ / ____/ /___/ /_/ / / / / /_/ / \n");
  print (" \\__, /_/   /_____/\\____/ /_/  \\____/  \n");
  print ("/____/                                \n");
  print ("******************************************\n");
  print ("Version: %s  \n", gPLUTO_VERSION);
  
  print ("\n> System:\n\n");

  if ( (fp = fopen("sysconf.out","r")) != NULL){

    while (fscanf (fp, "%s %s %s\n", str1, str2, str3) != EOF) {
      if (!strcmp(str1,"USER")) 
        print ("  %s:             %s\n",str1, str3);
      else if (!strcmp(str1,"WORKING_DIR"))
        print ("  %s:      %s\n",str1, str3);
      else if (!strcmp(str1,"SYSTEM_NAME"))
        print ("  %s:      %s\n",str1, str3);
      else if (!strcmp(str1,"NODE_NAME"))
        print ("  %s:        %s\n",str1, str3);
      else if (!strcmp(str1,"ARCH"))
        print ("  %s:             %s\n",str1, str3);
      else if (!strcmp(str1,"BYTE_ORDER"))
        print ("  %s:       %s\n\n",str1, str3);
    }
    fclose(fp);

  }else{
    print ("  sysconf.out file not found... \n\n");
  }
 
  print ("  - Basic data types size (Bytes):\n");
  print ("    o sizeof (char)     = %d\n", sizeof(char));
  print ("    o sizeof (uchar)    = %d\n", sizeof(unsigned char));
  print ("    o sizeof (uint32_t) = %d\n", sizeof(uint32_t));
  print ("    o sizeof (short)    = %d\n", sizeof(short));
  print ("    o sizeof (ushort)   = %d\n", sizeof(unsigned short));
  print ("    o sizeof (int)      = %d\n", sizeof(int));
  print ("    o sizeof (long)     = %d\n", sizeof(long));
  print ("    o sizeof (*int)     = %d\n", sizeof(int_pnt));
  print ("    o sizeof (float)    = %d\n", sizeof(float));
  print ("    o sizeof (double)   = %d\n", sizeof(double));
  print ("    o sizeof (*double)  = %d\n", sizeof(dbl_pnt));

  print ("\n");
  print ("  - Derived data types size (Bytes):\n");
  print ("    o sizeof (Ary1D)     = %d\n", sizeof(Ary1D));
  print ("    o sizeof (Ary2D)     = %d\n", sizeof(Ary2D));
  print ("    o sizeof (Ary3D)     = %d\n", sizeof(Ary3D));
  print ("    o sizeof (Ary4D)     = %d\n", sizeof(Ary4D));
  
  print ("\n");
  print ("  - Structure data types size (Bytes):\n");
  print ("    o sizeof (cmdLine)    = %d\n", sizeof(cmdLine));
  print ("    o sizeof (Data)       = %d\n", sizeof(Data));
  #if DIVB_CONTROL == CONSTRAINED_TRANSPORT
  print ("    o sizeof (EMF)        = %d\n", sizeof(EMF));  
  #endif
  print ("    o sizeof (Grid)       = %d\n", sizeof(Grid));
  print ("    o sizeof (RBox)       = %d\n", sizeof(RBox));
  print ("    o sizeof (Restart)    = %d\n", sizeof(Restart));
  print ("    o sizeof (Sweep)      = %d\n", sizeof(Sweep));
  print ("    o sizeof (timeStep)   = %d\n", sizeof(timeStep));


  time(&time_now);
  print("\n> Local time:       %s\n",asctime(localtime(&time_now)));
      
// -- printLog command line arguments -- 

//  print ("> Cmd line args:    ");
//  for (n = 1; n < argc; n++) print ("%s ",argv[n]);
// print ("\n\n");

// -- printLog problem configuration -- 

  print ("> Header configuration:\n\n");

  if (PHYSICS == HD)        print ("  PHYSICS:          HD\n");
  if (PHYSICS == RHD)       print ("  PHYSICS:          RHD\n");
  if (PHYSICS == MHD)       print ("  PHYSICS:          MHD    [div.B: ");
  if (PHYSICS == RMHD)      print ("  PHYSICS:          RMHD   [div.B: ");
  if (PHYSICS == ResRMHD)   print ("  PHYSICS:          ResMHD [div.B: ");
  if (PHYSICS == GRMHD)     print ("  PHYSICS:          GRMHD   [div.B: ");
#if (PHYSICS == MHD) || (PHYSICS == RMHD) || (PHYSICS == ResRMHD) || (PHYSICS == GRMHD)
  #if DIVB_CONTROL == NO
  print ("None]\n");
  #elif DIVB_CONTROL == EIGHT_WAVES
    print ("Powell's 8wave]\n");
  #elif DIVB_CONTROL == DIV_CLEANING
    #if GLM_EXTENDED == NO 
    print ("Divergence Cleaning (GLM)]\n");
    #elif GLM_EXTENDED == YES
    print ("Divergence Cleaning (Extended GLM)]\n");
    #endif
    
  #elif DIVB_CONTROL == CONSTRAINED_TRANSPORT
    print ("CT/");
    #if CT_EMF_AVERAGE == ARITHMETIC
    print ("Ar. average]\n");
    #elif CT_EMF_AVERAGE == CT_CONTACT
    print ("CT_CONTACT]\n");
    #elif CT_EMF_AVERAGE == UCT0
    print ("UCT0]\n");
    #elif CT_EMF_AVERAGE == UCT_HLL
    print ("UCT_HLL]\n");
    #elif CT_EMF_AVERAGE == CT_FLUX
    print ("CT_FLUX]\n");
    #elif CT_EMF_AVERAGE == CT_MAXWELL
    print ("CT_MAXWELL]\n");
    #elif CT_EMF_AVERAGE == UCT_HLLD
    print ("UCT_HLLD]\n");
    #elif CT_EMF_AVERAGE == UCT_GFORCE
    print ("UCT_GFORCE]\n");
    #endif
  #endif
#endif

#if (PHYSICS == ResRMHD)
  print ("                    ResMHD [div.E: ");
  #if DIVE_CONTROL == NO
  print ("None]\n");
  #elif DIVE_CONTROL == EIGHT_WAVES
    print ("Powell's 8wave]\n");
  #elif DIVE_CONTROL == DIV_CLEANING
    #if GLM_EXTENDED == NO 
    print ("Divergence Cleaning (GLM)]\n");
    #elif GLM_EXTENDED == YES
    print ("Divergence Cleaning (Extended GLM)]\n");
    #endif
  #elif DIVE_CONTROL == CONSTRAINED_TRANSPORT
    print ("CT]\n");
  #endif
#endif

  print ("  DIMENSIONS:       %d\n", DIMENSIONS);

#if PHYSICS != GRMHD
  print ("  GEOMETRY:         ");
  if (GEOMETRY == CARTESIAN)    print ("Cartesian\n");
  if (GEOMETRY == CYLINDRICAL)  print ("Cylindrical\n");
  if (GEOMETRY == SPHERICAL)    print ("Spherical\n");
#else
  print ("  METRIC:           ");
  if (GR_METRIC == CARTESIAN)       print ("Cartesian\n");
  if (GR_METRIC == CYLINDRICAL)     print ("Cylindrical\n");
  if (GR_METRIC == SPHERICAL)       print ("Spherical\n");
  if (GR_METRIC == BOYER_LINDQUIST) print ("Boyer-Lindquist\n");
  if (GR_METRIC == KERR_SCHILD)     print ("Kerr-Schild\n");
  if (GR_METRIC == USERDEF_METRIC)  print ("User-defined\n");
#endif

  print ("  BODY_FORCE:       ");
  print (BODY_FORCE == NO ? "NO\n":"EXPLICIT\n");

#if COOLING == H2_COOL
  print ("  COOLING:          H2_COOL\n");
#elif COOLING == KROME
  print ("  COOLING:          KROME\n");
#elif COOLING == MINEq
  print ("  COOLING:          MINEq\n");
#elif COOLING == POWER_LAW
  print ("  COOLING:          POWER_LAW\n");
#elif COOLING == SYNCH
  print ("  COOLING:          SYNCH\n");
#elif COOLING == SNEq
  print ("  COOLING:          SNEq\n");
#elif COOLING == TABULATED
  print ("  COOLING:          Tabulated\n");
#endif

  print ("  RECONSTRUCTION:   ");
  if (RECONSTRUCTION == FLAT)      print ("Flat ");
  if (RECONSTRUCTION == LINEAR)    print ("Linear ");
  if (RECONSTRUCTION == LimO3)     print ("LimO3 ");
  if (RECONSTRUCTION == WENO3)     print ("WENO3 ");
  if (RECONSTRUCTION == PARABOLIC) print ("Parabolic ");
  if (RECONSTRUCTION == WENOZ)     print ("WENOZ ");
  if (RECONSTRUCTION == MP5)       print ("MP5 ");
  #ifdef CHAR_LIMITING
  if (CHAR_LIMITING == YES) print (" (Characteristic lim)\n");
  else                      print (" (Primitive lim)\n");
  #endif

  print ("  TIME STEPPING:    ");
  if (TIME_STEPPING == EULER)      print ("Euler\n");
  if (TIME_STEPPING == RK2)        print ("RK2\n");
  if (TIME_STEPPING == RK3)        print ("RK3\n");
  if (TIME_STEPPING == RK4)        print ("RK4\n");
  if (TIME_STEPPING == ARK4)       print ("ARK4\n");
  if (TIME_STEPPING == IERK45)     print ("IERK45\n");


  print ("  TRACERS:          %d\n", NTRACER);
  print ("  VARIABLES:        %d\n", NVAR);
  print ("  ENTROPY_SWITCH:   %s\n",(ENTROPY_SWITCH != NO ? "ENABLED":"NO"));
#if PHYSICS == MHD 
  print ("  BACKGROUND_FIELD: %s\n",(BACKGROUND_FIELD == YES ? "YES":"NO"));
#endif

  print ("  LOADED MODULES:\n");

  #ifdef HIGH_ORDER
  if      (HO_LAP_LIMITER == NO)             
    snprintf(str1, sizeof(str1), "%s", "NO");
  else if (HO_LAP_LIMITER == HO_JAMESON_LIM) 
    snprintf(str1, sizeof(str1), "%s", "HO_JAMESON_LIM");
  else if (HO_LAP_LIMITER == HO_RDER_LIM)    
    snprintf(str1, sizeof(str1), "%s", "HO_RDER_LIM");
  else snprintf (str1,sizeof(str1), "%s", "[undef]");

  if      (HO_ORDER_REDUCTION == NO)     snprintf(str2, sizeof(str2), "%s", "NO");
  else if (HO_ORDER_REDUCTION == LINEAR) snprintf(str2, sizeof(str2), "%s", "LINEAR");
  else if (HO_ORDER_REDUCTION == WENO3)  snprintf(str2, sizeof(str2), "%s", "WENO3");
  else snprintf(str2, sizeof(str2), "%s", "[undef]");

  print ("\n  o [HIGH_ORDER]\n");
  print ("     - HO_LAP_LIMITER:         %s\n", str1);
  print ("     - HO_ORDER_REDUCTION:     %s\n", str2);
  #endif

 #ifdef SHEARINGBOX
    print ("\n  o [SHEARINGBOX]\n");
    print ("     - Order:             %d\n", SB_ORDER);
    print ("     - Sym Hydro Flux:    %s\n", 
             (SB_SYMMETRIZE_HYDRO == YES ? "YES":"NO"));
    print ("     - Sym Ey:            %s\n", 
             (SB_SYMMETRIZE_EY == YES ? "YES":"NO"));
    print ("     - Sym Ez:            %s\n", 
             (SB_SYMMETRIZE_EZ == YES ? "YES":"NO"));
    print ("     - Force EMF periods: %s\n", 
             (SB_FORCE_EMF_PERIODS == YES ? "YES":"NO"));
  #endif
  #ifdef FARGO
  print ("\n  o [FARGO]\n");
  print ("     - Order:         %d\n", FARGO_ORDER);
  print ("     - Av. Frequency: %d\n", FARGO_NSTEP_AVERAGE);
  #endif
  print ("\n");

  print ("  ROTATION:         ");
  print(ROTATING_FRAME == YES ? "YES\n":"NO\n");

  print ("  EOS:              ");
  if      (EOS == IDEAL)        print ("Ideal\n");
  else if (EOS == PVTE_LAW)     print ("PVTE_LAW\n");
  else if (EOS == BAROTROPIC)   print ("Barotropic\n");
  else if (EOS == ISOTHERMAL)   print ("Isothermal\n");
  else if (EOS == TAUB)         print ("Taub - TM\n");
  else                          print ("None\n");

  #if PARABOLIC_FLUX != NO
   print ("  DIFFUSION TERMS:");
   #if (RESISTIVITY == EXPLICIT) 
    print ("  Resistivity  [EXPLICIT]\n");  
   #elif (RESISTIVITY == SUPER_TIME_STEPPING)
    print ("  Resistivity  [STS]\n");  
   #elif (RESISTIVITY == RK_LEGENDRE)
    print ("  Resistivity  [RKL]\n");
   #endif

   #if (THERMAL_CONDUCTION == EXPLICIT) 
    print ("  Thermal Conduction [EXPLICIT]\n");  
   #elif (THERMAL_CONDUCTION == SUPER_TIME_STEPPING)
    print ("  Thermal Conduction [STS]\n");  
   #elif (THERMAL_CONDUCTION == RK_LEGENDRE)
    print ("  Thermal Conduction [RKL]\n");
   #endif

   #if (VISCOSITY == EXPLICIT) 
    print ("  Viscosity [EXPLICIT]\n");  
   #elif (VISCOSITY == SUPER_TIME_STEPPING)
    print ("  Viscosity [STS]\n");  
   #elif (VISCOSITY == RK_LEGENDRE)
    print ("  Viscosity [RKL]\n");
   #endif
  #endif

  print ("  WARNING_MESSAGES: %d\n", WARNING_MESSAGES);
  print ("  unit density:     %1.6e\n", run_config.unit_density);
  print ("  unit length:      %1.6e\n", run_config.unit_length);
  print ("  unit velocity:    %1.6e\n", run_config.unit_velocity);
  print ("  unit time:        %1.6e\n", run_config.unit_time);
  print ("\n");

  // -----------------------------------------------------
  // Print runtime configuration info (definitions.hpp 
  // and from pluto.ini)
  // -----------------------------------------------------
/*   
  print ("> Header file configuration (definitions.hpp):\n\n");
  print ("  +----------------------------------------------------------\n");
  fp = fopen("definitions.hpp","r");
  while ( fgets(sline, 512, fp) != NULL ) {
    print ("  | %s",sline);
  }
  fclose(fp);
  print ("  +---------------------------------------------------------\n\n");
*/
  print ("> Runtime configuration (%s):\n\n", ini_file);
  print ("  +----------------------------------------------------------\n");
  fp = fopen(ini_file,"r");
  while ( fgets(sline, 512, fp) != NULL ) {
    print ("  | %s",sline);
  }
  fclose(fp);
  print ("  +---------------------------------------------------------\n");

}

////////////////////////////////////////////////////////////////////////
/*!
  Show units when cooling is enabled.
*/
////////////////////////////////////////////////////////////////////////
void ShowUnits (RunConfig &run_config)
{
  double ul = run_config.unit_length;
  double ud = run_config.unit_density;
  double uv = run_config.unit_velocity;

#if COOLING != NO
  print ("> Cooling Module:    ");
  if (COOLING == SNEq)      print (" SNEq\n");
  if (COOLING == MINEq)     print (" MINEq\n");
  if (COOLING == TABULATED) print (" Tabulated\n");
  if (COOLING == H2_COOL)   print (" H2_COOL \n");
  if (COOLING == KROME)     print (" KROME \n");
#endif

  print ("> Normalization Units:\n\n");
  print ("  [Density]:      %8.3e (gr/cm^3), %8.3e (1/cm^3)\n",
          ud, ud/CONST_mp);
  print ("  [Pressure]:     %8.3e (dyne/cm^2)\n", ud*uv*uv);
  print ("  [Velocity]:     %8.3e (cm/s)\n", uv);
  print ("  [Length]:       %8.3e (cm)\n", ul);
  print ("  [Temperature]:  %8.3e X (p/rho*mu) (K)\n",Kelvin(run_config));
  print ("  [Time]:         %8.3e (sec), %8.3e (yrs) \n",
            ul/uv, ul/uv/86400./365.);
#if PHYSICS == MHD || PHYSICS == RMHD || PHYSICS == GRMHD
  print ("  [Mag Field]:    %8.3e (Gauss)\n",
           uv*sqrt(4.0*CONST_PI*ud));
#endif

  print (" \n");
}