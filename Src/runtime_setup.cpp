////////////////////////////////////////////////////////////////////////
/*!
  \file
  \brief Read runtime information from pluto.ini.

  Parse and read runtime information from the initialization file
  pluto.ini (default) and sets value of the runtime structure.

  \authors A. Mignone (andrea.mignone@unito.it)
  \date    Apr 17, 2026
*/
////////////////////////////////////////////////////////////////////////
#include "gpluto.hpp"

#define NOPT      32              //  # of possible options in a menu 
#define NLEN      128             //  # default string length         

#define COMPARE(s1,s2,ii) \
        for ( (ii) = 1 ; (ii) < NOPT && !(strcmp ( (s1), (s2)) == 0); (ii)++);

////////////////////////////////////////////////////////////////////////
/*!
  Open and parse the runtime initialization file.
  Assign values to the runtime structure.
  
  \param [out]  runtime   pointer to a Runtime structure
  \param [in]   cmd_line  pointer to a cmdLine structure (useful, e.g.,
  to create the domain using the \c -xres option)
  \param [in]   ini_file  the name of the initialization file (default
  is "pluto.ini") specified with the \c -i option.
*/
////////////////////////////////////////////////////////////////////////
int RuntimeSetup (Runtime *runtime, cmdLine *cmd_line, char *ini_file)
{
  int    idim, ip, ipos, itype, nlines;
  int    include_dir[] = {INCLUDE_IDIR, INCLUDE_JDIR, INCLUDE_KDIR};
  char  bound_opt[NOPT][NLEN], str_var[512], *str;
  char* dummy; // dummy variables to suppress unused variable warning
  int idummy; 
  char  const *glabel[]     = {"X1-grid", "X2-grid","X3-grid"};
  char  const *bbeg_label[] = {"X1-beg", "X2-beg","X3-beg"};
  char  const *bend_label[] = {"X1-end", "X2-end","X3-end"};
  double dbl_var, rx;
  Output *output;
  FILE *fp;

  for (itype = 0; itype < NOPT; itype++) {
    strcpy(bound_opt[itype], "0000");
  }

  // ------------------------------------------------------
  // 0. Map strings into actual constant values
  // ------------------------------------------------------

  strcpy(bound_opt[OUTFLOW],      "outflow");
  strcpy(bound_opt[REFLECTIVE],   "reflective");
  strcpy(bound_opt[AXISYMMETRIC], "axisymmetric");
  strcpy(bound_opt[EQTSYMMETRIC], "eqtsymmetric");
  strcpy(bound_opt[PERIODIC],     "periodic");
  strcpy(bound_opt[SHEARING],     "shearingbox");
  strcpy(bound_opt[ROTATED],      "rotated");
  strcpy(bound_opt[USERDEF],      "userdef");
  strcpy(bound_opt[SINGULAR_AXIS], "singular_axis");

  runtime->log_freq = 1; // Default

  nlines = ParamFileRead(ini_file);

  // ------------------------------------------------------
  // 1. [Grid] Section
  // ------------------------------------------------------

  for (idim = 0; idim < 3; idim++){
    runtime->npoint[idim] = 1;

    ipos = 0;
    runtime->xbeg[idim]   = atof(ParamFileGet(glabel[idim], ++ipos));
    runtime->npoint[idim] = atoi(ParamFileGet(glabel[idim], ++ipos));

    runtime->xend[idim] = atof(ParamFileGet(glabel[idim], ++ipos));

    if ( include_dir[idim] == FALSE && runtime->npoint[idim] != 1 ) {
      printf ("! RuntimeSetup(): %d point(s) on dim. %d is NOT valid, resetting to 1\n",
              runtime->npoint[idim],idim+1);
      runtime->npoint[idim] = 1;
    }
  }

  // Check validity of domain 

  #if GEOMETRY == CYLINDRICAL
  if (fabs(runtime->xbeg[JDIR] - runtime->xend[JDIR]) > 2.00000001*CONST_PI){
    printf ("! RuntimeSetup(): it seems that Lphi > 2\\pi\n");
    QuitPLUTO(" ",1);
  }
  #endif

  #if GEOMETRY == SPHERICAL
  if (fabs(runtime->xbeg[JDIR] - runtime->xend[JDIR]) > 1.0000001*CONST_PI){
    printf ("! RuntimeSetup(): it seems that Lth > \\pi\n");
    QuitPLUTO(" ",1);
  }
  #endif

  // ------------------------------------------------------
  // 1a. Change the resolution if cmd_line->xres has 
  //     been given
  // ------------------------------------------------------

  if (cmd_line->xres > 1) {
    rx =  (double)cmd_line->xres/(double)runtime->npoint[IDIR];
    for (idim = 0; idim < 3; idim++){

      if (!include_dir[idim]) continue;

      dbl_var = (double)runtime->npoint[idim];
      runtime->npoint[idim] = MAX( (int)(dbl_var*rx), 1);
    }
  }

  // ------------------------------------------------------
  // 2. [Time] Section
  // ------------------------------------------------------

  runtime->cfl         = atof(ParamFileGet("CFL", 1));

  if (ParamExist ("CFL_par")) runtime->cfl_par = atof(ParamFileGet("CFL_par", 1));
  else {
    int dimensions = INCLUDE_IDIR + INCLUDE_JDIR + INCLUDE_KDIR;
    runtime->cfl_par = 0.8/(double)dimensions;
  }

  if (ParamExist ("rmax_par")) runtime->rmax_par = atof(ParamFileGet("rmax_par", 1));
  else                         runtime->rmax_par = 100.0;

  runtime->cfl_max_var = atof(ParamFileGet("CFL_max_var", 1));

  runtime->tstop       = atof(ParamFileGet("tstop", 1));
  if (ParamExist ("tfreeze")) runtime->tfreeze = atof(ParamFileGet("tfreeze", 1));
  else                        runtime->tfreeze = runtime->tstop+1;
  runtime->first_dt    = atof(ParamFileGet("first_dt", 1));

  if (ParamExist ("rad_first_dt"))    runtime->rad_first_dt    = atof(ParamFileGet("rad_first_dt", 1));
  else                                runtime->rad_first_dt    = runtime->first_dt / 100.0;
  if (ParamExist ("rad_CFL"))         runtime->rad_cfl         = atof(ParamFileGet("rad_CFL", 1));
  else                                runtime->rad_cfl         = runtime->cfl;
  if (ParamExist ("rad_CFL_max_var")) runtime->rad_cfl_max_var = atof(ParamFileGet("rad_CFL_max_var", 1));
  else                                runtime->rad_cfl_max_var = runtime->cfl_max_var;
  if (ParamExist ("rad_max_iter"))    runtime->rad_max_iter    = atoi(ParamFileGet("rad_max_iter", 1));
  else                                runtime->rad_max_iter    = 100;
  if (ParamExist ("rad_tolerance"))   runtime->rad_tolerance   = atof(ParamFileGet("rad_tolerance", 1));
  else                                runtime->rad_tolerance   = 1e-7;
  if (ParamExist ("small_enr"))    runtime->small_enr    = atof(ParamFileGet("small_enr", 1));
  else                                runtime->small_enr    = 1e-16;

  // ------------------------------------------------------
  // 3. [Solver] Section
  // ------------------------------------------------------

  strcpy (runtime->solv_type, ParamFileGet("Solver",1));
  strcpy (runtime->rad_solv_type,"hll");  // Default 
  if (ParamExist ("RadSolver")) {
    strcpy (runtime->rad_solv_type, ParamFileGet("RadSolver",1));
  }

  // ------------------------------------------------------
  // 4. [Boundary] Section
  // ------------------------------------------------------

  for (idim = 0; idim < 3; idim++){

    str = ParamFileGet(bbeg_label[idim], 1);
    COMPARE (str, bound_opt[itype], itype);
    if (itype == NOPT) {
      printf ("! RuntimeSetup(): unknown left boundary '%s'  \n", str);
      QuitPLUTO(" ", 1);
    }
    runtime->left_bound[idim] = itype;
  }

  for (idim = 0; idim < 3; idim++){

    str = ParamFileGet(bend_label[idim], 1);
    COMPARE (str, bound_opt[itype], itype);
    if (itype == NOPT) {
      printf ("! RuntimeSetup(): unknown right boundary '%s'  \n", str);
      QuitPLUTO(" ", 1);
    }
    runtime->right_bound[idim] = itype;
  }

  // ------------------------------------------------------
  // 5. [Static Grid Output] Section
  // ------------------------------------------------------

  runtime->user_var = 0;
  if (ParamExist("uservar")){
    runtime->user_var = atoi(ParamFileGet("uservar", 1));
    for (ip = 0; ip < runtime->user_var; ip++){
      if ( (str = ParamFileGet("uservar", 2 + ip)) != NULL){
        strcpy (runtime->user_var_name[ip], str);
      }else{
        printf ("! RuntimeSetup(): missing name after user var name '%s'\n",
                runtime->user_var_name[ip-1]);
      QuitPLUTO(" ", 1);
      }
    }
  }

  // Set output directory

  strcpy (runtime->output_dir, "./");  // default is current directory 
  if (ParamExist("output_dir")){
    str = ParamFileGet("output_dir",1);
    strcpy (runtime->output_dir, str);
  }

  // Check if we have write access and if the directory exists

  snprintf (str_var,NLEN,"%s/tmp09123.txt",runtime->output_dir); // test file
  fp = fopen(str_var,"w");  // -- open test file for writing 
  if (fp == NULL){
    printf ("! Setup(): cannot access directory '%s'.\n", runtime->output_dir);
    printf ("!          Please check that the directory exists\n");
    printf ("!          and you have write permission.\n");
    QuitPLUTO(" ", 1);
  }else{
    fclose(fp);
    remove(str_var);  // -- remove file -- 
  }

  // 5a. Set output defaults 

  output = runtime->output; 
  for (int n=0; n < MAX_OUTPUT_TYPES; n++){
    output[n].user_var = runtime->user_var;
    output[n].type   = n;   // Make sure the array index corresponds to its type
    output[n].cgs    =  NO;
    output[n].dt     = -1.0;
    output[n].dn     = -1;
    output[n].dclock = -1.0;
  }

  // dbl output 

  output = &(runtime->output[DBL_OUTPUT]);
  GetOutputFrequency(output, "dbl");

  strcpy (output->mode, ParamFileGet("dbl",3));
  if (   strcmp(output->mode,"single_file")
      && strcmp(output->mode,"multiple_files")){
    printf (
     "! RuntimeSetup(): expecting 'single_file' or 'multiple_files' in dbl output\n");
    QuitPLUTO(" ", 1);
  }

  // flt output

  if (ParamExist("flt")){
    output = &(runtime->output[FLT_OUTPUT]);
    GetOutputFrequency(output, "flt");

    strcpy (output->mode, ParamFileGet("flt",3));
    if (    strcmp(output->mode,"single_file")
         && strcmp(output->mode,"multiple_files")){
       printf (
       "! RuntimeSetup(): expecting 'single_file' or 'multiple_files' in flt output\n");
      QuitPLUTO(" ", 1);
    }
    if (ParamFileHasBoth ("flt","cgs")) output->cgs = YES;
  }

  // hdf5 output

  if (ParamExist("dbl.h5")){
    output = &(runtime->output[DBL_H5_OUTPUT]);
    GetOutputFrequency(output, "dbl.h5");
  }
  if (ParamExist("flt.h5")){
    output = &(runtime->output[FLT_H5_OUTPUT]);
    GetOutputFrequency(output, "flt.h5");
  }

  // vtk

  if (ParamExist ("vtk")){
    output = &(runtime->output[VTK_OUTPUT]);
    GetOutputFrequency(output, "vtk");

    if (ParamFileGet("vtk",3) == NULL){
      printf ("! RuntimeSetup(): extra field missing in vtk output\n");
      QuitPLUTO(" ", 1);
    }
    strcpy (output->mode, ParamFileGet("vtk",3));
    if (   strcmp(output->mode,"single_file")
        && strcmp(output->mode,"multiple_files")){
       printf ("! RuntimeSetup(): expecting 'single_file' or 'multiple_files' in\n");
       printf ("                  vtk output\n");
      QuitPLUTO(" ", 1);
    }
    if (ParamFileHasBoth ("vtk","cgs")) output->cgs = 1;
  }

  // vtk err output

  if (ParamExist ("vtk_err")){
    output = &(runtime->output[VTK_ERR_OUTPUT]);
    GetOutputFrequency(output, "vtk_err");

    if (ParamFileGet("vtk_err",3) == NULL){
      printf ("! RuntimeSetup(): extra field missing in vtk_err output\n");
      QuitPLUTO(" ", 1);
    }
    strcpy (output->mode, ParamFileGet("vtk_err",3));
    if (   strcmp(output->mode,"single_file")
        && strcmp(output->mode,"multiple_files")){
       printf ("! RuntimeSetup(): expecting 'single_file' or 'multiple_files' in\n");
       printf ("                  vtk_err output\n");
      QuitPLUTO(" ", 1);
    }
  }

  // tab output

  if (ParamExist ("tab")){
    output = &(runtime->output[TAB_OUTPUT]);
    GetOutputFrequency(output, "tab");
    if (ParamFileHasBoth ("tab","cgs")) output->cgs = 1;
  }

  // ppm output

  if (ParamExist ("ppm")){
    output = &(runtime->output[PPM_OUTPUT]);
    GetOutputFrequency(output, "ppm");
  }

  // png output

  if (ParamExist ("png")){
    output = &(runtime->output[PNG_OUTPUT]);
    GetOutputFrequency(output, "png");
  }

  // log frequency

  strcpy (runtime->log_dir, runtime->output_dir);
  if (ParamExist ("log_dir")){
    str = ParamFileGet("log_dir",1);
    strcpy (runtime->log_dir, str);
  }

  runtime->log_freq = atoi(ParamFileGet("log", 1));
  runtime->log_freq = MAX(runtime->log_freq, 1);

  // Analysis

  if (ParamExist ("analysis")){
    runtime->anl_dt = atof(ParamFileGet("analysis", 1));
    runtime->anl_dn = atoi(ParamFileGet("analysis", 2));
  }else{
    runtime->anl_dt = -1.0;   // Default
    runtime->anl_dn = -1;
  }

  // ------------------------------------------------------
  // 6. [Particles] Section
  // ------------------------------------------------------

  #if PARTICLES != NO
//  runtime->Nparticles_cell = atoi(ParamFileGet("Nparticles", 1));

// Split the string nparticles/cell
{
int n,m;
if (sscanf(ParamFileGet("Nparticles", 1), "%d/%d", &n, &m) == 2) {
    // success: n = 12, m = 34
}
runtime->Nparticles_per_block   = n;  
runtime->Ncells_per_block       = m;  // cell size (in unit of base grid), 
                                     //    or  grid coarsening factor
printf ("n = %d; m = %d\n",n,m);
}

  runtime->particles_tstart = 0.0;
  if (ParamExist ("particles_tstart")) {
    runtime->particles_tstart = atof(ParamFileGet("particles_tstart", 1));
  }

  // Particles dbl output

  if (ParamExist("particles_dbl")){
    output = &(runtime->output[PARTICLES_DBL_OUTPUT]);
    GetOutputFrequency(output, "particles_dbl");
  }

  // Particles flt output

  if (ParamExist("particles_flt")){
    output = &(runtime->output[PARTICLES_FLT_OUTPUT]);
    GetOutputFrequency(output, "particles_flt");
  }

  // Particles vtk output

  if (ParamExist("particles_vtk")){
    output = &(runtime->output[PARTICLES_VTK_OUTPUT]);
    GetOutputFrequency(output, "particles_vtk");
  }

  // Particles tab output

  if (ParamExist("particles_tab")){
    output = &(runtime->output[PARTICLES_TAB_OUTPUT]);
    GetOutputFrequency(output, "particles_tab");
  }

  // Particles h5part output

  if (ParamExist("particles_hdf5")){
    output = &(runtime->output[PARTICLES_HDF5_OUTPUT]);
    GetOutputFrequency(output, "particles_hdf5");
  }

  // Particles analysis

  if (ParamExist ("particles_analysis")){
    runtime->Panl_dt = atof(ParamFileGet("particles_analysis", 1));
    runtime->Panl_dn = atoi(ParamFileGet("particles_analysis", 2));
  }else{
    runtime->Panl_dt = -1.0;   // -- defaults -- 
    runtime->Panl_dn = -1;
  }
  #endif

  // ------------------------------------------------------
  // 7. [Parameters] Section
  // ------------------------------------------------------

  fp = fopen(ini_file,"r");

  // Find position at "[Parameters" 

  for (ipos = 0; ipos <= nlines; ipos++){
    dummy = fgets(str_var, 512, fp);

    if (strlen(str_var) > 0) {
      str = strtok (str_var,"]");
      if (strcmp(str,"[Parameters") == 0) break;
    }
  }

  dummy = fgets(str_var, 512, fp);

  for (ip = 0; ip < USER_DEF_PARAMETERS; ip++){
    idummy= fscanf (fp,"%s \n", str_var);
    // if (strlen(str_var) <= 1 ) continue;
    dbl_var = atof(ParamFileGet(str_var,1));
    runtime->aux[ip] = dbl_var;
    dummy = fgets(str_var, sizeof(str_var), fp); // use fgets to advance to next line 
  }
  fclose(fp);

  return(0);
}
#undef COMPARE
#undef NOPT
#undef NLEN

////////////////////////////////////////////////////////////////////////
void GetOutputFrequency(Output *output, const char *output_format)
/*!
  Set the intervals between output files.
  This can be done in three different ways:
  
  - dt: time interval in code units
  - dn: step interval
  - dclock: actual clock time (in hours)
  
  However, dn and dclock are mutually exclusive.
*/
////////////////////////////////////////////////////////////////////////
{
  char *str;
  int len, nhrs, nmin, nsec;

// -- time interval in code units (dt) -- 

  output->dt = atof(ParamFileGet(output_format, 1));

// -- check the 2nd field and decide to set "dn" or "dclock" -- 

  str = ParamFileGet(output_format,2);
  len = strlen(str);
  if (str[len-1] == 'h'){
    output->dclock = atof(str);    // clock interval in hours 
    nhrs = (int)output->dclock;    // integer part 
    nmin = (int)((output->dclock - nhrs)*100.0); // remainder in minutes 
    if (nmin >= 60){
      printf ("! OutputFrequency: number of minutes exceeds 60 in %s output\n",
               output_format);
      QuitPLUTO(" ", 1);
    }
    output->dclock = nhrs*3600.0 + nmin*60;  // convert to seconds 
    output->dn     = -1;
  }else if (str[len-1] == 'm'){
    output->dclock = atof(str);      // clock interval in minutes 
    nmin = (int)output->dclock;      // integer part 
    nsec = (int)((output->dclock - nmin)*100.0); // remainder in seconds 
    if (nsec >= 60){
      printf ("! OutputFrequency: number of seconds exceeds 60 in %s output\n",
               output_format);
      QuitPLUTO(" ", 1);
    }
    output->dclock = nmin*60.0 + nsec;
    output->dn     = -1;
  }else if (str[len-1] == 's'){
    output->dclock = atof(str);           // clock interval in seconds 
    output->dn     = -1;
  }else{
    output->dclock = -1.0;
    output->dn     = atoi(ParamFileGet(output_format, 2));
  }
}

////////////////////////////////////////////////////////////////////////
static Runtime q;
void RuntimeSet(Runtime *runtime)
/*!
  Store a static copy of the runtime structure for later access.
*/
////////////////////////////////////////////////////////////////////////
{
  q = *runtime;
}

////////////////////////////////////////////////////////////////////////
Runtime *RuntimeGet(void)
/*!
  Return a pointer to runtime structure.
*/
////////////////////////////////////////////////////////////////////////
{
  return &q;
}
