////////////////////////////////////////////////////////////////////////
/*!
  //to be commented by D.Crocco
  \file
  \brief Resistive MHD module header file.

  Contains prototypes for the resistive MHD module.
  
  \authors A. Mignone (andrea.mignone@unito.it)\n
           T. Matsakos \n
           D. Crocco

  \date    Apr 17, 2025
*/
////////////////////////////////////////////////////////////////////////

void Resistive_eta (double *, double, double, double, double *, RunConfig &, double *);
void ComputeStaggeredEta (const Data *, Grid *);