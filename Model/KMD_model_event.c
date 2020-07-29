/* KMD_model_event.c for R deSolve package
   ___________________________________________________

   Model File:  C:/Software/KMD_CL_Model/Model/KMD_model_event.MODEL

   Date:  Wed Jul 22 16:18:14 2020

   Created by:  "C:/Software/R_MCSIM/MCSIM_~1/MCSim/mod.exe v6.1.0"
    -- a model preprocessor by Don Maszle
   ___________________________________________________

   Copyright (c) 1993-2019 Free Software Foundation, Inc.

   Model calculations for compartmental model:

   13 States:
     apls = 0.0,
     aliv = 0.0,
     aspf = 0.0,
     arpf = 0.0,
     amet = 0.0,
     aumet = 0.0,
     totiv = 0.0,
     ivswtch = 0.0,
     totodose = 0.0,
     aoral = 0.0,
     auc_cprnt = 0.0,
     auc_cmet = 0.0,
     auc_ctot = 0.0,

   12 Outputs:
    "cpls",
    "cliv",
    "cspf",
    "crpf",
    "mbal",
    "ramet",
    "raumet",
    "cumet",
    "cmet",
    "vbal",
    "qbal",
    "riv",

   0 Inputs:

   27 Parameters:
     bw = 0.0,
     qc = 0.0,
     vliv = 0.0,
     vpls = 0.0,
     vrbc = 0.0,
     hct = 0.42,
     vspf = 0.0,
     vrpf = 0.0,
     rvurine = 0,
     vurinec = 0.0,
     vke1 = 0.0,
     vmaxu = 0.0,
     kmu = 0.0,
     qliv = 0.0,
     qspf = 0.0,
     qrpf = 0.0,
     pliv = 0.0,
     prpf = 0.0,
     pspf = 0.0,
     vmax = 0.0,
     km = 0.0,
     vkm1 = 0.0,
     MW = 0.0,
     ivdose = 0.0,
     boral = 0.0,
     ka = 0.0,
     fa = 0.0,
*/

#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <R_ext/Rdynload.h>

/* Model variables: States */
#define ID_apls 0x00000
#define ID_aliv 0x00001
#define ID_aspf 0x00002
#define ID_arpf 0x00003
#define ID_amet 0x00004
#define ID_aumet 0x00005
#define ID_totiv 0x00006
#define ID_ivswtch 0x00007
#define ID_totodose 0x00008
#define ID_aoral 0x00009
#define ID_auc_cprnt 0x0000a
#define ID_auc_cmet 0x0000b
#define ID_auc_ctot 0x0000c

/* Model variables: Outputs */
#define ID_cpls 0x00000
#define ID_cliv 0x00001
#define ID_cspf 0x00002
#define ID_crpf 0x00003
#define ID_mbal 0x00004
#define ID_ramet 0x00005
#define ID_raumet 0x00006
#define ID_cumet 0x00007
#define ID_cmet 0x00008
#define ID_vbal 0x00009
#define ID_qbal 0x0000a
#define ID_riv 0x0000b

/* Parameters */
static double parms[27];

#define bw parms[0]
#define qc parms[1]
#define vliv parms[2]
#define vpls parms[3]
#define vrbc parms[4]
#define hct parms[5]
#define vspf parms[6]
#define vrpf parms[7]
#define rvurine parms[8]
#define vurinec parms[9]
#define vke1 parms[10]
#define vmaxu parms[11]
#define kmu parms[12]
#define qliv parms[13]
#define qspf parms[14]
#define qrpf parms[15]
#define pliv parms[16]
#define prpf parms[17]
#define pspf parms[18]
#define vmax parms[19]
#define km parms[20]
#define vkm1 parms[21]
#define MW parms[22]
#define ivdose parms[23]
#define boral parms[24]
#define ka parms[25]
#define fa parms[26]

/* Forcing (Input) functions */
static double forc[0];


/* Function definitions for delay differential equations */

int Nout=1;
int nr[1]={0};
double ytau[1] = {0.0};

static double yini[13] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}; /*Array of initial state variables*/

void lagvalue(double T, int *nr, int N, double *ytau) {
  static void(*fun)(double, int*, int, double*) = NULL;
  if (fun == NULL)
    fun = (void(*)(double, int*, int, double*))R_GetCCallable("deSolve", "lagvalue");
  return fun(T, nr, N, ytau);
}

double CalcDelay(int hvar, double dTime, double delay) {
  double T = dTime-delay;
  if (dTime > delay){
    nr[0] = hvar;
    lagvalue( T, nr, Nout, ytau );
}
  else{
    ytau[0] = yini[hvar];
}
  return(ytau[0]);
}

/*----- Initializers */
void initmod (void (* odeparms)(int *, double *))
{
  int N=27;
  odeparms(&N, parms);
}

void initforc (void (* odeforcs)(int *, double *))
{
  int N=0;
  odeforcs(&N, forc);
}


/* Calling R code will ensure that input y has same
   dimension as yini */
void initState (double *y)
{
  int i;

  for (i = 0; i < sizeof(yini) / sizeof(yini[0]); i++)
  {
    yini[i] = y[i];
  }
}

void getParms (double *inParms, double *out, int *nout) {
/*----- Model scaling */

  int i;

  for (i = 0; i < *nout; i++) {
    parms[i] = inParms[i];
  }


  rvurine = vurinec * bw / 24.0 ;

  for (i = 0; i < *nout; i++) {
    out[i] = parms[i];
  }
  }
/*----- Dynamics section */

void derivs (int *neq, double *pdTime, double *y, double *ydot, double *yout, int *ip)
{
  /* local */ double raoral;
  /* local */ double ramet_sat;
  /* local */ double ramet_lin;
  /* local */ double cv;
  /* local */ double raumet_sat;
  /* local */ double raumet_lin;
  /* local */ double totexpo;
  /* local */ double totbody;
  /* local */ double totmetab;
  /* local */ double totclear;

  ydot[ID_ivswtch] = 0 ;

  yout[ID_riv] = y[ID_ivswtch] * ivdose * 1000 / MW ;

  ydot[ID_totodose] = 0 ;

  raoral = ka * fa * y[ID_aoral] ;

  ydot[ID_aoral] = - raoral ;

  yout[ID_cpls] = y[ID_apls] / vpls ;

  yout[ID_cliv] = y[ID_aliv] / vliv ;

  yout[ID_cspf] = y[ID_aspf] / vspf ;

  yout[ID_crpf] = y[ID_arpf] / vrpf ;

  yout[ID_cmet] = y[ID_amet] / vpls ;

  ramet_sat = ( vmax * yout[ID_cliv] / pliv ) / ( ( yout[ID_cliv] / pliv ) + km ) ;

  ramet_lin = vkm1 * yout[ID_cliv] / pliv ;

  yout[ID_ramet] = ramet_sat + ramet_lin ;

  ydot[ID_aliv] = qliv * ( yout[ID_cpls] - yout[ID_cliv] / pliv ) - yout[ID_ramet] + ka * raoral ;

  ydot[ID_arpf] = qrpf * ( yout[ID_cpls] - yout[ID_crpf] / prpf ) ;

  ydot[ID_aspf] = qspf * ( yout[ID_cpls] - yout[ID_cspf] / pspf ) ;

  cv = ( qliv * yout[ID_cliv] / pliv + qspf * yout[ID_cspf] / pspf + qrpf * yout[ID_crpf] / prpf ) / qc ;

  ydot[ID_apls] = qc * cv - qc * yout[ID_cpls] + yout[ID_riv] ;

  raumet_sat = ( vmaxu * yout[ID_cmet] ) / ( ( yout[ID_cmet] ) + kmu ) ;

  raumet_lin = vke1 * yout[ID_cmet] ;

  yout[ID_raumet] = raumet_sat + raumet_lin ;

  ydot[ID_aumet] = yout[ID_raumet] ;

  ydot[ID_amet] = yout[ID_ramet] - yout[ID_raumet] ;

  yout[ID_cumet] = yout[ID_raumet] / rvurine ;

  ydot[ID_totiv] = yout[ID_riv] ;

  ydot[ID_auc_cprnt] = yout[ID_cpls] ;

  ydot[ID_auc_cmet] = yout[ID_cmet] ;

  ydot[ID_auc_ctot] = yout[ID_cpls] + yout[ID_cmet] ;

  totexpo = y[ID_totiv] + y[ID_totodose] ;

  totbody = y[ID_apls] + y[ID_aliv] + y[ID_arpf] + y[ID_aspf] + y[ID_aoral] ;

  totmetab = y[ID_amet] ;

  totclear = y[ID_aumet] ;

  yout[ID_mbal] = totexpo - totbody - ( totmetab + totclear ) ;

  yout[ID_vbal] = bw - ( vpls + vrbc + vliv + vspf + vrpf ) ;

  yout[ID_qbal] = qc - ( qliv + qrpf + qspf ) ;

} /* derivs */


/*----- Jacobian calculations: */
void jac (int *neq, double *t, double *y, int *ml, int *mu, double *pd, int *nrowpd, double *yout, int *ip)
{

} /* jac */


/*----- Events calculations: */
void event (int *n, double *t, double *y)
{

  y[ID_aoral] = ( boral > 0 ? y[ID_aoral] + ( boral * bw * 1000 / MW ) : y[ID_aoral] ) ;
  y[ID_totodose] = ( boral > 0 ? y[ID_totodose] + ( boral * bw * 1000 / MW ) : y[ID_totodose] ) ;
  y[ID_ivswtch] = ( ivdose > 0 ? ( y[ID_ivswtch] == 0 ? 1 : 0 ) : 0 ) ;

} /* event */

/*----- Roots calculations: */
void root (int *neq, double *t, double *y, int *ng, double *gout, double *out, int *ip)
{

} /* root */

