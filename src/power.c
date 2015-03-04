/* ==========================================================================*/
/*   Copyright (c) 2015       Cullan Howlett & Marc Manera,                  */
/*                            Institute of Cosmology and Gravitation         */
/*                            (University of Portsmouth) & University        */
/*                            College London.                                */
/*                                                                           */
/*   This file is part of L-PICOLA.                                          */
/*                                                                           */
/*   L-PICOLA is free software: you can redistribute it and/or modify        */
/*   it under the terms of the GNU General Public License as published by    */
/*   the Free Software Foundation, either version 3 of the License, or       */
/*   (at your option) any later version.                                     */
/*                                                                           */
/*   L-PICOLA is distributed in the hope that it will be useful,             */
/*   but WITHOUT ANY WARRANTY; without even the implied warranty of          */
/*   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           */
/*   GNU General Public License for more details.                            */
/*                                                                           */
/*   You should have received a copy of the GNU General Public License       */
/*   along with L-PICOLA.  If not, see <http://www.gnu.org/licenses/>.       */
/* ==========================================================================*/

/* =============================================================================================*/
/* This file contains the routines for setting up the initial power spectrum of the simulation. */
/* =============================================================================================*/

#include "vars.h"
#include "proto.h"

static int NPowerTable;
static int NTransferTable;
static double R8;
static double Norm;
static double klower;
static double r_tophat;
static struct pow_table {
  double logk, logP;
} *PowerTable;
static struct trans_table {
  double logk, logT;
} *TransferTable;

// Set up the Transfer function (normalised to 1)
// ==============================================
void initialize_transferfunction(void) {

#ifndef GAUSSIAN
   // The scale factor at which to compute the nonlinear potential
   FnlTime = 1.0 / (1.0 + Fnl_Redshift);
   DstartFnl = growthD(FnlTime);
#endif

   if(WhichTransfer == 1) read_transfer_table();
}

// Read the Transfer function from an input file (specified in run parameters)
// ========================================================================
void read_transfer_table(void) {

  FILE *fd;
  char buf[500];
  int i;
  double k, t;
  double Tlower;

  sprintf(buf, FileWithInputTransfer);

  if(!(fd = fopen(buf, "r"))) {
    if (ThisTask == 0) printf("\nERROR: Can't read input transfer function  in file '%s'.\n\n", buf);
    FatalError((char *)"power.c", 76);
  }
  fflush(stdout);

  NTransferTable = 0;
  do {
    if(fscanf(fd, " %le %le ", &k, &t) == 2) {
      NTransferTable++;
    } else {
      break;
    }
  } while(1);

  fclose(fd);

  if(ThisTask == 0) {
    printf("Found %d pairs of values in input transfer function table...\n", NTransferTable);
    fflush(stdout);
  }

  TransferTable = (struct trans_table*)malloc(NTransferTable * sizeof(struct trans_table));

  sprintf(buf, FileWithInputTransfer);

  if(!(fd = fopen(buf, "r"))) {
    if (ThisTask == 0) printf("\nERROR: Can't read input transfer function in file '%s'.\n\n", buf);
    FatalError((char *)"power.c", 102);
  }
  fflush(stdout);

  NTransferTable = 0;
  do {
    if(fscanf(fd, " %le %le ", &k, &t) == 2) {
      TransferTable[NTransferTable].logk = log10(k);
      TransferTable[NTransferTable].logT = log10(t);
      NTransferTable++;
    } else {
      break;
    }
  } while(1);

  fclose(fd);

  // Sort in k
  qsort(TransferTable, NTransferTable, sizeof(struct trans_table), compare_transfer_logk); 

  klower = pow(10.0, TransferTable[0].logk);

  if(ThisTask == 0) printf("klower for transfer function is %lf...\n", klower);

  if(TransferTable[0].logk >= -4.6 ) {
    if (ThisTask == 0) {
      printf("\nWARNING: klower may be too large to normalize transfer function.\n");
      printf("         Values outside the input range will be taken to be zero.\n\n");
    }
  }
  if(TransferTable[NTransferTable].logk <= log10(500./8.)) {
    if (ThisTask == 0) {
      printf("\nWARNING: kmax may be too small to normalize transfer function.\n");
      printf("         Values outside the input range will be taken to be zero.\n\n");
    }
  }

  Tlower = TransferTable[0].logT;      
  for(i=0; i < NTransferTable ; i++ ) TransferTable[i].logT -= Tlower;
      
  return;
}

// Frees the memory allocated for the tabulated Transfer function (if necessary)
// =============================================================================
void free_transfertable(void) {
  if(WhichTransfer == 1 && TransferTable != NULL) free(TransferTable);
}

// Comparison routine for qsort to sort tabulated Transfer function by k
// =====================================================================
int compare_transfer_logk(const void *a, const void *b) {
  if(((struct trans_table *) a)->logk < (((struct trans_table *) b)->logk)) return -1;
  if(((struct trans_table *) a)->logk > (((struct trans_table *) b)->logk)) return 1;
  return 0;
}

// Calculate the transfer function at k based on the run parameters
// ================================================================
double TransferFunc(double k) {

  double transfer;

  switch (WhichTransfer) {
    case 1:
      transfer = TransferFunc_Tabulated(k);
      break;
    default:
      transfer = TransferFunc_EH(k); 
      break;
  }

  return transfer;
}

// Calculate the transfer function at k for a tabulated transfer function
// ======================================================================
double TransferFunc_Tabulated(double k) {

  int binlow, binhigh, binmid;
  double logk, logT, T, kold, u, dlogk; 

  kold = k;

  // convert to h/Mpc
  k *= (InputSpectrum_UnitLength_in_cm / UnitLength_in_cm);     

  logk = log10(k);

  if((logk < TransferTable[0].logk) || (logk > TransferTable[NTransferTable - 1].logk)) return 0;

  binlow = 0;
  binhigh = NTransferTable - 1;

  while(binhigh - binlow > 1) {
    binmid = (binhigh + binlow) / 2;
    if(logk < TransferTable[binmid].logk) {
      binhigh = binmid;
    } else {
      binlow = binmid;
    }
  }

  dlogk = TransferTable[binhigh].logk - TransferTable[binlow].logk;

  if(dlogk == 0) FatalError((char *)"power.c", 207);

  u = (logk - TransferTable[binlow].logk) / dlogk;

  logT = (1 - u) * TransferTable[binlow].logT + u * TransferTable[binhigh].logT;

  T = pow(10.0, logT);

  return T;
}

// Set up the power spectrum
// =========================
void initialize_powerspectrum(void) {

  double res;

  // R8 = 8 MPc/h
  R8 = 8 * (3.085678e24 / UnitLength_in_cm);

  // Read power spectrum from input file if requested
  if(WhichSpectrum == 1) read_power_table();

  Norm = 1.0;
  res = TopHatSigma2(R8);
  if(ThisTask == 0 && WhichSpectrum == 1) printf("Normalization of spectrum in file: Sigma8 = %lf...\n",sqrt(res));

  Norm = Sigma8 * Sigma8 / res;
  if(ThisTask == 0) printf("Normalization adjusted to Sigma8=%lf (Normfac=%lf)...\n",Sigma8,Norm);

  // for WhichSpectrum == 0 do not use power spectrum, only transfer function,
  // the file of which is set in the run parameters
  if (WhichSpectrum == 0) Anorm = Norm;

  return;
}

// Read the power spectrum from an input file (specified in run parameters)
// ========================================================================
void read_power_table(void) {

  FILE *fd;
  char buf[500];
  double k, p;
  double kmin,kmax;

  sprintf(buf, FileWithInputSpectrum);

  if(!(fd = fopen(buf, "r"))) {
    if (ThisTask == 0) printf("\nERROR: Can't read input power spectrum in file '%s'.\n\n", buf);
    FatalError((char *)"power.c", 257);
  }
  fflush(stdout);

  NPowerTable = 0;
  do {
    if(fscanf(fd, " %lg %lg ", &k, &p) == 2) {
      NPowerTable++;
    } else {
      break;
    }
  } while(1);

  fclose(fd);

  if(ThisTask == 0) {
    printf("Task %d Found %d pairs of values in input power spectrum table...\n", ThisTask, NPowerTable);
    fflush(stdout);
  }

  PowerTable = (struct pow_table *)malloc(NPowerTable * sizeof(struct pow_table));

  sprintf(buf, FileWithInputSpectrum);

  if(!(fd = fopen(buf, "r"))) {
    if (ThisTask == 0) printf("\nERROR: Can't read input power spectrum in file '%s'.\n\n", buf);
    FatalError((char *)"power.c", 283);
  }
  fflush(stdout);

  NPowerTable = 0;
  kmin = 2 * PI / 0.0001 ;  // 10000 h/Mpc
  kmax = 2 * PI / 10000.0;  // 0.0001 h/Mpc
  do {
    if(fscanf(fd, " %lg %lg ", &k, &p) == 2) {
      PowerTable[NPowerTable].logk = log10(k);
      PowerTable[NPowerTable].logP = log10(p);
      NPowerTable++;

      k /= (InputSpectrum_UnitLength_in_cm/3.085678e24); // convert to h/Mpc

      if (k < kmin) kmin = k;
      if (k > kmax) kmax = k;
    } else {
      break;
    }
  } while(1);

  fclose(fd);

  //check if there is sufficient k-coverage 
  double k_Nyquist = PI*Nsample/(Box*UnitLength_in_cm/3.085678e24);
  double k_fundamental = 2.0*PI/(Box*UnitLength_in_cm/3.085678e24);

  if((kmin > k_fundamental) || (kmax < k_Nyquist)) {
    if (ThisTask == 0) printf("\nERROR: [kmin, kmax] = [%lf,%lf] h/Mpc are not sufficient to cover [k_fundamental, k_nyquist] = [%lf,%lf] h/Mpc.\n\n",kmin,kmax,k_fundamental,k_Nyquist);
    FatalError((char *)"power.c", 313);
  }

  // Sort by k
  qsort(PowerTable, NPowerTable, sizeof(struct pow_table), compare_logk);

  klower = pow(10.0, PowerTable[0].logk) * 1.0000001;

  if(ThisTask == 0) printf("klower for power spectrum is %f...\n",klower);

  if(PowerTable[0].logk >= -4.6 ) {
    if (ThisTask == 0) {
      printf("\nWARNING: klower may be too large to normalize power.\n");
      printf("         Values outside the input range will be taken to be zero.\n\n");
    }
  }
  if(PowerTable[NTransferTable].logk <= log10(500./8.)) {
    if (ThisTask == 0) {
      printf("\nWARNING: kmax may be too small to normalize the power spectrum.\n");
      printf("         Values outside the input range will be taken to be zero.\n\n");
    }
  }

  return;
}

// Frees the memory allocated for the tabulated Power Spectrum (if necessary)
// =============================================================================
void free_powertable(void) {
  if(WhichSpectrum == 1 && PowerTable != NULL) free(PowerTable);
}

// Comparison routine for qsort to sort tabulated power spectrum by k
// ==================================================================
int compare_logk(const void *a, const void *b) {
  if(((struct pow_table *) a)->logk < (((struct pow_table *) b)->logk)) return -1;
  if(((struct pow_table *) a)->logk > (((struct pow_table *) b)->logk)) return +1;
  return 0;
}

// Calculate the power spectrum at k based on the run parameters
// =============================================================
double PowerSpec(double k) {

  double power;

  switch (WhichSpectrum) {
    case 0:
      power = Norm * pow(k, PrimordialIndex)  * TransferFunc(k) * TransferFunc(k);
      break;
    case 1:
      power = PowerSpec_Tabulated(k);
      break;
    default:
      power = PowerSpec_EH(k);
      break;
  }

  return power;
}

// Calculate the power spectrum at k for a tabulated power spectrum
// ================================================================
double PowerSpec_Tabulated(double k) {

  int binlow, binhigh, binmid;
  double logk, logP, P, kold, u, dlogk, Delta2;

  kold = k;

  // convert k to h/Mpc
  k *= (InputSpectrum_UnitLength_in_cm / UnitLength_in_cm);	

  logk = log10(k);

  if((logk < PowerTable[0].logk) || (logk > PowerTable[NPowerTable - 1].logk)) return 0;

  binlow = 0;
  binhigh = NPowerTable - 1;

  while(binhigh - binlow > 1) {
    binmid = (binhigh + binlow) / 2;
    if(logk < PowerTable[binmid].logk) {
      binhigh = binmid;
    } else {
      binlow = binmid;
    }
  }

  dlogk = PowerTable[binhigh].logk - PowerTable[binlow].logk;

  if(dlogk == 0) FatalError((char *)"power.c", 404);

  u = (logk - PowerTable[binlow].logk) / dlogk;

  logP = (1 - u) * PowerTable[binlow].logP + u * PowerTable[binhigh].logP;

  Delta2 = pow(10.0, logP);

  P = Norm * Delta2;

  return P;
}

// Calculate the power spectrum at k for a power spectrum with Eisenstein & Hu parameterisation
// ============================================================================================
double PowerSpec_EH(double k) {
  return Norm * pow(k, PrimordialIndex) * pow(TransferFunc_EH(k), 2);
}

// Fitted analytic expressions for Eisenstein & Hu Transfer function from Martin White
// ===================================================================================
double TransferFunc_EH(double k) {
  double q, theta, ommh2, a, s, gamma, L0, C0;
  double tmp;
  double omegam, ombh2, hubble;

  // other input parameters
  hubble = HubbleParam;

  omegam = Omega;
  ombh2 = OmegaBaryon * HubbleParam * HubbleParam;

  if(OmegaBaryon == 0)
    ombh2 = 0.04 * HubbleParam * HubbleParam;

  // convert k to h/Mpc
  k *= (3.085678e24 / UnitLength_in_cm);	

  theta = 2.728 / 2.7;
  ommh2 = omegam * hubble * hubble;
  s = 44.5 * log(9.83 / ommh2) / sqrt(1. + 10. * exp(0.75 * log(ombh2))) * hubble;
  a = 1. - 0.328 * log(431. * ommh2) * ombh2 / ommh2
    + 0.380 * log(22.3 * ommh2) * (ombh2 / ommh2) * (ombh2 / ommh2);
  gamma = a + (1. - a) / (1. + exp(4 * log(0.43 * k * s)));
  gamma *= omegam * hubble;
  q = k * theta * theta / gamma;
  L0 = log(2. * exp(1.) + 1.8 * q);
  C0 = 14.2 + 731. / (1. + 62.5 * q);
  tmp = L0 / (L0 + C0 * q * q);
  return (tmp);
}

// Return the integral over a Tophat profile
// ==========================================
double TopHatSigma2(double R) {
  double alpha=0.0;
  double result, error;

  r_tophat = R;     // NB: 500/R is chosen as the integration boundary (infinity)

  gsl_function F;
  gsl_integration_workspace * w = gsl_integration_workspace_alloc (100000);
     
  F.function = &sigma2_int;
  F.params = &alpha;
     
  gsl_integration_qag(&F,0,500.0/R,1e-6,1e-4,100000,GSL_INTEG_GAUSS41,w,&result,&error); 
     
  gsl_integration_workspace_free (w);
      
  return result/(2.0*PI*PI);	
}

// The Tophat profile
// ==================
double sigma2_int(double k, void * params) {

  double kr, kr3, kr2, w, x;

  kr = r_tophat * k;
  kr2 = kr * kr;
  kr3 = kr2 * kr;

  if(kr < 1e-8) return 0;

  w = 3 * (sin(kr) / kr3 - cos(kr) / kr2);
  x = k * k * w * w * PowerSpec(k);

  return x;
}
