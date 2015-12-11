/* ==========================================================================*/
/*   Version 1.2.             Cullan Howlett & Marc Manera,                  */
/*   Copyright (c) 2015       Institute of Cosmology and Gravitation         */
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

/* =====================================================================*/
/* This file contains all the prototypes for function used in the code. */
/* v1.2: Updated prototypes for Output and Output_Info corresponding to */
/*       z=9 bug fix in main.c                                          */
/* =====================================================================*/

// main.c
void Output_Info(double A, double Z);
void Output(double A, double Z, double Dv, double Dv2);
void Kick(double AI, double AF, double A, double Di, double Di2);
void Drift(double A, double AFF, double AF, double Di, double Di2);
#ifdef TIMING
void Output_Timing(void);
#endif

// cosmo.c
double u32(double a);
double u52(double a);
double QdD1da(double a);
double QdD2da(double a);
double Hubblea(double a);
double Qfactor(double a);
double growthD(double a);
double growthD2(double a);
double TimeCOLA(double a);
double growthDtemp(double a);
double growthD2temp(double a);
double QdTimeCOLAda(double a);
double KickStd(double ai,double af);
double DriftStd(double ai,double af);
double u32func(double a, void * params);
double u52func(double a, void * params);
double KickStdfunc(double a, void * params);
double DriftStdfunc(double a, void * params);
double DriftCOLAfunc(double a, void * params);
double KickCOLA(double ai,double af,double ac);
double DriftCOLA(double ai,double af,double ac);
double growthDtempfunc(double a, void * params);

// auxPM.c
void Forces(void);
void PtoMesh(void);
void MtoParticles(void);
void MoveParticles(void);
void GetDisplacements(void);
void FatalError(char * filename, int linenum);
#if (MEMORY_MODE || SINGLE_PRECISION)
float periodic_wrap(float x);
#else
double periodic_wrap(double x);
#endif
size_t my_fread(void *ptr, size_t size, size_t nmemb, FILE * stream);
size_t my_fwrite(void *ptr, size_t size, size_t nmemb, FILE * stream);

// read_param.c
void read_outputs(void);
void read_parameterfile(char * fname);
int sort_redshift(const void * Item1, const void * Item2);

// kernel.c
#ifdef GENERIC_FNL
void read_kernel_table(void);
#endif

// lightcone.c
#ifdef LIGHTCONE
void set_lightcone(void);
void Output_Info_Lightcone(void);
void Output_Lightcone(unsigned int * pc, unsigned int blockmaxlen, float * block);
void flag_replicates(double Rcomov_old, double Rcomov_new, double boundary);
void Drift_Lightcone(double A, double AFF, double AF, double Di, double Di2);
double nearest_dist(double px, double py, double ix, double iy, double jx, double jy, double boundary);
#endif

// 2LPT.c
void set_units(void);
void initialize_ffts(void);
void initialize_parts(void);
void displacement_fields(void);

// power.c
void print_spec(void);
void free_powertable(void);
void read_power_table(void);
void free_transfertable(void);
void read_transfer_table(void);
void initialize_powerspectrum(void);
void initialize_transferfunction(void);
int compare_logk(const void *a, const void *b);
int compare_transfer_logk(const void *a, const void *b);
double fnl(double x);
double TransferFunc(double k);
double PowerSpec(double kmag);
double TopHatSigma2(double R);
double PowerSpec_EH(double k);
double TransferFunc_EH(double k);
double PowerSpec_Tabulated(double k);
double PowerSpec_Efstathiou(double k);
double TransferFunc_Tabulated(double k);
double sigma2_int(double k, void * params);
