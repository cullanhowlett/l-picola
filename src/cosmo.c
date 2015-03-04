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

/* ======================================================*/
/* This file contains all the cosmology related routines.*/
/* ======================================================*/

#include "vars.h"
#include "proto.h"

// The Hubble expansion
// ====================
double Hubblea(double a) { 
  return sqrt(Omega/(a*a*a)+(1.0-Omega));
  //return sqrt(Omega/(a*a*a)+OmegaLambda+(1.0-Omega-OmegaLambda)/(a*a));
}

// The Q factor from Tassev et. al., 2013 (Q = a^3 H(a)/H0)
// ========================================================
double Qfactor(double a) { 
  return Hubblea(a)*a*a*a;
  //return sqrt(Omega/(a*a*a)+OmegaLambda+(1.0-Omega-OmegaLambda)/(a*a))*a*a*a;
}

// Normalised growth factor for LCDM
// =================================
double growthD(double a){ 
  return Hubblea(a)*growthDtemp(a)/growthDtemp(1.0);
}

// Growth factor integral
// ======================
double growthDtemp(double a) {

  double alpha = 0.0;
  double result, error;

  gsl_function F;
  gsl_integration_workspace * w = gsl_integration_workspace_alloc (5000);
     
  F.function = &growthDtempfunc;
  F.params = &alpha;
     
  gsl_integration_qag(&F,0.0,a,0,1e-5,5000,6,w,&result,&error); 
      
  gsl_integration_workspace_free (w);
     
  return result;
}

// Growth factor integrand
// =======================
double growthDtempfunc(double a, void * params) {       
  return 1.0/(pow(a*Hubblea(a),3.0));
}

// Normalised Second order growth factor
// =====================================
double growthD2(double a) { 
  return growthD2temp(a)/growthD2temp(1.0);
}

// Second order growth factor
// ==========================
double growthD2temp(double a) {
  double D = growthD(a);
  double Ha = Hubblea(a);
  double Oma = Omega/(a*a*a*Ha*Ha);
  return D*D*pow(Oma,-1.0/143.0);
}

// returns Q*d(D_{1})/da, where Q=Qfactor(a) 
// =========================================
double QdD1da(double a) {
  return (1.0/Hubblea(a))*((1.0/growthDtemp(1.0))-(((3.0*Omega)/(2.0*a))*growthD(a)));
  //return (1.0/Hubblea(a))*((1.0/growthDtemp(1.0))-((((3.0*Omega)/(2.0*a))+(1.0-Omega-OmegaLambda))*growthD(a)));
}

// returns Q*d(D_{2})/da, where Q=Qfactor(a) 
// =========================================
double QdD2da(double a) {
  double D1 = growthD(a); 
  double D2 = growthD2(a);
  double Dv = QdD1da(a);
  return (2.0*D2*Dv/D1)+((3.0*D2*(1.0-Omega)*a*a)/(143.0*Hubblea(a)));
}

// Functions for COLA modified time-stepping (used when DeltaA=0,1)
// ================================================================
double DriftCOLA(double ai,double af,double ac) {

  double alpha=0.0;
  double result, error;

  gsl_function F;
  gsl_integration_workspace * w = gsl_integration_workspace_alloc (5000);
     
  F.function = &DriftCOLAfunc;
  F.params = &alpha;
     
  gsl_integration_qag(&F,(double)ai,(double)af,0,1e-5,5000,6,w, &result,&error); 
     
  gsl_integration_workspace_free(w);
      
  return result/TimeCOLA(ac);
}
     
double DriftCOLAfunc(double a, void * params) {     
  return TimeCOLA(a)/Qfactor(a);
}

double KickCOLA(double ai,double af,double ac) {
  return ac*(TimeCOLA(af)-TimeCOLA(ai))/QdTimeCOLAda(ac);
}

// This is the function u(t) in Tassev et. al, 2013
// ================================================
double TimeCOLA(double a) { 
  return pow(a,nLPT);      
}

// This is the derivative of u(t) in Tassev et. al, 2013 (multiplied by Q(a))
// ==========================================================================
double QdTimeCOLAda(double a) { 
  return Qfactor(a)*nLPT*pow(a,nLPT-1.0);
}

// Functions for Quinn et al time-stepping (used when DeltaA=1,2)
// ==============================================================
double DriftStd(double ai,double af) {
       
  double alpha = 0.0;
  double result, error;

  gsl_function F;
  gsl_integration_workspace * w = gsl_integration_workspace_alloc (5000);
     
  F.function = &DriftStdfunc;
  F.params = &alpha;
     
  gsl_integration_qag(&F,ai,af,0,1e-5,5000,6,w,&result,&error); 
      
  gsl_integration_workspace_free (w);
  
  return result;
}

double DriftStdfunc(double a, void * params) {
  return 1.0/Qfactor(a);
}
     
double KickStd(double ai,double af) {

  double alpha = 0.0;
  double result, error;

  gsl_function F;
  gsl_integration_workspace * w = gsl_integration_workspace_alloc (5000);
     
  F.function = &KickStdfunc;
  F.params = &alpha;
     
  gsl_integration_qag (&F,ai,af,0,1e-5,5000,6,w,&result,&error); 
      
  gsl_integration_workspace_free (w);
     
  return result;
}

double KickStdfunc(double a, void * params) {       
  return a/Qfactor(a);
}