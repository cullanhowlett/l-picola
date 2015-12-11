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

/* ======================================================================================*/
/* This file contains most of the routines for calculating the ZA and 2LPT displacements.*/
/* v1.2: Edited print statements to output in displacements in correct units             */
/* ======================================================================================*/

#include "vars.h"
#include "proto.h"

// Set some units
// ==============
void set_units(void) {

  UnitTime_in_s = UnitLength_in_cm / UnitVelocity_in_cm_per_s;
  G = GRAVITY / pow(UnitLength_in_cm, 3) * UnitMass_in_g * pow(UnitTime_in_s, 2);
  Hubble = HUBBLE * UnitTime_in_s;

  return;
}

// Set up the sizes for arrays and fft routines
// ============================================
void initialize_ffts(void) {

  int i;
  int * Slab_to_task_local;

#ifdef SINGLE_PRECISION
  alloc_local = fftwf_mpi_local_size_3d(Nmesh, Nmesh, Nmesh/2+1, MPI_COMM_WORLD, &Local_nx, &Local_x_start);
#else
  alloc_local = fftw_mpi_local_size_3d(Nmesh, Nmesh, Nmesh/2+1, MPI_COMM_WORLD, &Local_nx, &Local_x_start);
#endif

  Local_nx_table = (int *)malloc(sizeof(int) * NTask);

  MPI_Allgather(&Local_nx, 1, MPI_INT, Local_nx_table, 1, MPI_INT, MPI_COMM_WORLD);

  if(ThisTask == 0) {
    printf("\nLocal nx\n---------------------\n");
    for(i = 0; i < NTask; i++) printf("Task = %d: Local_nx = %d\n", i, Local_nx_table[i]);
    printf("---------------------\n");
    fflush(stdout);
  }

  // Set the neighbouring tasks
  if (Local_nx == 0) {
    LeftTask = MPI_PROC_NULL;
    RightTask = MPI_PROC_NULL;
  } else {
    LeftTask = ThisTask;
    do {
      LeftTask--;
      if(LeftTask < 0) LeftTask = NTask - 1;
    } while(Local_nx_table[LeftTask] == 0);
      
    RightTask = ThisTask;
    do {
      RightTask++;
      if(RightTask >= NTask) RightTask = 0;
    } while(Local_nx_table[RightTask] == 0);
  }

  // Let each processor know which parts of the fourier grids they all have
  Slab_to_task       = (int *)malloc(sizeof(int) * Nmesh);
  Slab_to_task_local = (int *)malloc(sizeof(int) * Nmesh);

  for(i = 0; i < Nmesh; i++)    Slab_to_task_local[i] = 0;
  for(i = 0; i < Local_nx; i++) Slab_to_task_local[Local_x_start + i] = ThisTask;

  MPI_Allreduce(Slab_to_task_local, Slab_to_task, Nmesh, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  // Add an additional plane
  alloc_slice = Nmesh*(Nmesh/2+1);
  last_slice = Local_nx*alloc_slice;
  Total_size = alloc_local+alloc_slice;

  free(Slab_to_task_local);

  return;
}

// Work out which particles each Task will have
// ============================================
void initialize_parts(void) {

  int i, slab;
  int * Part_to_task_local;

  Local_np = 0;
  Local_p_start = Nsample;
  for (i = 0; i < Nsample; i++) {
    slab = (int)((double)(i*Nmesh)/(double)Nsample);
    if (Slab_to_task[slab] == ThisTask) {
      Local_np++;
      if (i < Local_p_start) Local_p_start = i;
    }
  }

  Local_np_table = (int *)malloc(sizeof(int) * NTask);

  MPI_Allgather(&Local_np, 1, MPI_INT, Local_np_table, 1, MPI_INT, MPI_COMM_WORLD);

  NumPart = Local_np*Nsample*Nsample;
  TotNumPart = ((unsigned long long) Nsample) * ((unsigned long long) Nsample) *  ((unsigned long long) Nsample);

  if(ThisTask == 0) {
    printf("\nParticles\n---------------------\n");
    for(i = 0; i < NTask; i++) printf("Task = %d: Particles = %u\n", i, (unsigned int)Local_np_table[i]*(unsigned int)Nsample*(unsigned int)Nsample);
    printf("\n----------------------\n");
    printf("Total number of particles = %llu\n\n", TotNumPart);
    fflush(stdout);
  }

  for(i = 0; i < NTask; i++) {
    if(Local_np_table[i] > 0) NTaskWithN++;
  }

  // Let each processor know which parts of the particle grids they all have
  Part_to_task       = (int *)malloc(sizeof(int) * Nsample);
  Part_to_task_local = (int *)malloc(sizeof(int) * Nsample);

  for(i = 0; i < Nsample; i++)  Part_to_task_local[i] = 0;
  for(i = 0; i < Local_np; i++) Part_to_task_local[Local_p_start + i] = ThisTask;

  MPI_Allreduce(Part_to_task_local, Part_to_task, Nsample, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  free(Part_to_task_local);

  return;
}

// This is the largest routine in the code and is used to generate the 2LPT initial conditions. 
// It is adapted from the 2LPTic code provided by Marc Manera. 
//
// NOTE: We ALWAYS compute the power spectrum and displacements at redshift 0.0, regardless of 
// whether we use COLA or PM methods. These are then modified accordingly when the particle data is initialised in main.c
// ======================================================================================================================
void displacement_fields(void) {

  // There are a LOT of parameters here:
  // ===================================
  gsl_rng * random_generator;
  int i, j, k, n, m, p, ii, jj, kk, axes;
  unsigned int * seedtable, q, coord, bytes;
  unsigned long long nmesh3; 
  float_kind *(digrad[6]);
  float_kind *(disp[3]), *(disp2[3]);
  double u, v, w;
  double phase, ampl;
  double kvec[3], kmag, kmag2;
  double sumdis[3], sumdis2[3];
  double f1, f2, f3, f4, f5, f6, f7, f8;
  double dis[3], dis2[3], maxdisp, max_disp_glob;
  complex_kind *(cdigrad[6]);
  complex_kind *(cdisp[3]), *(cdisp2[3]);
  plan_kind Forward_plan, Inverse_plan;

// General parameters for any gaussianity/non-gaussianity
// ======================================================
#ifdef GAUSSIAN
  double delta;
  double p_of_k;
#else
  float_kind *(pot);
  double t_of_k;
  double twb, phig, Beta;
  complex_kind *(cpot);
#ifdef TIMING
  double startcpu, endcpu;
  double startwall, endwall;
#endif
#endif

// Parameters for generic non-gaussianity
// ======================================                                                
#ifdef GENERIC_FNL
  int ikernel;
  double kmag0, kmagA, kmagB;
  double kerCoef, ker0, kerA, kerB;
  float_kind *(ppA), *(ppB), *(ppC);                                        
  complex_kind *(cppA), *(cppB), *(cppC);                                  
#endif

// Parameters for equilateral or orthogonal fnl non-gaussianity
// ============================================================
#if (EQUIL_FNL || ORTHO_FNL)
  float_kind *(partpot);
  float_kind *(p1p2p3sym), *(p1p2p3sca), *(p1p2p3nab), *(p1p2p3tre);                                                 
  complex_kind *(cpartpot);
  complex_kind *(cp1p2p3sym), *(cp1p2p3sca), *(cp1p2p3nab), *(cp1p2p3tre);
#endif                                               

  if(ThisTask == 0) {
    printf("Computing displacement fields...\n");
    fflush(stdout);
  }

  maxdisp = 0;

  random_generator = gsl_rng_alloc(gsl_rng_ranlxd1);

  gsl_rng_set(random_generator, Seed);

  if(!(seedtable = (unsigned int *)malloc(Nmesh * Nmesh * sizeof(unsigned int)))) FatalError((char *)"2LPT.c", 223);

  for(i = 0; i < Nmesh / 2; i++) {
    for(j = 0; j < i; j++)     seedtable[i * Nmesh + j] = (unsigned int)(0x7fffffff * gsl_rng_uniform(random_generator));
    for(j = 0; j < i + 1; j++) seedtable[j * Nmesh + i] = (unsigned int)(0x7fffffff * gsl_rng_uniform(random_generator));
    for(j = 0; j < i; j++)     seedtable[(Nmesh - 1 - i) * Nmesh + j] = (unsigned int)(0x7fffffff * gsl_rng_uniform(random_generator));
    for(j = 0; j < i + 1; j++) seedtable[(Nmesh - 1 - j) * Nmesh + i] = (unsigned int)(0x7fffffff * gsl_rng_uniform(random_generator));
    for(j = 0; j < i; j++)     seedtable[i * Nmesh + (Nmesh - 1 - j)] = (unsigned int)(0x7fffffff * gsl_rng_uniform(random_generator));
    for(j = 0; j < i + 1; j++) seedtable[j * Nmesh + (Nmesh - 1 - i)] = (unsigned int)(0x7fffffff * gsl_rng_uniform(random_generator));
    for(j = 0; j < i; j++)     seedtable[(Nmesh - 1 - i) * Nmesh + (Nmesh - 1 - j)] = (unsigned int)(0x7fffffff * gsl_rng_uniform(random_generator));
    for(j = 0; j < i + 1; j++) seedtable[(Nmesh - 1 - j) * Nmesh + (Nmesh - 1 - i)] = (unsigned int)(0x7fffffff * gsl_rng_uniform(random_generator));
  }

// Gaussian initial conditions
// ===========================
#ifdef GAUSSIAN

  for(axes=0,bytes=0; axes < 3; axes++) {
    cdisp[axes] = (complex_kind *) malloc(bytes += sizeof(complex_kind) * Total_size);
    disp[axes] = (float_kind *)cdisp[axes];
  }

  if(ThisTask == 0) {
    printf("Starting gaussian calculations...\n");
    fflush(stdout);
  }

  // First, clean the array
  for(i = 0; i < Local_nx; i++) {
    for(j = 0; j < Nmesh; j++)     {
      for(k = 0; k <= Nmesh / 2; k++) {
        for(axes = 0; axes < 3; axes++)  {
          coord = (i * Nmesh + j) * (Nmesh / 2 + 1) + k;
          cdisp[axes][coord][0] = 0.0;
          cdisp[axes][coord][1] = 0.0;
    }
      }
    }
  }

  for(i = 0; i < Nmesh; i++) {
    ii = Nmesh - i;
    if(ii == Nmesh) ii = 0;
    if((i >= Local_x_start && i < (Local_x_start + Local_nx)) ||
       (ii >= Local_x_start && ii < (Local_x_start + Local_nx))) {
  
      for(j = 0; j < Nmesh; j++) {
        gsl_rng_set(random_generator, seedtable[i * Nmesh + j]); 

        for(k = 0; k < Nmesh / 2; k++) {
          phase = gsl_rng_uniform(random_generator) * 2 * PI;
          do {
            ampl = gsl_rng_uniform(random_generator);
          } while(ampl == 0);
 
          if(i == Nmesh / 2 || j == Nmesh / 2 || k == Nmesh / 2) continue;
          if(i == 0 && j == 0 && k == 0) continue;
 
          if(i < Nmesh / 2) {
            kvec[0] = i * 2 * PI / Box;
          } else {
            kvec[0] = -(Nmesh - i) * 2 * PI / Box;
          }
 
          if(j < Nmesh / 2) {
            kvec[1] = j * 2 * PI / Box;
          } else {
            kvec[1] = -(Nmesh - j) * 2 * PI / Box;
          }
 
          if(k < Nmesh / 2) {
	        kvec[2] = k * 2 * PI / Box;
          } else {
	        kvec[2] = -(Nmesh - k) * 2 * PI / Box;
          }
  
          kmag2 = kvec[0] * kvec[0] + kvec[1] * kvec[1] + kvec[2] * kvec[2];
          kmag  = sqrt(kmag2);
 
          if(SphereMode == 1) {
            if(kmag * Box / (2 * PI) > Nsample / 2) continue; // select a sphere in k-space
          } else {
            if(fabs(kvec[0]) * Box / (2 * PI) > Nsample / 2) continue;
            if(fabs(kvec[1]) * Box / (2 * PI) > Nsample / 2) continue;
            if(fabs(kvec[2]) * Box / (2 * PI) > Nsample / 2) continue;
          }

          p_of_k  = PowerSpec(kmag);
          p_of_k *= -log(ampl);
 
          delta = pow(Box,-1.5) * sqrt(p_of_k);  // keep at redshift 0.0

          if(k > 0) {
            if(i >= Local_x_start && i < (Local_x_start + Local_nx)) {
              for(axes = 0; axes < 3; axes++) {
                coord = ((i - Local_x_start) * Nmesh + j) * (Nmesh / 2 + 1) + k;
                cdisp[axes][coord][0] = -kvec[axes] / kmag2 * delta * sin(phase);
                cdisp[axes][coord][1] =  kvec[axes] / kmag2 * delta * cos(phase);
              }
            }
          } else {       // k=0 plane needs special treatment
            if(i == 0) {
              if(j >= Nmesh / 2) {
                continue;
              } else {
                if(i >= Local_x_start && i < (Local_x_start + Local_nx)) {
                  jj = Nmesh - j;  // note: j!=0 surely holds at this point
                  for(axes = 0; axes < 3; axes++) {
                    coord = ((i - Local_x_start) * Nmesh + j) * (Nmesh / 2 + 1) + k;
                    cdisp[axes][coord][0]  = -kvec[axes] / kmag2 * delta * sin(phase);
                    cdisp[axes][coord][1]  =  kvec[axes] / kmag2 * delta * cos(phase);
                    coord = ((i - Local_x_start) * Nmesh + jj) * (Nmesh / 2 + 1) + k;
                    cdisp[axes][coord][0] = -kvec[axes] / kmag2 * delta * sin(phase);
                    cdisp[axes][coord][1] = -kvec[axes] / kmag2 * delta * cos(phase);
                  }
                }
              }
            } else {    // here comes i!=0 : conjugate can be on other processor!
              if(i >= Nmesh / 2) {
                continue;
              } else {
                ii = Nmesh - i;
                jj = Nmesh - j;
                if(ii == Nmesh) ii = 0;
                if(jj == Nmesh) jj = 0;
                if(i >= Local_x_start && i < (Local_x_start + Local_nx)) {
                for(axes = 0; axes < 3; axes++) {
                    coord = ((i - Local_x_start) * Nmesh + j) * (Nmesh / 2 + 1) + k;
                    cdisp[axes][coord][0] = -kvec[axes] / kmag2 * delta * sin(phase);
                    cdisp[axes][coord][1] =  kvec[axes] / kmag2 * delta * cos(phase);
                  }
                }		  
                if(ii >= Local_x_start && ii < (Local_x_start + Local_nx)) {
                  for(axes = 0; axes < 3; axes++) {
                    coord = ((ii - Local_x_start) * Nmesh + jj) * (Nmesh / 2 + 1) + k;
                    cdisp[axes][coord][0] = -kvec[axes] / kmag2 * delta * sin(phase);
                    cdisp[axes][coord][1] = -kvec[axes] / kmag2 * delta * cos(phase);
                  }
                }
              }
            }
          }
        }
      }
    }
  }

// Non-gaussian initial conditions
// ===============================
#else

#ifdef TIMING
  startcpu = (double)clock();
  startwall = MPI_Wtime();
#endif  

  if(ThisTask == 0) {
    printf("Starting non-gaussian calculations...\n");
    fflush(stdout);
  }

  // Initialize
  bytes=0;
  cpot = (complex_kind *) malloc(bytes += sizeof(complex_kind) * Total_size);
  pot = (float_kind *) cpot;

  // first, clean the cpot array
  for(i = 0; i < Local_nx; i++) {
    for(j = 0; j < Nmesh; j++)     {
      for(k = 0; k <= Nmesh / 2; k++) {
        coord = (i * Nmesh + j) * (Nmesh / 2 + 1) + k;
        cpot[coord][0] = 0;
        cpot[coord][1] = 0;
      }
    }
  }

  // Ho in units of h/Mpc and c=1, i.e., internal units so far
  // Beta = 3/2 H(z)^2 a^2 Om(a) = 3/2 Ho^2 Om0 / a 
  Beta = 1.5 * Omega / FnlTime / (2998. * 2998. );

  for(i = 0; i < Nmesh; i++) {
    ii = Nmesh - i;
    if(ii == Nmesh) ii = 0;
    if((i >= Local_x_start && i < (Local_x_start + Local_nx)) ||
      (ii >= Local_x_start && ii < (Local_x_start + Local_nx))) {

      for(j = 0; j < Nmesh; j++) {
        gsl_rng_set(random_generator, seedtable[i * Nmesh + j]);

        for(k = 0; k < Nmesh / 2; k++) {
          phase = gsl_rng_uniform(random_generator) * 2* PI;
          do {
            ampl = gsl_rng_uniform(random_generator);
          } while(ampl == 0);

          if(i == Nmesh / 2 || j == Nmesh / 2 || k == Nmesh / 2) continue; 
          if(i == 0 && j == 0 && k == 0) continue;

          if(i < Nmesh / 2) {
            kvec[0] = i * 2 * PI / Box;
          } else {
            kvec[0] = -(Nmesh - i) * 2 * PI / Box;
          }
          
          if(j < Nmesh / 2) {
            kvec[1] = j * 2 * PI / Box;
          } else {
            kvec[1] = -(Nmesh - j) * 2 * PI / Box;
          }

          if(k < Nmesh / 2) {
            kvec[2] = k * 2 * PI / Box;
          } else {
            kvec[2] = -(Nmesh - k) * 2 * PI / Box;
          }
           
          kmag2 = kvec[0] * kvec[0] + kvec[1] * kvec[1] + kvec[2] * kvec[2];
          kmag = sqrt(kmag2);

          if(SphereMode == 1) {
            if(kmag * Box / (2 * PI) > Nsample / 2) continue;  // select a sphere in k-space
          } else {
            if(fabs(kvec[0]) * Box / (2 * PI) > Nsample / 2) continue;
            if(fabs(kvec[1]) * Box / (2 * PI) > Nsample / 2) continue;
            if(fabs(kvec[2]) * Box / (2 * PI) > Nsample / 2) continue;
          }

          phig = -log(ampl) * Anorm * pow(kmag, PrimordialIndex);            // initial normalized power
          phig = sqrt(phig) * pow(Box,-1.5) * Beta * DstartFnl / kmag2;      // amplitude of the initial gaussian potential
               
          if(k > 0) {
            if(i >= Local_x_start && i < (Local_x_start + Local_nx)) {
              coord = ((i - Local_x_start) * Nmesh + j) * (Nmesh / 2 + 1) + k;
              cpot[coord][0] =  phig * sin(phase);
              cpot[coord][1] = -phig * cos(phase);
            }
          } else {  // k=0 plane needs special treatment
            if(i == 0) {
              if(j >= Nmesh / 2) {
                continue;
              } else {
                if(i >= Local_x_start && i < (Local_x_start + Local_nx)) {
                  jj = Nmesh - j;   // note: j!=0 surely holds at this point

                  coord = ((i - Local_x_start) * Nmesh + j) * (Nmesh / 2 + 1) + k;
                  cpot[coord][0] =  phig * sin(phase);
                  cpot[coord][1] = - phig * cos(phase);

                  coord = ((i - Local_x_start) * Nmesh + jj) * (Nmesh / 2 + 1) + k; 
                  cpot[coord][0] = phig * sin(phase);
                  cpot[coord][1] = phig * cos(phase);
                }
              }
            } else {  // here comes i!=0 : conjugate can be on other processor!
              if(i >= Nmesh / 2) {
                  continue;
              } else {
                ii = Nmesh - i;
                jj = Nmesh - j;
                if(ii == Nmesh) ii = 0;
                if(jj == Nmesh) jj = 0;
                if(i >= Local_x_start && i < (Local_x_start + Local_nx)) {
                  coord = ((i - Local_x_start) * Nmesh + j) * (Nmesh / 2 + 1) + k;
                  cpot[coord][0] =  phig * sin(phase);
                  cpot[coord][1] = -phig * cos(phase);
                }
                if(ii >= Local_x_start && ii < (Local_x_start + Local_nx)) {
                  coord = ((ii - Local_x_start) * Nmesh + jj) * (Nmesh / 2 + 1) + k;
                  cpot[coord][0] = phig * sin(phase);
                  cpot[coord][1] = phig * cos(phase);
                }
              }
            }
          }
        }
      }
    }
  }

 // For non-local models it is important to keep all factors of SQRT(-1) as done below
 // Notice also that there is a minus to convert from Bardeen to gravitational potential

// Generic non-gaussian potential
// ==============================
#ifdef GENERIC_FNL

  cppA= (complex_kind *) calloc(Total_size, sizeof(complex_kind));
  ppA = (float_kind *) cppA;

  cppB= (complex_kind *) calloc(Total_size, sizeof(complex_kind));
  ppB = (float_kind *) cppB;

  cppC= (complex_kind *) calloc(Total_size, sizeof(complex_kind));
  ppC = (float_kind *) cppC;

  read_kernel_table();

  for(ikernel=0; ikernel < NKernelTable; ikernel++) {
    kerCoef=KernelTable[ikernel].Coef;
    ker0=KernelTable[ikernel].ker0;
    kerA=KernelTable[ikernel].kerA;
    kerB=KernelTable[ikernel].kerB;
        
    if (ThisTask == 0) printf("\nKernel table line %d\n-------------------\n", ikernel);

    if(fabs(ker0+kerA+kerB) < 0.000000001) {
      if(ker0+kerA+kerB != 0.0) ker0 = - ( kerA + kerB );
      if(ThisTask == 0) printf("Adjusting ker0 = - (kerA + kerB), ker0=%f, kerA=%f, kerB=%f\n", ker0,kerA,kerB);
    } else {
      if(ThisTask == 0) printf("\nERROR: ker0 + kerA + kerB does not equal 0\n"); 
      FatalError((char *)"2LPT.c", 534);
    }

    if(ThisTask == 0) printf("Values: %lf %lf %lf %lf\n",kerCoef,ker0,kerA,kerB);

    for(ii = 0; ii < Local_nx; ii++) {
      for(j = 0; j < Nmesh; j++)       {
        for(k = 0; k <= Nmesh / 2 ; k++) {
          i = ii + Local_x_start;
          coord = (ii * Nmesh + j) * (Nmesh / 2 + 1) + k;

          cppA[coord][0] = 0.0;
          cppA[coord][1] = 0.0;
          cppB[coord][0] = 0.0;
          cppB[coord][1] = 0.0;

          if((i == 0) && (j == 0) && (k == 0)) continue;

          if(i < Nmesh / 2) {
            kvec[0] = i * 2 * PI / Box;
          } else {
            kvec[0] = -(Nmesh - i) * 2 * PI / Box;
          }
          
          if(j < Nmesh / 2) {
            kvec[1] = j * 2 * PI / Box;
          } else {
            kvec[1] = -(Nmesh - j) * 2 * PI / Box;
          }

          if(k < Nmesh / 2) {
            kvec[2] = k * 2 * PI / Box;
          } else {
            kvec[2] = -(Nmesh - k) * 2 * PI / Box;
          }

          kmag2 = kvec[0] * kvec[0] + kvec[1] * kvec[1] + kvec[2] * kvec[2];
          kmag = sqrt(kmag2);

          if(SphereMode == 1) {
            if(kmag * Box / (2 * PI) > Nsample / 2) continue;  // select a sphere in k-space
          } else {
            if(fabs(kvec[0]) * Box / (2 * PI) > Nsample / 2) continue;
            if(fabs(kvec[1]) * Box / (2 * PI) > Nsample / 2) continue;
            if(fabs(kvec[2]) * Box / (2 * PI) > Nsample / 2) continue;
          }

          // Eq A1 of Scoccimarro et al 1108.5512, only k dependence is relevant since 
          // normalization terms in the power cancel because ker0+kerA+kerB == 0
          kmagA = exp((PrimordialIndex - 4.0) * kerA * log(kmag2) * 0.5); 
          kmagB = exp((PrimordialIndex - 4.0) * kerB * log(kmag2) * 0.5);

          cppA[coord][0] = kmagA * cpot[coord][0];
          cppA[coord][1] = kmagA * cpot[coord][1];

          cppB[coord][0] = kmagB * cpot[coord][0];
          cppB[coord][1] = kmagB * cpot[coord][1]; 
        }
      }
    }

    if(ThisTask == 0) printf("Fourier transforming initial potential ppA to configuration...\n");
#ifdef SINGLE_PRECISION
    Inverse_plan = fftwf_mpi_plan_dft_c2r_3d(Nmesh,Nmesh,Nmesh,cppA,ppA,MPI_COMM_WORLD,FFTW_ESTIMATE); 
    fftwf_execute(Inverse_plan);
    fftwf_destroy_plan(Inverse_plan);
#else
    Inverse_plan = fftw_mpi_plan_dft_c2r_3d(Nmesh,Nmesh,Nmesh,cppA,ppA,MPI_COMM_WORLD,FFTW_ESTIMATE);
    fftw_execute(Inverse_plan);
    fftw_destroy_plan(Inverse_plan);
#endif

    if(ThisTask == 0) printf("Fourier transforming initial potential ppB to configuration...\n");
#ifdef SINGLE_PRECISION
    Inverse_plan = fftwf_mpi_plan_dft_c2r_3d(Nmesh,Nmesh,Nmesh,cppB,ppB,MPI_COMM_WORLD,FFTW_ESTIMATE); 
    fftwf_execute(Inverse_plan);
    fftwf_destroy_plan(Inverse_plan);
#else
    Inverse_plan = fftw_mpi_plan_dft_c2r_3d(Nmesh,Nmesh,Nmesh,cppB,ppB,MPI_COMM_WORLD,FFTW_ESTIMATE);
    fftw_execute(Inverse_plan);
    fftw_destroy_plan(Inverse_plan);
#endif

    // add the terms in real space
    for(i = 0; i < Local_nx; i++) {
      for(j = 0; j < Nmesh; j++)    {
        for(k = 0; k < Nmesh; k++)    {
          coord = (i * Nmesh + j) * (2 * (Nmesh / 2 + 1)) + k; 
          ppA[coord] = ppA[coord]*ppB[coord];                     // ppA simply accoumulates A*B
        }
      }
    } 

    if(ThisTask == 0) printf("Fourier transforming convolution to fourier space...\n");
#ifdef SINGLE_PRECISION
    Forward_plan = fftwf_mpi_plan_dft_r2c_3d(Nmesh,Nmesh,Nmesh,ppA,cppA,MPI_COMM_WORLD,FFTW_ESTIMATE); 
    fftwf_execute(Forward_plan);
    fftwf_destroy_plan(Forward_plan);
#else
    Forward_plan = fftw_mpi_plan_dft_r2c_3d(Nmesh,Nmesh,Nmesh,ppA,cppA,MPI_COMM_WORLD,FFTW_ESTIMATE);
    fftw_execute(Forward_plan);
    fftw_destroy_plan(Forward_plan);
#endif

    // apply ker0 to the convolution of A and B
    // remove the N^3 I got by forward fourier transforming and put zero to zero mode
    nmesh3 = ((unsigned long long) Nmesh) * ((unsigned long long) Nmesh ) * ((unsigned long long) Nmesh);    
    for(ii = 0; ii < Local_nx; ii++) {
      for(j = 0; j < Nmesh; j++)       {
        for(k = 0; k <= Nmesh / 2 ; k++) {
          i = ii + Local_x_start;
          coord = (ii * Nmesh + j) * (Nmesh / 2 + 1) + k;

          // coord=0 is the fundamental mode in k-space
          if(i == 0 && j == 0 && k == 0) {
            cppC[coord][0]=0.; 
            cppC[coord][1]=0.; 
            continue; 
          }   

          if(i < Nmesh / 2) {
            kvec[0] = i * 2 * PI / Box;
          } else {
            kvec[0] = -(Nmesh - i) * 2 * PI / Box;
          }
 
          if(j < Nmesh / 2) {
            kvec[1] = j * 2 * PI / Box;
          } else {
            kvec[1] = -(Nmesh - j) * 2 * PI / Box;
          }

          if(k < Nmesh / 2) {
            kvec[2] = k * 2 * PI / Box;
          } else {
            kvec[2] = -(Nmesh - k) * 2 * PI / Box;
          }

          kmag2 = kvec[0] * kvec[0] + kvec[1] * kvec[1] + kvec[2] * kvec[2];
          kmag = sqrt(kmag2);

          if(SphereMode == 1) {
            if(kmag * Box / (2 * PI) > Nsample / 2) continue;  // select a sphere in k-space
          } else {
            if(fabs(kvec[0]) * Box / (2 * PI) > Nsample / 2) continue;
            if(fabs(kvec[1]) * Box / (2 * PI) > Nsample / 2) continue;
            if(fabs(kvec[2]) * Box / (2 * PI) > Nsample / 2) continue;
          }
		   
          // Eq A1 of Scoccimarro et al 1108.5512, only the k dependence is relevant since 
          // normalization terms in the power cancel because ker0+kerA+kerB == 0

          // All ikernel terms are accomulated to cppC
          kmag0 = exp((PrimordialIndex - 4.0) * ker0 * log(kmag2) * 0.5); 

          cppC[coord][0] = cppC[coord][0] + kerCoef * kmag0 * cppA[coord][0] / (double) nmesh3;   //cppA accumulated A*B earlier
          cppC[coord][1] = cppC[coord][1] + kerCoef * kmag0 * cppA[coord][1] / (double) nmesh3;   //cppA accumulated A*B earlier  
        }
      }
    }
  }
  
  // add the linear (i.e, non quadratic) part
  for(ii = 0; ii < Local_nx; ii++) {
    for(j = 0; j < Nmesh; j++)       {
      for(k = 0; k <= Nmesh / 2; k++)  {
        i = ii + Local_x_start;
        coord = (ii * Nmesh + j) * (Nmesh / 2 + 1) + k;

        // cpot has linear part + cppC has quadratic part */
        cpot[coord][0] = cpot[coord][0] + Fnl * cppC[coord][0];   
        cpot[coord][1] = cpot[coord][1] + Fnl * cppC[coord][1];
      }
    }
  } 

  free(cppA);
  free(cppB);
  free(cppC);


// Specific type of non-gaussianity (local, equilateral or orthogonal) for n_s=1 only
// ==================================================================================
#else

// Local primordial potential
// ==========================
#ifdef LOCAL_FNL  

  if(ThisTask == 0) printf("Fourier transforming initial potential to configuration...\n");
#ifdef SINGLE_PRECISION
  Inverse_plan = fftwf_mpi_plan_dft_c2r_3d(Nmesh,Nmesh,Nmesh,cpot,pot,MPI_COMM_WORLD,FFTW_ESTIMATE); 
  fftwf_execute(Inverse_plan);
  fftwf_destroy_plan(Inverse_plan);
#else
  Inverse_plan = fftw_mpi_plan_dft_c2r_3d(Nmesh,Nmesh,Nmesh,cpot,pot,MPI_COMM_WORLD,FFTW_ESTIMATE);
  fftw_execute(Inverse_plan);
  fftw_destroy_plan(Inverse_plan);
#endif

   // square the potential in configuration space
   for(i = 0; i < Local_nx; i++) {
     for(j = 0; j < Nmesh; j++)    {
       for(k = 0; k < Nmesh; k++)    {
         coord = (i * Nmesh + j) * (2 * (Nmesh / 2 + 1)) + k; 
         pot[coord] = pot[coord] + Fnl * pot[coord]*pot[coord];
       }
     }
   }

   if(ThisTask == 0) printf("Fourier transforming squared potential ...\n");
#ifdef SINGLE_PRECISION
    Forward_plan = fftwf_mpi_plan_dft_r2c_3d(Nmesh,Nmesh,Nmesh,pot,cpot,MPI_COMM_WORLD,FFTW_ESTIMATE); 
    fftwf_execute(Forward_plan);
    fftwf_destroy_plan(Forward_plan);
#else
    Forward_plan = fftw_mpi_plan_dft_r2c_3d(Nmesh,Nmesh,Nmesh,pot,cpot,MPI_COMM_WORLD,FFTW_ESTIMATE);
    fftw_execute(Forward_plan);
    fftw_destroy_plan(Forward_plan);
#endif

    // remove the N^3 I got by forward fourier transforming and put zero to zero mode */
    nmesh3 = ((unsigned long long) Nmesh) * ((unsigned long long) Nmesh ) * ((unsigned long long) Nmesh);    
    for(i = 0; i < Local_nx; i++) {
      for(j = 0; j < Nmesh; j++)     {
        for(k = 0; k <= Nmesh / 2; k++) {
          coord = (i * Nmesh + j) * (Nmesh / 2 + 1) + k;
          cpot[coord][0] /= (double) nmesh3; 
          cpot[coord][1] /= (double) nmesh3; 
        }
      }
    }

    if(ThisTask == 0) {
      cpot[0][1] = 0.0;
      cpot[0][0] = 0.0; 
    }

// Non-local primordial potential
// ==========================
#else
    
  // allocate partpotential
  cpartpot = (complex_kind *) malloc(bytes = sizeof(complex_kind) * Total_size);
  partpot = (float_kind *) cpartpot;

  cp1p2p3sym = (complex_kind *) malloc(bytes = sizeof(complex_kind) * Total_size);
  p1p2p3sym = (float_kind *) cp1p2p3sym;

  cp1p2p3sca = (complex_kind *) malloc(bytes = sizeof(complex_kind) * Total_size);
  p1p2p3sca = (float_kind *) cp1p2p3sca;

  cp1p2p3nab = (complex_kind *) malloc(bytes = sizeof(complex_kind) * Total_size);
  p1p2p3nab = (float_kind *) cp1p2p3nab;

  cp1p2p3tre = (complex_kind *) malloc(bytes = sizeof(complex_kind) * Total_size);
  p1p2p3tre = (float_kind *) cp1p2p3tre;

  // first, clean the array
  for(i = 0; i < Local_nx; i++) {
    for(j = 0; j < Nmesh; j++)     {
      for(k = 0; k <= Nmesh / 2; k++) {
        coord = (i * Nmesh + j) * (Nmesh / 2 + 1) + k;
        cp1p2p3sym[coord][0] = 0;
        cp1p2p3sym[coord][1] = 0;
        cp1p2p3sca[coord][0] = 0;
        cp1p2p3sca[coord][1] = 0;
        cp1p2p3nab[coord][0] = 0;
        cp1p2p3nab[coord][1] = 0;
        cp1p2p3tre[coord][0] = 0;
        cp1p2p3tre[coord][1] = 0;
        cpartpot[coord][0] = 0;
        cpartpot[coord][1] = 0;
      }
    }
  }

  // multiply by k  
  for(ii = 0; ii < Local_nx; ii++) {
    for(j = 0; j < Nmesh; j++)       {
      for(k = 0; k <= Nmesh / 2; k++)  {

        i = ii + Local_x_start;
        coord = (ii * Nmesh + j) * (Nmesh / 2 + 1) + k;

        if(i < Nmesh / 2) {
          kvec[0] = i * 2 * PI / Box;
        } else {
          kvec[0] = -(Nmesh - i) * 2 * PI / Box;
        }
        
        if(j < Nmesh / 2) {
          kvec[1] = j * 2 * PI / Box;
        } else {
          kvec[1] = -(Nmesh - j) * 2 * PI / Box;
        }
 
        if(k < Nmesh / 2) {
          kvec[2] = k * 2 * PI / Box;
        } else {
          kvec[2] = -(Nmesh - k) * 2 * PI / Box;
        }

        kmag2 = kvec[0] * kvec[0] + kvec[1] * kvec[1] + kvec[2] * kvec[2];
        kmag = sqrt(kmag2);

        cpartpot[coord][0] = kmag * cpot[coord][0];
        cpartpot[coord][1] = kmag * cpot[coord][1];

        cp1p2p3nab[coord][0] = kmag2 * cpot[coord][0];
        cp1p2p3nab[coord][1] = kmag2 * cpot[coord][1]; 
      }
    }
  }

  // fourier transform back to real
  if(ThisTask == 0) printf("Fourier transforming initial potential to configuration...\n");
#ifdef SINGLE_PRECISION
  Inverse_plan = fftwf_mpi_plan_dft_c2r_3d(Nmesh,Nmesh,Nmesh,cpot,pot,MPI_COMM_WORLD,FFTW_ESTIMATE); 
  fftwf_execute(Inverse_plan);
  fftwf_destroy_plan(Inverse_plan);
#else
  Inverse_plan = fftw_mpi_plan_dft_c2r_3d(Nmesh,Nmesh,Nmesh,cpot,pot,MPI_COMM_WORLD,FFTW_ESTIMATE);
  fftw_execute(Inverse_plan);
  fftw_destroy_plan(Inverse_plan);
#endif

  if(ThisTask == 0) printf("Fourier transforming partpotential to configuration...\n");
#ifdef SINGLE_PRECISION
  Inverse_plan = fftwf_mpi_plan_dft_c2r_3d(Nmesh,Nmesh,Nmesh,cpartpot,partpot,MPI_COMM_WORLD,FFTW_ESTIMATE); 
  fftwf_execute(Inverse_plan);
  fftwf_destroy_plan(Inverse_plan);
#else
  Inverse_plan = fftw_mpi_plan_dft_c2r_3d(Nmesh,Nmesh,Nmesh,cpartpot,partpot,MPI_COMM_WORLD,FFTW_ESTIMATE);
  fftw_execute(Inverse_plan);
  fftw_destroy_plan(Inverse_plan);
#endif

  if(ThisTask == 0) printf("Fourier transforming nabpotential to configuration...\n");
#ifdef SINGLE_PRECISION
  Inverse_plan = fftwf_mpi_plan_dft_c2r_3d(Nmesh,Nmesh,Nmesh,cp1p2p3nab,p1p2p3nab,MPI_COMM_WORLD,FFTW_ESTIMATE); 
  fftwf_execute(Inverse_plan);
  fftwf_destroy_plan(Inverse_plan);
#else
  Inverse_plan = fftw_mpi_plan_dft_c2r_3d(Nmesh,Nmesh,Nmesh,cp1p2p3nab,p1p2p3nab,MPI_COMM_WORLD,FFTW_ESTIMATE);
  fftw_execute(Inverse_plan);
  fftw_destroy_plan(Inverse_plan);
#endif

  // multiplying terms in real space
  for(i = 0; i < Local_nx; i++) {
    for(j = 0; j < Nmesh; j++)    {
      for(k = 0; k < Nmesh; k++)    {
        coord = (i * Nmesh + j) * (2 * (Nmesh / 2 + 1)) + k; 
        p1p2p3sym[coord] = partpot[coord]*partpot[coord];
        p1p2p3sca[coord] = pot[coord]*partpot[coord];   
        p1p2p3nab[coord] = pot[coord]*p1p2p3nab[coord];
        p1p2p3tre[coord] = p1p2p3nab[coord]*partpot[coord];
        partpot[coord] = pot[coord]*pot[coord];                 // NOTE: now partpot is potential squared
      }
    }
  }

  if(ThisTask == 0) printf("Fourier transforming potential ...\n");
#ifdef SINGLE_PRECISION
  Forward_plan = fftwf_mpi_plan_dft_r2c_3d(Nmesh,Nmesh,Nmesh,pot,cpot,MPI_COMM_WORLD,FFTW_ESTIMATE); 
  fftwf_execute(Forward_plan);
  fftwf_destroy_plan(Forward_plan);
#else
  Forward_plan = fftw_mpi_plan_dft_r2c_3d(Nmesh,Nmesh,Nmesh,pot,cpot,MPI_COMM_WORLD,FFTW_ESTIMATE);
  fftw_execute(Forward_plan);
  fftw_destroy_plan(Forward_plan);
#endif 

  if(ThisTask == 0) printf("Fourier transforming squared potential ...\n");
#ifdef SINGLE_PRECISION
  Forward_plan = fftwf_mpi_plan_dft_r2c_3d(Nmesh,Nmesh,Nmesh,partpot,cpartpot,MPI_COMM_WORLD,FFTW_ESTIMATE); 
  fftwf_execute(Forward_plan);
  fftwf_destroy_plan(Forward_plan);
#else
  Forward_plan = fftw_mpi_plan_dft_r2c_3d(Nmesh,Nmesh,Nmesh,partpot,cpartpot,MPI_COMM_WORLD,FFTW_ESTIMATE);
  fftw_execute(Forward_plan);
  fftw_destroy_plan(Forward_plan);
#endif

  if(ThisTask == 0) printf("Fourier transforming p1p2p3sym potential ...\n");
#ifdef SINGLE_PRECISION
  Forward_plan = fftwf_mpi_plan_dft_r2c_3d(Nmesh,Nmesh,Nmesh,p1p2p3psym,cp1p2p3sym,MPI_COMM_WORLD,FFTW_ESTIMATE); 
  fftwf_execute(Forward_plan);
  fftwf_destroy_plan(Forward_plan);
#else
  Forward_plan = fftw_mpi_plan_dft_r2c_3d(Nmesh,Nmesh,Nmesh,p1p2p3sym,cp1p2p3sym,MPI_COMM_WORLD,FFTW_ESTIMATE);
  fftw_execute(Forward_plan);
  fftw_destroy_plan(Forward_plan);
#endif

  if(ThisTask == 0) printf("Fourier transforming p1p2p3sca potential ...\n");
#ifdef SINGLE_PRECISION
  Forward_plan = fftwf_mpi_plan_dft_r2c_3d(Nmesh,Nmesh,Nmesh,p1p2p3sca,cp1p2p3sca,MPI_COMM_WORLD,FFTW_ESTIMATE); 
  fftwf_execute(Forward_plan);
  fftwf_destroy_plan(Forward_plan);
#else
  Forward_plan = fftw_mpi_plan_dft_r2c_3d(Nmesh,Nmesh,Nmesh,p1p2p3sca,cp1p2p3sca,MPI_COMM_WORLD,FFTW_ESTIMATE);
  fftw_execute(Forward_plan);
  fftw_destroy_plan(Forward_plan);
#endif

  if(ThisTask == 0) printf("Fourier transforming p1p2p3nab potential ...\n");
#ifdef SINGLE_PRECISION
  Forward_plan = fftwf_mpi_plan_dft_r2c_3d(Nmesh,Nmesh,Nmesh,p1p2p3nab,cp1p2p3nab,MPI_COMM_WORLD,FFTW_ESTIMATE); 
  fftwf_execute(Forward_plan);
  fftwf_destroy_plan(Forward_plan);
#else
  Forward_plan = fftw_mpi_plan_dft_r2c_3d(Nmesh,Nmesh,Nmesh,p1p2p3nab,cp1p2p3nab,MPI_COMM_WORLD,FFTW_ESTIMATE);
  fftw_execute(Forward_plan);
  fftw_destroy_plan(Forward_plan);
#endif

  if(ThisTask == 0) printf("Fourier transforming p1p2p3tre potential ...\n");
#ifdef SINGLE_PRECISION
  Forward_plan = fftwf_mpi_plan_dft_r2c_3d(Nmesh,Nmesh,Nmesh,p1p2p3tre,cp1p2p3tre,MPI_COMM_WORLD,FFTW_ESTIMATE); 
  fftwf_execute(Forward_plan);
  fftwf_destroy_plan(Forward_plan);
#else
  Forward_plan = fftw_mpi_plan_dft_r2c_3d(Nmesh,Nmesh,Nmesh,p1p2p3tre,cp1p2p3tre,MPI_COMM_WORLD,FFTW_ESTIMATE);
  fftw_execute(Forward_plan);
  fftw_destroy_plan(Forward_plan);
#endif

  // divide by appropiate k's, sum terms according to non-local model */ 
  // remove the N^3 I got by forward fourier transforming and put zero to zero mode */
  nmesh3 = ((unsigned long long) Nmesh) * ((unsigned long long) Nmesh ) * ((unsigned long long) Nmesh);
  for(ii = 0; ii < Local_nx; ii++) {
    for(j = 0; j < Nmesh; j++)       {
      for(k = 0; k <= Nmesh / 2 ; k++) {
        i = ii + Local_x_start;
        coord = (ii * Nmesh + j) * (Nmesh / 2 + 1) + k;

        if(i < Nmesh / 2) {
          kvec[0] = i * 2 * PI / Box;
        } else {
          kvec[0] = -(Nmesh - i) * 2 * PI / Box;
        }

        if(j < Nmesh / 2) {
          kvec[1] = j * 2 * PI / Box;
        } else {
          kvec[1] = -(Nmesh - j) * 2 * PI / Box;
        }

        if(k < Nmesh / 2) {
          kvec[2] = k * 2 * PI / Box;
        } else {
          kvec[2] = -(Nmesh - k) * 2 * PI / Box;
        }

        kmag2 = kvec[0] * kvec[0] + kvec[1] * kvec[1] + kvec[2] * kvec[2];
        kmag = sqrt(kmag2);
   
        if(i == 0 && j == 0 && k == 0) {
          cpot[0][0]=0.;
          cpot[0][1]=0.;
          continue;
        }

// F_nl equilateral
// ================ 
#ifdef EQUIL_FNL
 
        cpot[coord][0] = cpot[coord][0] + Fnl * (-3*cpartpot[coord][0] - 2*cp1p2p3sym[coord][0] / kmag2 + 4*cp1p2p3sca[coord][0] /kmag + 2*cp1p2p3nab[coord][0] / kmag2);
        cpot[coord][1] = cpot[coord][1] + Fnl * (-3*cpartpot[coord][1] - 2*cp1p2p3sym[coord][1] / kmag2 + 4*cp1p2p3sca[coord][1] /kmag + 2*cp1p2p3nab[coord][1] / kmag2); 
        cpot[coord][0] /= (double) nmesh3; 
        cpot[coord][1] /= (double) nmesh3; 

#endif

// F_nl othogonal
// ==============  
#ifdef ORTHO_FNL 

        cpot[coord][0] = cpot[coord][0] + Fnl * (-9*cpartpot[coord][0] - 8*cp1p2p3sym[coord][0] / kmag2 + 10*cp1p2p3sca[coord][0] /kmag + 8*cp1p2p3nab[coord][0] / kmag2);
        cpot[coord][1] = cpot[coord][1] + Fnl * (-9*cpartpot[coord][1] - 8*cp1p2p3sym[coord][1] / kmag2 + 10*cp1p2p3sca[coord][1] /kmag + 8*cp1p2p3nab[coord][1] / kmag2);
        cpot[coord][0] /= (double) nmesh3; 
        cpot[coord][1] /= (double) nmesh3; 

#endif

        if(i == 0 && j == 0 && k == 0) {
          cpot[0][0]=0.;
          cpot[0][1]=0.;
          continue;
        }
      }
    }
  }

  free(cpartpot);
  free(cp1p2p3sym);
  free(cp1p2p3sca);
  free(cp1p2p3nab);
  free(cp1p2p3tre);
 
#endif

// Compute displacements using nongaussian potential
// ================================================
#endif 

  if(ThisTask==0) {
    printf("Computing displacement using non-gaussian potential...\n");
    fflush(stdout);
  }

  for(axes=0,bytes=0; axes < 3; axes++) {
    cdisp[axes] = (complex_kind *) malloc(bytes += sizeof(complex_kind) * Total_size);
    disp[axes] = (float_kind *) cdisp[axes];
  }

  if(ThisTask == 0) {
    printf("Starting axes = %d...\n", axes);
    fflush(stdout);
  }

  // first, clean the array
  for(i = 0; i < Local_nx; i++) {
    for(j = 0; j < Nmesh; j++)     {
      for(k = 0; k <= Nmesh / 2; k++) {
        for(axes = 0; axes < 3; axes++) {
          cdisp[axes][(i * Nmesh + j) * (Nmesh / 2 + 1) + k][0] = 0;
          cdisp[axes][(i * Nmesh + j) * (Nmesh / 2 + 1) + k][1] = 0;
        }
      }
    }
  } 

  for(ii = 0; ii < Local_nx; ii++) {
    for(j = 0; j < Nmesh; j++)       {
      for(k = 0; k <= Nmesh / 2 ; k++) {
        i = ii + Local_x_start;
        coord = (ii * Nmesh + j) * (Nmesh / 2 + 1) + k;

        if(i < Nmesh / 2) {
          kvec[0] = i * 2 * PI / Box;
        } else {
          kvec[0] = -(Nmesh - i) * 2 * PI / Box;
        }

        if(j < Nmesh / 2) {
          kvec[1] = j * 2 * PI / Box;
        } else {
          kvec[1] = -(Nmesh - j) * 2 * PI / Box;
        }
           
        if(k < Nmesh / 2) {
          kvec[2] = k * 2 * PI / Box;
        } else {
          kvec[2] = -(Nmesh - k) * 2 * PI / Box;
        }
           
        kmag2 = kvec[0] * kvec[0] + kvec[1] * kvec[1] + kvec[2] * kvec[2];
        kmag = sqrt(kmag2);

        t_of_k = TransferFunc(kmag);
        twb = t_of_k / (Beta * DstartFnl);   

        for(axes = 0; axes < 3; axes++) {
          cdisp[axes][coord][1] = kvec[axes] * twb * cpot[coord][0];
          cdisp[axes][coord][0] = - kvec[axes] * twb * cpot[coord][1];
        }
      }
    }
  } 

  free(cpot);

#ifdef TIMING
  endcpu = (double)clock();
  endwall = MPI_Wtime();
  CpuTime_2LPTng = (endcpu-startcpu)/(double)CLOCKS_PER_SEC;
  WallTime_2LPTng = endwall-startwall;
#endif

#endif

  // Compute the displacements (identical regardless of gaussianity)
  // ==============================================================
 
  // At this point, cdisp[axes] contains the complex Zeldovich displacement
  if(ThisTask == 0) {
    printf("Computing 2LPT displacements...\n");
    fflush(stdout);
  }
      
  // Compute displacement gradient
  for(i = 0; i < 6; i++) {
    cdigrad[i] = (complex_kind *)malloc(bytes = sizeof(complex_kind) * Total_size);
    digrad[i] = (float_kind *)cdigrad[i];
  }

  for(i = 0; i < Local_nx; i++) {
    for(j = 0; j < Nmesh; j++)     {
      for(k = 0; k <= Nmesh / 2; k++) {
        coord = (i * Nmesh + j) * (Nmesh / 2 + 1) + k;

        if((i + Local_x_start) < Nmesh / 2) {
          kvec[0] = (i + Local_x_start) * 2 * PI / Box;
        } else {
          kvec[0] = -(Nmesh - (i + Local_x_start)) * 2 * PI / Box;
        }  
  
        if(j < Nmesh / 2) {
          kvec[1] = j * 2 * PI / Box;
        } else {
          kvec[1] = -(Nmesh - j) * 2 * PI / Box;
        }

        if(k < Nmesh / 2) {
          kvec[2] = k * 2 * PI / Box;
        } else {
          kvec[2] = -(Nmesh - k) * 2 * PI / Box;
        }

        // Derivatives of ZA displacement
        // d(dis_i)/d(q_j)  -> sqrt(-1) k_j dis_i

        cdigrad[0][coord][0] = -cdisp[0][coord][1] * kvec[0]; // disp0,0
        cdigrad[0][coord][1] =  cdisp[0][coord][0] * kvec[0];

        cdigrad[1][coord][0] = -cdisp[0][coord][1] * kvec[1]; // disp0,1
        cdigrad[1][coord][1] =  cdisp[0][coord][0] * kvec[1];

        cdigrad[2][coord][0] = -cdisp[0][coord][1] * kvec[2]; // disp0,2
        cdigrad[2][coord][1] =  cdisp[0][coord][0] * kvec[2];
 
        cdigrad[3][coord][0] = -cdisp[1][coord][1] * kvec[1]; // disp1,1
        cdigrad[3][coord][1] =  cdisp[1][coord][0] * kvec[1];

        cdigrad[4][coord][0] = -cdisp[1][coord][1] * kvec[2]; // disp1,2
        cdigrad[4][coord][1] =  cdisp[1][coord][0] * kvec[2];
 
        cdigrad[5][coord][0] = -cdisp[2][coord][1] * kvec[2]; // disp2,2
        cdigrad[5][coord][1] =  cdisp[2][coord][0] * kvec[2];

      }
    }
  }

  if(ThisTask == 0) printf("Fourier transforming displacement gradient...\n");
  for(i = 0; i < 6; i++) {
#ifdef SINGLE_PRECISION
    Inverse_plan = fftwf_mpi_plan_dft_c2r_3d(Nmesh,Nmesh,Nmesh,cdigrad[i],digrad[i],MPI_COMM_WORLD,FFTW_ESTIMATE);
    fftwf_execute(Inverse_plan);
    fftwf_destroy_plan(Inverse_plan);
#else
    Inverse_plan = fftw_mpi_plan_dft_c2r_3d(Nmesh,Nmesh,Nmesh,cdigrad[i],digrad[i],MPI_COMM_WORLD,FFTW_ESTIMATE);
    fftw_execute(Inverse_plan);
    fftw_destroy_plan(Inverse_plan);
#endif
  }     

  // Compute second order source and store it in digrad[3]
  for(i = 0; i < Local_nx; i++) {
    for(j = 0; j < Nmesh; j++)    {
      for(k = 0; k < Nmesh; k++)    {
        coord = (i * Nmesh + j) * (2 * (Nmesh / 2 + 1)) + k;
        digrad[3][coord] = digrad[0][coord]*(digrad[3][coord]+digrad[5][coord])+digrad[3][coord]*digrad[5][coord]
                          -digrad[1][coord]*digrad[1][coord]-digrad[2][coord]*digrad[2][coord]-digrad[4][coord]*digrad[4][coord];
      }
    }
  }

  if(ThisTask == 0) printf("Fourier transforming second order source...\n");
#ifdef SINGLE_PRECISION
  Forward_plan = fftwf_mpi_plan_dft_r2c_3d(Nmesh,Nmesh,Nmesh,digrad[3],cdigrad[3], MPI_COMM_WORLD, FFTW_ESTIMATE);
  fftwf_execute(Forward_plan);
  fftwf_destroy_plan(Forward_plan);
#else
  Forward_plan = fftw_mpi_plan_dft_r2c_3d(Nmesh,Nmesh,Nmesh,digrad[3],cdigrad[3], MPI_COMM_WORLD, FFTW_ESTIMATE);
  fftw_execute(Forward_plan);
  fftw_destroy_plan(Forward_plan);
#endif

  // The memory allocated for cdisp2[0], [1], and [2] will be used for 2nd order displacements
  // Freeing the rest. cdigrad[3] still has 2nd order displacement source, free later
  for(axes = 0; axes < 3; axes++) {
    cdisp2[axes] = cdigrad[axes]; 
    disp2[axes] = (float_kind *) cdisp2[axes];
  }
  free(cdigrad[4]); 
  free(cdigrad[5]); 

  // Solve Poisson eq. and calculate 2nd order displacements
  for(i = 0; i < Local_nx; i++) {
    for(j = 0; j < Nmesh; j++)     {
      for(k = 0; k <= Nmesh / 2; k++) {
        coord = (i * Nmesh + j) * (Nmesh / 2 + 1) + k;
        if((i + Local_x_start) < Nmesh / 2) {
          kvec[0] = (i + Local_x_start) * 2 * PI / Box;
        } else {
          kvec[0] = -(Nmesh - (i + Local_x_start)) * 2 * PI / Box;
        }
   
        if(j < Nmesh / 2) {
          kvec[1] = j * 2 * PI / Box;
        } else {
          kvec[1] = -(Nmesh - j) * 2 * PI / Box;
        }

        if(k < Nmesh / 2) {
          kvec[2] = k * 2 * PI / Box;
        } else {
          kvec[2] = -(Nmesh - k) * 2 * PI / Box;
        }

        kmag2 = kvec[0] * kvec[0] + kvec[1] * kvec[1] + kvec[2] * kvec[2];

        // cdisp2 = source * k / (sqrt(-1) k^2) */
        for(axes = 0; axes < 3; axes++) {
          if(kmag2 > 0.0) {
            cdisp2[axes][coord][0] = cdigrad[3][coord][1] * kvec[axes] / kmag2;
            cdisp2[axes][coord][1] = -cdigrad[3][coord][0] * kvec[axes] / kmag2;
          } else {
            cdisp2[axes][coord][0] = cdisp2[axes][coord][1] = 0.0;
          }
        }
      }
    }
  }
  
  // Free cdigrad[3]
  free(cdigrad[3]);

  // Now, both cdisp, and cdisp2 have the ZA and 2nd order displacements
  for(axes = 0; axes < 3; axes++) {
    if(ThisTask == 0) printf("Fourier transforming displacements, axis %d\n",axes);
#ifdef SINGLE_PRECISION
    Inverse_plan = fftwf_mpi_plan_dft_c2r_3d(Nmesh,Nmesh,Nmesh,cdisp[axes],disp[axes], MPI_COMM_WORLD, FFTW_ESTIMATE);
    fftwf_execute(Inverse_plan);
    fftwf_destroy_plan(Inverse_plan);
    Inverse_plan = fftwf_mpi_plan_dft_c2r_3d(Nmesh,Nmesh,Nmesh,cdisp2[axes],disp2[axes], MPI_COMM_WORLD, FFTW_ESTIMATE);
    fftwf_execute(Inverse_plan);
    fftwf_destroy_plan(Inverse_plan);
#else
    Inverse_plan = fftw_mpi_plan_dft_c2r_3d(Nmesh,Nmesh,Nmesh,cdisp[axes],disp[axes], MPI_COMM_WORLD, FFTW_ESTIMATE);
    fftw_execute(Inverse_plan);
    fftw_destroy_plan(Inverse_plan);
    Inverse_plan = fftw_mpi_plan_dft_c2r_3d(Nmesh,Nmesh,Nmesh,cdisp2[axes],disp2[axes], MPI_COMM_WORLD, FFTW_ESTIMATE);
    fftw_execute(Inverse_plan);
    fftw_destroy_plan(Inverse_plan);
#endif

    // now get the plane on the right side from neighbour on the right and send the left plane

    // send ZA disp
    MPI_Sendrecv(&(disp[axes][0]),sizeof(float_kind)*2*alloc_slice,MPI_BYTE,LeftTask,10,
                 &(disp[axes][2*last_slice]),sizeof(float_kind)*2*alloc_slice,MPI_BYTE,RightTask,10,MPI_COMM_WORLD,&status);

    // send 2nd order disp
    MPI_Sendrecv(&(disp2[axes][0]),sizeof(float_kind)*2*alloc_slice,MPI_BYTE,LeftTask,10,
                 &(disp2[axes][2*last_slice]),sizeof(float_kind)*2*alloc_slice,MPI_BYTE,RightTask,10,MPI_COMM_WORLD,&status);     
  }

  // read-out displacements
  gsl_rng_free(random_generator);
  free(seedtable);

  nmesh3 = ((unsigned long long ) Nmesh ) * ((unsigned long long) Nmesh) *  ((unsigned long long) Nmesh);

  for(axes = 0; axes < 3; axes++) {
#ifdef MEMORY_MODE
    ZA[axes]  = (float *)malloc(NumPart*sizeof(float));
    LPT[axes] = (float *)malloc(NumPart*sizeof(float));
#else
    ZA[axes]  = (float_kind *)malloc(NumPart*sizeof(float_kind));
    LPT[axes] = (float_kind *)malloc(NumPart*sizeof(float_kind));
#endif
    sumdis[axes] = 0;
    sumdis2[axes] = 0;
  }

  for (n = 0; n < Local_np; n++) {
    for (m = 0; m < Nsample; m++) {
      for (p = 0; p < Nsample; p++) {
        coord = (n * Nsample + m) * (Nsample) + p;

        u = (double)((n+Local_p_start)*Nmesh)/(double)Nsample;
        v = (double)(m*Nmesh)/(double)Nsample;
        w = (double)(p*Nmesh)/(double)Nsample;

        i = (int) u;
        j = (int) v;
        k = (int) w;

        if(i == (Local_x_start + Local_nx)) i = (Local_x_start + Local_nx) - 1;
        if(i < Local_x_start)               i = Local_x_start;
        if(j == Nmesh)                      j = Nmesh - 1;
        if(k == Nmesh)                      k = Nmesh - 1;

        u -= i;
        v -= j;
        w -= k;
 
        i -= Local_x_start;
        ii = i + 1;
        jj = j + 1;
        kk = k + 1;
 
        if(jj >= Nmesh) jj -= Nmesh;
        if(kk >= Nmesh) kk -= Nmesh;
  
        f1 = (1 - u) * (1 - v) * (1 - w);
        f2 = (1 - u) * (1 - v) * (w);
        f3 = (1 - u) * (v) * (1 - w);
        f4 = (1 - u) * (v) * (w);
        f5 = (u) * (1 - v) * (1 - w);
        f6 = (u) * (1 - v) * (w); 
        f7 = (u) * (v) * (1 - w);
        f8 = (u) * (v) * (w);

        for(axes = 0; axes < 3; axes++) {
          dis[axes] = disp[axes][(i * Nmesh + j)   * (2 * (Nmesh / 2 + 1)) + k]  * f1 +
                      disp[axes][(i * Nmesh + j)   * (2 * (Nmesh / 2 + 1)) + kk] * f2 +
                      disp[axes][(i * Nmesh + jj)  * (2 * (Nmesh / 2 + 1)) + k]  * f3 +
                      disp[axes][(i * Nmesh + jj)  * (2 * (Nmesh / 2 + 1)) + kk] * f4 +
                      disp[axes][(ii * Nmesh + j)  * (2 * (Nmesh / 2 + 1)) + k]  * f5 +
                      disp[axes][(ii * Nmesh + j)  * (2 * (Nmesh / 2 + 1)) + kk] * f6 +
                      disp[axes][(ii * Nmesh + jj) * (2 * (Nmesh / 2 + 1)) + k]  * f7 +
                      disp[axes][(ii * Nmesh + jj) * (2 * (Nmesh / 2 + 1)) + kk] * f8;

          dis2[axes] = disp2[axes][(i * Nmesh + j)   * (2 * (Nmesh / 2 + 1)) + k]  * f1 +
                       disp2[axes][(i * Nmesh + j)   * (2 * (Nmesh / 2 + 1)) + kk] * f2 +
                       disp2[axes][(i * Nmesh + jj)  * (2 * (Nmesh / 2 + 1)) + k]  * f3 +
                       disp2[axes][(i * Nmesh + jj)  * (2 * (Nmesh / 2 + 1)) + kk] * f4 +
                       disp2[axes][(ii * Nmesh + j)  * (2 * (Nmesh / 2 + 1)) + k]  * f5 +
                       disp2[axes][(ii * Nmesh + j)  * (2 * (Nmesh / 2 + 1)) + kk] * f6 +
                       disp2[axes][(ii * Nmesh + jj) * (2 * (Nmesh / 2 + 1)) + k]  * f7 +
                       disp2[axes][(ii * Nmesh + jj) * (2 * (Nmesh / 2 + 1)) + kk] * f8;

          dis2[axes] /= (double) nmesh3; 

          ZA[axes][coord]   =  dis[axes];
          LPT[axes][coord]  =  -3./7.*dis2[axes];
          sumdis[axes]  += dis[axes];
          sumdis2[axes] += -3./7.*dis2[axes];

          if(fabs(dis[axes] - 3./7. * dis2[axes]) > maxdisp) maxdisp = fabs(dis[axes] - 3./7. * dis2[axes]);
        }
      }
    }
  }

  // Make sure the average of the displacements is zero.
  for(axes = 0; axes < 3; axes++) {
    ierr = MPI_Allreduce(MPI_IN_PLACE, &(sumdis[axes]), 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    ierr = MPI_Allreduce(MPI_IN_PLACE, &(sumdis2[axes]), 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    sumdis[axes]  /= (double)TotNumPart;
    sumdis2[axes] /= (double)TotNumPart;
  }
  for(axes = 0; axes < 3; axes++) free(cdisp[axes]);
  for(axes = 0; axes < 3; axes++) free(cdigrad[axes]);

  for(q = 0; q < NumPart; q++) {
    for(axes = 0; axes < 3; axes++) {
      ZA[axes][q]  += sumdis[axes];
      LPT[axes][q] += sumdis2[axes];
    }
  }

  MPI_Reduce(&maxdisp, &max_disp_glob, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

  if (ThisTask == 0) {
    printf("Calculated Zeldovich and 2LPT displacements...\n");
    printf("Maximum displacement = %lf kpc/h (%lf in units of the particle separation)...\n\n",max_disp_glob/(InputSpectrum_UnitLength_in_cm/UnitLength_in_cm), max_disp_glob / (Box / Nmesh));
    fflush(stdout);
  }

  return;
}
