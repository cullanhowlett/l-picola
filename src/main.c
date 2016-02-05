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

/* ========================================================*/
/* This file contains the main driver routine for L-PICOLA.*/
/* v1.1: New routines for calculating the second-order     */
/*       growth factor and associated derivatives          */
/* v1.2: Particle positions and velocities are now output  */
/*       in user-defined units as opposed to always being  */
/*       in Mpc/h. Fixed a small bug in the naming         */
/*       output files that named z=9 snapshots as z8p100   */
/*       due to floating-point problems with the value 0.1 */
/* ========================================================*/
     
#include "vars.h"
#include "proto.h"

int main(int argc, char **argv) {
   
  // Set up MPI
  // ==========
  ierr = MPI_Init(&argc, &argv);
  ierr = MPI_Comm_rank(MPI_COMM_WORLD, &ThisTask);
  ierr = MPI_Comm_size(MPI_COMM_WORLD, &NTask);
#ifdef SINGLE_PRECISION
  fftwf_mpi_init();
#else
  fftw_mpi_init();
#endif

  if(argc < 2) {
    if(ThisTask == 0) {
      fprintf(stdout, "Input parameters not found\n");
      fprintf(stdout, "Call with <ParameterFile>\n");
    }
    ierr = MPI_Finalize();
    exit(1);
  }

  // Read the run parameters and setup code
  // ======================================
  int i;
#ifdef TIMING
  double startcpu, endcpu;
  double startwall, endwall;
#ifndef LIGHTCONE
  double startcpu2, endcpu2;
  double startwall2, endwall2;
#endif
  startcpu = (double)clock();
  startwall = MPI_Wtime();
#endif

  if (ThisTask == 0) {
    printf("\nReading Input Parameters and setting up L-PICOLA\n");
    printf("================================================\n");
  }
  read_parameterfile(argv[1]);
  read_outputs();
  set_units();
#ifdef LIGHTCONE
  set_lightcone();
#endif

  if(ThisTask == 0) {
    printf("\nRun Parameters\n");
    printf("==============\n");
    printf("Cosmology:\n");
    printf("  Omega Matter(z=0) = %lf\n",Omega);
    printf("  Omega Lambda(z=0) = %lf\n",OmegaLambda);
    printf("  Omega Baryon(z=0) = %lf\n",OmegaBaryon);
    printf("  Hubble Parameter(z=0) = %lf\n",HubbleParam);
    printf("  Sigma8(z=0) = %lf\n",Sigma8);
#ifndef GAUSSIAN
    printf("  F_nl = %lf\n",Fnl);
#endif
    printf("  Primordial Index = %lf\n",PrimordialIndex);
    printf("  Initial Redshift  = %lf\n",Init_Redshift);
#ifndef GAUSSIAN
    printf("  F_nl Redshift  = %lf\n",Fnl_Redshift);
#endif
    printf("\nSimulation:\n");
    printf("  Nmesh = %d\n", Nmesh);
    printf("  Nsample = %d\n", Nsample);
    printf("  Boxsize = %lf\n", Box);
    printf("  Buffer Size = %lf\n", Buffer);
    switch(WhichSpectrum) {
      case 0:
        switch (WhichTransfer) {
          case 1:
            printf("  Using Tabulated Transfer Function\n");
            break;
          default:
            printf("  Using Eisenstein & Hu Transfer Function\n");
            break;
        }
        break;
      case 1:
        printf("  Using Tabulated Power Spectrum\n");
        break;   
      default:
        printf("  Using Eisenstein & Hu Power Spectrum\n");
        break;
    }      
    if (UseCOLA) {
      printf("  Using COLA method\n"); 
    } else {
      printf("  Using Standard PM method\n");
    }
    printf("\nTimestepping:\n");
    printf("  nLPT = %lf\n", nLPT);
    printf("  DeltaA = %d\n", DeltaA);
    printf("  StepDist = %d\n", StepDist);
#ifdef LIGHTCONE
    printf("\nLightcone:\n");
    printf("  Maximum Comoving Radius = %lf\n", Light/Hubble*KickStd(1.0/(1.0+OutputList[0].Redshift),1.0));
    printf("  Origin (x, y, z) = %lf, %lf, %lf\n", Origin_x, Origin_y, Origin_z);
    printf("  Nrep_min (x, y, z) = %d (%lf Mpc/h), %d (%lf Mpc/h), %d (%lf Mpc/h)\n", Nrep_neg_x, -Nrep_neg_x*Box-Origin_x, Nrep_neg_y, -Nrep_neg_y*Box-Origin_y, Nrep_neg_z, -Nrep_neg_z*Box-Origin_z);
    printf("  Nrep_max (x, y, z) = %d (%lf Mpc/h), %d (%lf Mpc/h), %d (%lf Mpc/h)\n", Nrep_pos_x, (Nrep_pos_x+1)*Box-Origin_x, Nrep_pos_y, (Nrep_pos_y+1)*Box-Origin_y, Nrep_pos_z, (Nrep_pos_z+1)*Box-Origin_z);
#endif
    printf("\nOutputs:\n");
    for (i=0; i<Noutputs; i++) printf("  Redshift = %lf, Nsteps = %d\n", OutputList[i].Redshift, OutputList[i].Nsteps);
    fflush(stdout);
  }   

  if (ThisTask == 0) {
    printf("\nInitialising Transfer Function/Power Spectrum\n");
    printf("=============================================\n");
  }
  initialize_transferfunction();
  initialize_powerspectrum();
  initialize_ffts();
  initialize_parts();

  if(ThisTask == 0) {
    printf("Creating initial conditions\n");
    printf("===========================\n");
    fflush(stdout);
  }

#ifdef TIMING
  endcpu = (double)clock();
  endwall = MPI_Wtime();
  CpuTime_Init = (endcpu-startcpu)/(double)CLOCKS_PER_SEC;
  WallTime_Init = endwall-startwall;
  startcpu = (double)clock();
  startwall = MPI_Wtime();
#endif

  // Create the calculate the Zeldovich and 2LPT displacements and create the initial conditions
  // ===========================================================================================
  int j, k, m;
  int NoutputStart = 0;
  int timeStep;
  unsigned int coord=0;
  double da=0;
  double A=1.0/(1.0+Init_Redshift);  // This is the scale factor which we'll be advancing below.
  double Di=growthD(A);              // initial growth factor
  double Di2=growthD2(A);            // initial 2nd order growth factor  
  double Dv=QdD1da(A);               // Q*dD_{za}/da
  double Dv2=QdD2da(A);              // Q*dD_{2}/da

  // A is the scale factor of the particle positions.
  // AI is the scale factor of the particle velocities.
  // AF is the scale factor to which we should kick the particle velocities.
  // AFF is the scale factor to which we should drift the particle positions.
  double AI=A,AF=A,AFF=A;  

  displacement_fields();
    
  P = (struct part_data *) malloc((int)(ceil(NumPart*Buffer))*sizeof(struct part_data));

  // Generate the initial particle positions and velocities
  // If UseCOLA = 0 (non-COLA), then velocity is ds/dy, which is simply the 2LPT IC.
  // Else set vel = 0 if we subtract LPT. This is the same as the action of the operator L_- from TZE, as initial velocities are in 2LPT.
  for(i=0; i<Local_np; i++) {
    for (j=0; j<Nsample; j++) {
      for (k=0; k<Nsample; k++) {
        coord = (i * Nsample + j) * Nsample + k;

#ifdef PARTICLE_ID          
        P[coord].ID = ((unsigned long long)((i + Local_p_start) * Nsample + j)) * (unsigned long long)Nsample + (unsigned long long)k + 1;
#endif

        for (m=0; m<3; m++) {
          P[coord].Dz[m] = ZA[m][coord];
          P[coord].D2[m] = LPT[m][coord];
          if (UseCOLA == 0) {
            P[coord].Vel[m] = P[coord].Dz[m]*Dv+P[coord].D2[m]*Dv2;
          } else {
            P[coord].Vel[m] = 0.0;
          }
        }

        P[coord].Pos[0] = periodic_wrap((i+Local_p_start)*(Box/(double)Nsample)+P[coord].Dz[0]*Di+P[coord].D2[0]*Di2);
        P[coord].Pos[1] = periodic_wrap(j*(Box/(double)Nsample)+P[coord].Dz[1]*Di+P[coord].D2[1]*Di2);
        P[coord].Pos[2] = periodic_wrap(k*(Box/(double)Nsample)+P[coord].Dz[2]*Di+P[coord].D2[2]*Di2);

      }
    }
  }

  for (i=0; i<3; i++) {
    free(ZA[i]);
    free(LPT[i]);
  }

  // If we want to output or start the lightcone at the initial redshift this is where we do it (it is tricky to compare
  // floating point numbers due to rounding errors so instead we see whether they are close)
  // ===================================================================================================================
  if (((Init_Redshift-OutputList[0].Redshift)/Init_Redshift <= 1.0E-6) || (Init_Redshift <= 1.0e-6)) {

#ifndef LIGHTCONE

    // Output particles.
    if (ThisTask == 0) {
      printf("Outputting Initial Conditions\n");
      printf("=============================\n\n");
    }

    sumx=0;
    sumy=0;
    sumz=0;

#ifdef TIMING
    startcpu2 = (double)clock();
    startwall2 = MPI_Wtime();
#endif
    Output(A,Init_Redshift,Dv,Dv2);
#ifdef TIMING
    endcpu2 = (double)clock();
    endwall2 = MPI_Wtime();
    CpuTime_2LPToutput = (endcpu-startcpu)/(double)CLOCKS_PER_SEC;
    WallTime_2LPToutput = endwall-startwall;
#endif

    // If this is the only output timestep then simply skip to the end
    if(Noutputs == 1) {
#ifdef TIMING
      endcpu = (double)clock();
      endwall = MPI_Wtime();
      CpuTime_2LPT = (endcpu-startcpu)/(double)CLOCKS_PER_SEC;
      CpuTime_Move    = (double *)calloc(1, sizeof(double));
      CpuTime_PtoMesh = (double *)calloc(1, sizeof(double));
      CpuTime_Forces       = (double *)calloc(1, sizeof(double));
      CpuTime_MtoParticles = (double *)calloc(1, sizeof(double));
      CpuTime_Kick         = (double *)calloc(1, sizeof(double));
      CpuTime_Drift  = (double *)calloc(1, sizeof(double));
      CpuTime_Output = (double *)calloc(1, sizeof(double));
      WallTime_2LPT = endwall-startwall;
      WallTime_Move    = (double *)calloc(1, sizeof(double));
      WallTime_PtoMesh = (double *)calloc(1, sizeof(double));
      WallTime_Forces       = (double *)calloc(1, sizeof(double));
      WallTime_MtoParticles = (double *)calloc(1, sizeof(double));
      WallTime_Kick         = (double *)calloc(1, sizeof(double));
      WallTime_Drift  = (double *)calloc(1, sizeof(double));
      WallTime_Output = (double *)calloc(1, sizeof(double));

#endif
      goto finalize;
    }

#endif

    NoutputStart++;
  }

  // Now, we get to the N-Body part where we evolve with time via the Kick-Drift-Kick Method
  // =======================================================================================

  // The density grid and force grids  and associated fftw plans
#ifndef MEMORY_MODE
  density = (float_kind *)malloc(2*Total_size*sizeof(float_kind));
  N11  = (float_kind *)malloc(2*Total_size*sizeof(float_kind));
  N12  = (float_kind *)malloc(2*Total_size*sizeof(float_kind));
  N13  = (float_kind *)malloc(2*Total_size*sizeof(float_kind));
  P3D  = (complex_kind*)density;
  FN11 = (complex_kind*)N11;
  FN12 = (complex_kind*)N12;
  FN13 = (complex_kind*)N13;
#ifdef SINGLE_PRECISION
  plan = fftwf_mpi_plan_dft_r2c_3d(Nmesh,Nmesh,Nmesh,density,P3D,MPI_COMM_WORLD,FFTW_ESTIMATE);
  p11  = fftwf_mpi_plan_dft_c2r_3d(Nmesh,Nmesh,Nmesh,FN11,N11,MPI_COMM_WORLD,FFTW_ESTIMATE);
  p12  = fftwf_mpi_plan_dft_c2r_3d(Nmesh,Nmesh,Nmesh,FN12,N12,MPI_COMM_WORLD,FFTW_ESTIMATE);
  p13  = fftwf_mpi_plan_dft_c2r_3d(Nmesh,Nmesh,Nmesh,FN13,N13,MPI_COMM_WORLD,FFTW_ESTIMATE);
#else
  plan = fftw_mpi_plan_dft_r2c_3d(Nmesh,Nmesh,Nmesh,density,P3D,MPI_COMM_WORLD,FFTW_ESTIMATE);
  p11  = fftw_mpi_plan_dft_c2r_3d(Nmesh,Nmesh,Nmesh,FN11,N11,MPI_COMM_WORLD,FFTW_ESTIMATE);
  p12  = fftw_mpi_plan_dft_c2r_3d(Nmesh,Nmesh,Nmesh,FN12,N12,MPI_COMM_WORLD,FFTW_ESTIMATE);
  p13  = fftw_mpi_plan_dft_c2r_3d(Nmesh,Nmesh,Nmesh,FN13,N13,MPI_COMM_WORLD,FFTW_ESTIMATE);
#endif
#endif

#ifdef TIMING
  endcpu = (double)clock();
  endwall = MPI_Wtime();
  CpuTime_2LPT = (endcpu-startcpu)/(double)CLOCKS_PER_SEC;
  WallTime_2LPT = endwall-startwall;
  int nstepstot = 0;
#ifdef LIGHTCONE
  for (i=NoutputStart;i<Noutputs;i++) {
#else
  for (i=NoutputStart;i<=Noutputs;i++) {
#endif
    if (i == Noutputs) {
      nstepstot += 1;
    } else {
      nstepstot += OutputList[i].Nsteps;
    }
  }
  CpuTime_Move    = (double *)calloc(nstepstot, sizeof(double));
  CpuTime_PtoMesh = (double *)calloc(nstepstot, sizeof(double));
  CpuTime_Forces       = (double *)calloc(nstepstot, sizeof(double));
  CpuTime_MtoParticles = (double *)calloc(nstepstot, sizeof(double));
  CpuTime_Kick         = (double *)calloc(nstepstot, sizeof(double));
  CpuTime_Drift  = (double *)calloc(nstepstot, sizeof(double));
  CpuTime_Output = (double *)calloc(nstepstot, sizeof(double));
  WallTime_Move    = (double *)calloc(nstepstot, sizeof(double));
  WallTime_PtoMesh = (double *)calloc(nstepstot, sizeof(double));
  WallTime_Forces       = (double *)calloc(nstepstot, sizeof(double));
  WallTime_MtoParticles = (double *)calloc(nstepstot, sizeof(double));
  WallTime_Kick         = (double *)calloc(nstepstot, sizeof(double));
  WallTime_Drift  = (double *)calloc(nstepstot, sizeof(double));
  WallTime_Output = (double *)calloc(nstepstot, sizeof(double));
#endif

  if(ThisTask == 0) {
    printf("Beginning timestepping\n");
    printf("======================\n");
    fflush(stdout);
  }   

  // Loop over all the timesteps in the timestep list
  // ================================================
  timeSteptot=0;
#ifdef LIGHTCONE
  for (i=NoutputStart;i<Noutputs;i++) {
#else
  for (i=NoutputStart;i<=Noutputs;i++) {
#endif

    int nsteps=0;
    double ao=0;
    if (i == Noutputs) {
      nsteps = 1;
    } else {
      nsteps = OutputList[i].Nsteps;
      ao  = 1.0/(1.0+OutputList[i].Redshift);
      if (StepDist == 0) da=(ao-A)/((double)nsteps);
      if (StepDist == 1) da=(log(ao)-log(A))/((double)nsteps);
    }

    // Perform the required number of timesteps between outputs
    // ========================================================
    for (timeStep=0;timeStep<nsteps;timeStep++) {

#ifdef LIGHTCONE
      // For a lightcone simulation we always want the velocity set to mid-point of interval.
      if (StepDist == 0) AF=A+da*0.5;
      if (StepDist == 1) AF=A*exp(da*0.5);
#else
      // Calculate the time to update to
      // Half timestep for kick at output redshift to bring the velocities and positions to the same time
      if ((timeStep == 0) && (i != NoutputStart)) {
        AF=A; 
      } else { 
        // Set to mid-point of interval. In the infinitesimal timestep limit, these choices are identical. 
        // How one chooses the mid-point when not in that limit is really an extra degree of freedom in the code 
        // but Tassev et al. report negligible effects from the different choices below. 
        // Hence, this is not exported as an extra switch at this point.
        if (StepDist == 0) AF=A+da*0.5;
        if (StepDist == 1) AF=A*exp(da*0.5);
      }
#endif
      if (StepDist == 0) AFF=A+da;
      if (StepDist == 1) AFF=A*exp(da);

      timeSteptot++;
      if (ThisTask == 0) {
        printf("Iteration = %d\n------------------\n",timeSteptot);
        if (i != Noutputs) {
          printf("a = %lf -> a = %lf\n", A, AFF);
          printf("z = %lf -> z = %lf\n", 1.0/A-1.0, 1.0/AFF-1.0);
          fflush(stdout);
        } else {
          printf("Final half timestep to update velocities...\n");
          fflush(stdout);
        }
      }

      // Calculate the particle accelerations for this timestep
      // ======================================================
      GetDisplacements();

      // Kick the particle velocities
      // ============================
      if (ThisTask == 0) {
        printf("Kicking the particles...\n");
        fflush(stdout);
      }

      /**********************************************************************************************/
      // If we wanted to interpolate the lightcone velocities we could put section currently in the // 
      // Drift subroutine here. This would allow us to have the velocities at AF and AFF which we   //
      // can use to invert the previous drift step and get the particle positions at AF and AFF     //                                                                                             //
      /**********************************************************************************************/

#ifdef TIMING
      startcpu = (double)clock();
      startwall = MPI_Wtime();
#endif
      Kick(AI,AF,A,Di,Di2);
#ifdef TIMING
      endcpu = (double)clock();
      endwall = MPI_Wtime();
      CpuTime_Kick[timeSteptot-1] = (endcpu-startcpu)/(double)CLOCKS_PER_SEC;
      WallTime_Kick[timeSteptot-1] = endwall-startwall;
#endif

#ifndef LIGHTCONE

      // If we are at an output timestep we modify the velocities and output the
      // particles. Then, if we are not yet at the end of the simulation, we update  
      // the velocity again up to the middle of the next timestep as per the usual KDK method.
      // =====================================================================================
      if ((timeStep == 0) && (i != NoutputStart)) {

        if (ThisTask == 0) {
          printf("Outputting the particles...\n");
          fflush(stdout);
        }

        // At the output timestep, add back LPT velocities if we had subtracted them. 
        // This corresponds to L_+ operator in Tassev et. al, 2013
        Dv  = QdD1da(A);  // Q*dD_{1}/da
        Dv2 = QdD2da(A);  // Q*dD_{2}/da

#ifdef TIMING
        startcpu = (double)clock();
        startwall = MPI_Wtime();
#endif
        Output(A,OutputList[i-1].Redshift,Dv,Dv2);
#ifdef TIMING
        endcpu = (double)clock();
        endwall = MPI_Wtime();
        CpuTime_Output[timeSteptot-1] = (endcpu-startcpu)/(double)CLOCKS_PER_SEC;
        WallTime_Output[timeSteptot-1] = endwall-startwall;
#endif

        // If we have reached the last output timestep we skip to the end
        if(i == Noutputs) {
          if (ThisTask == 0) {
            printf("Iteration %d finished\n------------------\n\n", timeSteptot);
            fflush(stdout);
          }
          goto finalize;
        }

        // Otherwise we simply update the velocity to the middle of the timestep, where it
        // would have been if we weren't outputting. This involves only recalculating and applying `dda'
        // as the acceleration remains the same as calculated earlier
        AI = A;
        if (StepDist == 0) AF=A+da*0.5;
        if (StepDist == 1) AF=A*exp(da*0.5);

        sumDx=0;
        sumDy=0;
        sumDz=0;

#ifdef TIMING
        startcpu = (double)clock();
        startwall = MPI_Wtime();
#endif
        Kick(AI,AF,A,Di,Di2);
#ifdef TIMING
        endcpu = (double)clock();
        endwall = MPI_Wtime();
        CpuTime_Kick[timeSteptot-1] += (endcpu-startcpu)/(double)CLOCKS_PER_SEC;
        WallTime_Kick[timeSteptot-1] += endwall-startwall;
#endif
      }

#endif

      for (j=0; j<3; j++) free(Disp[j]);

      // Drift the particle positions
      // ============================
      if (ThisTask == 0) {
        printf("Drifting the particles...\n");
        fflush(stdout);
      }

#ifdef TIMING
      startcpu = (double)clock();
      endcpu = MPI_Wtime();
#endif
#ifdef LIGHTCONE
      if (i > 0) {
        Drift_Lightcone(A,AFF,AF,Di,Di2);
      } else {
        Drift(A,AFF,AF,Di,Di2);
      }
#else
      Drift(A,AFF,AF,Di,Di2);
#endif
#ifdef TIMING
      endcpu = (double)clock();
      endwall = MPI_Wtime();
      CpuTime_Drift[timeSteptot-1] = (endcpu-startcpu)/(double)CLOCKS_PER_SEC;
      WallTime_Drift[timeSteptot-1] = endwall-startwall;
#endif

      // Step in time
      // ================
      A  = AFF;
      AI = AF; 

      Di = growthD(A);
      Di2 = growthD2(A);

      if (ThisTask == 0) {
        printf("Iteration %d finished\n------------------\n\n", timeSteptot);
        fflush(stdout);
      }
     
      ierr = MPI_Barrier(MPI_COMM_WORLD);
    }
  }

  // Here is the last little bit
  // ===========================
#ifndef LIGHTCONE
  finalize:
#endif

  if (ThisTask == 0) {
    printf("Finishing up\n");
    printf("============\n");
    fflush(stdout);
  }

#ifdef LIGHTCONE
#ifdef TIMING
  startcpu = (double)clock();
  startwall = MPI_Wtime();
#endif
  Output_Info_Lightcone();
#ifdef TIMING
  endcpu = (double)clock();
  endwall = MPI_Wtime();
  CpuTime_Output[timeSteptot-1] += (endcpu-startcpu)/(double)CLOCKS_PER_SEC;
  WallTime_Output[timeSteptot-1] += endwall-startwall;
#endif
#endif

#ifdef TIMING
  Output_Timing();
#endif

  free_powertable();
  free_transfertable();

  free(P);
  free(OutputList);
  free(Slab_to_task);
  free(Part_to_task);
  free(Local_nx_table);
  free(Local_np_table);
#ifdef GENERIC_FNL
  free(KernelTable);
#endif
#ifdef LIGHTCONE
  free(Noutput);
  free(repflag);
#endif
#ifndef MEMORY_MODE
  free(density);
  free(N11);
  free(N12);
  free(N13);  
#ifdef SINGLE_PRECISION
  fftwf_destroy_plan(plan);
  fftwf_destroy_plan(p11);
  fftwf_destroy_plan(p12);
  fftwf_destroy_plan(p13);
#else
  fftw_destroy_plan(plan);
  fftw_destroy_plan(p11);
  fftw_destroy_plan(p12);
  fftw_destroy_plan(p13);
#endif
#endif

#ifdef SINGLE_PRECISION
  fftwf_mpi_cleanup();
#else
  fftw_mpi_cleanup();
#endif

  if (ThisTask == 0) printf("Done :)\n");

  MPI_Finalize();   

  return 0;
}

// Kicking the particle velocities
// ===============================
void Kick(double AI, double AF, double A, double Di, double Di2) {
  unsigned int n;
  double dda;
  double q1,q2;
  double ax,ay,az;
    
  sumx=0;
  sumy=0;
  sumz=0; 

  if (DeltaA == 0) {
    dda=KickCOLA(AI,AF,A);
  } else if (DeltaA == 1) {
    dda=KickCOLA(AI,AF,A);
  } else if (DeltaA == 2) {
    dda=KickStd(AI,AF);
  } else {
    dda=(AF-AI)*A/Qfactor(A);
  }  
    
  q1=1.5*Omega*Di*A;                                // Q*d/da(Q*dD_{1}/da)
  q2=1.5*Omega*Di2*(1.0+(7.0*Di*Di)/(3.0*Di2))*A;   // Q*d/da(Q*dD_{2}/da) (we have to sort out the stray factors of -3/7 as one of these is absorbed into D2 itself)
  
    
  for(n=0; n<NumPart; n++) {

    Disp[0][n] -= sumDx;
    Disp[1][n] -= sumDy;
    Disp[2][n] -= sumDz;

    ax=-1.5*Omega*Disp[0][n]-UseCOLA*(P[n].Dz[0]*q1+P[n].D2[0]*q2)/A;
    ay=-1.5*Omega*Disp[1][n]-UseCOLA*(P[n].Dz[1]*q1+P[n].D2[1]*q2)/A;
    az=-1.5*Omega*Disp[2][n]-UseCOLA*(P[n].Dz[2]*q1+P[n].D2[2]*q2)/A;

    P[n].Vel[0] += ax*dda;
    P[n].Vel[1] += ay*dda;
    P[n].Vel[2] += az*dda;

    sumx += P[n].Vel[0];
    sumy += P[n].Vel[1];
    sumz += P[n].Vel[2];
  }

  // Make sumx, sumy and sumz global averages
  ierr = MPI_Allreduce(MPI_IN_PLACE,&sumx,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  ierr = MPI_Allreduce(MPI_IN_PLACE,&sumy,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  ierr = MPI_Allreduce(MPI_IN_PLACE,&sumz,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);  
    
  sumx /= (double)TotNumPart;  // We will subtract these to conserve momentum. 
  sumy /= (double)TotNumPart;  // Should be conserved, but just in case 3-linear interpolation makes a problem.
  sumz /= (double)TotNumPart;  // Never checked whether this makes a difference.
}

// Drifting the particle positions
// ===============================
void Drift(double A, double AFF, double AF, double Di, double Di2) {

  unsigned int n;
  double dyyy;
  double da1,da2;
    
  if (DeltaA == 0) {
    dyyy=DriftCOLA(A,AFF,AF);
  } else if (DeltaA == 1) {
    dyyy=DriftStd(A,AFF);
  } else if (DeltaA == 2) {
    dyyy=DriftStd(A,AFF);
  } else {
    dyyy=(AFF-A)/Qfactor(AF);
  } 

  da1=growthD(AFF)-Di;    // change in D
  da2=growthD2(AFF)-Di2;  // change in D_{2lpt}

  for(n=0; n<NumPart; n++) {
    P[n].Pos[0] += (P[n].Vel[0]-sumx)*dyyy;
    P[n].Pos[1] += (P[n].Vel[1]-sumy)*dyyy;
    P[n].Pos[2] += (P[n].Vel[2]-sumz)*dyyy;

    P[n].Pos[0] = periodic_wrap(P[n].Pos[0]+UseCOLA*(P[n].Dz[0]*da1+P[n].D2[0]*da2));
    P[n].Pos[1] = periodic_wrap(P[n].Pos[1]+UseCOLA*(P[n].Dz[1]*da1+P[n].D2[1]*da2));
    P[n].Pos[2] = periodic_wrap(P[n].Pos[2]+UseCOLA*(P[n].Dz[2]*da1+P[n].D2[2]*da2));
  }
}

// Output the data
// ===============
void Output(double A, double Z, double Dv, double Dv2) {

  FILE * fp; 
  char buf[300];
  int nprocgroup, groupTask, masterTask;
  unsigned int n;
  double fac = Hubble/pow(A,1.5);
  double lengthfac = 1.0;   // Keep positions in user-specified units (Originally converted positions to Mpc/h)
  double velfac    = 1.0;   // Keep velocities in user-specified units (Originally converted velocities to km/s)

#ifdef GADGET_STYLE
  size_t bytes;
  int k, pc, dummy, blockmaxlen;
  float * block;
#ifdef PARTICLE_ID
  unsigned long long * blockid;
#endif
#endif

  nprocgroup = NTask / NumFilesWrittenInParallel;
  if (NTask % NumFilesWrittenInParallel) nprocgroup++;
  masterTask = (ThisTask / nprocgroup) * nprocgroup;
  for(groupTask = 0; groupTask < nprocgroup; groupTask++) {
    if (ThisTask == (masterTask + groupTask)) {
      if(NumPart > 0) {
        sprintf(buf, "%s/%s_z%dp%03d.%d", OutputDir, FileBase, (int)Z, (int)rint((Z-(int)Z)*1000), ThisTask);
        if(!(fp = fopen(buf, "w"))) {
          printf("\nERROR: Can't write in file '%s'.\n\n", buf);
          FatalError((char *)"main.c", 735);
        }
        fflush(stdout);
#ifdef GADGET_STYLE
        // Gadget header stuff
        for(k = 0; k < 6; k++) {
          header.npart[k] = 0;
          header.npartTotal[k] = 0;
          header.mass[k] = 0;
        }
        header.npart[1] = NumPart;
        header.npartTotal[1] = TotNumPart;
        header.npartTotal[2] = (TotNumPart >> 32);
        header.mass[1] = (3.0*Omega*Hubble*Hubble*Box*Box*Box) / (8.0*PI*G*TotNumPart);
        header.time = A;
        header.redshift = Z;

        header.flag_sfr = 0;
        header.flag_feedback = 0;
        header.flag_cooling = 0;
        header.flag_stellarage = 0;
        header.flag_metals = 0;
        header.flag_stellarage = 0;
        header.flag_metals = 0;
        header.hashtabsize = 0;

        header.num_files = NTaskWithN;

        header.BoxSize = Box;
        header.Omega0 = Omega;
        header.OmegaLambda = OmegaLambda;
        header.HubbleParam = HubbleParam;

        dummy = sizeof(header);
        my_fwrite(&dummy, sizeof(dummy), 1, fp);
        my_fwrite(&header, sizeof(header), 1, fp);
        my_fwrite(&dummy, sizeof(dummy), 1, fp);

        // We may have some spare memory from deallocating the force grids so use that for outputting
        // If not we are a little more conservative
#ifdef MEMORY_MODE
        block = (float *)malloc(bytes = 6*Total_size*sizeof(float_kind));
#else
        block = (float *)malloc(bytes = 10*1024*1024);
#endif
        blockmaxlen = bytes / (3 * sizeof(float));

        // Remember to add the ZA and 2LPT velocities back on and convert to PTHalos velocity units

        // write coordinates
        dummy = sizeof(float) * 3 * NumPart;
        my_fwrite(&dummy, sizeof(dummy), 1, fp);
        for(n = 0, pc = 0; n < NumPart; n++) {
          for(k = 0; k < 3; k++) block[3 * pc + k] = (float)(lengthfac*P[n].Pos[k]);
          pc++;
          if(pc == blockmaxlen) {
            my_fwrite(block, sizeof(float), 3 * pc, fp);
	    pc = 0;
	  }
        }
        if(pc > 0) my_fwrite(block, sizeof(float), 3 * pc, fp);
        my_fwrite(&dummy, sizeof(dummy), 1, fp);

        // write velocities
        my_fwrite(&dummy, sizeof(dummy), 1, fp);
        for(n = 0, pc = 0; n < NumPart; n++) {
          for(k = 0; k < 3; k++) block[3 * pc + k] = (float)(velfac*fac*(P[n].Vel[k]-sumx+(P[n].Dz[k]*Dv+P[n].D2[k]*Dv2)*UseCOLA));
          pc++;
          if(pc == blockmaxlen) {
            my_fwrite(block, sizeof(float), 3 * pc, fp);
           pc = 0;
          }
        }
        if(pc > 0) my_fwrite(block, sizeof(float), 3 * pc, fp);
        my_fwrite(&dummy, sizeof(dummy), 1, fp);

#ifdef PARTICLE_ID
        blockid = (unsigned long long *)block;
        blockmaxlen = bytes / sizeof(unsigned long long);

        // write particle ID
        dummy = sizeof(unsigned long long) * NumPart;
        my_fwrite(&dummy, sizeof(dummy), 1, fp);
        for(n = 0, pc = 0; n < NumPart; n++) {
          blockid[pc] = P[n].ID;
          pc++;
          if(pc == blockmaxlen) {
            my_fwrite(blockid, sizeof(unsigned long long), pc, fp);
            pc = 0;
          }
        }
        if(pc > 0) my_fwrite(blockid, sizeof(unsigned long long), pc, fp);
        my_fwrite(&dummy, sizeof(dummy), 1, fp);
#endif

        free(block);   
#else
        for(n=0; n<NumPart; n++){
          double P_Vel[3];
          P_Vel[0] = fac*(P[n].Vel[0]-sumx+(P[n].Dz[0]*Dv+P[n].D2[0]*Dv2)*UseCOLA);
          P_Vel[1] = fac*(P[n].Vel[1]-sumy+(P[n].Dz[1]*Dv+P[n].D2[1]*Dv2)*UseCOLA);
          P_Vel[2] = fac*(P[n].Vel[2]-sumz+(P[n].Dz[2]*Dv+P[n].D2[2]*Dv2)*UseCOLA);
#ifdef PARTICLE_ID
          fprintf(fp,"%12llu %15.6f %15.6f %15.6f %15.6f %15.6f %15.6f\n",
                      P[n].ID, (float)(lengthfac*P[n].Pos[0]),(float)(lengthfac*P[n].Pos[1]),(float)(lengthfac*P[n].Pos[2]),(float)(velfac*P_Vel[0]),(float)(velfac*P_Vel[1]),(float)(velfac*P_Vel[2]));
#else
          fprintf(fp,"%15.6f %15.6f %15.6f %15.6f %15.6f %15.6f\n",
                      (float)(lengthfac*P[n].Pos[0]),(float)(lengthfac*P[n].Pos[1]),(float)(lengthfac*P[n].Pos[2]),(float)(velfac*P_Vel[0]),(float)(velfac*P_Vel[1]),(float)(velfac*P_Vel[2]));
#endif
        }
#endif
        fclose(fp);
      }
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }
 
  Output_Info(A, Z);

  return;
}


// Generate the info file which contains a list of all the output files, the 8 corners of the slices on those files and the number of particles in the slice
// =========================================================================================================================================================
void Output_Info(double A, double Z) {

  FILE * fp; 
  char buf[300];
  int i;

  int * Local_p_start_table = (int *)malloc(sizeof(int) * NTask);
  unsigned int * Noutput_table = (unsigned int *)malloc(sizeof(unsigned int) * NTask);
  MPI_Allgather(&Local_p_start, 1, MPI_INT, Local_p_start_table, 1, MPI_INT, MPI_COMM_WORLD);
  MPI_Allgather(&NumPart, 1, MPI_UNSIGNED, Noutput_table, 1, MPI_UNSIGNED, MPI_COMM_WORLD);

  if (ThisTask == 0) {
    sprintf(buf, "%s/%s_z%dp%03d.info", OutputDir, FileBase, (int)Z, (int)rint((Z-(int)Z)*1000));
    if(!(fp = fopen(buf, "w"))) {
      printf("\nERROR: Can't write in file '%s'.\n\n", buf);
      FatalError((char *)"main.c", 876);
    }
    fflush(stdout);
    fprintf(fp, "#    FILENUM      XMIN         YMIN        ZMIN         XMAX         YMAX         ZMAX         NPART    \n");
    double y0 = 0.0;
    double z0 = 0.0;
    double y1 = Box;
    double z1 = Box;
    for (i=0; i<NTask; i++) {
      double x0 = Local_p_start_table[i]*(Box/(double)Nsample); 
      double x1 = (Local_p_start_table[i]+Local_np_table[i])*(Box/(double)Nsample);
      fprintf(fp, "%12d %15.6lf %15.6lf %15.6lf %15.6lf %15.6lf %15.6lf %12u\n", i, x0, y0, z0, x1, y1, z1, Noutput_table[i]);
    }
    fclose(fp);
  }

  free(Noutput_table);
  free(Local_p_start_table);

  return;
}

#ifdef TIMING
void Output_Timing(void) {

  FILE * fp; 
  char buf[300];
  int i;

  if (ThisTask == 0) {
    sprintf(buf, "%s/%s_cputime_0.dat", OutputDir, FileBase);
    if(!(fp = fopen(buf, "w"))) {
      printf("\nERROR: Can't write in file '%s'.\n\n", buf);
      FatalError((char *)"main.c", 909);
    }
    fflush(stdout);
    fprintf(fp, "Initialisation Time:\n");
    fprintf(fp, "%12.6lf\n\n", CpuTime_Init);
    fprintf(fp, "2LPT Time        (Non_Gaussian Kernel        Output):\n");
    fprintf(fp, "%12.6lf     %12.6lf     %12.6lf\n\n", CpuTime_2LPT, CpuTime_2LPTng, CpuTime_2LPToutput);
    fprintf(fp, "Timestep Times:\n");
    fprintf(fp, "Move          PtoMesh          Forces          MtoParticles          Kick           Drift           Output\n");
    for (i=0; i<timeSteptot; i++) {
      fprintf(fp, "%12.6lf     %12.6lf     %12.6lf     %12.6lf     %12.6lf     %12.6lf     %12.6lf\n", 
                   CpuTime_Move[i], CpuTime_PtoMesh[i], CpuTime_Forces[i], CpuTime_MtoParticles[i], CpuTime_Kick[i], CpuTime_Drift[i], CpuTime_Output[i]);
    }
    fclose(fp);
    sprintf(buf, "%s/%s_walltime_0.dat", OutputDir, FileBase);
    if(!(fp = fopen(buf, "w"))) {
      printf("\nERROR: Can't write in file '%s'.\n\n", buf);
      FatalError((char *)"main.c", 926);
    }
    fflush(stdout);
    fprintf(fp, "Initialisation Time:\n");
    fprintf(fp, "%12.6lf\n\n", WallTime_Init);
    fprintf(fp, "2LPT Time        (Non_Gaussian Kernel        Output):\n");
    fprintf(fp, "%12.6lf     %12.6lf     %12.6lf\n\n", WallTime_2LPT, WallTime_2LPTng, WallTime_2LPToutput);
    fprintf(fp, "Timestep Times:\n");
    fprintf(fp, "Move          PtoMesh          Forces          MtoParticles          Kick           Drift           Output\n");
    for (i=0; i<timeSteptot; i++) {
      fprintf(fp, "%12.6lf     %12.6lf     %12.6lf     %12.6lf     %12.6lf     %12.6lf     %12.6lf\n", 
                   WallTime_Move[i], WallTime_PtoMesh[i], WallTime_Forces[i], WallTime_MtoParticles[i], WallTime_Kick[i], WallTime_Drift[i], WallTime_Output[i]);
    }
    fclose(fp);
    MPI_Reduce(MPI_IN_PLACE, &CpuTime_Init, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(MPI_IN_PLACE, &CpuTime_2LPT, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);     
    MPI_Reduce(MPI_IN_PLACE, &CpuTime_2LPTng, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);   
    MPI_Reduce(MPI_IN_PLACE, &CpuTime_2LPToutput, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);     
    MPI_Reduce(MPI_IN_PLACE, &(CpuTime_Move[0]), timeSteptot, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);     
    MPI_Reduce(MPI_IN_PLACE, &(CpuTime_PtoMesh[0]), timeSteptot, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(MPI_IN_PLACE, &(CpuTime_Forces[0]), timeSteptot, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);         
    MPI_Reduce(MPI_IN_PLACE, &(CpuTime_MtoParticles[0]), timeSteptot, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);          
    MPI_Reduce(MPI_IN_PLACE, &(CpuTime_Kick[0]), timeSteptot, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);          
    MPI_Reduce(MPI_IN_PLACE, &(CpuTime_Drift[0]), timeSteptot, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);          
    MPI_Reduce(MPI_IN_PLACE, &(CpuTime_Output[0]), timeSteptot, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(MPI_IN_PLACE, &WallTime_Init, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(MPI_IN_PLACE, &WallTime_2LPT, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);     
    MPI_Reduce(MPI_IN_PLACE, &WallTime_2LPTng, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);   
    MPI_Reduce(MPI_IN_PLACE, &WallTime_2LPToutput, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);     
    MPI_Reduce(MPI_IN_PLACE, &(WallTime_Move[0]), timeSteptot, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);     
    MPI_Reduce(MPI_IN_PLACE, &(WallTime_PtoMesh[0]), timeSteptot, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(MPI_IN_PLACE, &(WallTime_Forces[0]), timeSteptot, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);         
    MPI_Reduce(MPI_IN_PLACE, &(WallTime_MtoParticles[0]), timeSteptot, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);          
    MPI_Reduce(MPI_IN_PLACE, &(WallTime_Kick[0]), timeSteptot, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);          
    MPI_Reduce(MPI_IN_PLACE, &(WallTime_Drift[0]), timeSteptot, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);          
    MPI_Reduce(MPI_IN_PLACE, &(WallTime_Output[0]), timeSteptot, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    sprintf(buf, "%s/%s_cputime_all.dat", OutputDir, FileBase);
    if(!(fp = fopen(buf, "w"))) {
      printf("\nERROR: Can't write in file '%s'.\n\n", buf);
      FatalError((char *)"main.c", 965);
    }
    fflush(stdout);
    fprintf(fp, "Initialisation Time:\n");
    fprintf(fp, "%12.6lf\n\n", CpuTime_Init);
    fprintf(fp, "2LPT Time        (Non_Gaussian Kernel        Output):\n");
    fprintf(fp, "%12.6lf     %12.6lf     %12.6lf\n\n", CpuTime_2LPT, CpuTime_2LPTng, CpuTime_2LPToutput);
    fprintf(fp, "Timestep Times:\n");
    fprintf(fp, "Move          PtoMesh          Forces          MtoParticles          Kick           Drift           Output\n");
    for (i=0; i<timeSteptot; i++) {
      fprintf(fp, "%12.6lf     %12.6lf     %12.6lf     %12.6lf     %12.6lf     %12.6lf     %12.6lf\n", 
                   CpuTime_Move[i], CpuTime_PtoMesh[i], CpuTime_Forces[i], CpuTime_MtoParticles[i], CpuTime_Kick[i], CpuTime_Drift[i], CpuTime_Output[i]);
    }
    fclose(fp);
    sprintf(buf, "%s/%s_walltime_all.dat", OutputDir, FileBase);
    if(!(fp = fopen(buf, "w"))) {
      printf("\nERROR: Can't write in file '%s'.\n\n", buf);
      FatalError((char *)"main.c", 982);
    }
    fflush(stdout);
    fprintf(fp, "Initialisation Time:\n");
    fprintf(fp, "%12.6lf\n\n", WallTime_Init);
    fprintf(fp, "2LPT Time        (Non_Gaussian Kernel        Output):\n");
    fprintf(fp, "%12.6lf     %12.6lf     %12.6lf\n\n", WallTime_2LPT, WallTime_2LPTng, WallTime_2LPToutput);
    fprintf(fp, "Timestep Times:\n");
    fprintf(fp, "Move          PtoMesh          Forces          MtoParticles          Kick           Drift           Output\n");
    for (i=0; i<timeSteptot; i++) {
      fprintf(fp, "%12.6lf     %12.6lf     %12.6lf     %12.6lf     %12.6lf     %12.6lf     %12.6lf\n", 
                   WallTime_Move[i], WallTime_PtoMesh[i], WallTime_Forces[i], WallTime_MtoParticles[i], WallTime_Kick[i], WallTime_Drift[i], WallTime_Output[i]);
    }
    fclose(fp);
  } else {
    MPI_Reduce(&CpuTime_Init, &CpuTime_Init, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&CpuTime_2LPT, &CpuTime_2LPT, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);     
    MPI_Reduce(&CpuTime_2LPTng, &CpuTime_2LPTng, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);   
    MPI_Reduce(&CpuTime_2LPToutput, &CpuTime_2LPToutput, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);     
    MPI_Reduce(&(CpuTime_Move[0]), &(CpuTime_Move[0]), timeSteptot, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);     
    MPI_Reduce(&(CpuTime_PtoMesh[0]), &(CpuTime_PtoMesh[0]), timeSteptot, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&(CpuTime_Forces[0]), &(CpuTime_Forces[0]), timeSteptot, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);         
    MPI_Reduce(&(CpuTime_MtoParticles[0]), &(CpuTime_MtoParticles[0]), timeSteptot, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);          
    MPI_Reduce(&(CpuTime_Kick[0]), &(CpuTime_Kick[0]), timeSteptot, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);          
    MPI_Reduce(&(CpuTime_Drift[0]), &(CpuTime_Drift[0]), timeSteptot, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);          
    MPI_Reduce(&(CpuTime_Output[0]), &(CpuTime_Output[0]), timeSteptot, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&WallTime_Init, &WallTime_Init, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&WallTime_2LPT, &WallTime_2LPT, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);     
    MPI_Reduce(&WallTime_2LPTng, &WallTime_2LPTng, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);   
    MPI_Reduce(&WallTime_2LPToutput, &WallTime_2LPToutput, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);     
    MPI_Reduce(&(WallTime_Move[0]), &(WallTime_Move[0]), timeSteptot, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);     
    MPI_Reduce(&(WallTime_PtoMesh[0]), &(WallTime_PtoMesh[0]), timeSteptot, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&(WallTime_Forces[0]), &(WallTime_Forces[0]), timeSteptot, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);         
    MPI_Reduce(&(WallTime_MtoParticles[0]), &(WallTime_MtoParticles[0]), timeSteptot, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);          
    MPI_Reduce(&(WallTime_Kick[0]), &(WallTime_Kick[0]), timeSteptot, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);          
    MPI_Reduce(&(WallTime_Drift[0]), &(WallTime_Drift[0]), timeSteptot, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);          
    MPI_Reduce(&(WallTime_Output[0]), &(WallTime_Output[0]), timeSteptot, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);                         
  }

  free(CpuTime_Move);
  free(CpuTime_PtoMesh);
  free(CpuTime_Forces);
  free(CpuTime_MtoParticles);
  free(CpuTime_Kick);
  free(CpuTime_Drift);
  free(CpuTime_Output);
  free(WallTime_Move);
  free(WallTime_PtoMesh);
  free(WallTime_Forces);
  free(WallTime_MtoParticles);
  free(WallTime_Kick);
  free(WallTime_Drift);
  free(WallTime_Output);

}
#endif

