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

/* =========================================================================================*/
/* This file contains some additional routines (parallel and serial) needed for any PM code.*/
/* =========================================================================================*/

#include "vars.h"
#include "proto.h"

// A master routine called from main.c to calculate the acceleration
// =================================================================
void GetDisplacements(void) {

  int j;
#ifdef TIMING
  double startcpu, endcpu;
  double startwall, endwall;
#endif

  // First we check whether all the particles are on the correct processor after the last time step/
  // original 2LPT displacement and move them if not
  if (ThisTask == 0) printf("Moving particles across task boundaries...\n");
#ifdef TIMING
  startcpu = (double)clock();
  startwall = MPI_Wtime();
#endif
  MoveParticles();
#ifdef TIMING
  endcpu = (double)clock();
  endwall = MPI_Wtime();
  CpuTime_Move[timeSteptot-1] = (endcpu-startcpu)/(double)CLOCKS_PER_SEC;  
  WallTime_Move[timeSteptot-1] = endwall-startwall;
#endif

#ifdef MEMORY_MODE
  density = (float_kind *)malloc(2*Total_size*sizeof(float_kind));
  P3D = (complex_kind*)density;
#ifdef SINGLE_PRECISION
  plan = fftwf_mpi_plan_dft_r2c_3d(Nmesh,Nmesh,Nmesh,density,P3D,MPI_COMM_WORLD,FFTW_ESTIMATE);
#else
  plan = fftw_mpi_plan_dft_r2c_3d(Nmesh,Nmesh,Nmesh,density,P3D,MPI_COMM_WORLD,FFTW_ESTIMATE);
#endif
#endif

  // Then we do the Cloud-in-Cell assignment to get the density grid and FFT it.  
  if (ThisTask == 0) printf("Calculating density using Cloud-in-Cell...\n");
#ifdef TIMING
  startcpu = (double)clock();
  startwall = MPI_Wtime();
#endif
  PtoMesh();
#ifdef TIMING
  endcpu = (double)clock();
  endwall = MPI_Wtime();
  CpuTime_PtoMesh[timeSteptot-1] = (endcpu-startcpu)/(double)CLOCKS_PER_SEC;  
  WallTime_PtoMesh[timeSteptot-1] = endwall-startwall;
#endif

#ifdef MEMORY_MODE
  N11  = (float_kind *)malloc(2*Total_size*sizeof(float_kind));
  N12  = (float_kind *)malloc(2*Total_size*sizeof(float_kind));
  N13  = (float_kind *)malloc(2*Total_size*sizeof(float_kind));
  FN11 = (complex_kind*)N11;
  FN12 = (complex_kind*)N12;
  FN13 = (complex_kind*)N13;
#ifdef SINGLE_PRECISION
  p11  = fftwf_mpi_plan_dft_c2r_3d(Nmesh,Nmesh,Nmesh,FN11,N11,MPI_COMM_WORLD,FFTW_ESTIMATE);
  p12  = fftwf_mpi_plan_dft_c2r_3d(Nmesh,Nmesh,Nmesh,FN12,N12,MPI_COMM_WORLD,FFTW_ESTIMATE);
  p13  = fftwf_mpi_plan_dft_c2r_3d(Nmesh,Nmesh,Nmesh,FN13,N13,MPI_COMM_WORLD,FFTW_ESTIMATE);
#else
  p11  = fftw_mpi_plan_dft_c2r_3d(Nmesh,Nmesh,Nmesh,FN11,N11,MPI_COMM_WORLD,FFTW_ESTIMATE);
  p12  = fftw_mpi_plan_dft_c2r_3d(Nmesh,Nmesh,Nmesh,FN12,N12,MPI_COMM_WORLD,FFTW_ESTIMATE);
  p13  = fftw_mpi_plan_dft_c2r_3d(Nmesh,Nmesh,Nmesh,FN13,N13,MPI_COMM_WORLD,FFTW_ESTIMATE);
#endif
#endif
    
  // This returns N11,N12,N13 which hold the components of
  // the vector (grad grad^{-2} density) on a grid.
  if (ThisTask == 0) printf("Calculating forces...\n");
#ifdef TIMING
  startcpu = (double)clock();
  startwall = MPI_Wtime();
#endif
  Forces();
#ifdef TIMING
  endcpu = (double)clock();
  endwall = MPI_Wtime();
  CpuTime_Forces[timeSteptot-1] = (endcpu-startcpu)/(double)CLOCKS_PER_SEC;  
  WallTime_Forces[timeSteptot-1] = endwall-startwall;
#endif 

#ifdef MEMORY_MODE
  free(density);
  for (j=0; j<3; j++) Disp[j] = (float *)malloc(NumPart*sizeof(float));
#ifdef SINGLE_PRECISION
  fftwf_destroy_plan(plan);
#else 
  fftw_destroy_plan(plan);
#endif
#else
  for (j=0; j<3; j++) Disp[j] = (float_kind *)malloc(NumPart*sizeof(float_kind));
#endif
    
  // Now find the accelerations at the particle positions using 3-linear interpolation. 
  if (ThisTask == 0) printf("Calculating accelerations...\n");
#ifdef TIMING
  startcpu = (double)clock();
  startwall = MPI_Wtime();
#endif
  MtoParticles();
#ifdef TIMING
  endcpu = (double)clock();
  endwall = MPI_Wtime();
  CpuTime_MtoParticles[timeSteptot-1] = (endcpu-startcpu)/(double)CLOCKS_PER_SEC;  
  WallTime_MtoParticles[timeSteptot-1] = endwall-startwall;
#endif

#ifdef MEMORY_MODE
  free(N11);
  free(N12);
  free(N13);  
#ifdef SINGLE_PRECISION
  fftwf_destroy_plan(p11);
  fftwf_destroy_plan(p12);
  fftwf_destroy_plan(p13);
#else
  fftw_destroy_plan(p11);
  fftw_destroy_plan(p12);
  fftw_destroy_plan(p13);
#endif
#endif
}

// A routine to check whether all the particles are on the correct processor and move them if not.
// ==============================================================================================
void MoveParticles(void) {

  // Note that there are some subtleties in this routine that deal with the fact in some instances there
  // may be no particles on the last N tasks depending on how the work is partioned, hence we need to 
  // skip over these tasks and copy to the correct ones. We include subtleties that deal with the fact that
  // a task may have no particles by skipping over them from the other tasks perspective and
  // setting any sendrecv commands on these tasks to null
 
  int X;
  int j;
  int neighbour, neighbour_left, neighbour_right, neighbour_count = 0;
  int send_count_max = (int)(ceil(Local_np*Nsample*Nsample*(Buffer-1.0)));
  int send_count_left = 0, send_count_right = 0;
  int recv_count_left = 0, recv_count_right = 0;
  int procdiff_left, procdiff_right, procdiffmax = 1, procdiffmaxglob = 1;
  unsigned int i;
  double scaleBox=(double)Nmesh/Box;
 
  // We assume that at least one send is needed and calculate the true number of sends needed in the first iteration.
  // (Yes, i know we shouldn't really modify the iteration counter inside the loop but it creates a good algorithm both here and
  // when we assign the particles to be copied) 
  for (j=1;j<=procdiffmaxglob;j++) {

    // Allocate memory to hold the particles to be transfered. We assume a maximum of Local_np*Nsample*Nsample*(buffer-1.0).
    struct part_data * P_send_left  = (struct part_data *)malloc(send_count_max*sizeof(struct part_data));
    struct part_data * P_send_right = (struct part_data *)malloc(send_count_max*sizeof(struct part_data));

    // The main purpose here is to calculate how many sendrecvs we need to perform (i.e., the maximum number 
    // of tasks a particle has moved across). However, we also assume that at least one send is needed 
    // and so set up the particles to be transferred to the neighbouring tasks
    send_count_left = 0; send_count_right = 0;
    recv_count_left = 0; recv_count_right = 0;
    if (j <= procdiffmax) {
      for (i=0;i<NumPart;i++) {
        X=(int)(P[i].Pos[0]*scaleBox);
        procdiff_left=0; procdiff_right=0;
        if (Slab_to_task[X] != ThisTask) {
          neighbour = ThisTask;
          do {
            procdiff_left++;
            neighbour--;
            if (neighbour < 0) neighbour += NTask;
            if (Local_np_table[neighbour] == 0) procdiff_left--;
          } while(Slab_to_task[X] != neighbour);
          neighbour = ThisTask;
          do {
            procdiff_right++;
            neighbour++;
            if (neighbour >= NTask) neighbour -= NTask;
            if (Local_np_table[neighbour] == 0) procdiff_right--;
          } while(Slab_to_task[X] != neighbour);
          if ((procdiff_left != 0) || (procdiff_right != 0)) {
            if (procdiff_left <= procdiff_right) {
              if (j == 1) {
                if (procdiff_left > procdiffmax) procdiffmax = procdiff_left;
              }
              if (procdiff_left == j) {
                P_send_left[send_count_left] = P[i];
                P[i] = P[NumPart-1];
                i--; NumPart--;
                send_count_left++;
                if (send_count_left >= send_count_max) {
                  printf("\nERROR: Number of particles to be sent left on task %d is greater than send_count_max\n", ThisTask);
                  printf("       You must increase the size of the buffer region.\n\n");
                  FatalError((char *)"auxPM.c", 219);
                }
              }
            } else {
              if (j == 1) {
                if (procdiff_right > procdiffmax) procdiffmax = procdiff_right;
              }
              if (procdiff_right == j) {
                P_send_right[send_count_right] = P[i];
                P[i] = P[NumPart-1];
                i--; NumPart--;
                send_count_right++;
                if (send_count_right >= send_count_max) {
                  printf("\nERROR: Number of particles to be sent right on task %d is greater than send_count_max\n", ThisTask);
                  printf("       You must increase the size of the buffer region.\n\n");
                  FatalError((char *)"auxPM.c", 234);
                }
              }
            }
          }
        }
      }
    } 

    // If we have to send to non-adjoining tasks then we have to recompute the neighbour's task number. For adjoining tasks 
    // we have already got these in the variables LeftTask and RightTask which are also used elsewhere
    if (j == 1) {
      neighbour_left = LeftTask;
      neighbour_right = RightTask;      
      ierr = MPI_Allreduce(&procdiffmax, &procdiffmaxglob, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
      if (ThisTask == 0) printf("Need to transfer particles %d times...\n", procdiffmaxglob);
    } else {
      if (Local_np == 0) {
        neighbour_left = MPI_PROC_NULL;
        neighbour_right = MPI_PROC_NULL;
      } else {

        neighbour_count = 0;
        neighbour_left = ThisTask;
        do {
          neighbour_left--;
          neighbour_count++;
          if(neighbour_left < 0) neighbour_left += NTask;
          if(Local_np_table[neighbour_left] == 0) neighbour_count--;
        } while(neighbour_count != j);

        neighbour_count = 0;
        neighbour_right = ThisTask;
        do {
          neighbour_right++;
          neighbour_count++;
          if(neighbour_right >= NTask) neighbour_right -= NTask;
          if(Local_np_table[neighbour_right] == 0) neighbour_count--;
        } while(neighbour_count != j);
      }
    }

    ierr = MPI_Sendrecv(&send_count_left,1,MPI_INT,neighbour_left,0,&recv_count_right,1,MPI_INT,neighbour_right,0,MPI_COMM_WORLD,&status);
    ierr = MPI_Sendrecv(&send_count_right,1,MPI_INT,neighbour_right,0,&recv_count_left,1,MPI_INT,neighbour_left,0,MPI_COMM_WORLD,&status);

    if (NumPart+recv_count_left+recv_count_right > Local_np*Nsample*Nsample*Buffer) {
      printf("\nERROR: Number of particles to be recieved on task %d is greater than available space\n", ThisTask);
      printf("       You must increase the size of the buffer region.\n\n");
      FatalError((char *)"auxPM.c", 282);
    }

    // Copy across the new particles and store them at the end (of the memory). Then modify NumPart to include them.
    ierr = MPI_Sendrecv(&(P_send_left[0]),send_count_left*sizeof(struct part_data),MPI_BYTE,neighbour_left,0,
                        &(P[NumPart]),recv_count_right*sizeof(struct part_data),MPI_BYTE,neighbour_right,0,MPI_COMM_WORLD,&status);
    ierr = MPI_Sendrecv(&(P_send_right[0]),send_count_right*sizeof(struct part_data),MPI_BYTE,neighbour_right,0,
                        &(P[NumPart+recv_count_right]),recv_count_left*sizeof(struct part_data),MPI_BYTE,neighbour_left,0,MPI_COMM_WORLD,&status);

    NumPart += (recv_count_left+recv_count_right);

    free(P_send_left);
    free(P_send_right);
  }
  return;  
}

// Does Cloud-in-Cell assignment.
// ==============================
void PtoMesh(void) {
      
  unsigned int i;
  unsigned int IX,IY,IZ;
  unsigned int IXneigh,IYneigh,IZneigh;
  double X,Y,Z;
  double TX,TY,TZ;
  double DX,DY,DZ;
  double scaleBox=(double)Nmesh/Box;
  double WPAR=pow((double)Nmesh/(double)Nsample,3);

  for(i=0;i<2*Total_size;i++) density[i] = -1.0;
  for(i=0;i<NumPart;i++) {
     
    X=P[i].Pos[0]*scaleBox;
    Y=P[i].Pos[1]*scaleBox;
    Z=P[i].Pos[2]*scaleBox;

    IX=(unsigned int)X;
    IY=(unsigned int)Y;
    IZ=(unsigned int)Z;
    DX=X-(double)IX;
    DY=Y-(double)IY;
    DZ=Z-(double)IZ;
    TX=1.0-DX;
    TY=1.0-DY;
    TZ=1.0-DZ;

    DY *= WPAR;
    TY *= WPAR;
            
    IX -= Local_x_start;
    if(IY >= (unsigned int)Nmesh) IY=0;
    if(IZ >= (unsigned int)Nmesh) IZ=0;

    IXneigh=IX+1;
    IYneigh=IY+1;
    IZneigh=IZ+1;
    if(IYneigh >= (unsigned int)Nmesh) IYneigh=0;
    if(IZneigh >= (unsigned int)Nmesh) IZneigh=0;

    density[(IX*Nmesh+IY)*2*(Nmesh/2+1)+IZ]           += TX*TY*TZ;
    density[(IX*Nmesh+IY)*2*(Nmesh/2+1)+IZneigh]      += TX*TY*DZ;
    density[(IX*Nmesh+IYneigh)*2*(Nmesh/2+1)+IZ]      += TX*DY*TZ;
    density[(IX*Nmesh+IYneigh)*2*(Nmesh/2+1)+IZneigh] += TX*DY*DZ;

    density[(IXneigh*Nmesh+IY)*2*(Nmesh/2+1)+IZ]           += DX*TY*TZ;
    density[(IXneigh*Nmesh+IY)*2*(Nmesh/2+1)+IZneigh]      += DX*TY*DZ;
    density[(IXneigh*Nmesh+IYneigh)*2*(Nmesh/2+1)+IZ]      += DX*DY*TZ;
    density[(IXneigh*Nmesh+IYneigh)*2*(Nmesh/2+1)+IZneigh] += DX*DY*DZ;
  }

  // Copy across the extra slice from the task on the left and add it to the leftmost slice
  // of the task on the right. Skip over tasks without any slices.
  float_kind * temp_density = (float_kind *)calloc(2*alloc_slice,sizeof(float_kind));

  ierr = MPI_Sendrecv(&(density[2*last_slice]),2*alloc_slice*sizeof(float_kind),MPI_BYTE,RightTask,0,
                      &(temp_density[0]),2*alloc_slice*sizeof(float_kind),MPI_BYTE,LeftTask,0,MPI_COMM_WORLD,&status);

  if (NumPart != 0) {
    for (i=0;i<2*alloc_slice;i++) density[i] += (temp_density[i]+1.0);
  }

  free(temp_density);

  // FFT the density field
#ifdef SINGLE_PRECISION
  fftwf_execute(plan);
#else
  fftw_execute(plan);
#endif

  return;
}

// Calculate the force grids from the density.
// ===========================================
void Forces(void) {

  int di, dj, dk;
  int iglobal, jrev, kmin;
  unsigned int i, j, k;
  double Scale=(2.0*PI)/Box;
  double RK, KK, grid_corr;
  double sinc_x, sinc_y, sinc_z;
  complex_kind dens;

  // We need global values for i as opposed to local values
  // Same goes for anything that relies on i (such as RK). 
  for (i=0;i<Local_nx;i++) {
    iglobal = i+Local_x_start;
    for (j=0;j<(unsigned int)(Nmesh/2+1);j++) {
      kmin = 0;
      if ((iglobal == 0) && (j == 0)) {
        FN11[0][0] = 0.0; FN11[0][1] = 0.0;
        FN12[0][0] = 0.0; FN12[0][1] = 0.0;
        FN13[0][0] = 0.0; FN13[0][1] = 0.0;
        kmin = 1;
      }
      for (k=kmin;k<(unsigned int)(Nmesh/2+1);k++) {

        unsigned int ind = (i*Nmesh+j)*(Nmesh/2+1)+k;
        if (iglobal > Nmesh/2) {
          di = iglobal-Nmesh;
          RK    = (double)(k*k+(Nmesh-iglobal)*(Nmesh-iglobal)+j*j);
        } else {
          di = iglobal;
          RK    = (double)(k*k+iglobal*iglobal+j*j) ;
        }
        dj = j;
        dk = k;

        // Deconvolve the CIC window function twice (once for density, once for force interpolation)
        // and add gaussian smoothing if requested

        KK = -1.0/RK;

        sinc_x = sinc_y = sinc_z = 1.0;
        if (di != 0) sinc_x = sin((PI*di)/(double)Nmesh)/((PI*di)/(double)Nmesh);
        if (dj != 0) sinc_y = sin((PI*dj)/(double)Nmesh)/((PI*dj)/(double)Nmesh);
        if (dk != 0) sinc_z = sin((PI*dk)/(double)Nmesh)/((PI*dk)/(double)Nmesh);
        grid_corr = 1.0/(sinc_x*sinc_y*sinc_z);
        grid_corr = pow(grid_corr, 4.0);
        grid_corr = 1.0;

        dens[0] = (     P3D[ind][0]*KK*grid_corr)/pow((double)Nmesh,3);
        dens[1] = (-1.0*P3D[ind][1]*KK*grid_corr)/pow((double)Nmesh,3);

        // dens now holds the potential so we can solve for the force. 

        FN11[ind][0] = dens[1]*di/Scale;
        FN11[ind][1] = dens[0]*di/Scale;
        FN12[ind][0] = dens[1]*dj/Scale;
        FN12[ind][1] = dens[0]*dj/Scale;
        FN13[ind][0] = dens[1]*dk/Scale;
        FN13[ind][1] = dens[0]*dk/Scale;
            
        // Do the mirror force along the y axis
        if ((j != (unsigned int)(Nmesh/2)) && (j != 0)) {
          jrev=Nmesh-j;

          int ind = (i*Nmesh+jrev)*(Nmesh/2+1)+k;

          dj = -j;

          // Deconvolve the CIC window function twice (once for density, once for force interpolation)
          // and add gaussian smoothing if requested

          sinc_y = 1.0;
          if (dj != 0) sinc_y = sin((PI*dj)/(double)Nmesh)/((PI*dj)/(double)Nmesh);
          grid_corr = 1.0/(sinc_x*sinc_y*sinc_z);
          grid_corr = pow(grid_corr, 4.0);
          grid_corr = 1.0;

          dens[0] = (     P3D[ind][0]*KK*grid_corr)/pow((double)Nmesh,3);
          dens[1] = (-1.0*P3D[ind][1]*KK*grid_corr)/pow((double)Nmesh,3);

          // dens now holds the potential so we can solve for the force. 
          
          FN11[ind][0] = dens[1]*di/Scale;
          FN11[ind][1] = dens[0]*di/Scale;
          FN12[ind][0] = dens[1]*dj/Scale;
          FN12[ind][1] = dens[0]*dj/Scale;
          FN13[ind][0] = dens[1]*dk/Scale;
          FN13[ind][1] = dens[0]*dk/Scale;       
        }
      }
    }
  }
   
#ifdef SINGLE_PRECISION  
  fftwf_execute(p11);
  fftwf_execute(p12);
  fftwf_execute(p13);
#else
  fftw_execute(p11);
  fftw_execute(p12);
  fftw_execute(p13);
#endif

  // Copy across the extra slice from the process on the right and save it at the 
  // end of the force array. Skip over tasks without any slices.
  ierr = MPI_Sendrecv(&(N11[0]),2*alloc_slice*sizeof(float_kind),MPI_BYTE,LeftTask,0,
                      &(N11[2*last_slice]),2*alloc_slice*sizeof(float_kind),MPI_BYTE,RightTask,0,MPI_COMM_WORLD,&status);

  ierr = MPI_Sendrecv(&(N12[0]),2*alloc_slice*sizeof(float_kind),MPI_BYTE,LeftTask,0,
                      &(N12[2*last_slice]),2*alloc_slice*sizeof(float_kind),MPI_BYTE,RightTask,0,MPI_COMM_WORLD,&status);

  ierr = MPI_Sendrecv(&(N13[0]),2*alloc_slice*sizeof(float_kind),MPI_BYTE,LeftTask,0,
                      &(N13[2*last_slice]),2*alloc_slice*sizeof(float_kind),MPI_BYTE,RightTask,0,MPI_COMM_WORLD,&status);

  return;
}

// Does 3-linear interpolation
// ===========================
void MtoParticles(void) {

  unsigned int i;
  unsigned int IX,IY,IZ;
  unsigned int IXneigh,IYneigh,IZneigh;
  double X,Y,Z;
  double TX,TY,TZ;
  double DX,DY,DZ;
  double scaleBox=(double)Nmesh/Box;
  double WPAR=1;

  sumDx=0; 
  sumDy=0; 
  sumDz=0; 

  for(i=0; i<NumPart; i++) {

    X=P[i].Pos[0]*scaleBox;
    Y=P[i].Pos[1]*scaleBox;
    Z=P[i].Pos[2]*scaleBox;

    IX=(unsigned int)X;
    IY=(unsigned int)Y;
    IZ=(unsigned int)Z;
    DX=X-(double)IX;
    DY=Y-(double)IY;
    DZ=Z-(double)IZ;
    TX=1.0-DX;
    TY=1.0-DY;
    TZ=1.0-DZ;

    DY *= WPAR;
    TY *= WPAR;
            
    IX -= Local_x_start;
    if(IY >= (unsigned int)Nmesh) IY=0;
    if(IZ >= (unsigned int)Nmesh) IZ=0;

    IXneigh=IX+1;
    IYneigh=IY+1;
    IZneigh=IZ+1;
    if(IYneigh >= (unsigned int)Nmesh) IYneigh=0;
    if(IZneigh >= (unsigned int)Nmesh) IZneigh=0;

    Disp[0][i] = N11[(IX*Nmesh+IY)*2*(Nmesh/2+1)+IZ]          *TX*TY*TZ +
                 N11[(IX*Nmesh+IY)*2*(Nmesh/2+1)+IZneigh]     *TX*TY*DZ +
                 N11[(IX*Nmesh+IYneigh)*2*(Nmesh/2+1)+IZ]     *TX*DY*TZ +
                 N11[(IX*Nmesh+IYneigh)*2*(Nmesh/2+1)+IZneigh]*TX*DY*DZ +
                 N11[(IXneigh*Nmesh+IY)*2*(Nmesh/2+1)+IZ]          *DX*TY*TZ +
                 N11[(IXneigh*Nmesh+IY)*2*(Nmesh/2+1)+IZneigh]     *DX*TY*DZ +
                 N11[(IXneigh*Nmesh+IYneigh)*2*(Nmesh/2+1)+IZ]     *DX*DY*TZ +
                 N11[(IXneigh*Nmesh+IYneigh)*2*(Nmesh/2+1)+IZneigh]*DX*DY*DZ;

    Disp[1][i] = N12[(IX*Nmesh+IY)*2*(Nmesh/2+1)+IZ]          *TX*TY*TZ +
                 N12[(IX*Nmesh+IY)*2*(Nmesh/2+1)+IZneigh]     *TX*TY*DZ +
                 N12[(IX*Nmesh+IYneigh)*2*(Nmesh/2+1)+IZ]     *TX*DY*TZ +
                 N12[(IX*Nmesh+IYneigh)*2*(Nmesh/2+1)+IZneigh]*TX*DY*DZ +
                 N12[(IXneigh*Nmesh+IY)*2*(Nmesh/2+1)+IZ]          *DX*TY*TZ +
                 N12[(IXneigh*Nmesh+IY)*2*(Nmesh/2+1)+IZneigh]     *DX*TY*DZ +
                 N12[(IXneigh*Nmesh+IYneigh)*2*(Nmesh/2+1)+IZ]     *DX*DY*TZ +
                 N12[(IXneigh*Nmesh+IYneigh)*2*(Nmesh/2+1)+IZneigh]*DX*DY*DZ;

    Disp[2][i] = N13[(IX*Nmesh+IY)*2*(Nmesh/2+1)+IZ]          *TX*TY*TZ +
                 N13[(IX*Nmesh+IY)*2*(Nmesh/2+1)+IZneigh]     *TX*TY*DZ +
                 N13[(IX*Nmesh+IYneigh)*2*(Nmesh/2+1)+IZ]     *TX*DY*TZ +
                 N13[(IX*Nmesh+IYneigh)*2*(Nmesh/2+1)+IZneigh]*TX*DY*DZ +
                 N13[(IXneigh*Nmesh+IY)*2*(Nmesh/2+1)+IZ]          *DX*TY*TZ +
                 N13[(IXneigh*Nmesh+IY)*2*(Nmesh/2+1)+IZneigh]     *DX*TY*DZ +
                 N13[(IXneigh*Nmesh+IYneigh)*2*(Nmesh/2+1)+IZ]     *DX*DY*TZ +
                 N13[(IXneigh*Nmesh+IYneigh)*2*(Nmesh/2+1)+IZneigh]*DX*DY*DZ;

    sumDx += Disp[0][i];
    sumDy += Disp[1][i];
    sumDz += Disp[2][i];       
  }

  // Make sumDx, sumDy and sumDz global averages
  ierr = MPI_Allreduce(MPI_IN_PLACE,&sumDx,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  ierr = MPI_Allreduce(MPI_IN_PLACE,&sumDy,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  ierr = MPI_Allreduce(MPI_IN_PLACE,&sumDz,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);  
    
  sumDx /= (double)TotNumPart; // We will subtract these below to conserve momentum. 
  sumDy /= (double)TotNumPart;
  sumDz /= (double)TotNumPart;   
  
  return;      
}

// Wrap the particles periodically
// ===============================
#if (MEMORY_MODE || SINGLE_PRECISION)
float periodic_wrap(float x)
{
  while(x >= (float)Box) x -= (float)Box;
  while(x < 0) x += (float)Box;
  if (x == (float)Box) x = 0.0;
  return x;
}
#else
double periodic_wrap(double x)
{
  while(x >= Box) x -= Box;
  while(x < 0) x += Box;
  if (x == Box) x = 0.0;
  return x;
}
#endif

// Error message
// =============
void FatalError(char* filename, int linenum) {
  printf("Fatal Error at line %d in file %s\n", linenum, filename);
  fflush(stdout);
  free(OutputList);
  MPI_Abort(MPI_COMM_WORLD, 1);
  exit(1);
}

// This catches I/O errors occuring for fwrite(). In this case we better stop.
// ===========================================================================
size_t my_fwrite(void *ptr, size_t size, size_t nmemb, FILE * stream) {
  size_t nwritten;
  if((nwritten = fwrite(ptr, size, nmemb, stream)) != nmemb) {
    printf("\nERROR: I/O error (fwrite) on task=%d has occured.\n\n", ThisTask);
    fflush(stdout);
    FatalError((char *)"auxPM.c", 621);
  }
  return nwritten;
}
