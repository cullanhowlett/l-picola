/* ==========================================================================*/
/*   Version 1.3.             Cullan Howlett & Marc Manera,                  */
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

/* ============================================================================*/
/* This file contains all extra routines associated with lightcone simulations */
/* v1.3: Particle positions and velocities are now output in user-defined      */
/*       units as opposed to always being in Mpc/h. Hardcoded 'boundary' value */
/*       for checking necessary replicates is now also in user-defined units   */
/* ============================================================================*/

#include "vars.h"
#include "proto.h"

// Set up the lightcone parameters
// ===============================
void set_lightcone(void) {

  Light = LIGHT / UnitLength_in_cm * UnitTime_in_s;

  // Calculate the actual maximum number of replicate boxes we need in all six directions. This allows
  // us to remove replicates that are unnecessary
  double Rcomov_max = Light/Hubble*KickStd(1.0/(1.0+OutputList[0].Redshift),1.0);
  if (ThisTask == 0) {
    if (((int)(ceil(fabs((Origin_x - Rcomov_max)/Box))) > Nrep_neg_x) &&
        ((int)(ceil(fabs((Origin_y - Rcomov_max)/Box))) > Nrep_neg_y)  &&
        ((int)(ceil(fabs((Origin_z - Rcomov_max)/Box))) > Nrep_neg_z)  &&
        ((int)(ceil((Rcomov_max + Origin_x)/Box))-1 > Nrep_pos_x) &&
        ((int)(ceil((Rcomov_max + Origin_y)/Box))-1 > Nrep_pos_y) &&
        ((int)(ceil((Rcomov_max + Origin_z)/Box))-1 > Nrep_pos_z)) {
       printf("\nWARNING: Initial comoving radius of the lightcone is outside all replicated boxes. Is this intentional?\n");
       printf("         If so, reducing the first output redshift will start the lightcone later and prove more efficient\n");
       printf("         If not then please increase the number of replicates in the input parameters file.\n\n");
    }
  } 
  if ((int)(ceil(fabs((Origin_x - Rcomov_max)/Box))) < Nrep_neg_x) Nrep_neg_x = (int)(ceil(fabs((Origin_x - Rcomov_max)/Box)));
  if ((int)(ceil(fabs((Origin_y - Rcomov_max)/Box))) < Nrep_neg_y) Nrep_neg_y = (int)(ceil(fabs((Origin_y - Rcomov_max)/Box)));
  if ((int)(ceil(fabs((Origin_z - Rcomov_max)/Box))) < Nrep_neg_z) Nrep_neg_z = (int)(ceil(fabs((Origin_z - Rcomov_max)/Box)));
  if ((int)(ceil((Rcomov_max + Origin_x)/Box))-1 < Nrep_pos_x) Nrep_pos_x = (int)(ceil((Rcomov_max + Origin_x)/Box))-1;
  if ((int)(ceil((Rcomov_max + Origin_y)/Box))-1 < Nrep_pos_y) Nrep_pos_y = (int)(ceil((Rcomov_max + Origin_y)/Box))-1;
  if ((int)(ceil((Rcomov_max + Origin_z)/Box))-1 < Nrep_pos_z) Nrep_pos_z = (int)(ceil((Rcomov_max + Origin_z)/Box))-1;

  // Create an array of flags that for each replicate tell us whether or not we need to check the particles in it.
  // A flag of 0 means we check it, 1 means the box is completely inside the lightcone so no need to check it this step 
  // 2 means the box has completely left the lightcone so no need to ever check it again.
  repflag = (int *)calloc((Nrep_neg_x+Nrep_pos_x+1)*(Nrep_neg_y+Nrep_pos_y+1)*(Nrep_neg_z+Nrep_pos_z+1),sizeof(int));
  writeflag = (int *)calloc((Nrep_neg_x+Nrep_pos_x+1)*(Nrep_neg_y+Nrep_pos_y+1)*(Nrep_neg_z+Nrep_pos_z+1),sizeof(int));
  Noutput = (unsigned int *)calloc((Nrep_neg_x+Nrep_pos_x+1)*(Nrep_neg_y+Nrep_pos_y+1)*(Nrep_neg_z+Nrep_pos_z+1),sizeof(unsigned int));
  
  Nrep_neg_max[0] = Nrep_neg_x; Nrep_pos_max[0] = Nrep_pos_x;
  Nrep_neg_max[1] = Nrep_neg_y; Nrep_pos_max[1] = Nrep_pos_y;
  Nrep_neg_max[2] = Nrep_neg_z; Nrep_pos_max[2] = Nrep_pos_z;
  

  return;
}

// Flag replicates that don't need to be looped over this iteration, as they are either completely inside or outside the lightcone
// ===============================================================================================================================
void flag_replicates(double Rcomov_old, double Rcomov_new, double boundary) {

  int i, j, k, ii, jj, kk; 
  int repcount_low, coord;
  double Xvert, Yvert, Zvert, Rvert;
  double dist_face_old, dist_face_new;
  double Rcomov_old2 = Rcomov_old*Rcomov_old;
  double Rcomov_new2 = Rcomov_new*Rcomov_new;

  // Check the maximum number of replicates in each direction and update if necessary
  if ((int)(ceil(fabs((Origin_x - Rcomov_old)/Box))) < Nrep_neg_x) Nrep_neg_x = (int)(ceil(fabs((Origin_x - Rcomov_old)/Box)));
  if ((int)(ceil(fabs((Origin_y - Rcomov_old)/Box))) < Nrep_neg_y) Nrep_neg_y = (int)(ceil(fabs((Origin_y - Rcomov_old)/Box)));
  if ((int)(ceil(fabs((Origin_z - Rcomov_old)/Box))) < Nrep_neg_z) Nrep_neg_z = (int)(ceil(fabs((Origin_z - Rcomov_old)/Box)));
  if ((int)(ceil((Rcomov_old + Origin_x)/Box))-1 < Nrep_pos_x) Nrep_pos_x = (int)(ceil((Rcomov_old + Origin_x)/Box))-1;
  if ((int)(ceil((Rcomov_old + Origin_y)/Box))-1 < Nrep_pos_y) Nrep_pos_y = (int)(ceil((Rcomov_old + Origin_y)/Box))-1;
  if ((int)(ceil((Rcomov_old + Origin_z)/Box))-1 < Nrep_pos_z) Nrep_pos_z = (int)(ceil((Rcomov_old + Origin_z)/Box))-1;

  // Loop over all replicates and if any replicate is completely within the lightcone (i.e. all eight vertices are within Rcomov_new)
  // then flag it so that we don't have to loop over it. NOTE: this method will not work to see if replicates are completely outside the
  // lightcone, as there can easily be a case where the lightcone is inside all eight vertices of a replicated box, yet the lightcone still 
  // passes through it. In this case the routine is much more complicated and involves looping over all 6 faces of the replicate, and computing
  // the shortest distance to the finite plane of the face.
  for (i = -Nrep_neg_x; i<=Nrep_pos_x; i++) {
    for (j = -Nrep_neg_y; j<=Nrep_pos_y; j++) {
      for (k = -Nrep_neg_z; k<=Nrep_pos_z; k++) {

        coord = ((i+Nrep_neg_max[0])*(Nrep_neg_max[1]+Nrep_pos_max[1]+1)+(j+Nrep_neg_max[1]))*(Nrep_neg_max[2]+Nrep_pos_max[2]+1)+(k+Nrep_neg_max[2]);

        // Skip this replicate if we already know it is completely outside the lightcone
        // the particle positions first
        if (repflag[coord] == 2) continue;
 
        // Loop over all the vertices
        repflag[coord] = 0;
        repcount_low = 0;
        for (ii = 0; ii < 2; ii++) {
          for (jj = 0; jj < 2; jj++) {
            for (kk = 0; kk < 2; kk++) {

              // Include a buffer region (boundary) to account for the fact that the particle might move beyond the box boundaries, 20Mpc should be enough.
              Xvert = (i+((ii*Local_np+Local_p_start)/(double)Nsample))*Box - Origin_x + pow(-1, ii+1)*boundary;
              Yvert = (j+jj)*Box - Origin_y + pow(-1, jj+1)*boundary;
              Zvert = (k+kk)*Box - Origin_z + pow(-1, kk+1)*boundary;
              Rvert = Xvert*Xvert+Yvert*Yvert+Zvert*Zvert;
                  
              if (Rvert < Rcomov_new2) repcount_low++;
            }
          }
        }
            
        // If box is completely inside the lightcone, we flag it and continue into the next replicate
        if (repcount_low == 8) {
          repflag[coord] = 1;
          continue;
        }

        // Otherwise we check to see if the replicate is completely outside the lightcone. This is where it gets more complicated. 
        // For all the necessary faces of the replicate, we calculate the shortest possible distance from the origin to the face.
        // This is given by the sum of the squares of the distance to the projected origin and the 
        // distance from the projected origin to the nearest line segment, which itself is the sum of the squares 
        // of the distance from the projection of the origin onto the plane to the projection of this onto the 
        // nearest line segment and the distance along the line segment.
        if (repcount_low == 0) {

          dist_face_old = 1.0e30;
          dist_face_new = 1.0e30;

          double LeftX = i+(Local_p_start/(double)Nsample);
          double RightX = i+((Local_np+Local_p_start)/(double)Nsample);

          // Check the two faces perpendicular to the x-axis (don't forget these are different depending on the task). Don't forget to include the boundary region
          if (Origin_x < LeftX*Box - boundary) {
            dist_face_old = Origin_x - LeftX*Box + boundary;
            dist_face_old *= dist_face_old;
            dist_face_old += nearest_dist(Origin_y, Origin_z, j, k, j+1, k+1, boundary);
          } else if (Origin_x > RightX*Box + boundary) {
            dist_face_old = Origin_x - RightX*Box - boundary;
            dist_face_old *= dist_face_old;
            dist_face_old += nearest_dist(Origin_y, Origin_z, j, k, j+1, k+1, boundary);
          }

          // Check the two faces perpendicular to the y-axis.
          if (Origin_y < j*Box - boundary) {
            dist_face_new = Origin_y - j*Box + boundary;
            dist_face_new *= dist_face_new;
            dist_face_new += nearest_dist(Origin_x, Origin_z, LeftX, k, RightX, k+1, boundary);
          } else if (Origin_y > (j+1)*Box + boundary) {
            dist_face_new = Origin_y - (j+1)*Box - boundary;
            dist_face_new *= dist_face_new;
            dist_face_new += nearest_dist(Origin_x, Origin_z, LeftX, k, RightX, k+1, boundary);
          }
          if (dist_face_old > dist_face_new) dist_face_old = dist_face_new;

          // Check the two faces perpendicular to the z-axis.
          if (Origin_z < k*Box - boundary) {
            dist_face_new = Origin_z - k*Box + boundary;
            dist_face_new *= dist_face_new;
            dist_face_new += nearest_dist(Origin_x, Origin_y, LeftX, j, RightX, j+1, boundary);
          } else if (Origin_z > (k+1)*Box + boundary) {
            dist_face_new = Origin_z - (k+1)*Box - boundary;
            dist_face_new *= dist_face_new;
            dist_face_new += nearest_dist(Origin_x, Origin_y, LeftX, j, RightX, j+1, boundary);
          }
          if (dist_face_old > dist_face_new) dist_face_old = dist_face_new;
          
          // If dist_face_old is greater than 9.9e29 then the origin is WITHIN the current replicate so we definitely have to loop over it.
          // Otherwise dist_face_old now contains the shortest distance from ANY point on the replicate to the lightcone. Hence if this is greater than Rcomov_old2
          // we never have to loop over this replicate again
          if (dist_face_old < 9.9e29) {
            if (dist_face_old > Rcomov_old2) repflag[coord] = 2;
          }        
        }
      }
    }
  }

  return;
}


// For a given face on a replicate this routine returns the shortest distance between the face, bounded by (ix, iy) and (jx, jy), and the point (px, py)
double nearest_dist(double px, double py, double ix, double iy, double jx, double jy, double boundary) {

  double dist1, dist2;
  double dist_line_old, dist_line_new;

  dist_line_old = 1.0e30;
  dist_line_new = 1.0e30;

  // Check the 4 line segments of the face of the replicate that is contained on the plane. For each necessary line segment we compute 
  // the distance, 'dist1', of the projection of the projected origin on the plane onto an infinite line containing the line segment (2D -> 1D).
  // The distance along the line between the 1D projection and the end of the line segment is then dist2. NOTE: for projections that land
  // inside the face, .i.e. within all four line segments, we don't care about dist1 or dist2 as the shortest distance between the point and the
  // replicate is just the distance to the plane (computed previously). Also, we only have to do a maximum of 2 line segments then as this is the most 
  // that any point can see.
  if (px < ix*Box - boundary) {
    dist1 = px - ix*Box + boundary;
    dist2 = 0.0;
    if (py < iy*Box - boundary) {
      dist2 = py - iy*Box + boundary; 
    } else if (py > jy*Box + boundary) {
      dist2 = py - jy*Box - boundary;
    }          
    dist_line_old = dist1*dist1+dist2*dist2;
  } else if (px > jx*Box + boundary) {
    dist1 = px - jx*Box - boundary;
    dist2 = 0.0;
    if (py < iy*Box - boundary) {
      dist2 = py - iy*Box + boundary; 
    } else if (py > jy*Box + boundary) {
      dist2 = py - jy*Box - boundary;
    }          
    dist_line_old = dist1*dist1+dist2*dist2;
  }
 
  // The third and fourth line segments
  if (py < iy*Box - boundary) {
    dist1 = py - iy*Box + boundary;
    dist2 = 0.0;
    if (px < ix*Box - boundary) {
      dist2 = px - ix*Box + boundary; 
    } else if (px > jx*Box + boundary) {
      dist2 = px - jx*Box - boundary;
    }        
    dist_line_new = dist1*dist1+dist2*dist2;
  } else if (py > jy*Box + boundary) {
    dist1 = py - jy*Box - boundary;
    dist2 = 0.0;
    if (px < ix*Box - boundary) {
      dist2 = px - ix*Box + boundary; 
    } else if (px > jx*Box + boundary) {
      dist2 = px - jx*Box - boundary;
    }        
    dist_line_new = dist1*dist1+dist2*dist2;
  }

  // Find the shortest distance to any of the line segments
  if (dist_line_old > dist_line_new) dist_line_old = dist_line_new;
  
  if (dist_line_old > 9.9e29) {
    return 0.0;
  } else {
    return dist_line_old;
  }

}

// Drift and output the particles for lightcone simulations
// ========================================================
void Drift_Lightcone(double A, double AFF, double AF, double Di, double Di2) {

  // We'll flag to see if the particle has left the lightcone here and output if it has
  // We don't bother interpolating the velocity, only the position, as the implicit assumption
  // with KDK anyway is that the velocity at the halfway point in the timestep is constant
  // between the initial and final particle positions for that timestep. If we did want to interpolate
  // the velocity, we could move this whole section to the location marked above, then reverse the 
  // updating and periodic wrapping of the particle positions to allow interpolation.
  // To avoid having to store all the replicate of the particles, we also output the particles as soon 
  // as they are flagged, allowing us to loop over all the replicated boxes instead. Finally, because all
  // the lightcone stuff is done here, we don't need memory to store the particle flags.

  size_t bytes;
  int NTAB = 1000;                                       // The length of the particle exit time lookup tables (we spline anyway so not really important)
  int i, j, k, coord, flag, repcount;
  unsigned int n, * pc, blockmaxlen, blockmaxlenglob;
  unsigned int outputflag, NumPartMax;
  float * block;
  double dyyy, da1, da2, dv1, dv2;
  double dyyy_tmp, da1_tmp, da2_tmp, AL;
  double Delta_Pos[3];
  double boundary = 20.0*(3.085678e24/UnitLength_in_cm); // Constant 20 Mpc/h (should be large enough but might need more if particles move a large distance in a timestep)
  double fac = Hubble/AF;                                // This differs from snapshot 'fac' by sqrt(AF) as we don't know after runtime what AF is. 
  double lengthfac = 1.0;                                // Keep positions in user-specified units (Originally converted positions to Mpc/h)
  double velfac    = 1.0;                                // Keep velocities in user-specified units (Originally converted velocities to km/s)
  double Rcomov_old  = Light/Hubble*KickStd(A,1.0);
  double Rcomov_new  = Light/Hubble*KickStd(AFF,1.0);
  double Rcomov_old2 = Rcomov_old*Rcomov_old; 
  double Rcomov_new2 = Rcomov_new*Rcomov_new;
  double Xpart, Ypart, Zpart, Rpart_old, Rpart_new, Rpart_old2, Rpart_new2;   
  double * AL_tab, * da1_tab, * da2_tab, * dyyy_tab;
  gsl_spline * da1_spline, * da2_spline, * dyyy_spline;                              // Spline fits to the exit-time lookup tables
  gsl_interp_accel * da1_acc, * da2_acc, * dyyy_acc;
#ifdef TIMING
  double startcpu, endcpu;
  double startwall, endwall;
#endif 

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

  // Add back LPT velocities if we had subtracted them. 
  // This corresponds to L_+ operator in Tassev et. al, 2013.
  dv1 = QdD1da(AF);  // Q*dD_{1}/da
  dv2 = QdD2da(AF);  // Q*dD_{2}/da

  // Flag the replicates that don't need looping over
  flag_replicates(Rcomov_old, Rcomov_new, boundary);

  // Create lookup tables for interpolating the time when the particle left the lightcone 
  // and the corresponding growth factors. Using lookup tables (for the growth factors especially) makes a HUGE difference. 
  // (In my small tests it reduced the time by a factor of TEN!!).
  AL_tab = (double *)malloc(NTAB*sizeof(double));
  da1_tab = (double *)malloc(NTAB*sizeof(double));
  da2_tab = (double *)malloc(NTAB*sizeof(double));
  dyyy_tab = (double *)malloc(NTAB*sizeof(double));
  for (i=0; i<NTAB; i++) {
    AL = (i*(AFF-A))/(NTAB-1.0) + A;
    AL_tab[i] = AL;
    da1_tab[i] = growthD(AL)- Di;
    da2_tab[i] = growthD2(AL) - Di2;
    if (DeltaA == 0) {
      dyyy_tab[i]=DriftCOLA(A,AL,AF);
    } else if (DeltaA == 1) {
      dyyy_tab[i]=DriftStd(A,AL);
    } else if (DeltaA == 2) {
      dyyy_tab[i]=DriftStd(A,AL);
    } else {
      dyyy_tab[i]=(AL-A)/Qfactor(AF);
    } 
  }
  da1_acc = gsl_interp_accel_alloc();
  da2_acc = gsl_interp_accel_alloc();
  dyyy_acc = gsl_interp_accel_alloc();
  da1_spline = gsl_spline_alloc(gsl_interp_cspline, NTAB);
  da2_spline = gsl_spline_alloc(gsl_interp_cspline, NTAB);
  dyyy_spline = gsl_spline_alloc(gsl_interp_cspline, NTAB);
  gsl_spline_init(da1_spline, AL_tab, da1_tab, NTAB);
  gsl_spline_init(da2_spline, AL_tab, da2_tab, NTAB);
  gsl_spline_init(dyyy_spline, AL_tab, dyyy_tab, NTAB);

  free(AL_tab);
  free(da1_tab);
  free(da2_tab);
  free(dyyy_tab);

  // For the outputting we can only moderate how many tasks output at once if they all enter the output stage at the same time.
  // As such we output at set times rather than when actually necessary. This also alleviates the need for
  // any inter-task communications, which would be necessary if we were to only output once a task has filled
  // it's quota. Unfortunately, as the number of particles on each task will be different we have to make 
  // each task loop over the maximum number of particles on any task, but most tasks will not have to actually do anything.
  // This will result in more outputs than is truly necessary. We also must assume that all looped over replicates of a 
  // given particle are to be outputted even if they are not.

  // Calculate the global maximum number of particles on a task and the global minimum number of particles we can store.
  // Then allocate memory to store the particles that we are outputting. We may have some spare memory from 
  // deallocating the force grids on top of that from deallocating the displacement arrays. If repcount == 0 then
  // no replicates are inside the lightcone. This can happen if we have an initial lightcone redshift beyond
  // the maximum extent of the replicates
  ierr = MPI_Allreduce(&NumPart, &NumPartMax, 1, MPI_UNSIGNED, MPI_MAX, MPI_COMM_WORLD);
#ifdef MEMORY_MODE
  block = (float *)malloc(bytes = 6*Total_size*sizeof(float_kind)+3*NumPart*sizeof(float));
#else
  block = (float *)malloc(bytes = 3*NumPart*sizeof(float_kind));
#endif

  // Create an array containing the replicates we actually need to loop over. The order of this array will also denote
  // The order in which we store the output and the number of particles needed to output. To output the replicates to different files
  // we can just make sure we store particles from the same replicate consecutively as we already assume that we have to output all replicated
  // particles (and hence have enough space).
  // How many particles can we store assuming all tasks replicate the particles by the maximum amount?
  repcount = 0;
  pc = (unsigned int *)calloc((Nrep_neg_x+Nrep_pos_x+1)*(Nrep_neg_y+Nrep_pos_y+1)*(Nrep_neg_z+Nrep_pos_z+1),sizeof(unsigned int));
  for (i = -Nrep_neg_x; i<=Nrep_pos_x; i++) {
    for (j = -Nrep_neg_y; j<=Nrep_pos_y; j++) {
      for (k = -Nrep_neg_z; k<=Nrep_pos_z; k++) {
        coord = ((i+Nrep_neg_max[0])*(Nrep_neg_max[1]+Nrep_pos_max[1]+1)+(j+Nrep_neg_max[1]))*(Nrep_neg_max[2]+Nrep_pos_max[2]+1)+(k+Nrep_neg_max[2]);
        if (repflag[coord] == 0) repcount++;
      }
    }
  }
  if (repcount == 0) {
    blockmaxlen = NumPartMax;
  } else {
    blockmaxlen = (unsigned int)(bytes / (6 * sizeof(float) * repcount));
  }
  ierr = MPI_Allreduce(&blockmaxlen, &blockmaxlenglob, 1, MPI_UNSIGNED, MPI_MIN, MPI_COMM_WORLD);

  // Loop over all particles, modifying the position based on the current replicate
  outputflag = 0;
  for(n=0; n<NumPartMax; n++) {

    outputflag++;
    if (n < NumPart) {

      Delta_Pos[0] = (P[n].Vel[0]-sumx)*dyyy+UseCOLA*(P[n].Dz[0]*da1+P[n].D2[0]*da2);
      Delta_Pos[1] = (P[n].Vel[1]-sumy)*dyyy+UseCOLA*(P[n].Dz[1]*da1+P[n].D2[1]*da2);   
      Delta_Pos[2] = (P[n].Vel[2]-sumz)*dyyy+UseCOLA*(P[n].Dz[2]*da1+P[n].D2[2]*da2);     

      // Check that 100Mpc^2/h^2 boundaries is enough
      if((Delta_Pos[0] > boundary) || (Delta_Pos[1] > boundary) || (Delta_Pos[2] > boundary)) {
        printf("\nERROR: Particle displacement greater than boundary for lightcone replicate estimate.\n");
        printf("       increase boundary condition in lightcone.c (line 56)\n\n");
        FatalError((char *)"lightcone.c", 390);
      }

      // Loop over all replicates
      repcount=0;
      for (i = -Nrep_neg_x; i<=Nrep_pos_x; i++) {
        for (j = -Nrep_neg_y; j<=Nrep_pos_y; j++) {
          for (k = -Nrep_neg_z; k<=Nrep_pos_z; k++) {

            coord = ((i+Nrep_neg_max[0])*(Nrep_neg_max[1]+Nrep_pos_max[1]+1)+(j+Nrep_neg_max[1]))*(Nrep_neg_max[2]+Nrep_pos_max[2]+1)+(k+Nrep_neg_max[2]);
            if (repflag[coord] == 0) {

              // Did the particle start the timestep inside the lightcone?
              flag = 0;
              Xpart = P[n].Pos[0] - Origin_x + (i*Box);
              Ypart = P[n].Pos[1] - Origin_y + (j*Box);
              Zpart = P[n].Pos[2] - Origin_z + (k*Box);
              Rpart_old2 = Xpart*Xpart+Ypart*Ypart+Zpart*Zpart;
 
              if (Rpart_old2 <= Rcomov_old2) flag = 1;

              // Have any particles that started inside the lightcone now exited?
              if (flag) {
                Xpart += Delta_Pos[0];
                Ypart += Delta_Pos[1];
                Zpart += Delta_Pos[2];
                Rpart_new2 = Xpart*Xpart+Ypart*Ypart+Zpart*Zpart;

                if (Rpart_new2 > Rcomov_new2) {
  
                  // Interpolate the particle position. We do this by first calculating the exact time at which
                  // the particle exited the lightcone, then updating the position to there.
                  Rpart_old = sqrt(Rpart_old2);
                  Rpart_new = sqrt(Rpart_new2);
                  AL = A + (AFF-A)*((Rcomov_old-Rpart_old)/((Rpart_new-Rpart_old)-(Rcomov_new-Rcomov_old)));
                  da1_tmp = gsl_spline_eval(da1_spline, AL, da1_acc);
                  da2_tmp = gsl_spline_eval(da2_spline, AL, da2_acc);
                  dyyy_tmp = gsl_spline_eval(dyyy_spline, AL, dyyy_acc);
                                                  
                  // Store the interpolated particle position and velocity.
                  unsigned int ind = 6*(blockmaxlen*repcount+pc[repcount]);
                  block[ind]     = (float)(lengthfac*(P[n].Pos[0] + (P[n].Vel[0]-sumx)*dyyy_tmp+UseCOLA*(P[n].Dz[0]*da1_tmp+P[n].D2[0]*da2_tmp) + (i*Box)));
                  block[ind + 1] = (float)(lengthfac*(P[n].Pos[1] + (P[n].Vel[1]-sumy)*dyyy_tmp+UseCOLA*(P[n].Dz[1]*da1_tmp+P[n].D2[1]*da2_tmp) + (j*Box)));
                  block[ind + 2] = (float)(lengthfac*(P[n].Pos[2] + (P[n].Vel[2]-sumz)*dyyy_tmp+UseCOLA*(P[n].Dz[2]*da1_tmp+P[n].D2[2]*da2_tmp) + (k*Box)));
                  block[ind + 3] = (float)(velfac*fac*(P[n].Vel[0]-sumx+(P[n].Dz[0]*dv1+P[n].D2[0]*dv2)*UseCOLA));
                  block[ind + 4] = (float)(velfac*fac*(P[n].Vel[1]-sumy+(P[n].Dz[1]*dv1+P[n].D2[1]*dv2)*UseCOLA));
                  block[ind + 5] = (float)(velfac*fac*(P[n].Vel[2]-sumz+(P[n].Dz[2]*dv1+P[n].D2[2]*dv2)*UseCOLA));
                  pc[repcount]++;   
                  Noutput[coord]++;
                }
              }
              repcount++;
            }
          }
        }
      }
 
      // Update the particle's position
      P[n].Pos[0] = periodic_wrap(P[n].Pos[0]+Delta_Pos[0]);
      P[n].Pos[1] = periodic_wrap(P[n].Pos[1]+Delta_Pos[1]);
      P[n].Pos[2] = periodic_wrap(P[n].Pos[2]+Delta_Pos[2]); 
    }
    
    if (outputflag == blockmaxlenglob) {
#ifdef TIMING
      startcpu = (double)clock();
      startwall = MPI_Wtime();
#endif
      //Output_Lightcone(pc, blockmaxlen, block_id, block_al, block);
      Output_Lightcone(pc, blockmaxlen, block);
#ifdef TIMING
      endcpu = (double)clock();
      endwall = MPI_Wtime();
      CpuTime_Output[timeSteptot-1] += (endcpu-startcpu)/(double)CLOCKS_PER_SEC;
      WallTime_Output[timeSteptot-1] += endwall-startwall;
#endif
      outputflag = 0;
      for (i=0; i<(Nrep_neg_x+Nrep_pos_x+1)*(Nrep_neg_y+Nrep_pos_y+1)*(Nrep_neg_z+Nrep_pos_z+1); i++) pc[i] = 0;
    }
  }

  if (outputflag > 0) {
#ifdef TIMING
    startcpu = (double)clock();
    startwall = MPI_Wtime();
#endif
    Output_Lightcone(pc, blockmaxlen, block);
#ifdef TIMING
    endcpu = (double)clock();
    endwall = MPI_Wtime();
    CpuTime_Output[timeSteptot-1] += (endcpu-startcpu)/(double)CLOCKS_PER_SEC;
    WallTime_Output[timeSteptot-1] += endwall-startwall;
#endif
  }
  free(pc);
  free(block);

  gsl_spline_free(da1_spline);
  gsl_spline_free(da2_spline);
  gsl_spline_free(dyyy_spline);
  gsl_interp_accel_free(da1_acc);
  gsl_interp_accel_free(da2_acc);
  gsl_interp_accel_free(dyyy_acc);

  return;
}

// Output the lightcone data
// =========================
void Output_Lightcone(unsigned int * pc, unsigned int blockmaxlen, float * block) {

  FILE * fp; 
  char buf[300];
  int i, j, k;
  int nprocgroup, groupTask, masterTask, repcount;
  unsigned int chunk, coord;
#ifdef UNFORMATTED
  int dummy1, dummy2;
#else
  unsigned int n;
#endif

  nprocgroup = NTask / NumFilesWrittenInParallel;
  if (NTask % NumFilesWrittenInParallel) nprocgroup++;
  masterTask = (ThisTask / nprocgroup) * nprocgroup;
  for(groupTask = 0; groupTask < nprocgroup; groupTask++) {
    if (ThisTask == (masterTask + groupTask)) {

      // Loop over all replicates
      repcount=0;
      for (i = -Nrep_neg_x; i<=Nrep_pos_x; i++) {
        for (j = -Nrep_neg_y; j<=Nrep_pos_y; j++) {
          for (k = -Nrep_neg_z; k<=Nrep_pos_z; k++) {

            coord = ((i+Nrep_neg_max[0])*(Nrep_neg_max[1]+Nrep_pos_max[1]+1)+(j+Nrep_neg_max[1]))*(Nrep_neg_max[2]+Nrep_pos_max[2]+1)+(k+Nrep_neg_max[2]);
            if (repflag[coord] == 0) {

              if(pc[repcount] > 0) {
                sprintf(buf, "%s/%s_lightcone.%d", OutputDir, FileBase, coord*NTask+ThisTask);
                if (writeflag[coord] == 0) {
                  // Overwrite any pre-existing output files otherwise we'll append onto the end of them.
                  if(!(fp = fopen(buf, "w"))) {
                    printf("\nERROR: Can't write in file '%s'.\n\n", buf);
                    FatalError((char *)"lightcone.c", 533);
                  }
                  fflush(stdout);
                  writeflag[coord] = 1;
                } else {
                  if(!(fp = fopen(buf, "a"))) {
                    printf("\nERROR: Can't write in file '%s'.\n\n", buf);
                    FatalError((char *)"lightcone.c", 540);
                  }
                  fflush(stdout);
                }

#ifdef UNFORMATTED
                // write coordinates and velocities in unformatted binary
                chunk = 6*blockmaxlen*repcount;
                dummy1 = sizeof(pc[repcount]); 
                dummy2 = sizeof(float) * 6 * pc[repcount];
                my_fwrite(&dummy1, sizeof(dummy1), 1, fp);
                my_fwrite(&(pc[repcount]), sizeof(unsigned int), 1, fp);
                my_fwrite(&dummy1, sizeof(dummy1), 1, fp);
                my_fwrite(&dummy2, sizeof(dummy2), 1, fp);
                my_fwrite(&(block[chunk]), sizeof(float), 6 * pc[repcount], fp);
                my_fwrite(&dummy2, sizeof(dummy2), 1, fp);
#else
                // write coordinates and velocities in ASCII
                for(n=0; n<pc[repcount]; n++) {
                  chunk = 6*(blockmaxlen*repcount+n);
                  fprintf(fp,"%12.6f %12.6f %12.6f %12.6f %12.6f %12.6f\n",block[chunk],block[chunk+1],block[chunk+2],block[chunk+3],block[chunk+4],block[chunk+5]);
                }
#endif
                fclose(fp);
              }
              repcount++;
            }
          }
        }
      } 
    }
    ierr = MPI_Barrier(MPI_COMM_WORLD);
  }

  return;
}

// Generate the info file which contains a list of all the output files, the 8 corners of the slices on those files and the number of particles in the slice
// =========================================================================================================================================================
void Output_Info_Lightcone(void) {

  FILE * fp=NULL; 
  char buf[300];
  int i, j, k, m;

  int * Local_p_start_table = (int *)malloc(sizeof(int) * NTask);
  MPI_Allgather(&Local_p_start, 1, MPI_INT, Local_p_start_table, 1, MPI_INT, MPI_COMM_WORLD);
 
  if (ThisTask == 0) {
    sprintf(buf, "%s/%s_lightcone.info", OutputDir, FileBase);
    if(!(fp = fopen(buf, "w"))) {
      printf("\nERROR: Can't write in file '%s'.\n\n", buf);
      FatalError((char *)"lightcone.c", 592);
    }
    fflush(stdout);
    fprintf(fp, "#    FILENUM      XMIN         YMIN        ZMIN         XMAX         YMAX         ZMAX         NPART    \n");
  }

  for (i=-Nrep_neg_max[0]; i<=Nrep_pos_max[0]; i++) {
    for (j=-Nrep_neg_max[1]; j<=Nrep_pos_max[1]; j++) {
      for (k=-Nrep_neg_max[2]; k<=Nrep_pos_max[2]; k++) {

        int coord = ((i+Nrep_neg_max[0])*(Nrep_neg_max[1]+Nrep_pos_max[1]+1)+(j+Nrep_neg_max[1]))*(Nrep_neg_max[2]+Nrep_pos_max[2]+1)+(k+Nrep_neg_max[2]);

        unsigned int * Noutput_table = (unsigned int *)malloc(sizeof(unsigned int) * NTask);
        MPI_Allgather(&(Noutput[coord]), 1, MPI_UNSIGNED, Noutput_table, 1, MPI_UNSIGNED, MPI_COMM_WORLD);

        if (ThisTask == 0) {
          double y0 = j*Box;
          double z0 = k*Box;
          double y1 = (j+1)*Box;
          double z1 = (k+1)*Box;
          for (m=0; m<NTask; m++) {
            double x0 = (i+(Local_p_start_table[m]/(double)Nsample))*Box; 
            double x1 = (i+((Local_np_table[m]+Local_p_start_table[m])/(double)Nsample))*Box; 
            fprintf(fp, "%12d %12.6lf %12.6lf %12.6lf %12.6lf %12.6lf %12.6lf %12u\n", coord*NTask+m, x0, y0, z0, x1, y1, z1, Noutput_table[m]);
          }
        }
        free(Noutput_table);

      }
    }
  }

  if (ThisTask == 0) fclose(fp);

  free(Local_p_start_table);

  return;
}
