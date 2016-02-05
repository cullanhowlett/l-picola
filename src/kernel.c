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

/* ================================================================================================================*/
/* This file contains the routines required to read in and save the generic kernel for the non-gaussian potentials.*/
/* v1.3: Removed some unneccesary print statements.                                                                */
/* ================================================================================================================*/

#include "vars.h"
#include "proto.h"

void read_kernel_table(void) {

  FILE *fd;
  char buf[500];
  double kerCoef,ker0,kerA,kerB;

  sprintf(buf, FileWithInputKernel);

  if(!(fd = fopen(buf, "r"))) {
    if (ThisTask == 0) printf("\nERROR: Can't read input kernel values in file '%s'.\n", buf);
    FatalError((char *)"kernel.c", 39);
  }

  NKernelTable = 0;
  do {
    if(fscanf(fd, " %le %le %le %le", &kerCoef, &ker0, &kerA, &kerB) == 4) {
      NKernelTable++;
    } else {
      break;
    }
  } while(1);
       
  fclose(fd);

  if(ThisTask == 0) {
    printf("Found %d lines of values in the kernel table...\n", NKernelTable);
    fflush(stdout);
  }

  KernelTable = (struct kern_table *)malloc(NKernelTable * sizeof(struct kern_table));

  sprintf(buf, FileWithInputKernel);

  if(!(fd = fopen(buf, "r"))) {
    if (ThisTask == 0) printf("\nERROR: Can't read kernel input values in file '%s'.\n", buf);
    FatalError((char *)"kernel.c", 66);
  }

  NKernelTable = 0;
  do {
    if(fscanf(fd, " %le %le %le %le", &kerCoef, &ker0, &kerA, &kerB) == 4) {
      KernelTable[NKernelTable].Coef = kerCoef;
      KernelTable[NKernelTable].ker0 = ker0;
      KernelTable[NKernelTable].kerA = kerA;
      KernelTable[NKernelTable].kerB = kerB;
      NKernelTable++;
    } else {
      break;
    }
  } while(1);

  fclose(fd);

  return;
}
