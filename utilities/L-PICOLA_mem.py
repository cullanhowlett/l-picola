import math
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import matplotlib.ticker as tk

# ============================================================================================= #
# This is a simple Python code to calculate the memory requirements of an L-PICOLA run, and     #
# the main portions of the code, based on the compilation options and run parameters.           #
# The output is a plot of the memory per processor as a function of the number of processors,   # 
# and the minimum number of processors required given some amount of memory available per       #
# processor. There is no support for calculating lightcone simulations, however the additional  #
# memory required during the lightcone is allocated based on how much 'free' memory there is,   #
# and hence lightcone simulations use negligible extra memory compared to that calculated here. #
# ============================================================================================= #

# Plotting parameters
# ===================
mem_limit = 4.0       # The maximum amount of memory available per processor
minproc = 1         # The smallest number of processors to consider
maxproc = 256       # The largest number of processors to consider

# Compilation options (0 = off, 1 = on)
# =====================================
SINGLE_PRECISION = 0
MEMORY_MODE = 1
PARTICLE_ID = 0

# The required run parameters
# ===========================
Nsample = 1280.0
Nmesh = 1280.0
Buffer = 1.3

# Perform the calculation
# =======================
# Memory per grid cell
if (SINGLE_PRECISION):
    grid = 4.0
else:
    grid = 8.0

# Memory per piece of particle data
if ((MEMORY_MODE) or (SINGLE_PRECISION)):
    part = 4.0
else:
    part = 8.0

# Memory per particle
if (PARTICLE_ID):
    partfull = 12.0*part+8.0
else:
    partfull = 12.0*part 

nproc = np.linspace(minproc,maxproc,maxproc-minproc+1).astype(int)
MEM = np.zeros(len(nproc))
MEM_2LPT = np.zeros(len(nproc))
MEM_DISP = np.zeros(len(nproc))
MEM_INIT = np.zeros(len(nproc))
MEM_MOVE = np.zeros(len(nproc))
MEM_DENS = np.zeros(len(nproc))
MEM_NBODY = np.zeros(len(nproc))
for i in range(len(nproc)):

    # 2LPT field calculation
    MEM_2LPT[i]  = (9.0*grid*Nmesh*Nmesh*(Nmesh+2.0))/nproc[i]     # Mesh split over processors
    MEM_2LPT[i] += 9.0*grid*Nmesh*(Nmesh+2.0)                      # Each processor gets an extra slice
    MEM_2LPT[i] += 4.0*Nmesh*Nmesh                                 # Seed table on each processor

    # 2LPT displacement calculation
    MEM_DISP[i]  = (6.0*grid*Nmesh*Nmesh*(Nmesh+2.0))/nproc[i]     # We can deallocate three 2LPT meshes before here
    MEM_DISP[i] += 6.0*grid*Nmesh*(Nmesh+2.0)                      # Each processor gets an extra slice
    MEM_DISP[i] += (6.0*part*Nsample*Nsample*Nsample)/nproc[i]     # ZA and 2LPT displacements split over processors

    # Initializing the particles
    MEM_INIT[i] += (6.0*part*Nsample*Nsample*Nsample)/nproc[i]         # ZA and 2LPT displacements split over processors
    MEM_INIT[i] += (partfull*Buffer*Nsample*Nsample*Nsample)/nproc[i]  # Particles split over processors

    # PM algorithm. In memory mode we deallocate certain parts as we go along
    # Moving the particles
    MEM_MOVE[i]  = (partfull*Buffer*Nsample*Nsample*Nsample)/nproc[i]            # Particles split over processors
    MEM_MOVE[i] += (partfull*2.0*(Buffer-1.0)*Nsample*Nsample*Nsample)/nproc[i]  # Extra memory for moving particles

    # Calculating the density and forces
    MEM_DENS[i]  = (partfull*Buffer*Nsample*Nsample*Nsample)/nproc[i]   # Particles split over processors
    MEM_DENS[i] += (4.0*grid*Nmesh*Nmesh*(Nmesh+2.0))/nproc[i]          # Force and density meshes split over processors
    MEM_DENS[i] += 4.0*grid*Nmesh*(Nmesh+2.0)                           # Each processor gets an extra slice
    MEM_DENS[i] += grid*Nmesh*(Nmesh+2.0)                               # Extra slice for transferring density

    # Calculating the N-Body displacements
    MEM_NBODY[i]  = (partfull*Buffer*Nsample*Nsample*Nsample)/nproc[i]  # Particles split over processors
    MEM_NBODY[i] += (3.0*part*Nsample*Nsample*Nsample)/nproc[i]         # N-Body displacments
    if (MEMORY_MODE):
        MEM_NBODY[i] += (3.0*grid*Nmesh*Nmesh*(Nmesh+2.0))/nproc[i]         # Force meshes split over processors (density is deallocated)
        MEM_NBODY[i] += 3.0*grid*Nmesh*(Nmesh+2.0)                          # Each processor gets an extra slice   
    else:
        MEM_MOVE[i]  += (4.0*grid*Nmesh*Nmesh*(Nmesh+2.0))/nproc[i]         # Force and density meshes are allocated before moving particles in this case
        MEM_MOVE[i]  += 4.0*grid*Nmesh*(Nmesh+2.0)                          # Each processor gets an extra slice  
        MEM_NBODY[i] += (4.0*grid*Nmesh*Nmesh*(Nmesh+2.0))/nproc[i]         # Density mesh is not deallocated.
        MEM_NBODY[i] += 4.0*grid*Nmesh*(Nmesh+2.0)                          # Each processor gets an extra slice   

    # Find out which contribution is largest overall
    MEM[i] = max(MEM_2LPT[i],MEM_DISP[i],MEM_INIT[i],MEM_MOVE[i],MEM_DENS[i],MEM_NBODY[i])

# Convert to Gigabytes
MEM       /= (1024.0*1024.0*1024.0)
MEM_2LPT  /= (1024.0*1024.0*1024.0)
MEM_DISP  /= (1024.0*1024.0*1024.0)
MEM_MOVE /= (1024.0*1024.0*1024.0)
MEM_DENS /= (1024.0*1024.0*1024.0)
MEM_INIT  /= (1024.0*1024.0*1024.0)
MEM_NBODY /= (1024.0*1024.0*1024.0)

# Find the minimum number of processors
index = np.where(MEM < mem_limit)
if (len(index[0]) == 0):
    outputstring = "Minimum Number of Processors = None"
else:
    outputstring = "Minimum Number of Processors = "+str(index[0][0]+minproc)

# Produce the plot
# ================
fig = plt.figure(1)
ax1=fig.add_axes([0.11,0.11,0.84,0.84])
ax1.plot(nproc,MEM_2LPT,linewidth=1.75,color=(1,0,0),label='2LPT - Fields')
ax1.plot(nproc,MEM_DISP,linewidth=1.75,color=(0,0.7255,0),label='2LPT - Displacements')
ax1.plot(nproc,MEM_INIT,linewidth=1.75,color=(0,0,1),label='Initialising Particle Data')
ax1.plot(nproc,MEM_MOVE,linewidth=1.75,color=(1,0.8039,0),label='Moving Particles')
ax1.plot(nproc,MEM_DENS,linewidth=1.75,color=(1,0,1),label='Force and Density Calculation')
ax1.plot(nproc,MEM_NBODY,linewidth=1.75,color=(0,0,0),label='N-Body Displacements')
ax1.axhline(y=mem_limit, color='r', linewidth=1.75, linestyle='--')
ax1.set_xlim(minproc, maxproc)
ymin = min(MEM_2LPT[len(nproc)-1],MEM_DISP[len(nproc)-1],MEM_INIT[len(nproc)-1],MEM_DENS[len(nproc)-1],MEM_MOVE[len(nproc)-1],MEM_NBODY[len(nproc)-1])
ax1.set_ylim(ymin-0.2,2.0*mem_limit)
ax1.set_xlabel('Number of Processors',fontsize=17)
ax1.set_ylabel('Memory per Processor (GB)',fontsize=17)
ax1.xaxis.set_minor_locator(tk.MultipleLocator(20))
ax1.yaxis.set_minor_locator(tk.MultipleLocator(0.2))
ax1.tick_params('both',length=10, which='major')
ax1.tick_params('both',length=5, which='minor')
for tick in ax1.xaxis.get_ticklabels():
    tick.set_fontsize(14)
for tick in ax1.yaxis.get_ticklabels():
    tick.set_fontsize(14)
for axis in ['top','left','bottom','right']:
    ax1.spines[axis].set_linewidth(1.75)
plt.text(0.06, 0.08, outputstring, transform=ax1.transAxes, fontsize=14)
plt.legend()

plt.show()

