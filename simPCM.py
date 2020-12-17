"""
Sim PCM Pack
"""
# simPCM: Simulate parallel-connected-module packs (cells are connected in
# parallel to make modules; these modules are connected in series to make
# packs). There are no inputs or direct outputs from this script.
#
# The parameters for each cell may be different (e.g., capacity,
# resistance, etc.) 

import numpy as np
import numpy.matlib
from scipy import interpolate
import scipy.io as sio 
from models import DataModel, ModelDyn
from OCVfromSOCtemp import OCVfromSOCtemp
from pathlib import Path
import matplotlib.pyplot as plt

# read model DYN file, previously computed by dyn_model.py
modeldyn = ModelDyn.load(Path(f'./A123modeldyn.pickle'))

modeldyn.RCParam = np.array([[1.5300], [1.8172], [2.2275], [1.7639], [3.3646], [3.9327], [0.8060], [0.5807]])

# Initialize some pack configuration parameters...
Ns = 7    # Number of modules connected in series to make a pack
Np = 48   # Number of cells connected in parallel in each module

# Initialize some simulation configuration parameters...
maxtime = (1)*4800                    # Simulation run time in simulated seconds
t0 = 4800                             # Pack rests after time t0
storez = np.zeros((maxtime, Ns, Np))  # create storage for SOC
storei = np.zeros((maxtime, Ns, Np))  # create storage for current

# Initialize states for ESC cell model
nump = np.shape(modeldyn.RCParam)[1]
z = 0.25* np.ones((Ns, Np, 1))      # State Of Charge (SOC) for each cell
irc = np.zeros((Ns, Np, 1))         # R-C resistor currents for each cell
h  = np.zeros((Ns, Np, 1))          # dynamic hysteresis for each cell
s = np.zeros((Ns, Np))              # current sign for each cell

# Default initialization for cells within the pack
kvec = np.array([0, maxtime + 1]) # Iteration (time) vector for temp. profile
tvec = np.array([25, 25]) # Default temperature profile

theT, = np.where(np.array(modeldyn.temps) == tvec[0]) # Temperature profile index

q = modeldyn.QParam[theT] * np.ones((Ns, Np, 1))   # Initialize Capacity [Ah]
rc = np.exp(-1 / abs(modeldyn.RCParam[theT][0]))*np.ones((Ns,Np,nump))  # The R-C time constant parameter RjCj [s]
r = np.transpose(modeldyn.RParam[theT])         # Resistance Rj of the R-C parameter 
m = modeldyn.MParam[theT] * np.ones((Ns,Np, 1)) # Hysteresis M parameter [V] 
g = modeldyn.GParam[theT] * np.ones((Ns,Np, 1))    # Hysteresis "gamma" parameter [unitless]
r0 = modeldyn.R0Param[theT] * np.ones((Ns,Np, 1))  # Series resistance parameter R0 [Ohm]
rt = 0.000125 # 125 microOhm resistance for each tab

# Eliminate model hysteresis for rough simulation: makes results 
# easier to interpret. Then, can put hysteresis back via "m = m" 
# for more realistic results.
m = 1*m 

# Modified initialization for cell variability...
# ------------------------------------------------------------------------
# Set individual random "initial SOC" values
if True: # set to "if True" to execute, or "if 0," to skip this code  
  zmin = 0.3
  zmax = 0.7
  z = zmin + (zmax-zmin) * np.random.rand(Ns, Np, 1) # rand. init. SOC for ea. cell. a + (b-a) * random()  

# Set individual random cell-capacity values
if True: # set to "if 1," to execute, or "if 0," to skip this code
  qmin = 4.5
  q = qmin + np.random.rand(Ns, Np, 1)      # random capacity for ea. cell


# Set individual random cell-resistance relationships
if True: # set to "if 1," to execute, or "if 0," to skip this code
  rmin = 0.005
  rmax = 0.025
  r0 = rmin + (rmax-rmin) * np.random.rand(Ns, Np, 1)

r0 = r0 + 2*rt # add tab resistance to cell resistance

# Load randomized matlab 
# fName = Path('C:/Users/rsani/Google Drive/2018-1 DOCUMENTACION MAESTRIA/Investigacion/06 Herramientas desarrolladas/02 Colorado University/PCMSCM/workspace_randomized.mat')
# mat = sio.loadmat(fName)
# q = mat['q']
# rc = mat['rc']
# r = mat['r']
# m = mat['m']
# g = mat['g']
# r0 = mat['r0']

# Add faults to pack: cells faulted open- and short-circuit
# ------------------------------------------------------------------------
# To delete a PCM (open-circuit fault), set a resistance to Inf
if False:
  r0[[1,1]] = np.inf # for example...

# To delete a cell from a PCM (short-circuit fault), set its SOC to NaN
if False:
  z[[0,1],[0,0]] = np.nan # for example, delete cell 1 in PCM 1 and cell 1 in PCM 2  

Rsc = 0.0025 # Resistance value to use for cell whose SOC < 0#

# Get ready to simulate... first compute pack capacity in Ah
totalCap = min(q.sum(axis=1)) # pack capacity [Ah] = minimum module capacity [Ah]
I = 10 * totalCap # C-rate, faster simulation for I=10*totalCap (10C) 

# Okay... now to simulate pack performance using ESC cell model.
for k in range(maxtime):
  T = np.interp(k, kvec, tvec)                            # cell temperature
  v = OCVfromSOCtemp(z, T, modeldyn).reshape(Ns, Np, 1)   # get OCV for each cell: Ns * Np matrix
  v = v + m * h - irc * r.T                               # add in capacitor voltages and hysteresis
  v = (v[:, :, 0]).reshape(Ns, Np, 1)                     # remove singleton dimension
  i,j = np.where(np.isnan(z[:,:,0]))                      # Indices where cells are short-circuited
  r0[i,j] = Rsc                                           # short-circuit fault has "short-circuit" resistance
  V = (np.sum(v/r0, axis = 1) - I) / np.sum(1/r0, axis=1)
  ik = (v - np.tile(V,(1, Np)).reshape(Ns, Np, 1)) / r0
  z = z - ((1/3600) * ik) / q                             # Update each cell SOC
  z[z<0] = np.nan                                         # set over-discharged cells to short-circuit fault
  irc = rc * irc + (1 - rc) * ik                          # Update capacitor voltages
  fac = np.exp(-abs(g * ik) / (3600 * q))
  h = fac * h + (1 - fac) * np.sign(ik)                   # Update hysteresis voltages
  minz = z.min()
  maxz = z.max()
  if minz < 0.05:
    I = -abs(I) # Stop discharging 

  if maxz > 0.95:
    I = abs(I) # Stop charging

  if k > t0: 
    I = 0 # rest 

  storez[k,:,:] = z[:,:,0]   # Store SOC for later plotting
  storei[k,:,:] = ik[:,:,0]  # Store current for later plotting 

# FIGURE 1.
# Plot some sample results. In Figure 1, plot the individual SOC vs. time
# for all cells in all series PCMs. There is one subplot for each PCM.
if True:
  nonocPCMs = ~np.isinf(np.sum(r0, axis=1)) # modules that are not open-circuit faulted
  fig = plt.figure(1) # Create plot 
  plt.clf()
  t = np.arange(0, len(storez[:,:,1])) / 60
  xplots = int(round(1.0 * np.ceil(np.sqrt(Ns))))
  yplots = int(np.ceil(Ns / xplots))
  means = np.zeros((Ns, maxtime))
  for k in range(Ns):
    zr = np.squeeze(100 * storez[:,k,:])
    ax = fig.add_subplot(yplots, xplots, k+1) # Plot
    ax.plot(t, zr, linewidth=0.01)   # Plot and customize plot, linewidth = 0.3 for interactive 
    plt.title(f'Cells in PCM {k+1}', fontweight='light')  #set subplot title
    ax.set(ylabel='SOC (%)', # set labels and axis limits
          xlabel='Time (min)',
          xlim=(0, int(np.ceil(maxtime/60))),
          ylim=(0, 100),
          xticks= np.arange(0, int(np.ceil(maxtime/60))+1, 20)
    )
                 
    plt.grid()
    i,j = np.where(np.isnan(zr))  # Indices where cells are short-circuited
    zr[i,j] = 0 # exclude dead cells (failed short) from mean
      
    if nonocPCMs[k]: # exclude from average if open-circuit!
      means[k] = np.mean(zr, axis=1) # ok<AGROW>
      
  plt.tight_layout() 
  plt.savefig('./Figures/SOC in each PCM Cell.pdf', bbox_inches='tight', dpi=1200) # Save plot
  plt.show(block=False) # Show plot 

# FIGURE 2
# plot the average SOC vs. time for all PCMs

if True:
  fig = plt.figure(2)               # Create plot
  bx = fig.add_subplot(111) # Plot
  bx.plot(t, means.T, linewidth=0.1, alpha=0.8)  # Plot and customize plot 
  plt.grid()
  bx.set(ylabel='SOC (%)',          # set labels and axis limits
         xlabel='Time (min)',
         xlim=(0, np.ceil(maxtime/60)),
         ylim=(10, 80),
           xticks = np.arange(0, np.ceil(maxtime/60), 10),
           yticks = np.arange(10, 81, 10),
           autoscale_on=True
        ) 
  plt.title('Average SOC for each PCM') #set subplot title
  legendstrings = []
  for k in range(Ns): 
    if nonocPCMs[k] == True: # exclude from average if open-circuit!
      legendstrings = legendstrings + [f'PCM {k}']  # ok<AGROW> 
  plt.legend(legendstrings)      
  plt.tight_layout() 
  plt.savefig('./Figures/Mean PCM SOC.pdf', bbox_inches='tight', dpi=1200) # Save plot
  plt.show(block=False) # Show plot 

# FIGURE 3

if True:
  fig = plt.figure(3) # Create plot
  cx = fig.add_subplot(111) # Plot
  cx.plot(t, (means.max(0) - means.min(0)), linewidth=0.1, alpha=0.8)
  plt.grid()
  cx.set(ylabel='Time (min)',
        xlabel='Difference in SOC (%)',
        # xlim=(),
        # ylim=(),
        # xticks=(),
        # yticks=()
        )
  plt.title('Maxium-average SOC minus minimum-average SOC')
  plt.tight_layout()
  plt.savefig('./Figures/Max-average minus min-average SOC.pdf', bbox_inches='tight', dpi=1200) # Save plot
  plt.show(block=False) # Show plot

# FIGURE 4

if True:
  fig = plt.figure(4) # Create plot 
  plt.clf()
  t = np.arange(0, len(storei[:,:,1])) / 60
  xplots = int(round(1.0 * np.ceil(np.sqrt(Ns))))
  yplots = int(np.ceil(Ns / xplots))
  means = np.zeros((Ns, maxtime))
  for k in range(Ns):
    zr = np.squeeze(storei[:,k,:])
    dx = fig.add_subplot(yplots, xplots, k+1) # Plot
    dx.plot(t, zr, linewidth=0.01)   # Plot and customize plot 
    plt.title(f'Cells in PCM {k+1}', fontweight='light')  #set subplot title
    dx.set(ylabel='Current (A)', # set labels and axis limits
          xlabel='Time (min)',
          xlim=(0, int(np.ceil(maxtime/60))),
          ylim=(0, 100),
          xticks= np.arange(0, int(np.ceil(maxtime/60))+1, 20)
    )
                 
    plt.grid()

  plt.tight_layout() 
  plt.savefig('./Figures/Current in each cell.pdf', bbox_inches='tight', dpi=1200) # Save plot
  plt.show(block=True) # Show plot 