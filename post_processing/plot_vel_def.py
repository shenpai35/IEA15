import matplotlib as mpl
import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
from cycler import cycler
from matplotlib import cm
from matplotlib.ticker import LinearLocator
from matplotlib.backends.backend_pdf import PdfPages
from mpl_toolkits.mplot3d import Axes3

mpl.rcParams['lines.linewidth'] = 4
mpl.rcParams['axes.titlesize'] = 34
mpl.rcParams['axes.labelsize'] = 34
mpl.rcParams['xtick.labelsize'] = 22
mpl.rcParams['ytick.labelsize'] = 22
mpl.rcParams['legend.fontsize'] = 15.0
#mpl.rcParams['figure.figsize'] = (6.328, 6.328)
mpl.rcParams["figure.figsize"] = [20, 7]
#plt.style.use('classic')

rho = 1.246
d = 240 

NPS_SST = 6
Nyz = 251 # number of points in y direction
Nyy = 151 # number of points in z direction

xmin = 0.0
xmax = 0.5

ymin = -1.04
ymax = 1.04

# SST
NInit_SST = 0
NFinal_SST = 2600
u_infty = 8.0

sp_sst = nc.Dataset('/home/shuang2/scratch/act_disk/two_turbines/ws08/0_yaw/post_processing/sampling00000.nc')

NTS_SST = NFinal_SST - NInit_SST

def get_mean_velocity(yz_plane):
    """Returns average velocity on a yz plane for a given velocity component"""
    #initialize some variables
    vel_mean = np.zeros((Nyz,Nyy)) #251 by 151
    vel = sp_sst['T0_Rotor']['velocityx'] #feed in velocityx
    
    #tstep_range = np.array(range(261)) #0-260 since save every 10 steps so that takes 2600 to 260 + last time step makes 261
    #for i in tstep_range: #loop over time steps
    vel_temp_tsi = vel[260,:].reshape(NPS_SST, Nyz, Nyy)
    vel_mean += vel_temp_tsi[yz_plane,:,:]
    
    #vel_mean = vel_mean / 261
    vel_mean = vel_mean / u_infty
    
    return vel_mean

def plot_velocity_deficit():

    nrows = 1
    ncols = 5 
 
    Xtot = [1, 2, 3, 4, 5, 6] 
    X = [1] 

    Y,Z = np.meshgrid(np.linspace(-250, 250, 251), np.linspace(-150, 150, 151))

    Y = Y/d 
    Z = Z/d 

    with PdfPages('velocity_deficit_z0_' + str(int(u_infty)) + '.pdf') as pfpgs:
         velocities = ['velocityx']
         plt.figure()
         fig, axs = plt.subplots(1, 1, sharex=True, sharey=True, squeeze=False)
         fig.add_subplot(111, frameon=False)
         for iplane in range(0, 1):
             ax = axs[0, iplane]
             index = Xtot.index(X[iplane])
             uavg_sst = get_mean_velocity(index)
             uavg_sst_at_z0 = uavg_sst[:,75]
             u_def_sst = 1.0 - uavg_sst_at_z0
             print(uavg_sst_at_z0)
             sst_plot, = ax.plot(u_def_sst, Y[0,:], label = 'SST', color='blue')
             ax.set_title('x/d = {}'.format(X[iplane]))
             ax.set_xlim([xmin, xmax])
             ax.set_ylim([ymin,ymax])
         plt.tick_params(labelcolor='none', which='both', top=False, bottom=False, left=False, right=False)
         plt.xlabel('$\Delta \overline{u} / u_\infty$')
         plt.ylabel('y/d')
         plt.tight_layout()
         pfpgs.savefig()
         plt.close(fig)

if __name__=="__main__":
     plot_velocity_deficit()


