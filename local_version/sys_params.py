import numpy as np
import pandas as pd

from functools import partial
from joblib import Parallel, delayed
from multiprocessing import Pool, cpu_count
import numba as nb

from scipy.fft import fft, fftfreq, rfft, fftshift, ifft, ifftshift
from scipy.linalg import null_space, expm, eig, eigh
from scipy.integrate import  complex_ode, solve_ivp, quad, quad_vec, dblquad, simpson, trapz, nquad, quadrature
import mpmath as mp
import scipy.special as spl
import matplotlib.pyplot as plt
import plotly.graph_objects as go
from matplotlib import cm
from IPython.core.display import Image
from mpl_toolkits.mplot3d import Axes3D  # For 3D plotting
import matplotlib.animation as animation
from matplotlib.animation import FuncAnimation
import statistics as st

from scipy.interpolate import PchipInterpolator
import matplotlib
from scipy.ndimage import gaussian_filter1d

from scipy.signal import fftconvolve, wiener

import plotly.graph_objs as go
from matplotlib.patches import Circle, FancyArrow
import matplotlib.colors as mcolors

from scipy.ndimage import gaussian_filter1d
from concurrent.futures import ThreadPoolExecutor
from matplotlib.lines import Line2D 
# matplotlib.use('Agg')  
plt.rcParams['text.usetex'] = True
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif"
})
import time
from scipy.signal import decimate

config_name = "45"
mode = 'polaron'
pump = 'sym'

mp.dps = 20

############################################################
'''Set up energies and couplings'''
'''homo dimer    : E_q1 = E_q2'''
'''hetero dimer  : E_q1 != E_q2'''
############################################################

detun = 0.0
E_q1 = 1.8 * (1.0 + detun)
E_q2 = 1.8

############################################################
'''Temperatures'''
############################################################
   

 
T_opt = 1.0 * 5800.0  # Temperature
T_vib = 1.0 * 300.0 #* np.inf



class Constants:
  """
  A class to store physical constants with their units.

  This class provides attributes for commonly used physical constants and allows
  for easy retrieval of their values and units. Derived constants can be
  calculated and stored during initialization to avoid redundant calculations
  and potential errors.
  """

  def __init__(self):
    self.eV = 1.6e-19  # J
    self.hbar = 1.0545718e-34  # J s
    self.k_B = 1.380649e-23  # J/K
    self.c = 2.99792458e8  # m/s
    self.eps_0 = 8.854187817e-12  # C^2 J^(-1) m^(-1)
    self.nm = 1.0e-9 # m
    self.debye =  (1.0e-21 / (3.0e8)) # C m

    # Derived constants with units
    self.hbar_ps = self.hbar * 1.0e12  # J-ps
    self.k_B_ps = self.k_B / self.hbar_ps  # 1 / ps-K
    self.c_ps = self.c * 1.0e-12  # m/ps







############################################################
'''Dipole directions and strengths'''
############################################################
class DipoleParameters:

    def __init__(self, config_name, data_dir="files/dipole_positions/"):
        if config_name not in ["H", "J", "ortho", "oppo", "45"]:
            raise ValueError(f"Invalid configuration name: {config_name}")
        self.config_name = config_name
        self.data_dir = data_dir
        self.config_files = {
            "H": "dipole_pos_dir_H.csv",
            "J": "dipole_pos_dir_J.csv",
            "ortho": "dipole_pos_dir_ortho.csv",
            "oppo": "dipole_pos_dir_oppo.csv",
            "45": "dipole_pos_dir_45.csv"
        }

    def read_positions_from_csv(self):
        # Read positions from CSV file
        file_path = self.data_dir + self.config_files[self.config_name]
        positions_df = pd.read_csv(file_path)
        return positions_df.to_numpy()

    def calculate_dipole_properties(self, pos1, pos2, dir1, dir2):
        # Extract positions
        pos1 = pos1[:3]
        pos2 = pos2[:3]

        # Calculate distance vector
        position_vector = pos2 - pos1

        # Calculate unit vector along the direction of the dipole
        dipole_distance = np.linalg.norm(position_vector)
        dipole_vector = position_vector / dipole_distance

        # Extract direction vectors
        dir1 = np.array([dir1[0], dir1[1], dir2[2]])
        dir2 = np.array([dir2[0], dir2[1], dir2[2]])

        return dipole_vector, dir1, dir2

    def get_dipole_parameters(self):
        """
        Gets dipole parameters (r_hat, d1_hat, d2_hat) for the configured dipole configuration.

        Returns:
            tuple: A tuple containing three NumPy arrays (r_hat, d1_hat, d2_hat).
        """

        # Read dipole positions from CSV
        dipole_positions = self.read_positions_from_csv()

        r_hat = np.array([])
        d1_hat = np.array([])
        d2_hat = np.array([])

        # Calculate properties for each pair of dipoles
        for i in range(0, len(dipole_positions), 2):
            position_dipole1 = dipole_positions[i]
            position_dipole2 = dipole_positions[i + 1]

            r_h, d1_h, d2_h = self.calculate_dipole_properties(position_dipole1, position_dipole2, position_dipole1[3:], position_dipole2[3:])

            # Append values to lists
            r_hat = np.append(r_hat, r_h)
            d1_hat = np.append(d1_hat, d1_h)
            d2_hat = np.append(d2_hat, d2_h)

        return r_hat, d1_hat, d2_hat 

dipole_params = DipoleParameters(config_name)
r_hat, d1_hat, d2_hat = dipole_params.get_dipole_parameters()
# print(d1_hat)











class SystemParameters(Constants):
    """
    A class to store and calculate system parameters for a light-matter interaction.

    This class inherits from the Constants class to access physical constants
    and provides methods to calculate derived parameters like frequencies, wavevectors,
    coupling strengths, lifetimes, and decay rates. It also depends on a separate
    `DipoleParameters` class to obtain dipole unit vectors.
    """

    def __init__(self, dipole_params, E_q1 = None, E_q2 = None):
        super().__init__()  # Call the base class constructor
        self.E_q1 = E_q1  # J
        self.E_q2 = E_q2  # J
        self.dipole_params = dipole_params

    def calculate_parameters(self):
        """
        Calculates derived system parameters based on provided energies and retrieves
        dipole unit vectors from the DipoleParameters class.

        This method should be called after setting E_q1 and E_q2.
        """
        self.omega_q1 = self.E_q1 / self.hbar_ps  # ps^-1
        self.k_q1 = self.omega_q1 / self.c_ps  # m^-1
        self.omega_q2 = self.E_q2 / self.hbar_ps  # ps^-1
        self.k_q2 = self.omega_q2 / self.c_ps  # m^-1

        self.r_hat, self.d1_hat, self.d2_hat = self.dipole_params.get_dipole_parameters()

        self.d_1 = 10 * self.debye  # Dipole strength (C m)
        self.d_2 = 10 * self.debye  # Dipole strength (C m)

        self.d_1_vec = self.d_1 * self.d1_hat  # Dipole vector (C m)
        self.d_2_vec = self.d_2 * self.d2_hat  # Dipole vector (C m)

        self.r = 2 * self.nm                   # m  
        self.r_vec = self.r * self.r_hat       # distance (vector)  

        self.kappa_0 = np.dot(self.d1_hat, self.d2_hat) - 3 * (np.dot(self.d1_hat, self.r_hat)) * (np.dot(self.d2_hat, self.r_hat))
        self.J = 1.0 * (self.kappa_0 * np.abs(self.d_1) * np.abs(self.d_2) / (4 * np.pi * self.eps_0 * np.linalg.norm(self.r)**3)) / self.hbar_ps

        self.tau_L = (3 * np.pi * self.eps_0 * self.hbar_ps * (self.c_ps**3)) / (((self.omega_q1)**3) * ((self.d_1)**2))  # ps
        self.gamma_q1 = self.gamma_q2 = 1 / self.tau_L  # ps^-1

            
        # self.T_opt = 0.0 * 5800.0              # Temperature of optical bath (K)
        # self.T_vib = 0.0 * 300.0               # Temperature of vibrational bath (K)


        self.gamma_p_1 = 1.0e0 * self.gamma_q1   # Incoherent Optical Pumping
        self.gamma_p_2 = 1.0e0 * self.gamma_q2   # Incoherent Optical Pumping

    def get_gamma(self, omega):
        """
        Calculates the radiative decay rate for a given frequency.
        """
        return (omega**3) / (8 * (np.pi**2) * self.eps_0 * self.hbar_ps * (self.c_ps**3))
    
    def gamma(self, omega):
        """
        Calculates the radiative decay rate for a given frequency.
        """
        return (self.d_1**2) * (omega**3) / (3 * np.pi * self.eps_0 * self.hbar_ps * (self.c_ps**3))





# a=np.arange(10)
# print(a[-2:])
# constants = Constants()
# dipole_params = DipoleParameters(config_name="J")  
# system_params = SystemParameters(E_q1=1.8*constants.eV, E_q2=1.8*constants.eV, dipole_params=dipole_params)
# system_params.calculate_parameters()



# ############################################################
# '''Propagation'''
# ############################################################

# # initial and final time
# ti,tf = 0, 5*system_params.tau_L

# # number of steps
# steps = 200

# # times
# times = np.linspace(0.0, 10 * system_params.tau_L, steps)   # Time scale through which the system evolves in ns
# times_ns = np.linspace(0.0, 10 * system_params.tau_L * 1.0e-3, steps)   # Time scale through which the system evolves in ns

# times_ss = np.linspace(0.0, 10000 * system_params.tau_L, steps)   # Time scale through which the system evolves in ns

# # finite difference
# dt = times[1]-times[0]

# constants = Constants()
# dipole_params = DipoleParameters(config_name)  
# system_params = SystemParameters(dipole_params=dipole_params, E_q1=E_q1*constants.eV, E_q2=E_q2*constants.eV)
# system_params.calculate_parameters()


# omega_c = 90e-3 * constants.eV / system_params.hbar_ps
# E_reorg = 5.0e-3 * constants.eV / system_params.hbar_ps
# k_vib_2 = E_reorg / (2 * omega_c**3)
# q = np.linspace(0, 2.0e3, 100)
# plt.plot(q, k_vib_2 * np.abs(q)**3 * np.exp(- (np.abs(q) / omega_c)))
# plt.show()