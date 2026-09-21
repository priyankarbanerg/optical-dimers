import numpy as np
import pandas as pd
from multiprocessing import Pool, cpu_count
import os
from scipy.fft import fft, fftfreq, rfft, fftshift, ifft, ifftshift
from scipy.linalg import null_space, expm, eig, eigh
from scipy.integrate import  complex_ode, solve_ivp, quad, quad_vec, dblquad, simpson, nquad


import plotly.graph_objects as go
from matplotlib import cm

from IPython.core.display import Image
from mpl_toolkits.mplot3d import Axes3D  # For 3D plotting
import matplotlib.animation as animation
from matplotlib.animation import FuncAnimation
import statistics as st

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
import matplotlib.pyplot as plt
plt.rcParams['text.usetex'] = True
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif"
})
import time
from scipy.signal import decimate






class Constants:
  """
  A class to store physical constants with their units.
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







class ConfigParameters:
    """
    Manages simulation configurations and extracts dipole geometries from data files.

    Parameters
    ----------
    config_name : str
        Name of the spatial configuration (e.g., 'H', 'J', 'ortho', 'oppo', '45').
    mode : str
        Coupling mode for the simulation (e.g., 'polaron', 'weak').
    pump : str
        Pumping configuration scheme (e.g., 'coh', 'sym', 'site', 'none').
    data_dir : str, default: 'files/dipole_positions/'
        Directory containing the dipole position CSV files.

    Attributes
    ----------
    config_name : str
        The selected spatial configuration.
    mode : str
        The selected coupling mode.
    pump : str
        The selected pumping scheme.
    file_path : str
        The full system path to the resolved configuration CSV file.
    """

    def __init__(self, config_name: str, mode: str, pump: str, integration_method: str = 'brute force', data_dir: str = "local_version/files/dipole_positions/"):
        valid_configs = {"H", "J", "ortho", "oppo", "45"}
        if config_name not in valid_configs:
            raise ValueError(f"Invalid config_name: '{config_name}'. Must be one of {valid_configs}")
        
        self.config_name = config_name
        self.mode = mode
        self.pump = pump
        self.integration_method = integration_method
        self.file_path = os.path.join(data_dir, f"dipole_pos_dir_{config_name}.csv")

    def calculate_dipole_properties(self, pos1: np.ndarray, pos2: np.ndarray, dir1: np.ndarray, dir2: np.ndarray) -> tuple:
        """
        Calculate the relative unit vector between two dipoles and extract their directions.

        Parameters
        ----------
        pos1 : ndarray
            The [x, y, z] coordinates of the first dipole.
        pos2 : ndarray
            The [x, y, z] coordinates of the second dipole.
        dir1 : ndarray
            The directional components of the first dipole.
        dir2 : ndarray
            The directional components of the second dipole.

        Returns
        -------
        dipole_vector : ndarray
            The normalized unit vector pointing from pos1 to pos2.
        dir1_vec : ndarray
            The extracted 3D direction vector for the first dipole.
        dir2_vec : ndarray
            The extracted 3D direction vector for the second dipole.
        """
        pos1_coords = pos1[:3]
        pos2_coords = pos2[:3]

        position_vector = pos2_coords - pos1_coords
        dipole_distance = np.linalg.norm(position_vector)
        dipole_vector = position_vector / dipole_distance

        dir1_vec = np.array([dir1[0], dir1[1], dir1[2]])
        dir2_vec = np.array([dir2[0], dir2[1], dir2[2]])

        return dipole_vector, dir1_vec, dir2_vec

    def get_dipole_parameters(self) -> tuple:
        """
        Compute the directional parameters for all dipoles.

        Returns
        -------
        r_hat : ndarray
            Array of unit distance vectors between each dipole pair.
        d1_hat : ndarray
            Array of direction unit vectors for the first dipole in each pair.
        d2_hat : ndarray
            Array of direction unit vectors for the second dipole in each pair.
        """
        dipole_positions = pd.read_csv(self.file_path).to_numpy()

        r_hat_list = []
        d1_hat_list = []
        d2_hat_list = []

        for i in range(0, len(dipole_positions), 2):
            pos1 = dipole_positions[i]
            pos2 = dipole_positions[i + 1]

            r_h, d1_h, d2_h = self.calculate_dipole_properties(
                pos1, pos2, pos1[3:], pos2[3:]
            )

            r_hat_list.append(r_h)
            d1_hat_list.append(d1_h)
            d2_hat_list.append(d2_h)

        if len(r_hat_list) != 1:
            raise ValueError(f"Expected one dipole pair, found {len(r_hat_list)} pairs")
        
        return r_hat_list[0], d1_hat_list[0], d2_hat_list[0]
        
        



class SystemParameters:
    """
    Store and calculate physical parameters for light-matter interaction.

    Parameters
    ----------
    constants : Constants
        Instantiated object containing universal physical constants.
    config_params : ConfigParameters
        Configuration object containing dipole geometry and directional vectors.
    E_q1_eV : float
        Transition energy of the first dipole in eV.
    E_q2_eV : float
        Transition energy of the second dipole in eV.
    d1_debye : float, default: 10.0
        Dipole moment magnitude of the first dipole in Debye.
    d2_debye : float, default: 10.0
        Dipole moment magnitude of the second dipole in Debye.
    gamma_p_1 : float, default: 0.0
        Incoherent pumping rate for the first dipole in ps^-1.
    gamma_p_2 : float, default: 0.0
        Incoherent pumping rate for the second dipole in ps^-1.
    Omega_1 : float, default: 0.0
        Coherent driving Rabi frequency for the first dipole in ps^-1.
    Omega_2 : float, default: 0.0
        Coherent driving Rabi frequency for the second dipole in ps^-1.
    r_nm : float, default: 2.0
        Separation distance between the dipoles in nanometers.
    T_opt : float, optional
        Temperature of the optical bath in Kelvin.
    T_vib : float, optional
        Temperature of the vibrational bath in Kelvin.
    """
    def __init__(self, constants, config_params, E_q1_eV: float, E_q2_eV: float, d1_debye: float = 10.0, d2_debye: float = 10.0, gamma_p_1: float = 1.0, gamma_p_2: float = 1.0, Omega_1: float = 0.0, Omega_2: float = 0.0, omega_L: float = 0.0, r_nm: float = 2.0, T_opt: float = 0.0, T_vib: float = 0.0):
        self.constants = constants
        self.config_params = config_params
        self.E_q1_eV = E_q1_eV
        self.E_q2_eV = E_q2_eV
        self.d1_debye = d1_debye
        self.d2_debye = d2_debye
        self.gamma_p_1 = gamma_p_1
        self.gamma_p_2 = gamma_p_2
        self.Omega_1 = Omega_1
        self.Omega_2 = Omega_2
        self.omega_L = omega_L
        self.r_nm = r_nm
        self.T_opt = T_opt
        self.T_vib = T_vib
        
    def calculate_parameters(self) -> None:
        """
        Compute derived quantum optical parameters from initial energies and spatial configuration.
        """
        # Convert input eV to Joules
        self.E_q1 = self.E_q1_eV * self.constants.eV
        self.E_q2 = self.E_q2_eV * self.constants.eV

        self.omega_q1 = self.E_q1 / self.constants.hbar_ps
        self.k_q1 = self.omega_q1 / self.constants.c_ps
        
        self.omega_q2 = self.E_q2 / self.constants.hbar_ps
        self.k_q2 = self.omega_q2 / self.constants.c_ps

        self.r_hat, self.d1_hat, self.d2_hat = self.config_params.get_dipole_parameters()

        # Convert Debye to Coulomb-meters
        self.d_1 = self.d1_debye * self.constants.debye
        self.d_2 = self.d2_debye * self.constants.debye

        self.d_1_vec = self.d_1 * self.d1_hat
        self.d_2_vec = self.d_2 * self.d2_hat

        self.r = self.r_nm * self.constants.nm
        self.r_vec = self.r * self.r_hat

        self.kappa_0 = np.dot(self.d1_hat, self.d2_hat) - 3.0 * np.dot(self.d1_hat, self.r_hat) * np.dot(self.d2_hat, self.r_hat)
        self.J = (self.kappa_0 * np.abs(self.d_1) * np.abs(self.d_2) / (4.0 * np.pi * self.constants.eps_0 * np.linalg.norm(self.r)**3)) / self.constants.hbar_ps
        
        self.tau_L1 = (3.0 * np.pi * self.constants.eps_0 * self.constants.hbar_ps * (self.constants.c_ps**3)) / ((self.omega_q1**3) * (self.d_1**2))
        self.tau_L2 = (3.0 * np.pi * self.constants.eps_0 * self.constants.hbar_ps * (self.constants.c_ps**3)) / ((self.omega_q2**3) * (self.d_2**2))
        
        self.gamma_q1 = 1.0 / self.tau_L1
        self.gamma_q2 = 1.0 / self.tau_L2


    def get_gamma(self, omega: float) -> float:
        """
        Calculate the bare radiative decay rate for a given frequency.

        Parameters
        ----------
        omega : float
            Angular frequency in ps^-1.

        Returns
        -------
        float
            Radiative decay rate.
        """
        return (omega**3) / (8.0 * (np.pi**2) * self.constants.eps_0 * self.constants.hbar_ps * (self.constants.c_ps**3))
    
    def gamma(self, omega: float, d_magnitude: float) -> float:
        """
        Calculate the dipole-weighted radiative decay rate for a given frequency.

        Parameters
        ----------
        omega : float
            Angular frequency in ps^-1.
        d_magnitude : float
            Dipole moment magnitude in Coulomb-meters.

        Returns
        -------
        float
            Dipole-weighted radiative decay rate.
        """
        return (d_magnitude**2) * (omega**3) / (3.0 * np.pi * self.constants.eps_0 * self.constants.hbar_ps * (self.constants.c_ps**3))