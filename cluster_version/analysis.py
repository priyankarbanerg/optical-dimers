from system_ops import *
from br_tensor import *
from scipy.linalg import expm
import os
 

class QuantumSystemAnalysis:
    def __init__(self, BR_calculator, constants, dipole_params, system_params, comb_rates, quantity, rho, final_time=None, steps=None, state=None):
        """
        Initialize the QuantumSystemAnalysis class.

        Parameters:
        - BR_calculator: BlochRedfieldCalculator instance for calculating BR tensors.
        - quantity: String specifying the quantity to compute ('intensity', 'population', 'g2', etc.).
        - rho: Initial density matrix of the system.
        - times: Array of time points for the calculation.
        - steps: Number of steps for the calculation.
        """
        self.BR_calculator = BR_calculator
        self.quantity = quantity
        self.H = BR_calculator.H
        self.dim = len(self.H)
        self.ekets = BR_calculator.ekets
        self.tf = final_time
        self.steps = steps
        self.rho = rho
        self.dt = self.tf / self.steps
        self.P = expm((self.BR_calculator.total_BR_tensor()) * self.dt)
        self.times = np.linspace(0, self.tf, steps)

        self.state = state
        self.constants = constants
        self.dipole_params = dipole_params
        self.system_params = system_params
        # Calculate spectral densities and rates
        self.comb_rates = comb_rates
        self.J_opt, self.J_vib = self.comb_rates.J_opt, self.comb_rates.J_opt
        self.r_hat, self.d1_hat, self.d2_hat = dipole_params.get_dipole_parameters()

        self.kappa = self.comb_rates.kappa()
        self._initialize_projections()

    def _initialize_projections(self):
        """
        Initialize the projection operators.
        """
        J_scaled = self.system_params.J * self.kappa**2 + 1.0e-8
        eta = np.sqrt(detun**2 + J_scaled**2)
        chi = np.arcsin(J_scaled / eta)

        self.psi_g1g2 = np.array([0, 0, 0, 1], dtype=complex)
        self.proj_g1g2 =  np.outer(self.psi_g1g2, self.psi_g1g2)  

        self.psi_e1g2 = np.array([0, 0, 1, 0], dtype=complex)
        self.proj_e1g2 = np.outer(self.psi_e1g2, self.psi_e1g2)

        self.psi_e2g1 = np.array([0, 1, 0, 0], dtype=complex)
        self.proj_e2g1 = np.outer(self.psi_e2g1, self.psi_e2g1)

        self.proj_e2g1e1g2 = np.outer(self.psi_e2g1, self.psi_e1g2)
        self.proj_e1g2e2g1 = np.outer(self.psi_e1g2, self.psi_e2g1)

        self.psi_e1e2 = np.array([1, 0, 0, 0], dtype=complex)
        self.proj_e1e2 = np.outer(self.psi_e1e2, self.psi_e1e2) 

        if self.system_params.J >= 0:
            self.psi_bright = np.cos(chi / 2) * self.psi_e1g2 + np.sin(chi / 2) * self.psi_e2g1 #self.ekets.conj().T @ self.ekets[:, 1]
            self.proj_bright = np.outer(self.psi_bright, self.psi_bright)
            
            self.psi_dark = -np.sin(chi / 2) * self.psi_e1g2 + np.cos(chi / 2) * self.psi_e2g1 #self.ekets.conj().T @ self.ekets[:, 2]
            self.proj_dark = np.outer(self.psi_dark, self.psi_dark) 


        else:
            self.psi_bright = -np.sin(chi / 2) * self.psi_e1g2 + np.cos(chi / 2) * self.psi_e2g1 #self.ekets.conj().T @ self.ekets[:, 2]
            self.proj_bright = np.outer(self.psi_bright, self.psi_bright)

            self.psi_dark = np.cos(chi / 2) * self.psi_e1g2 + np.sin(chi / 2) * self.psi_e2g1 #self.ekets.conj().T @ self.ekets[:, 1]
            self.proj_dark = np.outer(self.psi_dark, self.psi_dark) 


    def calculate_population(self):
        """
        Calculate population using the BR tensor.

        Returns:
        - population: Array of population values.
        """
        if self.state == 'ground':
            proj = self.proj_g1g2
        elif self.state == 'doubly excited':
            proj = self.proj_e1e2
        elif self.state == 'bright':
            proj = (self.proj_e1g2 + self.proj_e2g1 + self.proj_e1g2e2g1 + self.proj_e2g1e1g2) / 2
        elif self.state == 'dark':
            proj = (self.proj_e1g2 + self.proj_e2g1 - self.proj_e1g2e2g1 - self.proj_e2g1e1g2) / 2

        population = []
        for _ in range(self.steps):
            pop = np.real(np.trace(np.dot(proj, self.rho)))
            population.append(pop)
            
            # Propagate the density matrix
            rho_vec = np.reshape(self.rho, (self.dim**2, 1))
            rho_vec = self.P @ rho_vec
            self.rho = np.reshape(rho_vec, (self.dim, self.dim))
        return population

    def calculate_intensity(self):
        """
        Calculate intensity over time using the BR tensor.

        Returns:
        - norm_intensity: Array of normalized intensity values.
        """
        norm_intensity = []
        for _ in range(self.steps):
            int_e1g2 = np.real(np.trace(np.dot(self.proj_e1g2, self.rho)))
            int_e2g1 = np.real(np.trace(np.dot(self.proj_e2g1, self.rho)))
            int_e2g1e1g2 = np.real(np.trace(np.dot(self.proj_e2g1e1g2, self.rho))) * (self.kappa**2) * (self.d1_hat @ self.d2_hat)
            int_e1g2e2g1 = np.real(np.trace(np.dot(self.proj_e1g2e2g1, self.rho))) * (self.kappa**2) * (self.d1_hat @ self.d2_hat)
            int_e1e2 = np.real(np.trace(np.dot(self.proj_e1e2, self.rho))) + np.real(np.trace(np.dot(self.proj_e1e2, self.rho)))
            norm_intensity.append(int_e1g2 + int_e2g1 + int_e2g1e1g2 + int_e1g2e2g1 + int_e1e2)
            
            # Propagate the density matrix
            rho_vec = np.reshape(self.rho, (self.dim**2, 1))
            rho_vec = self.P @ rho_vec
            self.rho = np.reshape(rho_vec, (self.dim, self.dim))
        # print(np.round(np.array([int_e1g2, int_e2g1, int_e2g1e1g2, int_e1g2e2g1, int_e1e2]), 8))
        return norm_intensity

    # def calculate_intensity(self):
    #     """
    #     Calculate intensity over time using the BR tensor.

    #     Returns:
    #     - norm_intensity: Array of normalized intensity values.
    #     """
    #     norm_intensity = []
    #     gamma_p = self.system_params.gamma(self.system_params.omega_q1 + self.system_params.J) 
    #     gamma_m = self.system_params.gamma(self.system_params.omega_q2 - self.system_params.J)
    #     for _ in self.times:
    #         int_e1g2 = (gamma_p + gamma_m) * np.real(np.trace(np.dot(self.proj_e1g2, self.rho)))
    #         int_e2g1 = (gamma_p + gamma_m) * np.real(np.trace(np.dot(self.proj_e2g1, self.rho)))
    #         int_e2g1e1g2 = ((np.sqrt(gamma_p * gamma_p) + np.sqrt(gamma_m * gamma_m)) * (self.kappa**2) * (self.d1_hat @ self.d2_hat) + 0.5 * (gamma_p - gamma_m) - 0.5 * (gamma_p - gamma_m)) * np.real(np.trace(np.dot(self.proj_e2g1e1g2, self.rho))) 
    #         int_e1g2e2g1 = ((np.sqrt(gamma_p * gamma_p) + np.sqrt(gamma_m * gamma_m)) * (self.kappa**2) * (self.d1_hat @ self.d2_hat) + 0.5 * (gamma_p - gamma_m) - 0.5 * (gamma_p - gamma_m)) * np.real(np.trace(np.dot(self.proj_e1g2e2g1, self.rho))) 
    #         int_e1e2 = ((gamma_p + gamma_p) + (gamma_m + gamma_m)) * np.real(np.trace(np.dot(self.proj_e1e2, self.rho)))
    #         norm_intensity.append(0.5 * (int_e1g2 + int_e2g1 + int_e2g1e1g2 + int_e1g2e2g1 + int_e1e2) / self.system_params.gamma_q1) #(gamma_p + gamma_m)) #
        
    #         # Propagate the density matrix
    #         rho_vec = np.reshape(self.rho, (self.dim**2, 1))
    #         rho_vec = self.P @ rho_vec
    #         self.rho = np.reshape(rho_vec, (self.dim, self.dim))
    #     return norm_intensity


    def calculate_g2(self):
        """
        Calculate the second-order correlation function g2.

        Returns:
        - g2: g2 value.
        """
        BR_tensor = self.BR_calculator.total_BR_tensor()
        g2 = (np.trace(BR_tensor) ** 2).real  # Placeholder
        return g2


    def compute(self):
        """
        Compute the specified quantity.

        Returns:
        - Result of the calculation.
        """
        if self.quantity == 'population':
            return self.calculate_population()
        elif self.quantity == 'intensity':
            return self.calculate_intensity()
        elif self.quantity == 'g2':
            return self.calculate_g2()
        else:
            raise ValueError("Unknown quantity. Available options are: 'intensity', 'population', 'g2'.")



class DistributionAnalysis:
    def __init__(self, BR_calculator, constants, dipole_params, system_params, comb_rates, quantity, rho, final_time=None, steps=None):
        """
        Initialize the DistributionAnalysis class.

        Parameters:
        - BR_calculator: BlochRedfieldCalculator instance for calculating BR tensors.
        - quantity: String specifying the quantity to compute ('intensity distribution', 'g2 distribution').
        - rho: Initial density matrix of the system.
        - final_time: Final time for the simulation.
        - steps: Number of steps for the calculation.
        """
        self.BR_calculator = BR_calculator
        self.quantity = quantity
        self.H = BR_calculator.H
        self.dim = len(self.H)
        self.ekets = BR_calculator.ekets
        self.tf = final_time
        self.steps = steps
        self.rho = rho
        self.dt = self.tf / self.steps
        self.P = expm((self.BR_calculator.total_BR_tensor()) * self.dt)
        # self.P = expm((self.BR_calculator.Unitary_Dynamics() + self.BR_calculator.BR_tensor_opt() + self.BR_calculator.Liouvillian()) * self.dt)
        # self.P = expm((self.BR_calculator.Unitary_Dynamics() + self.BR_calculator.BR_tensor_opt() + self.BR_calculator.BR_tensor_vibcoup() + self.BR_calculator.Liouvillian()) * self.dt)
        self.times = np.linspace(0, self.tf, steps)
        self.times_ss = np.linspace(0, 1.0e3 * self.tf, steps)

        self.constants = constants
        self.dipole_params = dipole_params
        self.system_params = system_params
        # Calculate spectral densities and rates
        self.comb_rates = comb_rates
        self.J_opt, self.J_vib = self.comb_rates.J_opt, self.comb_rates.J_opt
        self.r_hat, self.d1_hat, self.d2_hat = dipole_params.get_dipole_parameters()

        self.kappa = self.comb_rates.kappa()
        self._initialize_projections()

    def _initialize_projections(self):
        """
        Initialize the projection operators.
        """
        J_scaled = self.system_params.J * self.kappa**2 + 1.0e-8
        eta = np.sqrt(detun**2 + J_scaled**2)
        chi = np.arcsin(J_scaled / eta)

        self.psi_g1g2 = np.array([0, 0, 0, 1], dtype=complex)
        self.proj_g1g2 =  np.outer(self.psi_g1g2, self.psi_g1g2)  

        self.psi_e1g2 = np.array([0, 0, 1, 0], dtype=complex)
        self.proj_e1g2 = np.outer(self.psi_e1g2, self.psi_e1g2)

        self.psi_e2g1 = np.array([0, 1, 0, 0], dtype=complex)
        self.proj_e2g1 = np.outer(self.psi_e2g1, self.psi_e2g1)

        self.proj_e2g1e1g2 = np.outer(self.psi_e2g1, self.psi_e1g2)
        self.proj_e1g2e2g1 = np.outer(self.psi_e1g2, self.psi_e2g1)

        self.psi_e1e2 = np.array([1, 0, 0, 0], dtype=complex)
        self.proj_e1e2 = np.outer(self.psi_e1e2, self.psi_e1e2) 

        if self.system_params.J >= 0:
            self.psi_bright = np.cos(chi / 2) * self.psi_e1g2 + np.sin(chi / 2) * self.psi_e2g1 #self.ekets.conj().T @ self.ekets[:, 1]
            self.proj_bright = np.outer(self.psi_bright, self.psi_bright)
            
            self.psi_dark = -np.sin(chi / 2) * self.psi_e1g2 + np.cos(chi / 2) * self.psi_e2g1 #self.ekets.conj().T @ self.ekets[:, 2]
            self.proj_dark = np.outer(self.psi_dark, self.psi_dark) 


        else:
            self.psi_bright = -np.sin(chi / 2) * self.psi_e1g2 + np.cos(chi / 2) * self.psi_e2g1 #self.ekets.conj().T @ self.ekets[:, 2]
            self.proj_bright = np.outer(self.psi_bright, self.psi_bright)

            self.psi_dark = np.cos(chi / 2) * self.psi_e1g2 + np.sin(chi / 2) * self.psi_e2g1 #self.ekets.conj().T @ self.ekets[:, 1]
            self.proj_dark = np.outer(self.psi_dark, self.psi_dark) 


    def intensity_distribution(self, theta, phi):
        """
        Calculate intensity distribution at steady state.

        Parameters:
        - theta: Angle in radians.
        - phi: Angle in radians.

        Returns:
        - intensity_dist: Intensity distribution value.
        """
        intensity_dist = []

        e_1_vec, e_2_vec = orthogonal_vectors(theta, phi)
        e_vec = np.array([e_1_vec, e_2_vec])

        # Evolve system
        # for _ in self.times:
        intensity = 0.0
            # Calculate intensity distribution
        for i in range(2):
            d1_pol = self.system_params.d_1_vec @ e_vec[i]
            d2_pol = self.system_params.d_2_vec @ e_vec[i]

            # if d1_pol == d2_pol == 0.0:
            #     psi_k = np.zeros(self.dim, dtype=complex)
            # else:
            #     Norm = np.sqrt(d1_pol**2 + d2_pol**2)
            psi_k_1 = ((d1_pol * np.exp(-1j * self.system_params.k_q1 * self.system_params.r * np.cos(theta) / 2)) * self.psi_e1g2 +
                        (d2_pol * np.exp(1j * self.system_params.k_q1 * self.system_params.r * np.cos(theta) / 2)) * self.psi_e2g1) 
            psi_k_2 = ((d1_pol * np.exp(1j * self.system_params.k_q1 * self.system_params.r * np.cos(theta) / 2)) * self.psi_e2g1 +
                        (d2_pol * np.exp(-1j * self.system_params.k_q1 * self.system_params.r * np.cos(theta) / 2)) * self.psi_e1g2) 
            
            sk_m = self.ekets.conj().T @ (np.outer(self.psi_g1g2, psi_k_1.conj()) + np.outer(psi_k_2, self.psi_e1e2)) @ self.ekets
            sk_p = sk_m.conj().T

            intensity += np.real(np.trace(np.dot(sk_p @ sk_m, self.rho))) * (self.system_params.get_gamma(self.system_params.omega_q1) / self.system_params.gamma_q1)
    
            
            # intensity_dist.append(intensity)
            # rho_vec = np.reshape(self.rho, (self.dim**2, 1))
            # rho_vec = self.P @ rho_vec
            # self.rho = np.reshape(rho_vec, (self.dim, self.dim))
        return intensity #intensity_dist


    # def calculate_steady_state_g2_correlations(self, theta, phi, theta_p, phi_p):
    #     """
    #     Calculate the g2 correlations at steady state.

    #     Parameters:
    #     - H: Hamiltonian of the system.
    #     - rho0: Initial density matrix.
    #     - Gamma_opt: Optical decay rates.
    #     - Gamma_vib: Vibrational decay rates.
    #     - theta: Angle in radians for first polarization.
    #     - phi: Angle in radians for first polarization.
    #     - theta_p: Angle in radians for second polarization.
    #     - phi_p: Angle in radians for second polarization.
    #     - k_q1: Wave vector component.

    #     Returns:
    #     - g2_corr: g2 correlation values.
    #     """
    #     g2_corr = 0.0
    #     # Initialize state
    #     rho_ss = self.rho
    #     rho_qr = np.zeros_like(self.rho, dtype = complex) # empty_dm

    #     # Propagate for a_ops
    #     for _ in self.times_ss:
    #         rho_vec_ss = np.reshape(rho_ss, (self.dim**2, 1))
    #         rho_vec_ss = self.P @ rho_vec_ss
    #         rho_ss = np.reshape(rho_vec_ss, (self.dim, self.dim))

    #     e_1_vec_0, e_2_vec_0 = orthogonal_vectors(theta, phi)
    #     e_vec_0 = np.array([e_1_vec_0, e_2_vec_0])

    #     e_1_vec_1, e_2_vec_1 = orthogonal_vectors(theta_p, phi_p)
    #     e_vec_1 = np.array([e_1_vec_1, e_2_vec_1])


    #     for i in range(2):
    #         d1_pol_0 = self.system_params.d_1_vec @ e_vec_0[i]
    #         d2_pol_0 = self.system_params.d_2_vec @ e_vec_0[i]


    #         psi_k1_0 = ((d1_pol_0 * np.exp(-1j * self.system_params.k_q1 * self.system_params.r * np.cos(theta) / 2)) * self.psi_e1g2 +
    #                    (d2_pol_0 * np.exp(1j * self.system_params.k_q1 * self.system_params.r * np.cos(theta) / 2)) * self.psi_e2g1) 
            
    #         psi_k2_0 = ((d1_pol_0 * np.exp(-1j * self.system_params.k_q1 * self.system_params.r * np.cos(theta) / 2)) * self.psi_e2g1 +
    #                    (d2_pol_0 * np.exp(1j * self.system_params.k_q1 * self.system_params.r * np.cos(theta) / 2)) * self.psi_e1g2) 

    #         sk_m_0 = self.ekets.conj().T @ (np.outer(self.psi_g1g2, psi_k1_0.conj()) + np.outer(psi_k2_0, self.psi_e1e2)) @ self.ekets
    #         sk_p_0 = sk_m_0.conj().T   

    #         rho_qr_0 = sk_m_0 @ rho_ss @ sk_p_0
    #         rho_qr += rho_qr_0


    #         for j in range(2):
    #             d1_pol_1 = self.system_params.d_1_vec @ e_vec_1[j]
    #             d2_pol_1 = self.system_params.d_2_vec @ e_vec_1[j]


    #             psi_k1_1 = ((d1_pol_1 * np.exp(-1j * self.system_params.k_q1 * self.system_params.r * np.cos(theta_p) / 2)) * self.psi_e1g2 +
    #                     (d2_pol_1 * np.exp(1j * self.system_params.k_q1 * self.system_params.r * np.cos(theta_p) / 2)) * self.psi_e2g1) 
                
    #             psi_k2_1 = ((d1_pol_1 * np.exp(-1j * self.system_params.k_q1 * self.system_params.r * np.cos(theta_p) / 2)) * self.psi_e2g1 +
    #                    (d2_pol_1 * np.exp(1j * self.system_params.k_q1 * self.system_params.r * np.cos(theta_p) / 2)) * self.psi_e1g2) 

    #             sk_m_1 = self.ekets.conj().T @ (np.outer(self.psi_g1g2, psi_k1_1.conj()) + np.outer(psi_k2_1, self.psi_e1e2)) @ self.ekets
    #             sk_p_1 = sk_m_1.conj().T

                
    #             g2_corr += np.real(np.trace(np.dot(sk_p_1 @ sk_m_1, rho_qr))) * (3 / (8 * np.pi))**2 #* (self.system_params.get_gamma(self.system_params.omega_q1) / self.system_params.gamma_q1)**2
                

    #             # g2 = self.ekets.conj().T @ np.outer(self.psi_e1e2, psi_k_0.conj()) @ np.outer(psi_k_1, psi_k_1.conj()) @ np.outer(psi_k_0, self.psi_e1e2) @ self.ekets 
    #             # g2_corr += np.real(np.trace(np.dot(g2, rho_ss)))* (3 / (8 * np.pi))**2 #* (self.system_params.get_gamma(self.system_params.omega_q1) / self.system_params.gamma_q1)**2
    
    #     return np.array(g2_corr) #/ g2_corr[-1]





    def calculate_steady_state_g2_correlations(self, theta, phi, theta_p, phi_p):
        """
        Calculate the g2 correlations at steady state.

        Parameters:
        - H: Hamiltonian of the system.
        - rho0: Initial density matrix.
        - Gamma_opt: Optical decay rates.
        - Gamma_vib: Vibrational decay rates.
        - theta: Angle in radians for first polarization.
        - phi: Angle in radians for first polarization.
        - theta_p: Angle in radians for second polarization.
        - phi_p: Angle in radians for second polarization.
        - k_q1: Wave vector component.

        Returns:
        - g2_corr: g2 correlation values.
        """
        g2_corr = []
        # Initialize state
        rho_ss = self.rho
        rho_qr = np.zeros_like(self.rho, dtype = complex) # empty_dm

        I0 = 0

        # Propagate for a_ops
        for _ in self.times_ss:
            rho_vec_ss = np.reshape(rho_ss, (self.dim**2, 1))
            rho_vec_ss = self.P @ rho_vec_ss
            rho_ss = np.reshape(rho_vec_ss, (self.dim, self.dim))
            # print(rho_ss.real)

        k_vec_0 = self.system_params.k_q1 * np.array([np.cos(phi) * np.sin(theta), np.sin(phi) * np.sin(theta), np.cos(theta)])
        e_1_vec_0, e_2_vec_0 = orthogonal_vectors(theta, phi)
        e_vec_0 = np.array([e_1_vec_0, e_2_vec_0])

        k_vec_1 = self.system_params.k_q1 * np.array([np.cos(phi_p) * np.sin(theta_p), np.sin(phi_p) * np.sin(theta_p), np.cos(theta_p)])
        e_1_vec_1, e_2_vec_1 = orthogonal_vectors(theta_p, phi_p)
        e_vec_1 = np.array([e_1_vec_1, e_2_vec_1])


        for i in range(2):
            d1_pol_0 = self.system_params.d1_hat @ e_vec_0[i]
            d2_pol_0 = self.system_params.d2_hat @ e_vec_0[i]
            
            if d1_pol_0 == d2_pol_0 == 0.0:
                psi_k1_0 = np.zeros(self.dim, dtype=complex)
                psi_k2_0 = np.zeros(self.dim, dtype=complex)

            else:
                Norm_0 = np.sqrt(d1_pol_0**2 + d2_pol_0**2)

                psi_k1_0 = ((d1_pol_0 * np.exp(-1j * (k_vec_0 @ self.system_params.r_vec) / 2)) * self.psi_e1g2 +
                            (d2_pol_0 * np.exp(1j * (k_vec_0 @ self.system_params.r_vec) / 2)) * self.psi_e2g1) 
                
                psi_k2_0 = ((d1_pol_0 * np.exp(1j * (k_vec_0 @ self.system_params.r_vec) / 2)) * self.psi_e2g1 +
                            (d2_pol_0 * np.exp(-1j * (k_vec_0 @ self.system_params.r_vec) / 2)) * self.psi_e1g2) 

                # print(k_vec_0, self.system_params.r_vec, (k_vec_0 @ self.system_params.r_vec))
            sk_m_0 = self.ekets.conj().T @ (np.outer(self.psi_g1g2, psi_k1_0.conj()) + np.outer(psi_k2_0, self.psi_e1e2)) @ self.ekets
            sk_p_0 = sk_m_0.conj().T   

            rho_qr_0 = sk_m_0 @ rho_ss @ sk_p_0 #* (self.system_params.get_gamma(self.system_params.omega_q1) / self.system_params.gamma_q1) 
            rho_qr += rho_qr_0

            I0 += np.real(np.trace(np.dot(sk_p_0 @ sk_m_0, rho_ss))) 

        for _ in np.linspace(0, 5.0 * self.system_params.tau_L, self.steps):
            g2 = 0.0

            for j in range(2):
                d1_pol_1 = self.system_params.d1_hat @ e_vec_1[j]
                d2_pol_1 = self.system_params.d2_hat @ e_vec_1[j]

                if d1_pol_1 == d2_pol_1 == 0.0:
                    psi_k1_1 = np.zeros(self.dim, dtype=complex)
                    psi_k2_1 = np.zeros(self.dim, dtype=complex)

                else:
                    Norm_1 = np.sqrt(d1_pol_1**2 + d2_pol_1**2)


                    psi_k1_1 = ((d1_pol_1 * np.exp(-1j * (k_vec_1 @ self.system_params.r_vec) / 2)) * self.psi_e1g2 +
                            (d2_pol_1 * np.exp(1j * (k_vec_1 @ self.system_params.r_vec) / 2)) * self.psi_e2g1) 
                    
                    psi_k2_1 = ((d1_pol_1 * np.exp(1j * (k_vec_1 @ self.system_params.r_vec) / 2)) * self.psi_e2g1 +
                            (d2_pol_1 * np.exp(-1j * (k_vec_1 @ self.system_params.r_vec) / 2)) * self.psi_e1g2) 

                sk_m_1 = self.ekets.conj().T @ (np.outer(self.psi_g1g2, psi_k1_1.conj()) + np.outer(psi_k2_1, self.psi_e1e2)) @ self.ekets
                sk_p_1 = sk_m_1.conj().T   

                sp_sm = sk_p_1 @ sk_m_1 #* (self.system_params.get_gamma(self.system_params.omega_q1) / self.system_params.gamma_q1)
                
                g2 += np.real(np.trace(np.dot(sp_sm, rho_qr))) 
               
    # #             g2 = self.ekets.conj().T @ np.outer(self.psi_e1e2, psi_k_0.conj()) @ np.outer(psi_k_1, psi_k_1.conj()) @ np.outer(psi_k_0, self.psi_e1e2) @ self.ekets 
    # #             g2_corr += n.conj()p.real(np.trace(np.dot(g2, rho_ss)))* (3 / (8 * np.pi))**2 #* (self.system_params.get_gamma(self.system_params.omega_q1) / self.system_params.gamma_q1)**2

            g2_corr.append(g2)
            rho_vec = np.reshape(rho_qr, (self.dim**2, 1))
            rho_vec = self.P @ rho_vec
            rho_qr = np.reshape(rho_vec, (self.dim, self.dim))
        
        print('next')
        return ((np.array(g2_corr)))  #/ ((g2_corr[-1]))#[0]





    def compute(self, theta=None, phi=None, theta_p=None, phi_p=None):
        """
        Compute the specified quantity.

        Parameters:
        - theta: Angle in radians (required for intensity distribution).
        - phi: Angle in radians (required for intensity distribution).
        - theta_p: Angle in radians for g2 distribution (optional).
        - phi_p: Angle in radians for g2 distribution (optional).
        - k_q1: Wave vector component for g2 distribution (optional).

        Returns:
        - Result of the calculation.
        """
        if self.quantity == 'intensity distribution':
            if theta is None or phi is None:
                raise ValueError("Theta and Phi must be provided for intensity distribution calculation.")
            return self.intensity_distribution(theta, phi)
        elif self.quantity == 'g2 distribution':
            if theta is None or phi is None or theta_p is None or phi_p is None:
                raise ValueError("Theta, Phi, Theta_p, Phi_p, and k_q1 must be provided for g2 distribution calculation.")
            return self.calculate_steady_state_g2_correlations(theta, phi, theta_p, phi_p)
        else:
            raise ValueError("Unknown quantity. Available options are: 'intensity distribution', 'g2 distribution'.")










    # def intensity_distribution(self, theta, phi):
    #     """
    #     Calculate intensity distribution at steady state.

    #     Parameters:
    #     - theta: Angle in radians.
    #     - phi: Angle in radians.

    #     Returns:
    #     - intensity_dist: Intensity distribution value.
    #     """
    #     intensity_dist = []

    #     e_1_vec, e_2_vec = orthogonal_vectors(theta, phi)
    #     e_vec = np.array([e_1_vec, e_2_vec])


    #     intensity = 0.0
    #     # Calculate intensity distribution
    #     for i in range(2):
    #         d1_pol = self.system_params.d1_hat @ e_vec[i]
    #         d2_pol = self.system_params.d2_hat @ e_vec[i]

    #         if d1_pol == d2_pol == 0.0:
    #             psi_k = np.zeros(self.dim, dtype=complex)
    #         else:
    #             Norm = np.sqrt(d1_pol**2 + d2_pol**2)
    #             psi_k = ((d1_pol * np.exp(-1j * self.system_params.k_q1 * self.system_params.r * np.cos(theta) / 2)) * self.psi_e1g2 +
    #                     (d2_pol * np.exp(1j * self.system_params.k_q1 * self.system_params.r * np.cos(theta) / 2)) * self.psi_e2g1) #/ Norm

    #         sk_m = self.ekets.conj().T @ (np.outer(self.psi_g1g2, psi_k) + np.outer(psi_k.conj(), self.psi_e1e2)) @ self.ekets
    #         sk_p = sk_m.conj().T

    #         intensity += np.real(np.trace(np.dot(sk_p @ sk_m, self.rho))) * (self.system_params.get_gamma(self.system_params.omega_q1) / self.system_params.gamma_q1) #* (Norm**2)  #
    
    #     return intensity




    # def calculate_steady_state_g2_correlations(self, theta, phi, theta_p, phi_p):
    #     """
    #     Calculate the g2 correlations at steady state.

    #     Parameters:
    #     - H: Hamiltonian of the system.
    #     - rho0: Initial density matrix.
    #     - Gamma_opt: Optical decay rates.
    #     - Gamma_vib: Vibrational decay rates.
    #     - theta: Angle in radians for first polarization.
    #     - phi: Angle in radians for first polarization.
    #     - theta_p: Angle in radians for second polarization.
    #     - phi_p: Angle in radians for second polarization.
    #     - k_q1: Wave vector component.

    #     Returns:
    #     - g2_corr: g2 correlation values.
    #     """
    #     g2_corr = 0.0
    #     # Initialize state
    #     rho_ss = self.rho
    #     rho_qr = np.zeros_like(self.rho, dtype = complex) # empty_dm

    #     # Propagate for a_ops
    #     for _ in self.times_ss:
    #         rho_vec_ss = np.reshape(rho_ss, (self.dim**2, 1))
    #         rho_vec_ss = self.P @ rho_vec_ss
    #         rho_ss = np.reshape(rho_vec_ss, (self.dim, self.dim))

    #     e_1_vec_0, e_2_vec_0 = orthogonal_vectors(theta, phi)
    #     e_vec_0 = np.array([e_1_vec_0, e_2_vec_0])

    #     e_1_vec_1, e_2_vec_1 = orthogonal_vectors(theta_p, phi_p)
    #     e_vec_1 = np.array([e_1_vec_1, e_2_vec_1])


    #     for i in range(2):
    #         d1_pol_0 = self.system_params.d1_hat @ e_vec_0[i]
    #         d2_pol_0 = self.system_params.d2_hat @ e_vec_0[i]


    #         psi_k1_0 = ((d1_pol_0 * np.exp(-1j * self.system_params.k_q1 * self.system_params.r * np.cos(theta) / 2)) * self.psi_e1g2 +
    #                    (d2_pol_0 * np.exp(1j * self.system_params.k_q1 * self.system_params.r * np.cos(theta) / 2)) * self.psi_e2g1) 
            
    #         psi_k2_0 = ((d1_pol_0 * np.exp(-1j * self.system_params.k_q1 * self.system_params.r * np.cos(theta) / 2)) * self.psi_e2g1 +
    #                    (d2_pol_0 * np.exp(1j * self.system_params.k_q1 * self.system_params.r * np.cos(theta) / 2)) * self.psi_e1g2) 

    #         sk_m_0 = self.ekets.conj().T @ (np.outer(self.psi_g1g2, psi_k1_0) + np.outer(psi_k2_0.conj(), self.psi_e1e2)) @ self.ekets
    #         sk_p_0 = sk_m_0.conj().T   

    #         rho_qr_0 = sk_m_0 @ rho_ss @ sk_p_0
    #         rho_qr += rho_qr_0


    #         for j in range(2):
    #             d1_pol_1 = self.system_params.d1_hat @ e_vec_1[j]
    #             d2_pol_1 = self.system_params.d2_hat @ e_vec_1[j]


    #             psi_k1_1 = ((d1_pol_1 * np.exp(-1j * self.system_params.k_q1 * self.system_params.r * np.cos(theta_p) / 2)) * self.psi_e1g2 +
    #                     (d2_pol_1 * np.exp(1j * self.system_params.k_q1 * self.system_params.r * np.cos(theta_p) / 2)) * self.psi_e2g1) 
                
    #             psi_k2_1 = ((d1_pol_1 * np.exp(-1j * self.system_params.k_q1 * self.system_params.r * np.cos(theta_p) / 2)) * self.psi_e2g1 +
    #                    (d2_pol_1 * np.exp(1j * self.system_params.k_q1 * self.system_params.r * np.cos(theta_p) / 2)) * self.psi_e1g2) 

    #             sk_m_1 = self.ekets.conj().T @ (np.outer(self.psi_g1g2, psi_k1_1) + np.outer(psi_k2_1.conj(), self.psi_e1e2)) @ self.ekets
    #             sk_p_1 = sk_m_1.conj().T

                
    #             g2_corr += np.real(np.trace(np.dot(sk_p_1 @ sk_m_1, rho_qr))) * (3 / (8 * np.pi))**2 #* (self.system_params.get_gamma(self.system_params.omega_q1) / self.system_params.gamma_q1)**2
                

    #             # g2 = self.ekets.conj().T @ np.outer(self.psi_e1e2, psi_k_0.conj()) @ np.outer(psi_k_1, psi_k_1.conj()) @ np.outer(psi_k_0, self.psi_e1e2) @ self.ekets 
    #             # g2_corr += np.real(np.trace(np.dot(g2, rho_ss)))* (3 / (8 * np.pi))**2 #* (self.system_params.get_gamma(self.system_params.omega_q1) / self.system_params.gamma_q1)**2
    
    #     return np.array(g2_corr) #/ g2_corr[-1]





# def calculate_intensity(H, rho0, Gamma_coup, Gamma_opt, P_vib):
    
#     sx, sz, sp, sm = pauli_matrices()
#     sx1, sx2, sz1, sz2, sp1, sp2, sm1, sm2 = two_qubit_pauli_matrices()
    
#     # Integrate over times
#     norm_intensity = []
    
#     # Coupling operator
#     a_ops_coup = [sp1 @ sm2, sm1 @ sp2, sp2 @ sm1, sm2 @ sp1]
#     a_ops_opt = [sp1, sm1, sp2, sm2]
    
#     # Bloch-Redfield tensor
#     R_coup, ekets = BR_tensor_coup(H, a_ops_coup, Gamma_coup)

#     R_opt, ekets = BR_tensor_opt(H, a_ops_opt, Gamma_opt, P_vib)

#     # calculate propagator
#     P = expm((R_opt + R_coup) * dt) 

#     # initialize state
#     rho = rho0
#     dim = len(rho[0])
    
#     psi_g1g2 = np.array([0, 0, 0, 1])
#     proj_g1g2 = ekets.conj().T @ np.outer(psi_g1g2, psi_g1g2) @ ekets

#     psi_e1g2 = np.array([0, 0, 1, 0])
#     proj_e1g2 = ekets.conj().T @ np.outer(psi_e1g2, psi_e1g2) @ ekets
    
#     psi_e2g1 = np.array([0, 1, 0, 0])
#     proj_e2g1 = ekets.conj().T @ np.outer(psi_e2g1, psi_e2g1) @ ekets
    
#     proj_e2g1e1g2 = ekets.conj().T @ np.outer(psi_e2g1, psi_e1g2) @ ekets
#     proj_e1g2e2g1 = ekets.conj().T @ np.outer(psi_e1g2, psi_e2g1) @ ekets

#     psi_e1e2 = np.array([1, 0, 0, 0])
#     proj_e1e2 = ekets.conj().T @ np.outer(psi_e1e2, psi_e1e2) @ ekets
    
#     print((proj_e1g2 + proj_e2g1 - proj_e2g1e1g2 - proj_e1g2e2g1)/2)
    
#     # propagate for a_ops
#     for _ in range(steps):
#         int_e1g2 = np.real(np.trace(np.dot(proj_e1g2, rho)))
#         int_e2g1 = np.real(np.trace(np.dot(proj_e2g1, rho)))
#         int_e2g1e1g2 = np.real(np.trace(np.dot(proj_e2g1e1g2, rho))) * (d1_hat @ d2_hat)
#         int_e1g2e2g1 = np.real(np.trace(np.dot(proj_e1g2e2g1, rho))) * (d1_hat @ d2_hat)
#         int_e1e2 = np.real(np.trace(np.dot(proj_e1e2, rho))) + np.real(np.trace(np.dot(proj_e1e2, rho)))
#         norm_intensity.append((int_e1g2 + int_e2g1 + int_e2g1e1g2 + int_e1g2e2g1 + int_e1e2))
#         # norm_intensity.append((int_e1g2 + int_e2g1 - int_e2g1e1g2 - int_e1g2e2g1)/2)
#         rho_vec = np.reshape(rho, (dim**2, 1))
#         rho_vec = P @ rho_vec
#         rho = np.reshape(rho_vec, (dim, dim))

   

#     return norm_intensity





































    # def projectors(self):

    #     '''Gotta make this more streamlined, suitable for single TLS or projections into bright and dark states.'''
    #     psi_g1g2 = self.ekets[:,0]
    #     proj_g1g2 = self.ekets.conj().T @ np.outer(psi_g1g2, psi_g1g2) @ self.ekets

    #     psi_e1g2 = self.ekets[:,2]
    #     proj_e1g2 = self.ekets.conj().T @ np.outer(psi_e1g2, psi_e1g2) @ self.ekets
        
    #     psi_e2g1 = self.ekets[:,1] 
    #     proj_e2g1 = self.ekets.conj().T @ np.outer(psi_e2g1, psi_e2g1) @ self.ekets
        
    #     proj_e2g1e1g2 = np.outer(psi_e2g1, psi_e1g2)
    #     proj_e1g2e2g1 = np.outer(psi_e1g2, psi_e2g1)

    #     psi_e1e2 = self.ekets[:,3]
    #     proj_e1e2 = self.ekets.conj().T @ np.outer(psi_e1e2, psi_e1e2) @ self.ekets

    #     return proj_g1g2, proj_e1g2, proj_e2g1, proj_e2g1e1g2, proj_e1g2e2g1, proj_e1e2





















    # def calculate_steady_state_g2_correlations(self, theta, phi, theta_p, phi_p):
    #     """
    #     Calculate the g2 correlations at steady state.

    #     Parameters:
    #     - H: Hamiltonian of the system.
    #     - rho0: Initial density matrix.
    #     - Gamma_opt: Optical decay rates.
    #     - Gamma_vib: Vibrational decay rates.
    #     - theta: Angle in radians for first polarization.
    #     - phi: Angle in radians for first polarization.
    #     - theta_p: Angle in radians for second polarization.
    #     - phi_p: Angle in radians for second polarization.
    #     - k_q1: Wave vector component.

    #     Returns:
    #     - g2_corr: g2 correlation values.
    #     """
    #     g2_corr = 0 #[]
    #     # Initialize state
    #     rho_ss = self.rho
    #     rho_qr = np.zeros_like(self.rho, dtype = complex) # empty_dm

    #     # Propagate for a_ops
    #     for _ in self.times_ss:
    #         rho_vec_ss = np.reshape(rho_ss, (self.dim**2, 1))
    #         rho_vec_ss = self.P @ rho_vec_ss
    #         rho_ss = np.reshape(rho_vec_ss, (self.dim, self.dim))

    #     e_1_vec_0, e_2_vec_0 = orthogonal_vectors(theta, phi)
    #     e_vec_0 = np.array([e_1_vec_0, e_2_vec_0])

    #     e_1_vec_1, e_2_vec_1 = orthogonal_vectors(theta_p, phi_p)
    #     e_vec_1 = np.array([e_1_vec_1, e_2_vec_1])


    #     for i in range(2):
    #         d1_pol_0 = self.system_params.d_1_vec @ e_vec_0[i]
    #         d2_pol_0 = self.system_params.d_2_vec @ e_vec_0[i]

    #         psi_k_0 = ((d1_pol_0 * np.exp(-1j * self.system_params.k_q1 * self.system_params.r * np.cos(theta) / 2)) * self.psi_e1g2 +
    #                    (d2_pol_0 * np.exp(1j * self.system_params.k_q1 * self.system_params.r * np.cos(theta) / 2)) * self.psi_e2g1) 

    #         sk_m_0 = self.ekets.conj().T @ (np.outer(self.psi_g1g2, psi_k_0) + np.outer(psi_k_0.conj(), self.psi_e1e2)) @ self.ekets
    #         sk_p_0 = sk_m_0.conj().T #self.ekets.conj().T @ (np.outer(psi_k_0.conj(), self.psi_g1g2) + np.outer(self.psi_e1e2, psi_k_0)) @ self.ekets

    #         rho_qr_0 = sk_m_0 @ rho_ss @ sk_p_0
    #         rho_qr += rho_qr_0

    #     # for _ in range(self.steps):
    #     # g2 = 0.0

    #         for j in range(2):
    #             d1_pol_1 = self.system_params.d_1_vec @ e_vec_1[j]
    #             d2_pol_1 = self.system_params.d_2_vec @ e_vec_1[j]

    #             psi_k_1 = ((d1_pol_1 * np.exp(-1j * self.system_params.k_q1 * self.system_params.r * np.cos(theta_p) / 2)) * self.psi_e1g2 +
    #                         (d2_pol_1 * np.exp(1j * self.system_params.k_q1 * self.system_params.r * np.cos(theta_p) / 2)) * self.psi_e2g1)

    #             sk_m_1 = self.ekets.conj().T @ (np.outer(self.psi_g1g2, psi_k_1) + np.outer(psi_k_1.conj(), self.psi_e1e2)) @ self.ekets
    #             sk_p_1 = sk_m_1.conj().T # self.ekets.conj().T @ (np.outer(psi_k_1.conj(), self.psi_g1g2) + np.outer(self.psi_e1e2, psi_k_1)) @ self.ekets

    #             # g2 += np.real(np.trace(np.dot(sk_p_1 @ sk_m_1, rho_qr))) * (self.system_params.get_gamma(self.system_params.omega_q1) / self.system_params.gamma_q1)**2

    #             g2 = self.ekets.conj().T @ np.outer(self.psi_e1e2, psi_k_0.conj()) @ np.outer(psi_k_1, psi_k_1.conj()) @ np.outer(psi_k_0, self.psi_e1e2) @ self.ekets 
    #             g2_corr += np.real(np.trace(np.dot(g2, rho_ss))) * (self.system_params.get_gamma(self.system_params.omega_q1) / self.system_params.gamma_q1)**2

    #         # g2_corr.append(g2)
    #         # rho_vec = np.reshape(rho_qr, (self.dim**2, 1))
    #         # rho_vec = self.P @ rho_vec
    #         # rho_qr = np.reshape(rho_vec, (self.dim, self.dim))
        
    #     return np.array(g2_corr) # / g2_corr[-1]









    # def intensity_distribution(self, theta, phi):
    #     """
    #     Calculate intensity distribution at steady state.

    #     Parameters:
    #     - theta: Angle in radians.
    #     - phi: Angle in radians.

    #     Returns:
    #     - intensity_dist: Intensity distribution value.
    #     """
    #     intensity_dist = []

    #     e_1_vec, e_2_vec = orthogonal_vectors(theta, phi)
    #     e_vec = np.array([e_1_vec, e_2_vec])


    #     intensity = 0.0
    #     # Calculate intensity distribution
    #     for i in range(2):
    #         d1_pol = self.system_params.d1_hat @ e_vec[i]
    #         d2_pol = self.system_params.d2_hat @ e_vec[i]

    #         # if d1_pol == d2_pol == 0.0:
    #         #     psi_k = np.zeros(self.dim, dtype=complex)
    #         # else:
    #         #     Norm = np.sqrt(d1_pol**2 + d2_pol**2)
    #         psi_k = ((d1_pol * np.exp(-1j * self.system_params.k_q1 * self.system_params.r * np.cos(theta) / 2)) * self.psi_e1g2 +
    #                     (d2_pol * np.exp(1j * self.system_params.k_q1 * self.system_params.r * np.cos(theta) / 2)) * self.psi_e2g1) #/ Norm

    #         sk_m = self.ekets.conj().T @ (np.outer(self.psi_g1g2, psi_k.conj()) + np.outer(psi_k, self.psi_e1e2)) @ self.ekets
    #         sk_p = sk_m.conj().T

    #         intensity += np.real(np.trace(np.dot(sk_p @ sk_m, self.rho))) * (3 / (8 * np.pi)) #* (self.system_params.get_gamma(self.system_params.omega_q1) / self.system_params.gamma_q1) #* (Norm**2)  #
    
    #     return intensity


    # def calculate_steady_state_g2_correlations(self, theta, phi, theta_p, phi_p):
    #     """
    #     Calculate the g2 correlations at steady state.

    #     Parameters:
    #     - H: Hamiltonian of the system.
    #     - rho0: Initial density matrix.
    #     - Gamma_opt: Optical decay rates.
    #     - Gamma_vib: Vibrational decay rates.
    #     - theta: Angle in radians for first polarization.
    #     - phi: Angle in radians for first polarization.
    #     - theta_p: Angle in radians for second polarization.
    #     - phi_p: Angle in radians for second polarization.
    #     - k_q1: Wave vector component.

    #     Returns:
    #     - g2_corr: g2 correlation values.
    #     """
    #     g2_corr = []
    #     # Initialize state
    #     rho_ss = self.rho
    #     rho_qr = np.zeros_like(self.rho, dtype = complex) # empty_dm

    #     # Propagate for a_ops
    #     for _ in self.times_ss:
    #         rho_vec_ss = np.reshape(rho_ss, (self.dim**2, 1))
    #         rho_vec_ss = self.P @ rho_vec_ss
    #         rho_ss = np.reshape(rho_vec_ss, (self.dim, self.dim))

    #     e_1_vec_0, e_2_vec_0 = orthogonal_vectors(theta, phi)
    #     e_vec_0 = np.array([e_1_vec_0, e_2_vec_0])

    #     e_1_vec_1, e_2_vec_1 = orthogonal_vectors(theta_p, phi_p)
    #     e_vec_1 = np.array([e_1_vec_1, e_2_vec_1])


    #     for i in range(2):
    #         d1_pol_0 = self.system_params.d1_hat @ e_vec_0[i]
    #         d2_pol_0 = self.system_params.d2_hat @ e_vec_0[i]

    #         if d1_pol_0 == d2_pol_0 == 0.0:
    #             psi_k1_0 = np.zeros(self.dim, dtype=complex)
    #             psi_k2_0 = np.zeros(self.dim, dtype=complex)
    #         else:
    #             Norm_0 = np.sqrt(d1_pol_0**2 + d2_pol_0**2)
    #             psi_k1_0 = ((d1_pol_0 * np.exp(-1j * self.system_params.k_q1 * self.system_params.r * np.cos(theta) / 2)) * self.psi_e1g2 +
    #                     (d2_pol_0 * np.exp(1j * self.system_params.k_q1 * self.system_params.r * np.cos(theta) / 2)) * self.psi_e2g1) #/ Norm_0
                
    #             psi_k2_0 = ((d1_pol_0 * np.exp(1j * self.system_params.k_q1 * self.system_params.r * np.cos(theta) / 2)) * self.psi_e2g1 +
    #                     (d2_pol_0 * np.exp(-1j * self.system_params.k_q1 * self.system_params.r * np.cos(theta) / 2)) * self.psi_e1g2) #/ Norm_0

    #         sk_m_0 = self.ekets.conj().T @ (np.outer(self.psi_g1g2, psi_k1_0.conj()) + np.outer(psi_k2_0, self.psi_e1e2)) @ self.ekets
    #         sk_p_0 = sk_m_0.T   

    #         rho_qr_0 = sk_m_0 @ rho_ss @ sk_p_0
    #         rho_qr += rho_qr_0

    #     for _ in np.linspace(0, 5.0 * self.system_params.tau_L, self.steps):
    #         g2 = 0.0

    #         for j in range(2):
    #             d1_pol_1 = self.system_params.d1_hat @ e_vec_1[j]
    #             d2_pol_1 = self.system_params.d2_hat @ e_vec_1[j]

    #             if d1_pol_1 == d2_pol_1 == 0.0:
    #                 psi_k1_1 = np.zeros(self.dim, dtype=complex)
    #                 psi_k2_1 = np.zeros(self.dim, dtype=complex)
    #             else:
    #                 Norm_1 = np.sqrt(d1_pol_1**2 + d2_pol_1**2)
    #                 psi_k1_1 = ((d1_pol_1 * np.exp(-1j * self.system_params.k_q1 * self.system_params.r * np.cos(theta_p) / 2)) * self.psi_e1g2 +
    #                         (d2_pol_1 * np.exp(1j * self.system_params.k_q1 * self.system_params.r * np.cos(theta_p) / 2)) * self.psi_e2g1) #/ Norm_1
                    
    #                 psi_k2_1 = ((d1_pol_1 * np.exp(1j * self.system_params.k_q1 * self.system_params.r * np.cos(theta_p) / 2)) * self.psi_e2g1 +
    #                     (d2_pol_1 * np.exp(-1j * self.system_params.k_q1 * self.system_params.r * np.cos(theta_p) / 2)) * self.psi_e1g2) #/ Norm_1

    #             sk_m_1 = self.ekets.conj().T @ (np.outer(self.psi_g1g2, psi_k1_1.conj()) + np.outer(psi_k2_1, self.psi_e1e2)) @ self.ekets
    #             sk_p_1 = sk_m_1.T

                
    #             g2 += np.real(np.trace(np.dot(sk_p_1 @ sk_m_1, rho_qr))) * (3 / (8 * np.pi))**2 #* (self.system_params.get_gamma(self.system_params.omega_q1) / self.system_params.gamma_q1)**2
               
    # # #             g2 = self.ekets.conj().T @ np.outer(self.psi_e1e2, psi_k_0.conj()) @ np.outer(psi_k_1, psi_k_1.conj()) @ np.outer(psi_k_0, self.psi_e1e2) @ self.ekets 
    # # #             g2_corr += n.conj()p.real(np.trace(np.dot(g2, rho_ss)))* (3 / (8 * np.pi))**2 #* (self.system_params.get_gamma(self.system_params.omega_q1) / self.system_params.gamma_q1)**2

    #         g2_corr.append(g2)
    #         rho_vec = np.reshape(rho_qr, (self.dim**2, 1))
    #         rho_vec = self.P @ rho_vec
    #         rho_qr = np.reshape(rho_vec, (self.dim, self.dim))
        
    #     return np.array(g2_corr) / g2_corr[-1]
