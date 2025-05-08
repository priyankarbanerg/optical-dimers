
from sys_params import *
from system_ops import *
from br_tensor import *
from analysis import *
from matplotlib.ticker import MaxNLocator
from vonmises_fisher import *
# from propagate import *


t0 = time.time()


sx1, sx2, sz1, sz2, sp1, sp2, sm1, sm2 = two_qubit_pauli_matrices()

a_ops_vib = [sz1, sz2]
a_ops_coup = [sp1 @ sm2, sm1 @ sp2]
a_ops_opt = [sp1, sm1, sp2, sm2]

# # initial state in the system basis
psi0 = np.array([[1, 0, 0, 0]], dtype = complex)  # 
# psi0 = np.array([[0, 1, 1, 0] / np.sqrt(2)], dtype = complex) # Dimer

# initial state as density operator
rho_sz = psi2rho(psi0)



class G2CorrelationSimulator:
    def __init__(self, constants, dipole_params, system_params, rates, config_name="H", mode="polaron", pump="sym", faulty_detector="yes", tumbling="slow", jiggling="yes", kappa_detec=0, kappa_0=10, 
                 num_perp_samples=10000, steps=10000, instrument_response=50, 
                 save_path='files/g2_con/', g2_save_path='files/g2/'):
        self.config_name = config_name
        self.mode = mode
        self.pump = pump
        self.faulty_detector = faulty_detector
        self.tumbling = tumbling
        self.jiggling = jiggling
        self.kappa_detec = kappa_detec
        self.kappa_0 = kappa_0
        self.num_perp_samples = num_perp_samples
        self.steps = steps
        self.sigma = instrument_response
        self.save_path = save_path
        self.g2_save_path = g2_save_path

        # Make sure directories exist
        os.makedirs(save_path, exist_ok=True)
        os.makedirs(g2_save_path, exist_ok=True)

        self.constants = constants
        self.dipole_params = dipole_params
        self.system_params = system_params
        self.rates = rates
        self.mean_direction = np.array([0.0, 1.0, 0.0]) if config_name in ['H', 'J'] else np.array([0.0, 0.0, 1.0])

        # Define the mean direction for the von Mises-Fisher distribution
        if self.config_name == 'H':
            self.mu = np.array([1.0, 0.0, 0.0])
        elif self.config_name == 'J':
            self.mu = np.array([0.0, 0.0, 1.0])
        elif self.config_name == 'ortho':
            self.mu = np.array([0.0, 1.0, 0.0])

    def generate_perp_samples(self):
        if self.tumbling == "fast":
            return rand_von_mises_fisher(mu=self.mean_direction, kappa=self.kappa_detec, N=self.num_perp_samples)
        else:
            return np.full((1, 3), self.mean_direction)

    def generate_orient_samples(self):
        if self.jiggling == "yes":
            # Sample random orientations (new_dir) using von Mises-Fisher
            return rand_von_mises_fisher(self.mu, kappa=self.kappa_0, N=10000)
        else:
            return np.full((1, 3), self.mu)
        
    def compute_g2_for_orientation_disorder(self, new_dir, perp):
        
        # Initialize constants and system parameters as before
        constants = Constants()
        dipole_params = DipoleParameters(self.config_name, new_dir=new_dir, average=True)
        system_params = SystemParameters(dipole_params=dipole_params, E_q1=E_q1 * constants.eV, E_q2=E_q2 * constants.eV)
        system_params.calculate_parameters()
        # perp = self.mean_direction

        # Calculate rates and Hamiltonian
        Gamma_opt = self.rates.Opt_rate()
        Gamma_vib = self.rates.Vib_rate()
        kappa = self.rates.kappa()
        Gamma_coup = self.rates.Coup_rate()
        P_vib = self.rates.P_weight()
        
        # Define the Hamiltonian and initial state
        ham_calc = Hamiltonian(system_params, self.rates, self.mode, kappa)
        H = ham_calc.get_hamiltonian()
        H_diag, rho0 = initial_state(H, psi0)
        
        # Convert cartesian to spherical for both directions
        _, theta_k, phi_k = cart2sph(self.mean_direction)
        _, theta_k_prime, phi_k_prime = cart2sph(perp)
        
        # Create Bloch-Redfield calculator instance
        BR_calculator = BlochRedfieldCalculator(H, constants, dipole_params, system_params, 
                                                a_ops_opt, Gamma_opt, self.rates, 
                                                a_ops_vib=a_ops_vib, Gamma_vib=Gamma_vib, 
                                                a_ops_coup=a_ops_coup, Gamma_coup=Gamma_coup, 
                                                P_vib=P_vib, mode=self.mode, pump=self.pump)
        
        # Compute g2 correlation with paired directions (theta_k, phi_k) and (theta_k_prime, phi_k_prime)
        quantum_analysis_g2_dist = DistributionAnalysis(BR_calculator, constants, dipole_params, system_params, self.rates, 
                                                        quantity='g2 distribution', rho=rho0, 
                                                        final_time=5 * system_params.tau_L, steps=self.steps)
        
        return quantum_analysis_g2_dist.compute(theta_k, phi_k, theta_k_prime, phi_k_prime)


    def run_simulation(self):
        """
        Simulates g2 correlation for multiple orientations and perturbations.
        """
        # Generate `perp` samples based on tumbling type
        perp_samples = self.generate_perp_samples()
        points = self.generate_orient_samples()
        # Initialize a list to store futures
        futures = []
        
        # # Parallel execution with ProcessPoolExecutor
        # with ProcessPoolExecutor(max_workers=20) as executor:
        #     # Submit tasks for all combinations of `points` and `perp_samples`
        #     for perp in perp_samples:
        #         for new_dir in points:
        #             futures.append(
        #                 executor.submit(self.compute_g2_for_orientation_disorder, new_dir, perp)
        #             )

        num_pairings = 10000

        # Generate random pairings
        indices = np.random.randint(0, len(perp_samples), size=(num_pairings, 2))

        # Initialize a list to store futures
        futures = []

        # Parallel execution with ProcessPoolExecutor
        with ProcessPoolExecutor(max_workers=5) as executor:
            # Submit tasks for 10,000 random pairings
            for i in range(num_pairings):
                random_perp = perp_samples[indices[i, 0]]
                random_point = points[indices[i, 1]]
                
                futures.append(
                    executor.submit(self.compute_g2_for_orientation_disorder, random_point, random_perp)
                )
            

            # Collect results as they complete
            g2_results = []
            for future in as_completed(futures):
                g2_results.append(future.result())

        # # Reshape results into a 2D array (len(perp_samples), len(points)) for further averaging
        # g2_results_array = np.array(g2_results).reshape(len(perp_samples), len(points), -1)
        g2_results_array = np.array(g2_results).reshape(num_pairings, -1)

        # print(np.shape(g2_results_array))

        # Average along specific axes
        # g2_corr = np.mean(g2_results_array, axis=(0, 1))
        g2_corr = np.mean(g2_results_array, axis=0)
        self.save_g2_corr(g2_corr, normalized=False)

        # Normalize g2 correlation
        normalized_g2_corr = g2_corr / g2_corr[-1]
        self.save_g2_corr(normalized_g2_corr, normalized=True)

        # plt.plot(g2_corr)
        # plt.show()

        # Apply Gaussian convolution for the final result
        g2_corr_appended, times_appended = self.apply_gaussian_convolution(normalized_g2_corr)

        # Plot and save results
        self.save_and_plot(g2_corr_appended, times_appended)

        return normalized_g2_corr
        
    def apply_gaussian_convolution(self, g2_corr):
        tf = 5 * self.system_params.tau_L
        times = np.linspace(0, tf, self.steps)
        times_ns = np.linspace(0, 1.0e-3 * tf, self.steps)

        dt = (5 * self.system_params.tau_L - (-5 * self.system_params.tau_L)) / self.steps
        t_kernel = np.arange(-10 * self.sigma, 10 * self.sigma, dt)
        gaussian_kernel = np.exp(-t_kernel**2 / (2 * self.sigma**2))
        gaussian_kernel /= np.sum(gaussian_kernel)
        
        g2_corr_appended = fftconvolve(np.append(np.flip(g2_corr[1:]), g2_corr), gaussian_kernel, mode="same")
        times_appended = np.append(-np.flip((times_ns[1:] - times_ns[0])), times_ns)
        
        return g2_corr_appended, times_appended

    def save_g2_corr(self, g2_corr, normalized):
        # Generate filename suffix based on the specified format
        suffix = "normalized" if normalized else "unnormalized"
        mode_str = "weakC" if self.mode == "weak" else "strongC"
        config_name_str = self.config_name  # "H" or "J"
        tumbling_str = "slow_tumb" if self.tumbling == "slow" else "fast_tumb"
        jiggling_str = "jigg" if self.jiggling == "yes" else "nojigg"
        
        # New filename format
        file_name = f"{config_name_str}_{mode_str}_2nm_symm_{tumbling_str}_{jiggling_str}_{suffix}.npy"
        file_path = os.path.join(self.g2_save_path, file_name)
        np.save(file_path, g2_corr)

    def save_and_plot(self, g2_corr_appended, times_appended):
        # Generate filename based on the specified format
        mode_str = "weakC" if self.mode == "weak" else "strongC"
        config_name_str = self.config_name  # "H" or "J"
        tumbling_str = "slow_tumb" if self.tumbling == "slow" else "fast_tumb"
        jiggling_str = "jigg" if self.jiggling == "yes" else "nojigg"

        # New filename format
        file_name = f"{config_name_str}_{mode_str}_2nm_symm_{tumbling_str}__{jiggling_str}_{self.sigma}.npy"
        file_path = os.path.join(self.save_path, file_name)
        np.save(file_path, g2_corr_appended)
        
        # # Plotting the results
        # plt.plot(times_appended, g2_corr_appended)
        # plt.ylim([0.0, 1.25])
        # plt.grid(True)
        # plt.show()

    def run(self):
        g2_corr = self.run_simulation()
        g2_corr_appended, times_appended = self.apply_gaussian_convolution(g2_corr)
        self.save_and_plot(g2_corr_appended, times_appended)

configurations = [
    {"config_name": "H", "tumbling": "slow", "jiggling": "no", "mode": "polaron"},
    {"config_name": "H", "tumbling": "fast", "jiggling": "no", "mode": "polaron"},
    {"config_name": "H", "tumbling": "fast", "jiggling": "yes", "mode": "polaron"},
    {"config_name": "J", "tumbling": "slow", "jiggling": "no", "mode": "polaron"},
    {"config_name": "J", "tumbling": "fast", "jiggling": "no", "mode": "polaron"},
    {"config_name": "J", "tumbling": "fast", "jiggling": "yes", "mode": "polaron"}
]

class G2InstrumentResponseSimulator:
    def __init__(self, configurations, instrument_responses=np.linspace(1, 250, 100, dtype=int), save_path='files/instrument_response'):
        self.configurations = configurations
        self.instrument_responses = instrument_responses
        self.save_path = save_path

        # Ensure the save directory exists
        os.makedirs(save_path, exist_ok=True)

    def run_simulation_for_responses(self):
        T_opt = 1.0 * 5800.0              # Temperature
        T_vib = 1.0 * 300.0

        for config in self.configurations:
            # Generate and store g2_corr only once per configuration
            constants, dipole_params, system_params = self.setup_system_params(config)
            rates = CombinedRates(constants, dipole_params, system_params, T_opt, T_vib, mode="polaron")
            simulator = G2CorrelationSimulator(constants=constants, dipole_params=dipole_params,
                                               system_params=system_params, rates=rates, **config)
            
            # Generate the base g2 correlation without instrument response
            g2_corr = simulator.run_simulation()
            
            # plt.plot(g2_corr)
            # plt.show()

            # Store results for zero-delay convolved g2 values for the current configuration
            zero_delay_convolved_array = []
            
            for instrument_response in self.instrument_responses:
                # Set the simulator's instrument response and apply convolution
                simulator.sigma = int(instrument_response)
                g2_corr_appended, _ = simulator.apply_gaussian_convolution(g2_corr)

                # plt.plot(g2_corr_appended)
                # plt.show()

                # Append the zero-delay value for this convolved g2
                zero_delay_convolved_array.append(g2_corr_appended[simulator.steps])  # Adjusted for zero-delay index
                
            # Save the zero-delay convolved array for each configuration
            self.save_zero_delay_array(zero_delay_convolved_array, config)

            plt.plot(self.instrument_responses, zero_delay_convolved_array)
            plt.show()

    def setup_system_params(self, config):
        """ Set up and return constants, dipole_params, and system_params for each configuration. """
        constants = Constants()
        if config["config_name"] == "H":
            dipole_params = DipoleParameters(config_name="H", new_dir=None, average=False)
        elif config["config_name"] == "J":
            dipole_params = DipoleParameters(config_name="J", new_dir=None, average=False)
        system_params = SystemParameters(dipole_params=dipole_params, E_q1=E_q1 * constants.eV, E_q2=E_q2 * constants.eV)
        system_params.calculate_parameters()
        
        return constants, dipole_params, system_params

    def save_zero_delay_array(self, zero_delay_convolved_array, config):
        """Save zero-delay values for each instrument response as a .npy file."""
        zero_delay_convolved_array = np.array(zero_delay_convolved_array)
        
        # Generate filename based on configuration details
        mode_str = "weakC" if config.get("mode", "weak") == "weak" else "strongC"
        config_name_str = config["config_name"]
        tumbling_str = "slow_tumb" if config["tumbling"] == "slow" else "fast_tumb"
        jiggling_str = "jigg" if config["jiggling"] == "yes" else "nojigg"
        
        # Save array for the current configuration
        file_name = f"{config_name_str}_{mode_str}_2nm_symm_{tumbling_str}_{jiggling_str}_convolved.npy"
        file_path = os.path.join(self.save_path, file_name)
        np.save(file_path, zero_delay_convolved_array)

# Instantiate and run the instrument response simulator
instrument_simulator = G2InstrumentResponseSimulator(configurations=configurations)
instrument_simulator.run_simulation_for_responses()


