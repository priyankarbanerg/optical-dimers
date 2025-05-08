
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

def compute_g2(T_vib, config_name):
    constants = Constants()
    dipole_params = DipoleParameters(config_name, new_dir=None, average=False)
    system_params = SystemParameters(dipole_params=dipole_params, E_q1=E_q1 * constants.eV, E_q2=E_q2 * constants.eV)
    system_params.calculate_parameters()
    
    mean_direction = np.array([0.0, 1.0, 0.0]) if config_name in ['H', 'J'] else np.array([0.0, 0.0, 1.0])
    
    rates = CombinedRates(constants, dipole_params, system_params, T_opt, T_vib, mode=mode)
    Gamma_opt = rates.Opt_rate()
    Gamma_vib = rates.Vib_rate()
    kappa = rates.kappa()
    Gamma_coup = rates.Coup_rate()
    P_vib = rates.P_weight()

    ham_calc = Hamiltonian(system_params, rates, mode, kappa)
    H = ham_calc.get_hamiltonian()
    H_diag, rho0 = initial_state(H, psi0)

    _, theta_k, phi_k = cart2sph(mean_direction)
    _, theta_k_prime, phi_k_prime = cart2sph(mean_direction)

    BR_calculator = BlochRedfieldCalculator(
        H, constants, dipole_params, system_params, 
        a_ops_opt, Gamma_opt, rates, 
        a_ops_vib=a_ops_vib, Gamma_vib=Gamma_vib, 
        a_ops_coup=a_ops_coup, Gamma_coup=Gamma_coup, 
        P_vib=P_vib, mode=mode, pump=pump
    )

    quantum_analysis_g2_dist = DistributionAnalysis(
        BR_calculator, constants, dipole_params, system_params, rates, 
        quantity='g2 distribution', rho=rho0, 
        final_time=5 * system_params.tau_L, steps=steps
    )

    g2_corr = quantum_analysis_g2_dist.compute(theta_k, phi_k, theta_k_prime, phi_k_prime)
    g2_corr = g2_corr / g2_corr[-1]
    return g2_corr[0]

# Function to compute zero-delay g2 value for a given T_vib
def compute_g2_wrapper(args):
    """Wrapper function to unpack arguments since ProcessPoolExecutor can't use lambdas."""
    T_vib, config_name = args
    return compute_g2(T_vib, config_name)

configurations = ['H', 'J']
T_opt = 5800
steps = 10000
Temp_vib = np.linspace(1, 300, 500)
mode = "polaron"
pump = "sym"

start_time = time.time()

for config_name in configurations:
    print(f"Running for config: {config_name}")
    
    # Create argument tuples for parallel execution
    args_list = [(T, config_name) for T in Temp_vib]

    with ProcessPoolExecutor(max_workers=20) as executor:
        zero_delay_array = list(executor.map(compute_g2_wrapper, args_list))

    # Save results
    output_file = f'g2_results_{config_name}.npy'
    np.save(output_file, zero_delay_array)
    print(f"Results saved to {output_file}")

end_time = time.time()
print(f"Total execution time: {end_time - start_time} seconds")