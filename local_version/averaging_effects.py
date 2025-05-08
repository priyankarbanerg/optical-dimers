from sys_params import *
from system_ops import *
from br_tensor import *
# from g2_correlations import *
from spectral_dens import *
from analysis import *


t0 = time.time()


# Initialize constants and system parameters
constants = Constants()
dipole_params = DipoleParameters(config_name)


# Define parameters for g2 plot
n_theta = 20
n_phi = 20

# Define parameters
theta_vals = np.linspace(0, 0.01 * np.pi, n_theta)
phi_vals = np.linspace(0,  0.01 * 2*np.pi, n_phi)

for theta in theta_vals:
    for phi in phi_vals:
        system_params = SystemParameters(dipole_params=dipole_params, E_q1=E_q1 * constants.eV, E_q2=E_q2 * constants.eV)
        system_params.calculate_parameters()

# Calculate spectral densities and rates
J_opt, J_vib = calculate_spectral_densities()
rates = CombinedRates(T_opt, T_vib, J_opt, J_vib, mode)

Gamma_opt = rates.Opt_rate()
Gamma_vib = rates.Vib_rate()
kappa = rates.kappa()
Gamma_coup = rates.Coup_rate()
P_vib = rates.P_weight()

# Hamiltonian
H = hamiltonian(mode, kappa)

sx1, sx2, sz1, sz2, sp1, sp2, sm1, sm2 = two_qubit_pauli_matrices()

a_ops_vib = [sz1, sz2]
a_ops_coup = [sp1 @ sm2, sm1 @ sp2, sp2 @ sm1, sm2 @ sp1]
a_ops_opt = [sp1, sm1, sp2, sm2]

# # Create a Bloch-Redfield calculator instance
# BR_calculator = BlochRedfieldCalculator(H, a_ops_opt, Gamma_opt, a_ops_vib = None, Gamma_vib = Gamma_vib, a_ops_coup = a_ops_coup, Gamma_coup = Gamma_coup, c_ops_S = None, P_vib = P_vib, mode=mode)

# initial state in the system basis
psi0 = np.array([[1, 0, 0, 0]], dtype = complex)  # np.array([[0, 1, 1, 0] / np.sqrt(2)], dtype = complex) # Dimer

# initial state as density operator
rho_sz = psi2rho(psi0)

H_diag, rho0 = initial_state(H, psi0)

# Define time points and steps for the simulation
tf = 10 * system_params.tau_L
steps = 200 
times = np.linspace(0, tf, steps)
times_ns = np.linspace(0, 1.0e-3 * tf, steps)

# # Initialize the QuantumSystemAnalysis class for calculating intensity
# quantum_analysis_intensity = QuantumSystemAnalysis(BR_calculator, quantity='intensity', rho=rho0, final_time=tf, steps=steps)

# # Compute intensity
# intensity_results = quantum_analysis_intensity.compute()





# # Initialize the QuantumSystemAnalysis class for calculating population
# state = 'bright'
# quantum_analysis_population = QuantumSystemAnalysis(BR_calculator, quantity='population', rho=rho0, final_time=tf, steps=steps, state=state)

# Compute population
# pop_results = quantum_analysis_population.compute()



# Define parameters for g2 plot
n_theta_p = 100
n_phi_p = 100

# Define parameters
theta_vals_p = np.linspace(0, np.pi, n_theta_p)
phi_vals_p = np.linspace(0,  2*np.pi, n_phi_p)
r_p = 2  # radius of sphere

# Generate x, y, z coordinates of points on the sphere
x = r_p * np.outer(np.cos(phi_vals_p), np.sin(theta_vals_p))
y = r_p * np.outer(np.sin(phi_vals_p), np.sin(theta_vals_p))
z = r_p * np.outer(np.ones(np.size(phi_vals_p)), np.cos(theta_vals_p))


# # Initialize the DistributionAnalysis class for calculating intensity
# quantum_analysis_intensity_dist = DistributionAnalysis(BR_calculator, quantity='intensity distribution', rho=rho0, final_time=tf, steps=steps)

# # Compute intensity distribution
# norm_int = np.zeros((n_phi_p, n_theta_p, steps))
# for l, phi_p in enumerate(phi_vals_p):
#     for k, theta_p in enumerate(theta_vals_p):
#         norm_int[l, k, :] = quantum_analysis_intensity_dist.compute(theta_p, phi_p)





perp_H = [0, 0, 1]
r, theta, phi = cart2sph(perp_H) #cart2sph(d2_hat) 
# r_p, theta_p, phi_p = r, theta, phi # cart2sph(d2_hat) # 

# perp_J = [1, 0, 0]
# r, theta, phi = cart2sph(perp_J)
# r_p, theta_p, phi_p = r, theta, phi

# perp_ortho = [0.0, 0, 1.0]
# r, theta, phi = cart2sph(perp_ortho)

# perp_oppo = [0, 0, 1]
# r, theta, phi = cart2sph(perp_oppo)

# collapse operators (for initial incoherent pumping)
c_ops_S = [np.sqrt(system_params.gamma_p_1) * sp1, np.sqrt(system_params.gamma_p_2) * sp2]
# Create a Bloch-Redfield calculator instance
BR_calculator = BlochRedfieldCalculator(H, a_ops_opt, Gamma_opt, a_ops_vib = a_ops_vib, Gamma_vib = Gamma_vib, a_ops_coup = a_ops_coup, Gamma_coup = Gamma_coup, c_ops_S = c_ops_S, P_vib = P_vib, mode=mode)

# # Initialize the DistributionAnalysis class for calculating g2
# quantum_analysis_g2_dist = DistributionAnalysis(BR_calculator, quantity='g2 distribution', rho=rho0, final_time=tf, steps=steps)

# # Compute g2 distribution
# g2_corr = np.zeros((n_phi_p, n_theta_p, steps))
# for l, phi_p in enumerate(phi_vals_p):
#     for k, theta_p in enumerate(theta_vals_p):
#         g2_corr[l, k, :] = quantum_analysis_g2_dist.compute(theta, phi, theta_p, phi_p)

# # Normalize the g2 values
# g2_corr = g2_corr #/ g2_corr.max()



'''Averaging over quantities'''

#######################################################
'''Averaging over k, k_prime (same) == The tumbling effect is slower that detection rate'''
#######################################################

# Initialize the DistributionAnalysis class for calculating g2
quantum_analysis_g2_dist = DistributionAnalysis(BR_calculator, quantity='g2 distribution', rho=rho0, final_time=tf, steps=steps)

# Compute intensity distribution
g2_corr = np.zeros((n_phi_p, n_theta_p, steps))
for l, phi_p in enumerate(phi_vals_p):
    for k, theta_p in enumerate(theta_vals_p):
        theta_p, phi_p = theta, phi 
        g2_corr[l, k, :] = quantum_analysis_g2_dist.compute(theta, phi, theta_p, phi_p)

g2_corr = np.mean(g2_corr, axis=(0, 1))



########################################################################
'''Averaging over relative orientation, phi or phi_prime for dipoles'''
########################################################################

# Initialize the DistributionAnalysis class for calculating g2
quantum_analysis_g2_dist = DistributionAnalysis(BR_calculator, quantity='g2 distribution', rho=rho0, final_time=tf, steps=steps)

# Compute intensity distribution
g2_corr = np.zeros((n_phi_p, n_theta_p, steps))
for l, phi_p in enumerate(phi_vals_p):
    for k, theta_p in enumerate(theta_vals_p):
        theta_p, phi_p = theta, phi 
        g2_corr[l, k, :] = quantum_analysis_g2_dist.compute(theta, phi, theta_p, phi_p)

g2_corr = np.mean(g2_corr, axis=(0, 1))






# sx1, sx2, sz1, sz2, sp1, sp2, sm1, sm2 = two_qubit_pauli_matrices()

# # initial state in the system basis
# psi0 = np.array([[1, 0, 0, 0]], dtype = complex)  # Dimer

# # initial state as density operator
# rho_sz = psi2rho(psi0)

# rho0 = initial_state(H, psi0)

# # Define spectral density functions G_1(w) and G_2(w)
# Gamma_coup = Coup_rate(kappa)
# Gamma_opt = Opt_rate(T_opt)
# P_vib = P_weight(kappa)



# # Define parameters for intensity plot
# n_theta = 10
# n_phi = 10
# n_theta_p = 10
# n_phi_p = 10

# # Define parameters
# theta_vals = np.linspace(0, np.pi, n_theta)
# phi_vals = np.linspace(0,  2*np.pi, n_phi)
# theta_vals_p = np.linspace(0, np.pi, n_theta)
# phi_vals_p = np.linspace(0,  2*np.pi, n_phi)
# r = 1  # radius of sphere

# # Generate x, y, z coordinates of points on the sphere
# x = r * np.outer(np.cos(phi_vals), np.sin(theta_vals))
# y = r * np.outer(np.sin(phi_vals), np.sin(theta_vals))
# z = r * np.outer(np.ones(np.size(phi_vals)), np.cos(theta_vals))

# # Calculate normalized intensity for each point on the sphere
# rho_qr = np.zeros((n_phi_p, n_theta_p))

# for j in range(n_phi):
#     for i in range(n_theta):
#         rho_qr[j, i] = calculate_rho_qr(H, rho0, Gamma_coup, Gamma_opt, P_vib, theta_vals[i], phi_vals[j], k_q1)

# print(rho_qr)