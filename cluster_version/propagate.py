
from sys_params import *
from system_ops import *
from br_tensor import *
from analysis import *
from matplotlib.ticker import MaxNLocator
from vonmises_fisher import *

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


config_name = "J"
mode = 'polaron'
pump = 'sym'
faulty_detector = 'no'
tumbling = 'slow'


'''Averaging over quantities'''

#######################################################
'''Averaging over k, k_prime (same) == The tumbling effect is slower than detection rate'''
#######################################################


# Constants
kappa_detec = 0
kappa_0 = 0
# num_dir = 1000
num_perp_samples = 100  # Number of samples for varying `perp`
steps = 1000
theta_max = np.pi / 6  # 7.5 degrees

# Define the configuration's mean direction
if config_name == 'H':
    mean_direction = np.array([0.0, 1.0, 0.0])
elif config_name == 'J':
    mean_direction = np.array([0.0, 1.0, 0.0])
elif config_name == 'ortho':
    mean_direction = np.array([0.0, 0.0, 1.0])

# # Sample multiple `perp` directions around the mean direction
perp_samples = rand_von_mises_fisher(mu=mean_direction, kappa=kappa_0, N=num_perp_samples) # vonmises_fisher.rvs(mu=mean_direction, kappa=kappa_0, size=num_perp_samples)

def sample_direction_pairs(perp, theta_max, num_samples):
    """
    Generate `num_samples` pairs of directions within a cone of half-angle `theta_max` around `perp`.
    """
    if faulty_detector == 'yes' and tumbling == 'slow':
        points_k = rand_von_mises_fisher(mu=mean_direction, kappa=kappa_detec, N=num_perp_samples)
        points_k_prime = rand_von_mises_fisher(mu=mean_direction, kappa=kappa_detec, N=num_perp_samples)

        angles_k = np.arccos(np.clip(np.dot(points_k, mean_direction), -1.0, 1.0))
        angles_k_prime = np.arccos(np.clip(np.dot(points_k_prime, mean_direction), -1.0, 1.0))
        points_k = points_k[angles_k <= theta_max]
        points_k_prime = points_k_prime[angles_k_prime <= theta_max]

    elif faulty_detector == 'yes' and tumbling == 'fast':
        points_k = rand_von_mises_fisher(mu=mean_direction, kappa=kappa_detec, N=num_perp_samples)
        points_k_prime = rand_von_mises_fisher(mu=perp, kappa=kappa_detec, N=num_perp_samples)

        angles_k = np.arccos(np.clip(np.dot(points_k, mean_direction), -1.0, 1.0))
        angles_k_prime = np.arccos(np.clip(np.dot(points_k_prime, perp), -1.0, 1.0))
        points_k = points_k[angles_k <= theta_max]
        points_k_prime = points_k_prime[angles_k_prime <= theta_max]

    elif faulty_detector == 'no' and tumbling == 'slow':
        points_k = np.full((1, 3), mean_direction)
        points_k_prime = np.full((1, 3), mean_direction)

        angles_k = np.arccos(np.clip(np.dot(points_k, mean_direction), -1.0, 1.0))
        angles_k_prime = np.arccos(np.clip(np.dot(points_k_prime, mean_direction), -1.0, 1.0))
        points_k = points_k[angles_k <= theta_max]
        points_k_prime = points_k_prime[angles_k_prime <= theta_max]

    elif faulty_detector == 'no' and tumbling == 'fast':
        points_k = np.full((1, 3), mean_direction)
        points_k_prime = np.full((1, 3), perp)

        angles_k = np.arccos(np.clip(np.dot(points_k, mean_direction), -1.0, 1.0))
        angles_k_prime = np.arccos(np.clip(np.dot(points_k_prime, perp), -1.0, 1.0))
        points_k = points_k[angles_k <= theta_max]
        points_k_prime = points_k_prime[angles_k_prime <= theta_max]

    
    # Check for empty arrays after filtering
    if len(points_k) == 0 or len(points_k_prime) == 0:
        return []  # Return an empty list if there are no valid points

    return list(zip(points_k[:num_samples], points_k_prime[:num_samples]))  # Return paired points

def compute_g2_for_pair(pair):
    samp_dir_k, samp_dir_k_prime = pair
    # Initialize constants and system parameters as before
    constants = Constants()
    dipole_params = DipoleParameters(config_name, new_dir=None, average=False)
    system_params = SystemParameters(dipole_params=dipole_params, E_q1=E_q1 * constants.eV, E_q2=E_q2 * constants.eV)
    system_params.calculate_parameters()
    
    # Calculate rates and Hamiltonian
    rates = CombinedRates(constants, dipole_params, system_params, T_opt, T_vib, mode)
    Gamma_opt = rates.Opt_rate()
    Gamma_vib = rates.Vib_rate()
    kappa = rates.kappa()
    Gamma_coup = rates.Coup_rate()
    P_vib = rates.P_weight()
    
    # Define the Hamiltonian and initial state
    ham_calc = Hamiltonian(system_params, rates, mode, kappa)
    H = ham_calc.get_hamiltonian()
    H_diag, rho0 = initial_state(H, psi0)
    
    # Convert cartesian to spherical for both directions
    _, theta_k, phi_k = cart2sph(samp_dir_k)
    _, theta_k_prime, phi_k_prime = cart2sph(samp_dir_k_prime)
    
    # Create Bloch-Redfield calculator instance
    BR_calculator = BlochRedfieldCalculator(H, constants, dipole_params, system_params, 
                                            a_ops_opt, Gamma_opt, 
                                            a_ops_vib=a_ops_vib, Gamma_vib=Gamma_vib, 
                                            a_ops_coup=a_ops_coup, Gamma_coup=Gamma_coup, 
                                            P_vib=P_vib, mode=mode, pump=pump)
    
    # Compute g2 correlation with paired directions (theta_k, phi_k) and (theta_k_prime, phi_k_prime)
    quantum_analysis_g2_dist = DistributionAnalysis(BR_calculator, constants, dipole_params, system_params, rates, 
                                                    quantity='g2 distribution', rho=rho0, 
                                                    final_time=5 * system_params.tau_L, steps=steps)
    
    return quantum_analysis_g2_dist.compute(theta_k, phi_k, theta_k_prime, phi_k_prime)

def compute_g2_for_perp(perp):
    # Generate direction pairs for the given `perp`
    direction_pairs = sample_direction_pairs(perp, theta_max, num_perp_samples)
    with ProcessPoolExecutor(max_workers=4) as executor:
        results = list(executor.map(compute_g2_for_pair, direction_pairs))  # Using compute_g2_for_pair here
    return np.mean(results, axis=0)  # Average g2 results for the sampled `perp`


# Compute g2 for all sampled `perp` orientations in parallel
with ProcessPoolExecutor(max_workers=20) as executor:
    if tumbling == 'fast':
        g2_results = list(executor.map(compute_g2_for_perp, perp_samples))
    elif tumbling == 'slow':
        perp_samples = np.full((1, 3), mean_direction) 
        g2_results = list(executor.map(compute_g2_for_perp, perp_samples))


# Average g2 results across all `perp` orientations
g2_corr = np.mean(g2_results, axis=0)
g2_corr = g2_corr / g2_corr[-1]  # Normalize

# Plotting the results
# Initialize constants and system parameters
constants = Constants()
dipole_params = DipoleParameters(config_name, new_dir=None, average=False)
system_params = SystemParameters(dipole_params=dipole_params, E_q1=E_q1 * constants.eV, E_q2=E_q2 * constants.eV)
system_params.calculate_parameters()

tf = 5 * system_params.tau_L
times = np.linspace(0, tf, steps)
times_ns = np.linspace(0, 1.0e-3 * tf, steps)


# np.save('files/g2/crystal/H_dimer/H_strongC_2nm_symm_fd.npy', g2_corr)  # nfd = not faulty detector, fd = faulty detector

# g2_corr = np.load('files/g2/crystal/H_dimer/H_strongC_2nm_symm_fd.npy')

# np.save('files/g2/J_strongC_2nm_symm_2.npy', g2_corr)

# sigma = 50  # in picoseconds

# dt = (tf - (-tf)) / steps
# t_kernel = np.arange(-10 * sigma, 10 * sigma, dt)
# gaussian_kernel = np.exp(-t_kernel**2 / (2 * sigma**2))
# gaussian_kernel /= np.sum(gaussian_kernel)  

# g2_corr_appended = fftconvolve(np.append(np.flip(g2_corr[1:]), g2_corr), gaussian_kernel, mode = "same") #fftconvolve(np.append(np.flip(np.load('files/g2/H_dimer/H_strongC_T0.npy')[1:]), np.load('files/g2/H_dimer/H_strongC_T0.npy')), gaussian_kernel, mode = "same")
# times_appended = np.append(-np.flip((times_ns[1:] - times_ns[0])), times_ns)

# # np.save('files/g2_con/crystal/H_dimer/H_strongC_2nm_symm_fd.npy', g2_corr_appended)         #pbs = pumping bright state

# Plot the g2 correlation
plt.plot(np.append(-np.flip((times_ns[1:] - times_ns[0])), times_ns), np.append(np.flip(g2_corr[1:]), g2_corr))
# plt.plot(np.append(-np.flip((times_ns[1:] - times_ns[0])), times_ns), g2_corr_appended)
# plt.ylim([0.0, 1.25])
plt.grid(True)
plt.show()





# #######################################################
# '''Averaging over mutual orientation == The jiggling effect'''
# #######################################################


# kappa = 10
# num_dir = 1000
# steps = 100000


# def compute_g2(new_dir):
    
#     # Initialize constants and system parameters
#     constants = Constants()
    
#     # Pass the tuple (theta, phi) as disorder angles
#     dipole_params = DipoleParameters(config_name, new_dir=new_dir, average=True)
    
#     system_params = SystemParameters(dipole_params=dipole_params, E_q1=E_q1 * constants.eV, E_q2=E_q2 * constants.eV)
#     system_params.calculate_parameters()

#     # Calculate spectral densities and rates
#     rates = CombinedRates(constants, dipole_params, system_params, T_opt, T_vib, mode)

#     Gamma_opt = rates.Opt_rate()
#     Gamma_vib = rates.Vib_rate()
#     kappa = rates.kappa()
#     Gamma_coup = rates.Coup_rate()
#     P_vib = rates.P_weight()

#     # Hamiltonian
#     ham_calc = Hamiltonian(system_params, rates, mode, kappa)
#     H = ham_calc.get_hamiltonian()

#     H_diag, rho0 = initial_state(H, psi0)

#     # Define time points and steps for the simulation
#     tf = 5 * system_params.tau_L

#     # Convert cartesian to spherical coordinates for initial angles
#     if config_name == 'H':
#         perp = [0.0, 1.0, 0.0]
#         # perp = [0.0, 0.0, 1.0]
#     elif config_name == 'J':
#         perp = [0.0, 1.0, 0.0]
#     elif config_name == 'ortho':
#         perp = [0.0, 0.0, 1.0]
#         # perp = [0.5, 0.5, 0.0]
#     r, theta_init, phi_init = cart2sph(perp)
    
#     r_p, theta_p, phi_p = r, theta_init, phi_init

#     # Collapse operators (for initial incoherent pumping)
#     c_ops_S = [np.sqrt(system_params.gamma_p_1) * sp1, np.sqrt(system_params.gamma_p_2) * sp2]

#     # Create a Bloch-Redfield calculator instance
#     BR_calculator = BlochRedfieldCalculator(H, constants, dipole_params, system_params, 
#                                             a_ops_opt, Gamma_opt, 
#                                             a_ops_vib=a_ops_vib, Gamma_vib=Gamma_vib, 
#                                             a_ops_coup=a_ops_coup, Gamma_coup=Gamma_coup, 
#                                             P_vib=P_vib, mode=mode, pump = pump)

#     # Initialize the DistributionAnalysis class for calculating g2
#     quantum_analysis_g2_dist = DistributionAnalysis(BR_calculator, constants, dipole_params, system_params, rates, 
#                                                     quantity='g2 distribution', rho=rho0, 
#                                                     final_time=tf, steps=steps)

#     # Use the provided (theta, phi) for the tilt angles
#     return quantum_analysis_g2_dist.compute(theta_init, phi_init, theta_p, phi_p)

# if config_name == 'H':
#     mu = np.array([1.0, 0.0, 0.0])
# elif config_name == 'J':
#     mu = np.array([0.0, 0.0, 1.0])
# elif config_name == 'ortho':
#     mu = np.array([0.0, 1.0, 0.0]) # The second dipole is along y direction

# # Generate points from von Mises-Fisher distribution
# points = rand_von_mises_fisher(mu, kappa=kappa, N=num_dir)

# # print(points)

# # # Compute intensity distribution
# # g2_corr = np.zeros((num_dir, steps))
# # for i, new_dir in enumerate(points):
# #     g2_corr[i, :] = compute_g2(new_dir)


# # # Create an array for all random dipole directions
# new_dir = [new_dir for new_dir in points]

# # Initialize an array for storing the results
# g2_corr = np.zeros((num_dir, steps))

# # Use ProcessPoolExecutor to parallelize the computation
# with ProcessPoolExecutor() as executor:
#     results = list(executor.map(compute_g2, new_dir))

# # Reshape the results back into a 3D array
# g2_corr = np.array(results).reshape(num_dir, steps)

# # Compute the mean across all angles
# g2_corr = np.mean(g2_corr, axis=(0))

# # print(g2_corr[0])

# # Plotting the results
# # Initialize constants and system parameters
# constants = Constants()
# dipole_params = DipoleParameters(config_name, new_dir=None, average=False)
# system_params = SystemParameters(dipole_params=dipole_params, E_q1=E_q1 * constants.eV, E_q2=E_q2 * constants.eV)
# system_params.calculate_parameters()

# tf = 5 * system_params.tau_L
# times = np.linspace(0, tf, steps)
# times_ns = np.linspace(0, 1.0e-3 * tf, steps)

# # sigma = 40  # in picoseconds

# # dt = (tf - (-tf)) / steps
# # t_kernel = np.arange(-10 * sigma, 10 * sigma, dt)
# # gaussian_kernel = np.exp(-t_kernel**2 / (2 * sigma**2))
# # gaussian_kernel /= np.sum(gaussian_kernel)  

# # g2_corr_appended = fftconvolve(np.append(np.flip(np.load('files/g2/H_dimer/H_strongC_2nm_mav_k10.npy')[1:]), np.load('files/g2/H_dimer/H_strongC_2nm_mav_k10.npy')), gaussian_kernel, mode = "same") #np.append(np.flip(np.load('files/g2/H_dimer/H_strongC_T0.npy')[1:]), np.load('files/g2/H_dimer/H_strongC_T0.npy')), gaussian_kernel, mode = "same")
# times_appended = np.append(-np.flip((times_ns[1:] - times_ns[0])), times_ns)

# # np.save('files/g2/J_strongC_2nm_symm.npy', g2_corr)

# # np.save('files/g2_con/H_dimer/H_strongC_2nm_mav_k10.npy', g2_corr_appended)

# # # Plot the g2 correlation
# # plt.plot(np.append(-np.flip((times_ns[1:] - times_ns[0])), times_ns), np.append(np.flip(g2_corr[1:]), g2_corr))
# # # plt.plot(times_appended, g2_corr_appended)
# # plt.grid(True)
# # plt.show()



# # # Create the main figure and axis
# # fig, ax = plt.subplots(figsize=(14, 8))

# # # # Plot the main data
# # # time_data = np.append(-np.flip(np.linspace(0, 1.0e-3 * tf, 20000)[1:] - np.linspace(0, 1.0e-3 * tf, 20000)[0]), np.linspace(0, 1.0e-3 * tf, 20000))
# # # line0, = ax.plot(time_data, np.append(np.flip(np.load('files/g2/H_dimer/H_strongC_T0_2nm_oav_k1.npy')[1:]), np.load('files/g2/H_dimer/H_strongC_T0_2nm_oav_k1.npy')), 'blue')
# # # line1, = ax.plot(time_data, np.append(np.flip(np.load('files/g2/H_dimer/H_strongC_T0_2nm_oav_k10.npy')[1:]), np.load('files/g2/H_dimer/H_strongC_T0_2nm_oav_k10.npy')), 'mediumblue')
# # # line2, = ax.plot(time_data, np.append(np.flip(np.load('files/g2/H_dimer/H_strongC_T0_2nm_oav_k100.npy')[1:]), np.load('files/g2/H_dimer/H_strongC_T0_2nm_oav_k100.npy')), 'midnightblue')

# # time_data = times_appended 
# # line0, = ax.plot(time_data, np.load('files/g2_con/H_dimer/H_strongC_2nm_mav_k10.npy'), 'b', linewidth = 3)
# # line1, = ax.plot(time_data, np.load('files/g2_con/J_dimer/J_strongC_2nm_mav_k10.npy'), 'r', linewidth = 3)

# # # Set axis labels and limits for the main plot
# # ax.set_xlabel(r'$\tau$ (ns)', fontsize=35)
# # ax.set_ylabel(r'$g^{(2)}(\infty, \tau)$', fontsize=35)
# # ax.set_xlim([-10, 10])  # Example range
# # ax.tick_params(axis='x', labelsize=30, direction='in', length=6)
# # ax.tick_params(axis='y', labelsize=30, direction='in', length=6)

# # # Add a legend to the main plot
# # # ax.legend(handles=[line0, line1, line2], labels=[r"$\kappa = 1$", r"$\kappa = 10$", r"$\kappa = 100$"], loc='best', fontsize=30, frameon=False)
# # ax.legend(handles=[line0, line1], labels=[r"H-dimer", r"J-dimer"], loc='best', fontsize=30, frameon=False)

# # # Create the inset plot
# # ax_inset = plt.axes([0.63, 0.3, 0.25, 0.25])

# # # # Plot the same data in the inset with zoomed-in x-range
# # # ax_inset.plot(time_data, np.append(np.flip(np.load('files/g2/H_dimer/H_strongC_T0_2nm_oav_k1.npy')[1:]), np.load('files/g2/H_dimer/H_strongC_T0_2nm_oav_k1.npy')), 'blue')
# # # ax_inset.plot(time_data, np.append(np.flip(np.load('files/g2/H_dimer/H_strongC_T0_2nm_oav_k10.npy')[1:]), np.load('files/g2/H_dimer/H_strongC_T0_2nm_oav_k10.npy')), 'mediumblue')
# # # ax_inset.plot(time_data, np.append(np.flip(np.load('files/g2/H_dimer/H_strongC_T0_2nm_oav_k100.npy')[1:]), np.load('files/g2/H_dimer/H_strongC_T0_2nm_oav_k100.npy')), 'midnightblue')

# # line0, = ax_inset.plot(time_data, np.load('files/g2_con/H_dimer/H_strongC_2nm_mav_k10.npy'), 'b', linewidth = 3)
# # line1, = ax_inset.plot(time_data, np.load('files/g2_con/J_dimer/J_strongC_2nm_mav_k10.npy'), 'r', linewidth = 3)


# # # Set axis limits for the inset (zoom in on a specific region)
# # ax_inset.set_xlim([-0.3, 0.3])
# # ax_inset.set_ylim([0.27, 0.62])

# # ax_inset.set_ylabel(r'$g^{(2)}(\infty, \tau)$', fontsize=27)
# # ax_inset.tick_params(axis='x', labelsize=25, direction='in', length=6)
# # ax_inset.tick_params(axis='y', labelsize=25, direction='in', length=6)
# # ax_inset.yaxis.set_major_locator(MaxNLocator(nbins=4))

# # # Save and show the plot with the inset
# # plt.savefig('orient_g2_con_strongC_2nm_mav', dpi=500)
# # plt.show()



# # def plot_points_with_axes(points, mu):
# #     fig = plt.figure(figsize=(6, 6))
# #     ax = fig.add_subplot(111, projection='3d')

# #     # Plot the points sampled from the von Mises-Fisher distribution
# #     ax.scatter(points[:, 0], points[:, 1], points[:, 2], c='b', s=30, alpha=0.6)

# #     axis_length = 1.5  # Length of the axes
# #     arrow_head_ratio = 0.05  # Block arrow head by reducing arrow length ratio

# #     # Move the origin to between the two red arrows (same origin for x, y, z)
# #     origin_x = 0
# #     origin_y = 0
# #     origin_z = -1  # Center between the two red arrows (-1.5 * mu[0] and 0)

# #     # X axis
# #     ax.quiver(origin_x, origin_y, origin_z, axis_length, 0, 0, color="black", arrow_length_ratio=arrow_head_ratio, linewidth=2)
# #     # Y axis
# #     ax.quiver(origin_x, origin_y, origin_z, 0, axis_length, 0, color="black", arrow_length_ratio=arrow_head_ratio, linewidth=2)
# #     # Z axis - Both sides (up and down)
# #     ax.quiver(origin_x, origin_y, origin_z, 0, 0, axis_length, color="black", arrow_length_ratio=arrow_head_ratio, linewidth=2)  # Z upwards
# #     ax.quiver(origin_x, origin_y, origin_z, 0, 0, -axis_length, color="black", arrow_length_ratio=arrow_head_ratio, linewidth=2)  # Z downwards


# #     # Plot the two red arrows for the mean direction along the z direction
# #     ax.quiver(0, 0, -2.25, 0, 0, 0.75, color="red", arrow_length_ratio=0.1, linewidth=6)  # From -1.5 to -0.5
# #     ax.quiver(0, 0, -0.5, 0, 0, 0.75, color="red", arrow_length_ratio=0.1, linewidth=6)   # From 0.5 to 1.5

# #     # Set equal aspect ratio and limits
# #     ax.set_box_aspect([1, 1, 1])
# #     ax.set_xlim([-1.5, 1.5])
# #     ax.set_ylim([-1.5, 1.5])
# #     ax.set_zlim([-1.5, 1.5])

# #     # Remove axes, grid, and background
# #     ax.grid(False)
# #     ax.axis('off')  # Remove the axes

# #     # Show the figure
# #     plt.show()

# # # Save the plot
# # plot_points_with_axes(points, mu)






# '''distance dependence on g(2) - get back Moritz's results on being secular'''

# # def compute_g2(r):
    
# #     # Initialize constants and system parameters
# #     constants = Constants()
    
# #     # Pass the tuple (theta, phi) as disorder angles
# #     dipole_params = DipoleParameters(config_name, new_dir=None, average=False)
    
# #     system_params = SystemParameters(dipole_params=dipole_params, r = r, E_q1=E_q1 * constants.eV, E_q2=E_q2 * constants.eV)
# #     system_params.calculate_parameters()

# #     # Calculate spectral densities and rates
# #     rates = CombinedRates(constants, dipole_params, system_params, T_opt, T_vib, mode)

# #     Gamma_opt = rates.Opt_rate()
# #     Gamma_vib = rates.Vib_rate()
# #     kappa = rates.kappa()
# #     Gamma_coup = rates.Coup_rate()
# #     P_vib = rates.P_weight()

# #     # Hamiltonian
# #     ham_calc = Hamiltonian(system_params, rates, mode, kappa)
# #     H = ham_calc.get_hamiltonian()

# #     H_diag, rho0 = initial_state(H, psi0)

# #     # Define time points and steps for the simulation
# #     tf = 5 * system_params.tau_L

# #     # Convert cartesian to spherical coordinates for initial angles
# #     if config_name == 'H':
# #         perp = [0.0, 1.0, 0.0]
# #         # perp = [0.0, 0.0, 1.0]
# #     elif config_name == 'J':
# #         perp = [0.0, 1.0, 0.0]
# #     elif config_name == 'ortho':
# #         perp = [0.0, 0.0, 1.0]
# #         # perp = [0.5, 0.5, 0.0]
# #     r, theta_init, phi_init = cart2sph(perp)
    
# #     r_p, theta_p, phi_p = r, theta_init, phi_init

# #     # Collapse operators (for initial incoherent pumping)
# #     c_ops_S = [np.sqrt(system_params.gamma_p_1) * sp1, np.sqrt(system_params.gamma_p_2) * sp2]

# #     # Create a Bloch-Redfield calculator instance
# #     BR_calculator = BlochRedfieldCalculator(H, constants, dipole_params, system_params, 
# #                                             a_ops_opt, Gamma_opt, 
# #                                             a_ops_vib=a_ops_vib, Gamma_vib=Gamma_vib, 
# #                                             a_ops_coup=a_ops_coup, Gamma_coup=Gamma_coup, 
# #                                             P_vib=P_vib, mode=mode, pump = pump)

# #     # Initialize the DistributionAnalysis class for calculating g2
# #     quantum_analysis_g2_dist = DistributionAnalysis(BR_calculator, constants, dipole_params, system_params, rates, 
# #                                                     quantity='g2 distribution', rho=rho0, 
# #                                                     final_time=tf, steps=steps)

# #     # Use the provided (theta, phi) for the tilt angles
# #     return quantum_analysis_g2_dist.compute(theta_init, phi_init, theta_p, phi_p)

# # # g2_corr = []
# # distance = np.linspace(2, 1000, num_dir)


# # # Initialize an array for storing the results
# # g2_corr = np.zeros((num_dir, steps))

# # # Use ProcessPoolExecutor to parallelize the computation
# # with ProcessPoolExecutor() as executor:
# #     results = list(executor.map(compute_g2, distance))

# # # Reshape the results back into a 3D array
# # g2_corr = np.array(results).reshape(num_dir, steps)


# # # for r in distance:
# # #     g2_corr.append(compute_g2(r))

# # # np.save('g2_corr_H_distance_dep_seq', g2_corr)
# # plt.plot(distance, g2_corr[:,1], 'b', label = 'H-dimer')
# # # plt.plot(distance, np.load('g2_corr_H_distance_dep_seq.npy')[:,1], 'b', label = 'H-dimer')
# # # plt.plot(distance, np.load('g2_corr_J_distance_dep_seq.npy')[:,1], '--r', label = 'J-dimer')
# # plt.ylabel(r'Zero-delay $g^{(2)}(\infty, 0)$')
# # plt.xlabel(r'distance r (nm)')
# # plt.ylim(0, 2.0)
# # plt.legend()
# # # plt.savefig('g2_corr_distance_dep_seq.jpg', dpi = 500)
# # plt.show()



# # # Initialize constants and system parameters
# # constants = Constants()
# # dipole_params = DipoleParameters(config_name, disorder_angles=None, average=False)
# # system_params = SystemParameters(dipole_params=dipole_params, E_q1=E_q1 * constants.eV, E_q2=E_q2 * constants.eV)
# # system_params.calculate_parameters()

# # # Calculate spectral densities and rates
# # rates = CombinedRates(constants, dipole_params, system_params, T_opt, T_vib, mode)

# # Gamma_opt = rates.Opt_rate()
# # Gamma_vib = rates.Vib_rate()
# # kappa = rates.kappa()
# # Gamma_coup = rates.Coup_rate()
# # P_vib = rates.P_weight()
# # print(kappa)

# # # Hamiltonian
# # ham_calc = Hamiltonian(system_params, rates, mode, kappa)
# # H = ham_calc.get_hamiltonian()

# # H_diag, rho0 = initial_state(H, psi0)

# # # Define time points and steps for the simulation
# # steps = 2000
# # tf = 10 * system_params.tau_L
# # times = np.linspace(0, tf, steps)
# # times_ns = np.linspace(0, 1.0e-3 * tf, steps)



# # # Create a Bloch-Redfield calculator instance
# # BR_calculator = BlochRedfieldCalculator(H, constants, dipole_params, system_params, a_ops_opt, Gamma_opt, a_ops_vib=a_ops_vib, Gamma_vib=Gamma_vib, a_ops_coup=a_ops_coup, Gamma_coup=Gamma_coup, c_ops_S=None, P_vib=P_vib, mode=mode)

# # # Initialize the QuantumSystemAnalysis class for calculating intensity
# # quantum_analysis_intensity = QuantumSystemAnalysis(BR_calculator, constants, dipole_params, system_params, rates, quantity='intensity', rho=rho0, final_time=tf, steps=steps)

# # # Compute intensity
# # intensity_results = quantum_analysis_intensity.compute()

# # # np.save('H_noC_T0.npy', intensity_results)

# # # # plt.plot(system_params.gamma_q1*times, np.load('H_weakC_T0.npy') - np.load('H_noC_T0.npy'), label = 'H dimer')
# # # plt.plot(system_params.gamma_q1*times, np.load('J_strongerC.npy') - np.load('J_noC_T0.npy'), label = 'J dimer')
# # # plt.legend()
# # # # plt.savefig('Int.png', dpi = 500)
# # plt.plot(intensity_results)
# # plt.show()


# # # Create a Bloch-Redfield calculator instance
# # BR_calculator = BlochRedfieldCalculator(H, constants, dipole_params, system_params, a_ops_opt, Gamma_opt, a_ops_vib=a_ops_vib, Gamma_vib=Gamma_vib, a_ops_coup=a_ops_coup, Gamma_coup=Gamma_coup, c_ops_S=None, P_vib=P_vib, mode=mode)

# # # Initialize the QuantumSystemAnalysis class for calculating population
# # state = 'doubly excited'
# # quantum_analysis_population = QuantumSystemAnalysis(BR_calculator, constants, dipole_params, system_params, rates, quantity='population', rho=rho0, final_time=tf, steps=steps, state=state)

# # pop_results_de = quantum_analysis_population.compute()

# # # Initialize the QuantumSystemAnalysis class for calculating population
# # state = 'bright'
# # quantum_analysis_population = QuantumSystemAnalysis(BR_calculator, constants, dipole_params, system_params, rates, quantity='population', rho=rho0, final_time=tf, steps=steps, state=state)

# # pop_results_br = quantum_analysis_population.compute()

# # # Initialize the QuantumSystemAnalysis class for calculating population
# # state = 'dark'
# # quantum_analysis_population = QuantumSystemAnalysis(BR_calculator, constants, dipole_params, system_params, rates, quantity='population', rho=rho0, final_time=tf, steps=steps, state=state)

# # pop_results_d = quantum_analysis_population.compute()

# # # Initialize the QuantumSystemAnalysis class for calculating population
# # state = 'ground'
# # quantum_analysis_population = QuantumSystemAnalysis(BR_calculator, constants, dipole_params, system_params, rates, quantity='population', rho=rho0, final_time=tf, steps=steps, state=state)

# # pop_results_g = quantum_analysis_population.compute()



# # plt.plot(pop_results_br, label = 'bright')
# # plt.plot(pop_results_d, label = 'dark')
# # plt.plot(pop_results_de, label = 'doubly excited')
# # plt.plot(pop_results_g, label = 'ground')
# # # plt.ylim([0, 1])
# # plt.legend()
# # plt.show()


# # # Initialize the QuantumSystemAnalysis class for calculating population
# # state = 'bright'
# # quantum_analysis_population = QuantumSystemAnalysis(BR_calculator, constants, dipole_params, system_params, rates, quantity='population', rho=rho0, final_time=tf, steps=steps, state=state)

# # pop_results_br = quantum_analysis_population.compute()

# # # Compute population
# # pop_results_br = quantum_analysis_population.compute()
# # np.save('J_strongerC_bright.npy', pop_results_br)
# # # np.save('files/datafiles/ortho_weakC_bright.npy', pop_results_br)

# # # Initialize the QuantumSystemAnalysis class for calculating population
# # state = 'dark'
# # quantum_analysis_population = QuantumSystemAnalysis(BR_calculator, constants, dipole_params, system_params, rates, quantity='population', rho=rho0, final_time=tf, steps=steps, state=state)

# # # Compute population
# # pop_results_d = quantum_analysis_population.compute()
# # np.save('J_strongerC_dark.npy', pop_results_d)
# # # # np.save('files/datafiles/ortho_strongC_dark.npy', pop_results_d)


# # # plt.plot(system_params.gamma_q1 * times, np.load('H_strongC_bright.npy'), label = 'H_strongC_bright')
# # plt.plot(system_params.gamma_q1 * times, np.load('H_strongC_dark.npy'), label = 'H_strongC_dark')
# # # plt.plot(system_params.gamma_q1 * times, np.load('J_strongerC_bright.npy'), label = 'J_strongerC_bright')
# # plt.plot(system_params.gamma_q1 * times, np.load('J_strongerC_dark.npy'), label = 'J_strongerC_dark')
# # # plt.savefig('pop_results.png')
# # plt.legend(loc='best')
# # # plt.savefig('pop_strongC.png', dpi = 500)
# # # plt.xlim([0, 5])
# # plt.show()



# # plt.plot(system_params.gamma_q1 * times, np.load('H_weakC_bright.npy'), label = 'H_weakC_bright')
# # plt.plot(system_params.gamma_q1 * times, np.load('H_weakC_dark.npy'), label = 'H_weakC_dark')
# # plt.plot(system_params.gamma_q1 * times, np.load('J_weakC_bright.npy'), label = 'J_weakC_bright')
# # plt.plot(system_params.gamma_q1 * times, np.load('J_weakC_dark.npy'), label = 'J_weakC_dark')
# # plt.legend(loc='best')
# # # plt.savefig('pop_weakC.png', dpi = 500)
# # # plt.xlim([0, 5])
# # plt.show()


# # plt.plot(system_params.gamma_q1 * times, np.load('files/datafiles/H_strongC_bright.npy'), label = 'H_strongC_bright')
# # plt.plot(system_params.gamma_q1 * times, np.load('files/datafiles/H_strongC_dark.npy'), label = 'H_strongC_dark')
# # plt.plot(system_params.gamma_q1 * times, np.load('files/datafiles/J_strongerC_bright.npy'), label = 'J_strongerC_bright')
# # plt.plot(system_params.gamma_q1 * times, np.load('files/datafiles/J_strongerC_dark.npy'), label = 'J_strongerC_dark')
# # # plt.savefig('pop_results.png')
# # plt.legend(loc='best')
# # plt.savefig('pop_strongerC.png', dpi = 500)
# # plt.xlim([0, 5])
# # plt.show()


# # plt.plot(pop_results_br)
# # plt.plot(pop_results_d)
# # # plt.ylim([0, 1])
# # plt.show()

# # perp = [0.0, 0.0, 1.0]
# # r, theta, phi = cart2sph(perp)
# # r_p, theta_p, phi_p = r, theta, phi


# # # Create a Bloch-Redfield calculator instance
# # BR_calculator = BlochRedfieldCalculator(H, constants, dipole_params, system_params, a_ops_opt, Gamma_opt, a_ops_vib=a_ops_vib, Gamma_vib=Gamma_vib, a_ops_coup=a_ops_coup, Gamma_coup=Gamma_coup, c_ops_S=None, P_vib=P_vib, mode=mode)

# # # Initialize the DistributionAnalysis class for calculating intensity
# # quantum_analysis_intensity_dist = DistributionAnalysis(BR_calculator, constants, dipole_params, system_params, rates, quantity='intensity distribution', rho=rho0, final_time=tf, steps=steps)

# # norm_int = quantum_analysis_intensity_dist.compute(theta_p, phi_p)

# # plt.plot(system_params.gamma_q1 * times, norm_int)
# # plt.show()




# # # Define parameters for g2 plot
# # n_theta_p = 20
# # n_phi_p = 20

# # # Define parameters
# # theta_vals_p = np.linspace(0, np.pi, n_theta_p)
# # phi_vals_p = np.linspace(0,  2*np.pi, n_phi_p)
# # r_p = 1  # radius of sphere

# # # Generate x, y, z coordinates of points on the sphere
# # x = r_p * np.outer(np.cos(phi_vals_p), np.sin(theta_vals_p))
# # y = r_p * np.outer(np.sin(phi_vals_p), np.sin(theta_vals_p))
# # z = r_p * np.outer(np.ones(np.size(phi_vals_p)), np.cos(theta_vals_p))


# # # Create a Bloch-Redfield calculator instance
# # BR_calculator = BlochRedfieldCalculator(H, constants, dipole_params, system_params, a_ops_opt, Gamma_opt, a_ops_vib=a_ops_vib, Gamma_vib=Gamma_vib, a_ops_coup=a_ops_coup, Gamma_coup=Gamma_coup, c_ops_S=None, P_vib=P_vib, mode=mode)

# # # Initialize the DistributionAnalysis class for calculating intensity
# # quantum_analysis_intensity_dist = DistributionAnalysis(BR_calculator, constants, dipole_params, system_params, rates, quantity='intensity distribution', rho=rho0, final_time=tf, steps=steps)

# # # Compute intensity distribution
# # norm_int = np.zeros((n_phi_p, n_theta_p))
# # for l, phi_p in enumerate(phi_vals_p):
# #     for k, theta_p in enumerate(theta_vals_p):
# #         # theta_p = np.pi / 2
# #         norm_int[l, k] = quantum_analysis_intensity_dist.compute(theta_p, phi_p)



# # # Generate R_vals and Theta_vals for the polar plot
# # R_vals, Theta_vals = np.meshgrid(np.sin(theta_vals_p), phi_vals_p)

# # # Create the polar heatmap
# # fig = plt.figure(layout='constrained')
# # ax1 = fig.add_subplot(projection='polar')
# # c = ax1.pcolormesh(Theta_vals, R_vals, norm_int, cmap='coolwarm', shading='auto', vmin=0, vmax=1)

# # ax1.set_xlabel(r'$\theta$ (degrees)', fontsize=18, labelpad=8) 
# # ax1.set_ylabel(r'$I(\theta)$', fontsize=18, labelpad=28)
# # ax1.tick_params(axis='x', labelsize=16, direction='in', length=6)
# # ax1.tick_params(axis='y', labelsize=16, direction='in', length=6)
# # ax1.set_yticklabels([])
# # # # Set radial limits
# # # ax1.set_ylim(0, 1)
# # # Add colorbar
# # # cbar = plt.colorbar(c)
# # # # Add arrows
# # # # Arrow facing theta = 0 (along the x-axis)
# # # ax1.annotate('', xy=(0, 0.2), xytext=(np.pi, 0.2),
# # #              arrowprops=dict(facecolor='blue', shrink=0, width=2, headwidth=10))

# # # # Second arrow (also stretching from -0.25 to 0.25 along the x-axis)
# # # ax1.annotate('', xy=(0, 0.3), xytext=(np.pi, 0.25),
# # #              arrowprops=dict(facecolor='blue', shrink=0, width=2, headwidth=10))

# # # # Save the figure
# # plt.savefig('dist_int_J.png')

# # plt.show()

# # perp_H = [0, 0, 1]
# # r, theta, phi = cart2sph(perp_H) # cart2sph(d1_hat)#
# # r_p, theta_p, phi_p = r, theta, phi 

# # perp_J = [1, 0, 0]
# # r, theta, phi = cart2sph(perp_J)
# # r_p, theta_p, phi_p = r, theta, phi

# # perp_ortho = [0.0, 1.0, 0.0]
# # perp_ortho = [1.0, 0.0, 0.0]
# # perp_ortho = [0.0, 0.0, 1.0]
# perp_ortho = [0.5, 0.5, 0.0]
# r, theta, phi = cart2sph(perp_ortho)
# r_p, theta_p, phi_p = r, theta, phi

# # perp_oppo = [0, 0, 1]
# # r, theta, phi = cart2sph(perp_oppo)


# # # Initialize the DistributionAnalysis class for calculating intensity
# # quantum_analysis_intensity_dist = DistributionAnalysis(BR_calculator, quantity='intensity distribution', rho=rho0, final_time=tf, steps=steps)

# # # Compute intensity distribution
# # intensity_results = quantum_analysis_intensity_dist.compute(theta_p, phi_p)




# # tf = 5 * system_params.tau_L
# # steps = 200000
# # times = np.linspace(0, tf, steps)
# # times_ns = np.linspace(0, 1.0e-3 * tf, steps)

# # # collapse operators (for initial incoherent pumping)
# # c_ops_S = [np.sqrt(system_params.gamma_p_1) * sp1, np.sqrt(system_params.gamma_p_2) * sp2]
# # # c_ops_S = [np.sqrt(system_params.gamma_p_1) * sp1, np.sqrt(system_params.gamma_p_2) * sp2, np.sqrt(system_params.gamma_q1) * sp1@sm1, np.sqrt(system_params.gamma_q1) * sp2@sm2]

# # # Create a Bloch-Redfield calculator instance
# # BR_calculator = BlochRedfieldCalculator(H, a_ops_opt, Gamma_opt, a_ops_vib = a_ops_vib, Gamma_vib = Gamma_vib, a_ops_coup = a_ops_coup, Gamma_coup = Gamma_coup, c_ops_S = c_ops_S, P_vib = P_vib, mode=mode)

# # # Initialize the DistributionAnalysis class for calculating g2
# # quantum_analysis_g2_dist = DistributionAnalysis(BR_calculator, quantity='g2 distribution', rho=rho0, final_time=tf, steps=steps)

# # g2_corr = quantum_analysis_g2_dist.compute(theta, phi, theta_p, phi_p)

# # # np.save('files/g2/H_dimer/H_weakC.npy', g2_corr)

# # # # # # print(g2_corr)
# # # # np.save('g2_par', g2_corr)
# # plt.plot(np.append(- np.flip((times_ns[1:] - times_ns[0])), (times_ns)), np.append(np.flip(g2_corr[1:]), g2_corr))
# # plt.grid(True)
# # # plt.xlim([-5, 5])
# # plt.show()







# # # # print(norm_int_dir)
# # # print(g2_corr[0])
# # plt.plot(np.append(- np.flip((times_ns[1:] - times_ns[0])), (times_ns)), np.append(np.flip(g2_corr[1:]), g2_corr[:]))
# # plt.grid(True)
# # # plt.ylim([0.0,1.1])
# # plt.show()

# # plt.plot(np.append(- np.flip((times_ns[1:] - times_ns[0])), (times_ns)), np.append(np.flip(g2_corr[0, 0, 1:]), g2_corr[0, 0, :]))
# # plt.show()
























# t1 = time.time()

# total = t1-t0
# print(total)

















































# # # Load the data
# # H_dimer_weak = np.load('files/intensity/H_dimer/H_weakC.npy')
# # J_dimer_weak = np.load('files/intensity/H_dimer/H_weakC.npy')
# # H_dimer_strong = np.load('files/intensity/H_dimer/H_strongC.npy')
# # J_dimer_strong = np.load('files/intensity/H_dimer/H_strongC.npy')
# # H_dimer_stronger = np.load('files/intensity/H_dimer/H_strongC.npy')
# # J_dimer_stronger = np.load('files/intensity/H_dimer/H_strongC.npy')
# # ortho_dimer_weak = np.load('files/intensity/ortho_dimer/ortho_weakC.npy')
# # ortho_dimer_strong = np.load('files/intensity/ortho_dimer/ortho_strongC.npy')
# # ortho_dimer_stronger = np.load('files/intensity/ortho_dimer/ortho_strongerC.npy')
# # anti_parallel = np.load('files/intensity/oppo_dimer/oppo_strongerC.npy')
# # H_dimer_noC = np.load('files/intensity/H_dimer/H_noC.npy')
# # J_dimer_noC = np.load('files/intensity/H_dimer/H_noC.npy')
# # ortho_dimer_noC = np.load('files/intensity/ortho_dimer/ortho_noC.npy')


# # # Create main plot
# # fig, ax = plt.subplots(figsize=(12, 8))

# # # Plot the data with appropriate labels
# # line2, = ax.plot(system_params.gamma_q1 * times, H_dimer_weak, color='blue', linewidth=0.6)
# # line3, = ax.plot(system_params.gamma_q1 * times, J_dimer_weak, color='tomato', linewidth=0.6)
# # line4, = ax.plot(system_params.gamma_q1 * times, H_dimer_strong, color='mediumblue', linewidth=0.6)
# # line5, = ax.plot(system_params.gamma_q1 * times, J_dimer_strong, color='red', linewidth=0.6)
# # line6, = ax.plot(system_params.gamma_q1 * times, H_dimer_stronger, color='darkblue', linewidth=0.6)
# # line7, = ax.plot(system_params.gamma_q1 * times, J_dimer_stronger, color='darkred', linewidth=0.6)
# # line8, = ax.plot(system_params.gamma_q1 * times, ortho_dimer_weak, '--g')

# # # Set labels and title
# # ax.set_xlabel(r'$\gamma t$', fontsize=30)
# # ax.set_ylabel(r'$I(t) / I_0$', fontsize=30)
# # ax.set_ylim([0, 2.1])
# # ax.set_xlim([system_params.gamma_q1 * times[0], system_params.gamma_q1 * times[-1]])

# # # Set tick parameters
# # ax.tick_params(axis='x', labelsize=20, direction='in', length=6)
# # ax.tick_params(axis='y', labelsize=20, direction='in', length=6)
# # ax.text(0.1, 2, '(a)', size=20)

# # # Change number of ticks on y-axis
# # ax.yaxis.set_major_locator(MaxNLocator(nbins=5))

# # # Add legend and grid
# # ax.legend(handles=[line2, line3, line4, line5, line6, line7, line8],
# #           labels=['H dimer-weak', 'J dimer-weak', 'H dimer-strong', 'J dimer-strong', 'H dimer-stronger', 'J dimer-stronger', 'Orthogonal'],
# #           loc='best',
# #           fontsize=16,
# #           frameon=False)

# # # Inset plot for H-dimer vs J-dimer differences
# # ax_inset_1 = plt.axes([0.6, 0.2, 0.25, 0.35])  # [left, bottom, width, height]
# # diff_data_1 = H_dimer_weak - J_dimer_weak
# # diff_data_2 = H_dimer_strong - J_dimer_strong
# # diff_data_3 = H_dimer_stronger - J_dimer_stronger
# # ax_inset_1.plot(system_params.gamma_q1 * times, diff_data_1, label=r'$\lambda = 5 \mu eV$, $\omega_{c} = 90 \mu eV$')
# # ax_inset_1.plot(system_params.gamma_q1 * times, diff_data_2, label=r'$\lambda = 5 meV$, $\omega_{c} = 90 meV$')
# # ax_inset_1.plot(system_params.gamma_q1 * times, diff_data_3, label=r'$\lambda = 20 meV$, $\omega_{c} = 90 meV$')
# # ax_inset_1.set_xlabel(r'$\gamma t$', fontsize=15)
# # ax_inset_1.set_ylabel(r'$\Delta I(t) / I_0$', fontsize=15)
# # ax_inset_1.tick_params(axis='x', labelsize=10)
# # ax_inset_1.tick_params(axis='y', labelsize=10)
# # ax_inset_1.set_xlim([system_params.gamma_q1 * times[0], system_params.gamma_q1 * times[-1]])
# # ax_inset_1.set_title(r'H-dimer vs J-dimer (difference)')
# # ax_inset_1.legend(fontsize=10, frameon=False)


# # # Save the main plot
# # plt.savefig("orient_intensity_strongC.png", dpi=500)
# # plt.show()






# # # Independent plot for J-dimer vs NoC differences
# # fig2, ax2 = plt.subplots(figsize=(10, 6))
# # diff_data_4 = J_dimer_weak - J_dimer_noC 
# # diff_data_5 = J_dimer_strong - J_dimer_noC 
# # diff_data_6 = J_dimer_stronger - J_dimer_noC
# # ax2.plot(system_params.gamma_q1 * times, diff_data_4, label=r'$\lambda = 5 \mu eV$, $\omega_{c} = 90 \mu eV$')
# # ax2.plot(system_params.gamma_q1 * times, diff_data_5, label=r'$\lambda = 5 meV$, $\omega_{c} = 90 meV$')
# # ax2.plot(system_params.gamma_q1 * times, diff_data_6, label=r'$\lambda = 20 meV$, $\omega_{c} = 90 meV$')
# # ax2.set_xlabel(r'$\gamma t$', fontsize=25)
# # ax2.set_ylabel(r'$\Delta I(t) / I_0$', fontsize=25)
# # ax2.tick_params(axis='x', labelsize=15)
# # ax2.tick_params(axis='y', labelsize=15)
# # ax2.set_xlim([0, 5])
# # ax2.legend(fontsize=15, frameon=False)
# # plt.savefig("j_dimer_vs_noc_diff.png", dpi=500)
# # plt.show()


# # # New plot for H-dimer vs NoC differences
# # fig3, ax3 = plt.subplots(figsize=(10, 6))
# # diff_data_7 = H_dimer_weak - H_dimer_noC 
# # diff_data_8 = H_dimer_strong - H_dimer_noC 
# # diff_data_9 = H_dimer_stronger - H_dimer_noC
# # ax3.plot(system_params.gamma_q1 * times, diff_data_7, label=r'$\lambda = 5 \mu eV$, $\omega_{c} = 90 \mu eV$')
# # ax3.plot(system_params.gamma_q1 * times, diff_data_8, label=r'$\lambda = 5 meV$, $\omega_{c} = 90 meV$')
# # ax3.plot(system_params.gamma_q1 * times, diff_data_9, label=r'$\lambda = 20 meV$, $\omega_{c} = 90 meV$')
# # ax3.set_xlabel(r'$\gamma t$', fontsize=25)
# # ax3.set_ylabel(r'$\Delta I(t) / I_0$', fontsize=25)
# # ax3.tick_params(axis='x', labelsize=15)
# # ax3.tick_params(axis='y', labelsize=15)
# # ax3.set_xlim([0, 5])
# # ax3.legend(fontsize=15, frameon=False)
# # plt.savefig("h_dimer_vs_noc_diff.png", dpi=500)
# # plt.show()





# # # Independent plot for ortho-dimer different couplings
# # fig4, ax4 = plt.subplots(figsize=(10, 6))
# # ax4.plot(system_params.gamma_q1 * times, ortho_dimer_weak - ortho_dimer_noC, label=r'$\lambda = 5 \mu eV$, $\omega_{c} = 90 \mu eV$')
# # ax4.plot(system_params.gamma_q1 * times, ortho_dimer_strong - ortho_dimer_noC, label=r'$\lambda = 5 meV$, $\omega_{c} = 90 meV$')
# # ax4.plot(system_params.gamma_q1 * times, ortho_dimer_stronger - ortho_dimer_noC, label=r'$\lambda = 20 meV$, $\omega_{c} = 90 meV$')
# # ax4.set_xlabel(r'$\gamma t$', fontsize=25)
# # ax4.set_ylabel(r'$\Delta I(t) / I_0$', fontsize=25)
# # ax4.tick_params(axis='x', labelsize=15)
# # ax4.tick_params(axis='y', labelsize=15)
# # ax4.set_xlim([system_params.gamma_q1 * times[0], system_params.gamma_q1 * times[-1]])
# # ax4.legend(fontsize=15, frameon=False)
# # plt.savefig("ortho_dimer_vs_noc_diff.png", dpi=500)
# # plt.show()








# # # Load the data
# # H_dimer_weak = np.load('files/intensity/H_dimer/H_weakC.npy')
# # J_dimer_weak = np.load('files/intensity/H_dimer/H_weakC.npy')
# # H_dimer_strong = np.load('files/intensity/H_dimer/H_strongC.npy')
# # J_dimer_strong = np.load('files/intensity/H_dimer/H_strongC.npy')
# # H_dimer_stronger = np.load('files/intensity/H_dimer/H_strongC.npy')
# # J_dimer_stronger = np.load('files/intensity/H_dimer/H_strongC.npy')
# # ortho_dimer_weak = np.load('files/intensity/ortho_dimer/ortho_weakC.npy')
# # ortho_dimer_strong = np.load('files/intensity/ortho_dimer/ortho_strongC.npy')
# # ortho_dimer_stronger = np.load('files/intensity/ortho_dimer/ortho_strongerC.npy')
# # H_dimer_noC = np.load('files/intensity/H_dimer/H_noC.npy')
# # J_dimer_noC = np.load('files/intensity/H_dimer/H_noC.npy')
# # ortho_dimer_noC = np.load('files/intensity/ortho_dimer/ortho_noC.npy')


# # # Create main plot
# # fig, ax = plt.subplots(figsize=(12, 8))

# # # Plot the data with appropriate labels
# # line2, = ax.plot(system_params.gamma_q1 * times, H_dimer_weak, color='blue', linewidth=0.6)
# # line3, = ax.plot(system_params.gamma_q1 * times, J_dimer_weak, color='tomato', linewidth=0.6)
# # line4, = ax.plot(system_params.gamma_q1 * times, H_dimer_strong, color='mediumblue', linewidth=0.6)
# # line5, = ax.plot(system_params.gamma_q1 * times, J_dimer_strong, color='red', linewidth=0.6)
# # line6, = ax.plot(system_params.gamma_q1 * times, H_dimer_stronger, color='darkblue', linewidth=0.6)
# # line7, = ax.plot(system_params.gamma_q1 * times, J_dimer_stronger, color='darkred', linewidth=0.6)
# # line8, = ax.plot(system_params.gamma_q1 * times, ortho_dimer_weak, '--g')

# # # Set labels and title
# # ax.set_xlabel(r'$\gamma t$', fontsize=30)
# # ax.set_ylabel(r'$I(t) / I_0$', fontsize=30)
# # ax.set_ylim([0, 2.1])
# # ax.set_xlim([system_params.gamma_q1 * times[0], system_params.gamma_q1 * times[-1]])

# # # Set tick parameters
# # ax.tick_params(axis='x', labelsize=20, direction='in', length=6)
# # ax.tick_params(axis='y', labelsize=20, direction='in', length=6)
# # ax.text(0.1, 2, '(a)', size=20)

# # # Change number of ticks on y-axis
# # ax.yaxis.set_major_locator(MaxNLocator(nbins=5))

# # # Add legend and grid
# # ax.legend(handles=[line2, line3, line4, line5, line6, line7, line8],
# #           labels=['H dimer-weak', 'J dimer-weak', 'H dimer-strong', 'J dimer-strong', 'H dimer-stronger', 'J dimer-stronger', 'Orthogonal'],
# #           loc='best',
# #           fontsize=16,
# #           frameon=False)

# # # Inset plot for H-dimer vs J-dimer differences
# # ax_inset_1 = plt.axes([0.6, 0.2, 0.25, 0.35])  # [left, bottom, width, height]
# # diff_data_1 = H_dimer_weak - J_dimer_weak
# # diff_data_2 = H_dimer_strong - J_dimer_strong
# # diff_data_3 = H_dimer_stronger - J_dimer_stronger
# # ax_inset_1.plot(system_params.gamma_q1 * times, diff_data_1, label=r'$\lambda = 5 \mu eV$, $\omega_{c} = 90 \mu eV$')
# # ax_inset_1.plot(system_params.gamma_q1 * times, diff_data_2, label=r'$\lambda = 5 meV$, $\omega_{c} = 90 meV$')
# # ax_inset_1.plot(system_params.gamma_q1 * times, diff_data_3, label=r'$\lambda = 20 meV$, $\omega_{c} = 90 meV$')
# # ax_inset_1.set_xlabel(r'$\gamma t$', fontsize=15)
# # ax_inset_1.set_ylabel(r'$\Delta I(t) / I_0$', fontsize=15)
# # ax_inset_1.tick_params(axis='x', labelsize=10)
# # ax_inset_1.tick_params(axis='y', labelsize=10)
# # ax_inset_1.set_xlim([system_params.gamma_q1 * times[0], system_params.gamma_q1 * times[-1]])
# # ax_inset_1.set_title(r'H-dimer vs J-dimer (difference)')
# # ax_inset_1.legend(fontsize=10, frameon=False)


# # # Save the main plot
# # plt.savefig("orient_intensity_strongC.png", dpi=500)
# # plt.show()






# # # Independent plot for J-dimer vs NoC differences
# # fig2, ax2 = plt.subplots(figsize=(10, 6))
# # diff_data_4 = H_dimer_weak - H_dimer_noC 
# # diff_data_5 = J_dimer_weak - J_dimer_noC 
# # diff_data_6 = ortho_dimer_weak - ortho_dimer_noC
# # ax2.plot(system_params.gamma_q1 * times, diff_data_4, label=r'H-dimer')
# # ax2.plot(system_params.gamma_q1 * times, diff_data_5, label=r'J-dimer')
# # ax2.plot(system_params.gamma_q1 * times, diff_data_6, label=r'Ortho-dimer')
# # ax2.set_xlabel(r'$\gamma t$', fontsize=25)
# # ax2.set_ylabel(r'$\Delta I(t) / I_0$', fontsize=25)
# # ax2.tick_params(axis='x', labelsize=15)
# # ax2.tick_params(axis='y', labelsize=15)
# # ax2.set_xlim([0, 5])
# # ax2.legend(fontsize=15, frameon=False)
# # ax2.set_title(r'Weak - $\lambda = 5 \mu eV$, $\omega_{c} = 90 \mu eV$')
# # plt.savefig("weakC.png", dpi=500)
# # plt.show()


# # # New plot for H-dimer vs NoC differences
# # fig3, ax3 = plt.subplots(figsize=(10, 6))
# # diff_data_7 = H_dimer_strong - H_dimer_noC 
# # diff_data_8 = J_dimer_strong - J_dimer_noC 
# # diff_data_9 = ortho_dimer_stronger - ortho_dimer_noC
# # ax3.plot(system_params.gamma_q1 * times, diff_data_7, label=r'H-dimer')
# # ax3.plot(system_params.gamma_q1 * times, diff_data_8, label=r'J-dimer')
# # ax3.plot(system_params.gamma_q1 * times, diff_data_9, label=r'Ortho-dimer')
# # ax3.set_xlabel(r'$\gamma t$', fontsize=25)
# # ax3.set_ylabel(r'$\Delta I(t) / I_0$', fontsize=25)
# # ax3.tick_params(axis='x', labelsize=15)
# # ax3.tick_params(axis='y', labelsize=15)
# # ax3.set_xlim([0, 5])
# # ax3.legend(fontsize=15, frameon=False)
# # ax3.set_title(r'Strong - $\lambda = 5 meV$, $\omega_{c} = 90 meV$')
# # plt.savefig("strongC.png", dpi=500)
# # plt.show()





# # # Independent plot for ortho-dimer different couplings
# # fig4, ax4 = plt.subplots(figsize=(10, 6))
# # ax4.plot(system_params.gamma_q1 * times, H_dimer_stronger - H_dimer_noC, label=r'H-dimer')
# # ax4.plot(system_params.gamma_q1 * times, J_dimer_stronger - J_dimer_noC, label=r'J-dimer')
# # ax4.plot(system_params.gamma_q1 * times, ortho_dimer_stronger - ortho_dimer_noC, label=r'Ortho-dimer')
# # ax4.set_xlabel(r'$\gamma t$', fontsize=25)
# # ax4.set_ylabel(r'$\Delta I(t) / I_0$', fontsize=25)
# # ax4.tick_params(axis='x', labelsize=15)
# # ax4.tick_params(axis='y', labelsize=15)
# # ax4.set_xlim([system_params.gamma_q1 * times[0], system_params.gamma_q1 * times[-1]])
# # ax4.legend(fontsize=15, frameon=False)
# # ax4.set_title(r'$\lambda = 20 meV$, $\omega_{c} = 90 meV$')
# # plt.savefig("strongerC.png", dpi=500)
# # plt.show()




# # plt.plot(J_dimer_noC)
# # plt.plot(J_dimer_weak)
# # plt.show()




























# '''Wierd spikes instead of lobes'''


# # # Define parameters for g2 plot
# # n_theta_p = 200
# # n_phi_p = 200

# # # Define parameters
# # theta_vals_p = np.linspace(0, np.pi, n_theta_p)
# # phi_vals_p = np.linspace(0,  2*np.pi, n_phi_p)
# # r_p = 1  # radius of sphere



# # theta_values, phi_values = np.meshgrid(theta_vals_p, phi_vals_p)
# # # Create a Bloch-Redfield calculator instance
# # BR_calculator = BlochRedfieldCalculator(H, constants, dipole_params, system_params, a_ops_opt, Gamma_opt, a_ops_vib=a_ops_vib, Gamma_vib=Gamma_vib, a_ops_coup=a_ops_coup, Gamma_coup=Gamma_coup, c_ops_S=None, P_vib=P_vib, mode=mode)

# # # Initialize the DistributionAnalysis class for calculating intensity
# # quantum_analysis_intensity_dist = DistributionAnalysis(BR_calculator, constants, dipole_params, system_params, rates, quantity='intensity distribution', rho=rho0, final_time=tf, steps=steps)

# # # Compute intensity distribution
# # norm_int = np.zeros((n_phi_p, n_theta_p))
# # for l, phi_p in enumerate(range(phi_values.shape[0])):
# #     for k, theta_p in enumerate(range(theta_values.shape[1])):
# #         # theta_p = np.pi / 2
# #         norm_int[l, k] = quantum_analysis_intensity_dist.compute(theta_p, phi_p)


# # r_p = norm_int / np.max(norm_int)
# # # Generate x, y, z coordinates of points on the sphere
# # x = r_p * np.outer(np.cos(phi_vals_p), np.sin(theta_vals_p))
# # y = r_p * np.outer(np.sin(phi_vals_p), np.sin(theta_vals_p))
# # z = r_p * np.outer(np.ones(np.size(phi_vals_p)), np.cos(theta_vals_p))

# # fig = plt.figure()
# # ax = fig.add_subplot(111,projection="3d")
# # ax.plot_surface(x, y, z)

# # plt.show()









































# # kappa = 10
# # num_dir = 100
# # steps = 2000

# # mode = "tumbling"

# # def compute_g2(mode, orientation=None, sample_direction=None):
# #     # Initialize constants and system parameters
# #     constants = Constants()
    
# #     # Adjust disorder angles based on the mode
# #     if mode == "tumbling":
# #         # Tumbling effect - angles (theta, phi) are None, direction is sampled
# #         dipole_params = DipoleParameters(config_name, new_dir=None, average=False)
# #     elif mode == "jiggling":
# #         # Jiggling effect - disorder angles provided
# #         dipole_params = DipoleParameters(config_name, new_dir=orientation, average=True)
    
# #     system_params = SystemParameters(dipole_params=dipole_params, E_q1=E_q1 * constants.eV, E_q2=E_q2 * constants.eV)
# #     system_params.calculate_parameters()

# #     # Calculate spectral densities and rates
# #     rates = CombinedRates(constants, dipole_params, system_params, T_opt, T_vib, mode)
# #     Gamma_opt = rates.Opt_rate()
# #     Gamma_vib = rates.Vib_rate()
# #     kappa = rates.kappa()
# #     Gamma_coup = rates.Coup_rate()
# #     P_vib = rates.P_weight()

# #     # Hamiltonian
# #     ham_calc = Hamiltonian(system_params, rates, mode, kappa)
# #     H = ham_calc.get_hamiltonian()

# #     H_diag, rho0 = initial_state(H, psi0)

# #     # Define time points and steps for the simulation
# #     tf = 5 * system_params.tau_L

# #     # Convert cartesian to spherical coordinates for initial angles
# #     if config_name == 'H':
# #         perp = [0.0, 1.0, 0.0]
# #     elif config_name == 'J':
# #         perp = [0.0, 1.0, 0.0]
# #     elif config_name == 'ortho':
# #         perp = [0.0, 0.0, 1.0]

# #     r, theta_init, phi_init = cart2sph(perp)

# #     if mode == "tumbling":
# #         # Random direction for second photon in tumbling effect
# #         theta_p, phi_p = cart2sph(sample_direction)
# #     elif mode == "jiggling":
# #         # Fixed dipole orientation for jiggling effect
# #         theta_p, phi_p = theta_init, phi_init

# #     # Collapse operators (for initial incoherent pumping)
# #     c_ops_S = [np.sqrt(system_params.gamma_p_1) * sp1, np.sqrt(system_params.gamma_p_2) * sp2]

# #     # Create a Bloch-Redfield calculator instance
# #     BR_calculator = BlochRedfieldCalculator(H, constants, dipole_params, system_params, 
# #                                             a_ops_opt, Gamma_opt, 
# #                                             a_ops_vib=a_ops_vib, Gamma_vib=Gamma_vib, 
# #                                             a_ops_coup=a_ops_coup, Gamma_coup=Gamma_coup, 
# #                                             P_vib=P_vib, mode=mode, pump = pump)

# #     # Initialize the DistributionAnalysis class for calculating g2
# #     quantum_analysis_g2_dist = DistributionAnalysis(BR_calculator, constants, dipole_params, system_params, rates, 
# #                                                     quantity='g2 distribution', rho=rho0, 
# #                                                     final_time=tf, steps=steps)

# #     # Compute the g2 for the respective angles
# #     return quantum_analysis_g2_dist.compute(theta_init, phi_init, theta_p, phi_p)


# # g2_corr = np.zeros((num_dir, steps))

# # if mode == "tumbling":
# #     if config_name == 'H':
# #         mu = [1.0, 0.0, 0.0]
# #         perp = [0.0, 1.0, 0.0]
# #     elif config_name == 'J':
# #         mu = [0.0, 0.0, 1.0]
# #         perp = [0.0, 1.0, 0.0]
# #     elif config_name == 'ortho':
# #         mu = [1.0, 0.0, 0.0]
# #         perp = [0.0, 0.0, 1.0]

# #     # Generate random directions from von Mises-Fisher distribution
# #     points = vonmises_fisher.rvs(perp, kappa=kappa, size=num_dir)
# #     # Tumbling effect: loop over sampled directions
# #     samp_dir = [samp_dir for samp_dir in points]

# #     orientation, sample_direction = mu, samp_dir

# #     with ProcessPoolExecutor() as executor:
# #         results = list(executor.map(compute_g2, samp_dir))

# # elif mode == "jiggling":
# #     if config_name == 'H':
# #         mu = [1.0, 0.0, 0.0]
# #         perp = [0.0, 1.0, 0.0]
# #     elif config_name == 'J':
# #         mu = [0.0, 0.0, 1.0]
# #         perp = [0.0, 1.0, 0.0]
# #     elif config_name == 'ortho':
# #         mu = [1.0, 0.0, 0.0]
# #         perp = [0.0, 0.0, 1.0]

# #     # Generate random directions from von Mises-Fisher distribution
# #     points = vonmises_fisher.rvs(mu, kappa=kappa, size=num_dir)
# #     # Jiggling effect: loop over mutual orientations
# #     new_dir = [new_dir for new_dir in points]

# #     orientation, sample_direction = new_dir, perp

# #     with ProcessPoolExecutor() as executor:
# #         results = list(executor.map(compute_g2, new_dir))

# # # Reshape and compute mean across all angles
# # g2_corr = np.array(results).reshape(num_dir, steps)
# # g2_corr = np.mean(g2_corr, axis=0)

# # # Plotting the results
# # constants = Constants()
# # dipole_params = DipoleParameters(config_name, new_dir=None, average=False)
# # system_params = SystemParameters(dipole_params=dipole_params, E_q1=E_q1 * constants.eV, E_q2=E_q2 * constants.eV)
# # system_params.calculate_parameters()

# # tf = 5 * system_params.tau_L
# # times_ns = np.linspace(0, 1.0e-3 * tf, steps)


# # # sigma = 50  # in picoseconds

# # dt = (tf - (-tf)) / steps
# # # t_kernel = np.arange(-10 * sigma, 10 * sigma, dt)
# # # gaussian_kernel = np.exp(-t_kernel**2 / (2 * sigma**2))
# # # gaussian_kernel /= np.sum(gaussian_kernel)  

# # # g2_corr_appended = np.append(np.flip(np.load('files/g2/H_dimer/H_strongC_T0_2nm_oav_k1.npy')[1:]), np.load('files/g2/H_dimer/H_strongC_T0_2nm_oav_k1.npy')), gaussian_kernel, mode = "same") #np.append(np.flip(np.load('files/g2/H_dimer/H_strongC_T0.npy')[1:]), np.load('files/g2/H_dimer/H_strongC_T0.npy')), gaussian_kernel, mode = "same")
# # # times_appended = np.append(-np.flip((times_ns[1:] - times_ns[0])), times_ns)

# # # np.save('files/g2/J_dimer/J_strongC_2nm_mav_k10.npy', g2_corr)


# # plt.plot(np.append(-np.flip((times_ns[1:] - times_ns[0])), times_ns), np.append(np.flip(g2_corr[1:]), g2_corr))
# # # plt.plot(np.append(-np.flip((times_ns[1:] - times_ns[0])), times_ns), g2_corr_appended)
# # plt.grid(True)
# # plt.show()



