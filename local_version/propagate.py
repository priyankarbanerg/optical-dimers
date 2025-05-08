
from sys_params import *
from system_ops import *
from br_tensor import *
from spectral_dens import *
from analysis import *
from matplotlib.ticker import MaxNLocator

t0 = time.time()


# Initialize constants and system parameters
constants = Constants()
dipole_params = DipoleParameters(config_name)
system_params = SystemParameters(dipole_params=dipole_params, E_q1=E_q1 * constants.eV, E_q2=E_q2 * constants.eV)
system_params.calculate_parameters()


# # Calculate spectral densities and rates
# J_opt, J_vib = calculate_spectral_densities()
rates = CombinedRates(T_opt, T_vib, mode)

Gamma_opt = rates.Opt_rate()
Gamma_vib = rates.Vib_rate()
kappa = rates.kappa()
Gamma_coup = rates.Coup_rate()
# print(np.shape(Gamma_coup))
P_vib = rates.P_weight()
print(kappa)
# Hamiltonian
H = hamiltonian(mode, kappa)


# # Constants
# hbar = 1.0545718e-34  # J s
# hbar_ps = hbar * 1.0e12  # J-ps

# # Frequency in Hz
# f_hz = 1/(1.0e-12/3.5)  # Hz

# # Angular frequency in rad/s
# omega = 2 * 3.141592653589793 * f_hz  # rad/s

# # Convert omega to rad/ps
# omega_ps = omega * 1.0e-12  # rad/ps
# print(system_params.J*kappa**2)
# print(omega_ps)


sx1, sx2, sz1, sz2, sp1, sp2, sm1, sm2 = two_qubit_pauli_matrices()

a_ops_vib = [sz1, sz2]
a_ops_coup = [sp1 @ sm2, sm1 @ sp2]
a_ops_opt = [sp1, sm1, sp2, sm2]

# # initial state in the system basis
psi0 = np.array([[1, 0, 0, 0]], dtype = complex)  # 
# psi0 = np.array([[0, 1, 1, 0] / np.sqrt(2)], dtype = complex) # Dimer

# initial state as density operator
rho_sz = psi2rho(psi0)

H_diag, rho0 = initial_state(H, psi0)

# Define time points and steps for the simulation
tf = 5 * system_params.tau_L
steps = 10000
times = np.linspace(0, tf, steps)
times_ns = np.linspace(0, 1.0e-3 * tf, steps)


'''Population'''

# # Create a Bloch-Redfield calculator instance
# BR_calculator = BlochRedfieldCalculator(H, a_ops_opt, Gamma_opt, rates, a_ops_vib = a_ops_vib, Gamma_vib = Gamma_vib, a_ops_coup = a_ops_coup, Gamma_coup = Gamma_coup, P_vib = P_vib, mode=mode, pump = pump)

# # # Initialize the QuantumSystemAnalysis class for calculating population
# # state = 'doubly excited' 
# # quantum_analysis_population = QuantumSystemAnalysis(BR_calculator, rates,quantity='population', rho=rho0, final_time=tf, steps=steps, state=state)

# # # Compute population
# # pop_results_de = quantum_analysis_population.compute()
# # # np.save('J_weakC_bright_Tv0.npy', pop_results_br)

# # Initialize the QuantumSystemAnalysis class for calculating population
# state = 'bright' 
# quantum_analysis_population = QuantumSystemAnalysis(BR_calculator, rates,quantity='population', rho=rho0, final_time=tf, steps=steps, state=state)

# # Compute population
# pop_results_br = quantum_analysis_population.compute()
# # # np.save('J_weakC_bright_Tv0.npy', pop_results_br)

# # Initialize the QuantumSystemAnalysis class for calculating population
# state = 'dark' 
# quantum_analysis_population = QuantumSystemAnalysis(BR_calculator, rates,quantity='population', rho=rho0, final_time=tf, steps=steps, state=state)

# # Compute population
# pop_results_d = quantum_analysis_population.compute()
# # np.save('J_weakC_dark_Tv0.npy', pop_results_d)

# # # Initialize the QuantumSystemAnalysis class for calculating population
# # state = 'ground' 
# # quantum_analysis_population = QuantumSystemAnalysis(BR_calculator, rates,quantity='population', rho=rho0, final_time=tf, steps=steps, state=state)

# # # Compute population
# # pop_results_g = quantum_analysis_population.compute()

# # plt.plot(system_params.gamma_q1 * times, pop_results_de, 'b')
# plt.plot(system_params.gamma_q1 * times, pop_results_br, 'r')
# plt.plot(system_params.gamma_q1 * times, pop_results_d, 'g')
# # plt.plot(system_params.gamma_q1 * times, pop_results_g, 'k')
# # times = np.linspace(0, tf, 20000)
# # plt.xlim([0, 5])
# # plt.ylim([1.0e-6, 1])
# plt.show()



'''Spectra'''
# # Create a Bloch-Redfield calculator instance
# BR_calculator = BlochRedfieldCalculator(H, a_ops_opt, Gamma_opt, rates, a_ops_vib = a_ops_vib, Gamma_vib = Gamma_vib, a_ops_coup = a_ops_coup, Gamma_coup = Gamma_coup, P_vib = P_vib, mode=mode, pump = pump)

# # Initialize the QuantumSystemAnalysis class for calculating intensity
# quantum_analysis_spectra = QuantumSystemAnalysis(BR_calculator, rates, quantity='spectra', rho=rho0, final_time=tf, steps=steps)

# # Compute intensity
# spectra_results, omega = quantum_analysis_spectra.compute()
# # omega = omega * system_params.hbar_ps / constants.eV
# # spectra_results = quantum_analysis_spectra.compute()
# # np.save('files/spectra/H_strongC_2nm_wosb.npy', spectra_results)

# plt.semilogy(omega, spectra_results)
# # plt.semilogy(omega, np.load('files/spectra/H_strongC_2nm_wsb.npy'))
# # plt.ylim([1.0e-5, 1])
# plt.xlim([0, 2])
# plt.show()

'''Intensity evolution'''
# # Create a Bloch-Redfield calculator instance
# BR_calculator = BlochRedfieldCalculator(H, a_ops_opt, Gamma_opt, rates, a_ops_vib = a_ops_vib, Gamma_vib = Gamma_vib, a_ops_coup = a_ops_coup, Gamma_coup = Gamma_coup, P_vib = P_vib, mode=mode, pump = pump)

# # Initialize the QuantumSystemAnalysis class for calculating intensity
# quantum_analysis_intensity = QuantumSystemAnalysis(BR_calculator, rates,quantity='intensity', rho=rho0, final_time=tf, steps=steps)

# # Compute intensity
# intensity_results = quantum_analysis_intensity.compute()

# np.save('files/distinguishable_strongC.npy', intensity_results)

# plt.plot(times, intensity_results)
# plt.show()

'''Intensity distribution'''
# # Define parameters for g2 plot
# n_theta_p = 10
# n_phi_p = 10

# # Define parameters
# theta_vals_p = np.linspace(0, np.pi, n_theta_p)
# phi_vals_p = np.linspace(0,  2*np.pi, n_phi_p)
# r_p = 1  # radius of sphere

# # Generate x, y, z coordinates of points on the sphere
# x = r_p * np.outer(np.cos(phi_vals_p), np.sin(theta_vals_p))
# y = r_p * np.outer(np.sin(phi_vals_p), np.sin(theta_vals_p))
# z = r_p * np.outer(np.ones(np.size(phi_vals_p)), np.cos(theta_vals_p))


# # Create a Bloch-Redfield calculator instance
# BR_calculator = BlochRedfieldCalculator(H, a_ops_opt, Gamma_opt, rates, a_ops_vib = a_ops_vib, Gamma_vib = Gamma_vib, a_ops_coup = a_ops_coup, Gamma_coup = Gamma_coup, P_vib = P_vib, mode=mode, pump = pump)

# # Initialize the DistributionAnalysis class for calculating intensity
# quantum_analysis_intensity_dist = DistributionAnalysis(BR_calculator, rates,quantity='intensity distribution', rho=rho0, final_time=tf, steps=steps)

# # Compute intensity distribution
# norm_int = np.zeros((n_phi_p, n_theta_p))
# for l, phi_p in enumerate(phi_vals_p):
#     for k, theta_p in enumerate(theta_vals_p):
#         # theta_p = np.pi / 2
#         norm_int[l, k] = quantum_analysis_intensity_dist.compute(theta_p, phi_p)

# # print(norm_int.max())
# # # plt.plot(system_params.gamma_q1 * times, np.mean(norm_int, axis = (0, 1)))
# # # plt.show()


'''k-dependent intensity evolution'''

# if config_name == 'H':
#     perp = [0.0, 1.0, 0.0]
# elif config_name == 'J':
#     perp = [0.0, 1.0, 0.0]
# elif config_name == 'ortho':
#     perp = [0.0, 0.0, 1.0]
#     # perp = [1.0, 0.0, 0.0]
#     # perp = [-1.0, 1.0, 0.0]    

# r, theta, phi = cart2sph(perp)

# tf = 5 * system_params.tau_L
# steps = 10000
# times = np.linspace(0, tf, steps)
# times_ns = np.linspace(0, 1.0e-3 * tf, steps)


# # Create a Bloch-Redfield calculator instance
# BR_calculator = BlochRedfieldCalculator(H, a_ops_opt, Gamma_opt, rates, a_ops_vib = a_ops_vib, Gamma_vib = Gamma_vib, a_ops_coup = a_ops_coup, Gamma_coup = Gamma_coup, P_vib = P_vib, mode=mode, pump = pump)

# # Initialize the DistributionAnalysis class for calculating intensity
# quantum_analysis_intensity_dist = DistributionAnalysis(BR_calculator, rates,quantity='intensity distribution', rho=rho0, final_time=tf, steps=steps)

# norm_int = quantum_analysis_intensity_dist.compute(theta, phi)

# plt.plot(times_ns, norm_int)
# # plt.plot(times_ns, np.exp(-times*system_params.gamma_q1), '--')
# plt.show()


'''g2 evolution'''

# if config_name == 'H':
#     perp = [0.0, 1.0, 0.0]
# elif config_name == 'J':
#     perp = [0.0, 1.0, 0.0]
# elif config_name == 'ortho':
#     # perp = [0.0, 0.0, 1.0]
#     # perp = [1.0, 0.0, 0.0]
#     perp = [1.0/np.sqrt(2), -1.0/np.sqrt(2), 0.0]
# elif config_name == '45':
#     # perp = [1.0, 1.0e-7, 0.0]
#     # perp = [-1.0, 1.0, 0.0]    
#     perp = [0.0, 0.0, 1.0]
#     # perp = [1.0, -1.0, 0.0]    

# r, theta, phi = cart2sph(perp)
# # print(theta, phi)
# # r_p, theta_p, phi_p = cart2sph([0.25, 0.25, 0.0]) #r, theta, phi
# r_p, theta_p, phi_p = r, theta, phi


# tf = 5 * system_params.tau_L
# steps = 5000000#len(np.load('files/g2/45_strongC_2nm_symm_d1.npy')) #5000000
# times = np.linspace(0, tf, steps)
# times_ns = np.linspace(0, 1.0e-3 * tf, steps)


# # Create a Bloch-Redfield calculator instance
# BR_calculator = BlochRedfieldCalculator(H, a_ops_opt, Gamma_opt, rates, a_ops_vib = a_ops_vib, Gamma_vib = Gamma_vib, a_ops_coup = a_ops_coup, Gamma_coup = Gamma_coup, P_vib = P_vib, mode=mode, pump = pump)

# # Initialize the DistributionAnalysis class for calculating g2
# quantum_analysis_g2_dist = DistributionAnalysis(BR_calculator, rates, quantity='g2 distribution', rho=rho0, final_time=tf, steps=steps)

# g2_corr = quantum_analysis_g2_dist.compute(theta, phi, theta_p, phi_p)

# # # # # times_ns = np.linspace(0, 1.0e-3 * tf, int(len(g2_corr)))
# # # # print(g2_corr[0], g2_corr[-1])
# # np.save('files/g2/45_strongC_2nm_symm_z.npy', g2_corr)

# # dt = times_ns[1] - times_ns[0]
# # # Compute FFT with zero-padding
# # spec = 2 * fft(np.real(np.load('files/g2/45_strongC_2nm_symm.npy'))) * dt 
# # spec = fftshift(spec)

# # # Compute new frequency axis
# # wlist = 2 * np.pi * fftfreq(int(len(np.load('files/g2/45_strongC_2nm_symm.npy'))), dt)
# # wlist = fftshift(wlist)
        

# # plt.semilogy(wlist, spec)
# # # plt.semilogy(omega, np.load('files/spectra/H_strongC_2nm_wsb.npy'))
# # # plt.ylim([1.0e-5, 1])
# # # plt.xlim([0, 2])
# # plt.show()


# # np.array([])
# # # Create 100 points between [1.0, 1.0, 0.0] and [-1.0, 1.0, 0.0]
# # perp_start = np.array([1.0, 1.0, 0.0])
# # perp_end = np.array([-1.0, 1.0, 0.0])
# # perp_points = np.linspace(perp_start, perp_end, 100)


# # # Create directory to store plots
# # output_dir = "frames"
# # os.makedirs(output_dir, exist_ok=True)

# # # Function to compute and save a single frame
# # def compute_and_plot_frame(args):
# #     i, perp = args
# #     r, theta, phi = cart2sph(perp)
# #     r_p, theta_p, phi_p = r, theta, phi

# #     # Compute g2 correlation
# #     g2_corr = quantum_analysis_g2_dist.compute(theta, phi, theta_p, phi_p)

# #     # Plot g2 correlation
# #     fig, ax = plt.subplots(figsize=(8, 6))
# #     ax.plot(
# #         np.append(-np.flip((times_ns[1:] - times_ns[0])), times_ns),
# #         np.append(np.flip(g2_corr[1:]), g2_corr),
# #         label=f"Frame {i+1}"
# #     )

# #     ax.set_ylim(0.0, 0.27)
# #     ax.grid(True)
# #     ax.set_title(f"g2 Correlation for Frame {i+1}")
# #     ax.set_xlabel("Time (ns)")
# #     ax.set_ylabel("g2 Correlation")
# #     ax.legend()

# #     # Add polar plot inset
# #     phi_deg = np.rad2deg(phi)  # Convert phi to degrees
# #     theta_deg = np.rad2deg(theta)  # Convert theta to degrees

# #     inset_ax = fig.add_axes([0.6, 0.6, 0.25, 0.25], polar=True)  # Create polar inset
# #     inset_ax.set_thetamin(45)  # Start angle: π/4 (in degrees)
# #     inset_ax.set_thetamax(135)  # End angle: 3π/4 (in degrees)

# #     # Plot a dashed line and the red dot for phi and theta
# #     inset_ax.plot([np.deg2rad(45), phi], [0, theta], linestyle="--", color="gray")  # Convert 45° to radians
# #     inset_ax.plot(phi, theta, "ro", label=r"$\theta, \phi$")  # Use LaTeX for the label

# #     # Set radial limits for theta
# #     inset_ax.set_ylim(0, np.pi)  # Theta from 0 to π
# #     inset_ax.set_yticks([0, np.pi / 2, np.pi])  # Tick marks at 0, π/2, π
# #     inset_ax.set_yticklabels([r"$0$", r"$\pi/2$", r"$\pi$"], fontsize=10)  # Use LaTeX for labels

# #     # Convert theta ticks to degrees for display
# #     inset_ax.set_xticks(np.radians([45, 60, 75, 90, 105, 120, 135]))  # Ticks every 15° within the range
# #     inset_ax.set_xticklabels([r"$45^\circ$", r"$60^\circ$", r"$75^\circ$", r"$90^\circ$", 
# #                             r"$105^\circ$", r"$120^\circ$", r"$135^\circ$"], fontsize=10)

# #     # Add a legend
# #     inset_ax.legend(loc="lower left", fontsize="small")

# #     # Save the frame
# #     plt.savefig(f"{output_dir}/frame_{i:03d}.png")
# #     plt.close()

# # # Prepare arguments for parallel execution
# # args = [(i, perp) for i, perp in enumerate(perp_points)]

# # # Use multiprocessing to generate frames in parallel
# # if __name__ == "__main__":
# #     num_processes = min(cpu_count(), 10)  # Adjust the number of processes if needed
# #     with Pool(num_processes) as pool:
# #         pool.map(compute_and_plot_frame, args)

# #     # Use ffmpeg to create a video
# #     os.system("ffmpeg -r 10 -i frames/frame_%03d.png -vcodec libx264 -crf 25 -pix_fmt yuv420p video_sym.mp4")




# # np.save('files/g2/ortho_strongC_2nm_symm_intermediate.npy', g2_corr)

# # g2_corr = np.load('files/g2/ortho_strongC_2nm_symm.npy')

# # # sigma = 50  # in picoseconds

# # dt = (tf - (-tf)) / steps
# # t_kernel = np.arange(-10 * sigma, 10 * sigma, dt)
# # gaussian_kernel = np.exp(-t_kernel**2 / (2 * sigma**2))
# # gaussian_kernel /= np.sum(gaussian_kernel)  

# # g2_corr_appended = fftconvolve(np.append(np.flip(g2_corr[1:]), g2_corr), gaussian_kernel, mode = "same") 
# # times_appended = np.append(-np.flip((times_ns[1:] - times_ns[0])), times_ns)

# # # print(np.shape(g2_corr_appended), np.shape(times_appended))
# # np.save('files/g2_con/ortho_strongC_2nm_symm.npy', g2_corr_appended)         #pbs = pumping bright state

# # # # g2_corr_kkpsame = np.load('files/g2/H_strongC_2nm_symm_kkpsame.npy') / np.load('files/g2/H_strongC_2nm_symm_kkpsame.npy')[-1]
# # # # g2_corr_kkpdiff = np.load('files/g2/H_strongC_2nm_symm_kkpdiff.npy') / np.load('files/g2/H_strongC_2nm_symm_kkpdiff.npy')[-1]

# # # plt.plot(np.append(-np.flip((times_ns[1:] - times_ns[0])), times_ns), np.append(np.flip(g2_corr_kkpsame[1:]), g2_corr_kkpsame))
# # plt.plot(np.append(-np.flip((times_ns[1:] - times_ns[0])), times_ns), np.append(np.flip(g2_corr_kkpdiff[1:]), g2_corr_kkpdiff))
# plt.plot(np.append(-np.flip((times_ns[1:] - times_ns[0])), times_ns), np.append(np.flip(g2_corr[1:]), g2_corr))
# # plt.plot(np.append(-np.flip(np.linspace(0, tf, 2000000)[1:] - np.linspace(0, tf, 2000000)[0]), np.linspace(0, tf, 2000000)), fftconvolve(np.append(np.flip(np.load('files/g2/H_dimer/H_strongC_T0.npy')[1:]), np.load('files/g2/H_dimer/H_strongC_T0.npy')), gaussian_kernel, mode = "same"))
# # plt.plot(np.append(-np.flip((times_ns[1:] - times_ns[0])), times_ns), g2_corr_appended)
# # plt.plot(np.append(-np.flip((times_ns[1:] - times_ns[0])), times_ns),np.load('files/g2_con/H_strongC_2nm_symm_nfd_slow_tumb_150.npy'))
# plt.grid(True)
# plt.show()


# # [23.02848026713309, 0.011424406565876627, 9.557031225474333e-08] - 0
# # [23.02843221541343, 0.011424285922020591, 4.2447492431440256e-08] - 5
# # [23.028289467676586, 0.011424084654110595, 4.203550484285212e-08] - 10
# # [23.02805395320229, 0.011423843898617746, 1.3354398463955862e-07] - 15
# # [23.02772758137371, 0.011425020353186984, 1.7722378524510584e-06] - 20
# # [23.027312214484027, 0.011432758523788403, 1.0101657382161636e-05] - 25
# # [23.02680966909611, 0.011454544922860159, 3.2607710057288345e-05] - 30
# # [23.025550088262207, 0.011562555421096883, 0.00014244135584260407] - 40
# # [23.023962505586496, 0.011767075341648467, 0.00034929460858079556] - 50
# # [23.011539226732566 0.01377383946539112, 0.002375256121222725] - 100

# # # np.array([23.02848026713309, 23.02843221541343, 23.028289467676586, 23.02805395320229, 23.02772758137371, 23.027312214484027, 23.02680966909611, 23.025550088262207, 23.023962505586496, 23.011539226732566])
# # Data
# x = np.array([0, 5, 10, 15, 20, 25, 30, 40, 50, 100, 300, 500, 600, 1000, 3000, 5000, 10000])
# y = np.array([0.011424406565876627, 0.011424285922020591, 0.011424084654110595, 
#               0.011423843898617746, 0.011425020353186984, 0.011432758523788403, 
#               0.011454544922860159, 0.011562555421096883, 0.011767075341648467, 
#               0.01377383946539112, 0.025378997060922472, 0.03719530121034322,
#               0.042891400954189585, 0.06392489187567942, 0.13067351568370383,
#               0.15305141546166268, 0.12864983071925648])

# z = np.array([9.557031225474333e-08, 4.2447492431440256e-08, 4.203550484285212e-08, 
#               1.3354398463955862e-07, 1.7722378524510584e-06, 1.0101657382161636e-05,
#               3.2607710057288345e-05, 0.00014244135584260407, 0.00034929460858079556, 
#               0.002375256121222725, 0.01416371277386932, 0.026275399432181618,
#               0.03214278588458, 0.0539283524043761, 0.12414879162344039,
#               0.14887564131652054, 0.12729218749119972])

# # Create figure and set size
# plt.figure(figsize=(8, 6))

# # Plot data
# # plt.plot(x, z/y, marker='o', linestyle='-')
# plt.plot(x, y, marker='o', linestyle='--')
# plt.plot(x, z, marker='o', linestyle='--')

# # Labels with larger font size
# plt.xlabel('Temperature (K)', fontsize=14)
# plt.ylabel(r'Phonon Induced Decay $\gamma_{pn}$ $(ps^{-1})$', fontsize=14)

# # Limits
# # plt.xlim([0, 30])
# # plt.ylim([0.01142, 0.01145])

# # Ticks with larger font size
# plt.xticks(fontsize=12)
# plt.yticks(fontsize=12)

# # Adjust layout to avoid cutting off labels
# plt.tight_layout()

# # Save and show
# # plt.savefig('rate_phonon.png', dpi=300, bbox_inches='tight')
# plt.show()
# # # g2_corr_convolved = gaussian_filter1d(g2_corr_appended, sigma)

# # Plot the downsampled data
# fig, ax = plt.subplots(figsize=(14, 10))

# # # # Plot the data
# # # # line0, = ax.plot(np.append(-np.flip(np.linspace(0, 1.0e-3 * tf, 20000)[1:] - np.linspace(0, 1.0e-3 * tf, 20000)[0]), np.linspace(0, 1.0e-3 * tf, 20000)), np.append(np.flip(np.load('files/g2/J_dimer/J_weakC_T0_50nm_wcme.npy')[1:]), np.load('files/g2/J_dimer/J_weakC_T0_50nm_wcme.npy')), '-b')
# # # # line1, = ax.plot(np.append(-np.flip(np.linspace(0, 1.0e-3 * tf, 20000)[1:] - np.linspace(0, 1.0e-3 * tf, 20000)[0]), np.linspace(0, 1.0e-3 * tf, 20000)), np.append(np.flip(np.load('files/g2/J_dimer/J_strongerC_T0_50nm_wcme.npy')[1:]), np.load('files/g2/J_dimer/J_strongerC_T0_50nm_wcme.npy')), '-r')

# # line0, = ax.plot(np.append(-np.flip(np.linspace(0, 1.0e-3 * tf, 100000)[1:] - np.linspace(0, 1.0e-3 * tf, 100000)[0]), np.linspace(0, 1.0e-3 * tf, 100000)), np.append(np.flip(np.load('files/g2/H_strongC_2nm_symm.npy')[1:]), np.load('files/g2/H_strongC_2nm_symm.npy')), '-b', linewidth = 3)
# # line1, = ax.plot(np.append(-np.flip(np.linspace(0, 1.0e-3 * tf, 100000)[1:] - np.linspace(0, 1.0e-3 * tf, 100000)[0]), np.linspace(0, 1.0e-3 * tf, 100000)), np.append(np.flip(np.load('files/g2/J_strongC_2nm_symm.npy')[1:]), np.load('files/g2/J_strongC_2nm_symm.npy')), '-r', linewidth = 3) #ax.plot(np.append(-np.flip(np.linspace(0, 1.0e-3 * tf, 20000)[1:] - np.linspace(0, 1.0e-3 * tf, 20000)[0]), np.linspace(0, 1.0e-3 * tf, 20000)), np.append(np.flip(np.load('files/g2/H_dimer/H_strongerC_T0_50nm_wcme.npy')[1:]), np.load('files/g2/H_dimer/H_strongerC_T0_50nm_wcme.npy')), '-r')
# # # line2, = ax.plot(np.append(-np.flip(np.linspace(0, 1.0e-3 * tf, 100000)[1:] - np.linspace(0, 1.0e-3 * tf, 100000)[0]), np.linspace(0, 1.0e-3 * tf, 100000)), np.append(np.flip(np.load('files/g2_con/ortho_strongC_2nm_symm.npy')[1:]), np.load('files/g2_con/ortho_strongC_2nm_symm.npy')), '--g', linewidth = 3)

# # # line0, = ax.plot(np.linspace(0, 1.0e-3 * tf, 100000), np.load('files/g2/H_dimer/H_strongC_2nm_symm.npy'), '-b', linewidth = 3)
# # # line1, = ax.plot(np.linspace(0, 1.0e-3 * tf, 100000), np.load('files/g2/J_dimer/J_strongC_2nm_symm.npy'), '-r', linewidth = 3) #ax.plot(np.append(-np.flip(np.linspace(0, 1.0e-3 * tf, 20000)[1:] - np.linspace(0, 1.0e-3 * tf, 20000)[0]), np.linspace(0, 1.0e-3 * tf, 20000)), np.append(np.flip(np.load('files/g2/H_dimer/H_strongerC_T0_50nm_wcme.npy')[1:]), np.load('files/g2/H_dimer/H_strongerC_T0_50nm_wcme.npy')), '-r')
# # # line2, = ax.plot(np.linspace(0, 1.0e-3 * tf, 100000), np.load('files/g2/ortho_strongC_2nm_symm.npy'), '--g', linewidth = 3)

# line0, = ax.plot(np.append(-np.flip(np.linspace(0, 1.0e-3 * tf, 10000)[1:] - np.linspace(0, 1.0e-3 * tf, 10000)[0]), np.linspace(0, 1.0e-3 * tf, 10000)), np.append(np.flip(np.load('files/g2/ortho_strongC_2nm_symm_perp.npy')[1:]), np.load('files/g2/ortho_strongC_2nm_symm_perp.npy')), '-b', linewidth = 3)
# line1, = ax.plot(np.append(-np.flip(np.linspace(0, 1.0e-3 * tf, 10000)[1:] - np.linspace(0, 1.0e-3 * tf, 10000)[0]), np.linspace(0, 1.0e-3 * tf, 10000)), np.append(np.flip(np.load('files/g2/ortho_strongC_2nm_symm_par.npy')[1:]), np.load('files/g2/ortho_strongC_2nm_symm_par.npy')), '-r', linewidth = 3) #ax.plot(np.append(-np.flip(np.linspace(0, 1.0e-3 * tf, 20000)[1:] - np.linspace(0, 1.0e-3 * tf, 20000)[0]), np.linspace(0, 1.0e-3 * tf, 20000)), np.append(np.flip(np.load('files/g2/H_dimer/H_strongerC_T0_50nm_wcme.npy')[1:]), np.load('files/g2/H_dimer/H_strongerC_T0_50nm_wcme.npy')), '-r')
# line2, = ax.plot(np.append(-np.flip(np.linspace(0, 1.0e-3 * tf, 10000)[1:] - np.linspace(0, 1.0e-3 * tf, 10000)[0]), np.linspace(0, 1.0e-3 * tf, 10000)), np.append(np.flip(np.load('files/g2/ortho_strongC_2nm_symm_intermediate.npy')[1:]), np.load('files/g2/ortho_strongC_2nm_symm_intermediate.npy')), '-g', linewidth = 3)
# line3, = ax.plot(np.append(-np.flip(np.linspace(0, 1.0e-3 * tf, 10000)[1:] - np.linspace(0, 1.0e-3 * tf, 10000)[0]), np.linspace(0, 1.0e-3 * tf, 10000)), np.append(np.flip(np.load('files/g2/ortho_strongC_2nm_symm_intermediate_neg.npy')[1:]), np.load('files/g2/ortho_strongC_2nm_symm_intermediate_neg.npy')), 'purple', linewidth = 3)


# # # # line0, = ax.plot(np.append(-np.flip(np.linspace(0, 1.0e-3 * tf, 20000)[1:] - np.linspace(0, 1.0e-3 * tf, 20000)[0]), np.linspace(0, 1.0e-3 * tf, 20000)), np.append(np.flip(np.load('files/g2/ortho_dimer/ortho_weakC_T0_50nm_wcme.npy')[1:]), np.load('files/g2/ortho_dimer/ortho_weakC_T0_50nm_wcme.npy')), '-b')
# # # # line1, = ax.plot(np.append(-np.flip(np.linspace(0, 1.0e-3 * tf, 20000)[1:] - np.linspace(0, 1.0e-3 * tf, 20000)[0]), np.linspace(0, 1.0e-3 * tf, 20000)), np.append(np.flip(np.load('files/g2/ortho_dimer/ortho_strongerC_T0_50nm_wcme.npy')[1:]), np.load('files/g2/ortho_dimer/ortho_strongerC_T0_50nm_wcme.npy')), '-r')


# # # line0, = ax.plot(times_appended, np.load('files/g2_con/H_dimer/H_strongC_2nm_symm.npy'), '-b', linewidth = 3)
# # # line1, = ax.plot(times_appended, np.load('files/g2_con/J_dimer/J_strongC_2nm_symm.npy'), '-r', linewidth = 3) #ax.plot(np.append(-np.flip(np.linspace(0, 1.0e-3 * tf, 20000)[1:] - np.linspace(0, 1.0e-3 * tf, 20000)[0]), np.linspace(0, 1.0e-3 * tf, 20000)), np.append(np.flip(np.load('files/g2/H_dimer/H_strongerC_T0_50nm_wcme.npy')[1:]), np.load('files/g2/H_dimer/H_strongerC_T0_50nm_wcme.npy')), '-r')
# # # line2, = ax.plot(times_appended, np.load('files/g2_con/ortho_strongC_2nm_symm.npy'), '--g', linewidth = 3)


# # # line0, = ax.plot(times_ns, np.load('files/g2/crystal/H_dimer/H_weakC_2nm_symm_unnorm.npy'), '-b', linewidth = 3)
# # # line1, = ax.plot(times_ns, np.load('files/g2/crystal/H_dimer/H_strongC_2nm_symm_fd_unnorm.npy'), '-r', linewidth = 3)


# # # line0, = ax.plot(times_appended, np.load('files/g2_con/H_strongC_2nm_symm_nfd_slow_tumb_150.npy'), '-b', linewidth = 3)
# # # line1, = ax.plot(times_appended, np.load('files/g2_con/J_strongC_2nm_symm_nfd_slow_tumb_150.npy'), '-r', linewidth = 3)

# # Set axis labels and limits
# ax.set_xlabel(r'$\tau$ (ns)', fontsize=40)
# ax.set_ylabel(r'$g^{(2)}(\infty, \tau)$', fontsize=40)
# # ax.set_xlim([-1.0e-2, 1.0e-2])  
# ax.set_xlim([-40, 40]) 
# # ax.set_ylim([0.433, 0.436]) 
# ax.tick_params(axis='x', labelsize=35, direction='in', length=6)
# ax.tick_params(axis='y', labelsize=35, direction='in', length=6)

# # ax.set_xscale('log')
# # Add a legend with appropriate styling
# # ax.legend(handles=[line0, line1, line2, line3], labels=[r"$\theta = 0$", r"$\theta = \pi/2, \phi = n\pi/2$" , r"$\theta = \pi/2$, $\phi = \pi/4$", r"$\theta = \pi/2$, $\phi = 3\pi/4$"], loc='best', fontsize=30, frameon=False)
# # ax.legend(handles=[line0, line1, line2, line3], labels=[r"$\textbf{q} = (0, 0, \pm 1)$", r"$\textbf{q} = (\pm 1, 0, 0), (0, \pm 1, 0)$" , r"$\textbf{q} = (\pm 1, \pm 1, 0)$", r"$\textbf{q} = (\pm 1, \mp 1, 0)$"], loc='best', fontsize=28, frameon=False)
# # ax.legend(handles=[line0, line1], labels=[r"H dimer", r"J dimer"], fontsize=30, frameon=False)

# # # ax.text(-0.1, 1.05, r'(b)', size=30, transform=ax.transAxes)  # Adjusted position to upper-left

# # # Inset plot
# # ax_inset = plt.axes([0.63, 0.21, 0.24, 0.25])  # [left, bottom, width, height]

# # # # ax_inset.plot(np.append(-np.flip(np.linspace(0, 1.0e-3 * tf, 20000)[1:] - np.linspace(0, 1.0e-3 * tf, 20000)[0]), np.linspace(0, 1.0e-3 * tf, 20000)), np.append(np.flip(np.load('files/g2/J_dimer/J_weakC_T0_50nm_wcme.npy')[1:]), np.load('files/g2/J_dimer/J_weakC_T0_50nm_wcme.npy')), '-b')
# # # # ax_inset.plot(np.append(-np.flip(np.linspace(0, 1.0e-3 * tf, 20000)[1:] - np.linspace(0, 1.0e-3 * tf, 20000)[0]), np.linspace(0, 1.0e-3 * tf, 20000)), np.append(np.flip(np.load('files/g2/J_dimer/J_strongerC_T0_50nm_wcme.npy')[1:]), np.load('files/g2/J_dimer/J_strongerC_T0_50nm_wcme.npy')), '-r')

# # line0, = ax_inset.plot(np.append(-np.flip(np.linspace(0, 1.0e-3 * tf, 100000)[1:] - np.linspace(0, 1.0e-3 * tf, 100000)[0]), np.linspace(0, 1.0e-3 * tf, 100000)), np.append(np.flip(np.load('files/g2/H_strongC_2nm_symm.npy')[1:]), np.load('files/g2/H_strongC_2nm_symm.npy')), '-b', linewidth = 3)
# # line1, = ax_inset.plot(np.append(-np.flip(np.linspace(0, 1.0e-3 * tf, 100000)[1:] - np.linspace(0, 1.0e-3 * tf, 100000)[0]), np.linspace(0, 1.0e-3 * tf, 100000)), np.append(np.flip(np.load('files/g2/J_strongC_2nm_symm.npy')[1:]), np.load('files/g2/J_strongC_2nm_symm.npy')), '-r', linewidth = 3) #ax.plot(np.append(-np.flip(np.linspace(0, 1.0e-3 * tf, 20000)[1:] - np.linspace(0, 1.0e-3 * tf, 20000)[0]), np.linspace(0, 1.0e-3 * tf, 20000)), np.append(np.flip(np.load('files/g2/H_dimer/H_strongerC_T0_50nm_wcme.npy')[1:]), np.load('files/g2/H_dimer/H_strongerC_T0_50nm_wcme.npy')), '-r')
# # # line2, = ax_inset.plot(np.append(-np.flip(np.linspace(0, 1.0e-3 * tf, 100000)[1:] - np.linspace(0, 1.0e-3 * tf, 100000)[0]), np.linspace(0, 1.0e-3 * tf, 100000)), np.append(np.flip(np.load('files/g2_con/ortho_strongC_2nm_symm.npy')[1:]), np.load('files/g2_con/ortho_strongC_2nm_symm.npy')), '--g', linewidth = 3)

# # # # line0, = ax_inset.plot(times_appended, np.load('files/g2_con/H_dimer/H_strongC_2nm_symm.npy'), '-b', linewidth = 3)
# # # # line1, = ax_inset.plot(times_appended, np.load('files/g2_con/J_dimer/J_strongC_2nm_symm.npy'), '-r', linewidth = 3) #ax.plot(np.append(-np.flip(np.linspace(0, 1.0e-3 * tf, 20000)[1:] - np.linspace(0, 1.0e-3 * tf, 20000)[0]), np.linspace(0, 1.0e-3 * tf, 20000)), np.append(np.flip(np.load('files/g2/H_dimer/H_strongerC_T0_50nm_wcme.npy')[1:]), np.load('files/g2/H_dimer/H_strongerC_T0_50nm_wcme.npy')), '-r')
# # # # line2, = ax_inset.plot(times_appended, np.load('files/g2_con/ortho_strongC_2nm_symm.npy'), '--g', linewidth = 3)


# # # # line0, = ax_inset.plot(times_appended, np.load('files/g2_con/H_dimer/H_strongC_2nm.npy'), '-b', linewidth = 3)
# # # # line1, = ax_inset.plot(times_appended, np.load('files/g2_con/J_dimer/J_strongC_2nm.npy'), '-r', linewidth = 3) #ax.plot(np.append(-np.flip(np.linspace(0, 1.0e-3 * tf, 20000)[1:] - np.linspace(0, 1.0e-3 * tf, 20000)[0]), np.linspace(0, 1.0e-3 * tf, 20000)), np.append(np.flip(np.load('files/g2/H_dimer/H_strongerC_T0_50nm_wcme.npy')[1:]), np.load('files/g2/H_dimer/H_strongerC_T0_50nm_wcme.npy')), '-r')

 
# # ax_inset.set_xlim([-0.25, 0.25])
# # ax_inset.set_ylim([0.4, 1.25])
# # ax_inset.set_xlabel(r'$\tau$ (ns)', fontsize=27)
# # ax_inset.set_ylabel(r'$g^{(2)}(\infty, \tau)$', fontsize=27)
# # ax_inset.tick_params(axis='x', labelsize=25, direction='in', length=6)
# # ax_inset.tick_params(axis='y', labelsize=25, direction='in', length=6)
# # ax_inset.yaxis.set_major_locator(MaxNLocator(nbins=4))


# plt.savefig('orient_g2_ortho', dpi=300, bbox_inches='tight')
# plt.show()

g2_corr = np.load('files/g2/ortho_strongC_2nm_symm_fast_tumb_jigg_normalized_smart.npy')

zero_delay_convolved_array = []
step_con = 1.0e4

# Generate the instrument response
instrument_responses = np.linspace(0, 300, 10000)
for instrument_response in instrument_responses:
    sigma = instrument_response

    dt = (5 * system_params.tau_L - (-5 * system_params.tau_L)) / step_con
    g2_corr_appended = gaussian_filter1d(np.append(np.flip(g2_corr[1:]), g2_corr), sigma=sigma / dt, mode='reflect')

    # g2_corr_appended = fftconvolve(np.append(np.flip(g2_corr[1:]), g2_corr), gaussian_kernel, mode="same")
    times_appended = np.append(-np.flip((times_ns[1:] - times_ns[0])), times_ns)
    

    # Append the zero-delay value for this convolved g2
    zero_delay_convolved_array.append(np.max(g2_corr_appended[steps-10:steps+10]))  # Adjusted for zero-delay index
    # plt.plot(times_appended[steps-500:steps+500], g2_corr_appended[steps-500:steps+500])
    # plt.show()
np.save('files/instrument_responses/ortho_strongC_2nm_symm_fast_tumb_jigg_convolved_smart_k5.npy', zero_delay_convolved_array)
# zero_delay_convolved_array = np.load('files/instrument_responses/ortho_strongC_2nm_symm_fast_tumb_jigg_convolved_smart.npy')
plt.plot(instrument_responses, zero_delay_convolved_array)
plt.show()

# # # Define instrument responses
# # # instrument_responses = np.linspace(1, 250, 100)

# instrument_responses = np.linspace(1, 300, 10000)

# # Load data

# # Fast tumbling without any wiggling would give the same result as the static case because, the photon correlation doesn't change based on the direction of measurement for parallel dipoles.


# data_H_nfd_slow = np.load('files/instrument_responses/H_strongC_2nm_symm_slow_tumb_nojigg_convolved_smart.npy')
# data_H_nfd_slow_jigg_k5 = np.load('files/instrument_responses/H_strongC_2nm_symm_slow_tumb_jigg_convolved_smart_k5.npy')
# data_H_nfd_slow_jigg_k10 = np.load('files/instrument_responses/H_strongC_2nm_symm_slow_tumb_jigg_convolved_smart_k10.npy')
# data_H_nfd_fast_jigg_k5 = np.load('files/instrument_responses/H_strongC_2nm_symm_fast_tumb_jigg_convolved_smart_k5.npy')
# # data_H_fd_slow = np.load('files/H_strongC_2nm_symm_fd_slow_tumb_convolved.npy')

# data_J_nfd_slow = np.load('files/instrument_responses/J_strongC_2nm_symm_slow_tumb_nojigg_convolved_smart.npy')
# data_J_nfd_slow_jigg_k5 = np.load('files/instrument_responses/J_strongC_2nm_symm_slow_tumb_jigg_convolved_smart_k5.npy')
# data_J_nfd_slow_jigg_k10 = np.load('files/instrument_responses/J_strongC_2nm_symm_slow_tumb_jigg_convolved_smart_k10.npy')
# data_J_nfd_fast_jigg_k5 = np.load('files/instrument_responses/J_strongC_2nm_symm_fast_tumb_jigg_convolved_smart_k5.npy')
# # data_J_fd_slow = np.load('files/J_weakC_2nm_symm_fd_slow_tumb_convolved.npy')

# # Colors for H dimer (shades of blue)
# colors_H = ['lightsteelblue', 'royalblue', 'mediumblue', 'midnightblue']  # Variations of blue
# colors_J = ['rosybrown', 'lightcoral', 'red', 'maroon']  # Variations of red

# # Updated line styles for more differentiation
# line_styles = ['-', '--', (0, (3, 1, 1, 1)), ':']  # Solid, dashed, dash-dot, dotted


# # Load g2 data
# g2_corr_J = np.load('files/g2/J_strongC_2nm_symm_slow_tumb_jigg_normalized_smart_k5.npy')
# g2_corr_H = np.load('files/g2/H_strongC_2nm_symm_slow_tumb_jigg_normalized_smart_k5.npy')

# # Inset configurations
# insets_H = [
#     {'ir': 50, 'pos': [0.2, 0.72, 0.21, 0.24]}, 
#     {'ir': 100, 'pos': [0.09, 0.12, 0.21, 0.24]},
#     {'ir': 200, 'pos': [0.42, 0.03, 0.21, 0.24]}
# ]

# insets_J = [
#     {'ir': 50, 'pos': [0.18, 0.68, 0.21, 0.24]},
#     {'ir': 100, 'pos': [0.09, 0.08, 0.21, 0.24]},
#     {'ir': 200, 'pos': [0.42, 0.05, 0.21, 0.24]}
# ]


# def add_inset(ax, g2_corr, instrument_response, inset_pos, x_lim, y_lim, color=None):
#     step_con = 1.0e4
#     dt = (5 * system_params.tau_L - (-5 * system_params.tau_L)) / step_con
#     sigma = instrument_response

#     g2_corr_appended = gaussian_filter1d(
#         np.append(np.flip(g2_corr[1:]), g2_corr), 
#         sigma=sigma / dt, 
#         mode='reflect'
#     )
#     times_appended = np.append(-np.flip((times_ns[1:] - times_ns[0])), times_ns)

#     inset_ax = ax.inset_axes(inset_pos)
#     inset_ax.plot(times_appended, g2_corr_appended, color=color)
#     inset_ax.set_xlim(x_lim)
#     inset_ax.set_ylim(y_lim)

#     inset_ax.set_xlabel(r'$\tau$ (ns)', fontsize=18)
#     inset_ax.set_ylabel(r'$g^{(2)}(\infty, \tau)$', fontsize=18)
#     # Add y-ticks and two numerical labels
#     inset_ax.set_yticks([y_lim[0], y_lim[1]])  
#     inset_ax.set_yticklabels([f'{y_lim[0]:.2f}', f'{y_lim[1]:.2f}'], fontsize=18)

#     # Hide x-ticks
#     inset_ax.set_xticklabels([])
#     inset_ax.tick_params(axis='x', which='both', length=0, labelbottom=False)

#     inset_ax.text(0.5, 0.95, rf'$\Delta \tau = {instrument_response}\,\mathrm{{ps}}$', 
#               ha='center', va='top', transform=inset_ax.transAxes, fontsize=18)



# # Create a single figure with two subplots, stacked vertically, sharing the x-axis
# fig, (ax_H, ax_J) = plt.subplots(2, 1, figsize=(9.5, 13), sharex=True, gridspec_kw={'hspace': 0.05, 'height_ratios': [1, 1]})

# # Plot for H dimer
# ax_H.plot(instrument_responses, data_H_nfd_slow, color=colors_H[0], linestyle=line_styles[0], label=r'$[g^{(2)}(\infty, 0)]$', linewidth=3)
# ax_H.plot(instrument_responses, data_H_nfd_slow_jigg_k10, color=colors_H[1], linestyle=line_styles[1], label=r'$[g^{(2)}(\infty, 0)]_{O}^{\kappa = 10}$', linewidth=3)
# ax_H.plot(instrument_responses, data_H_nfd_slow_jigg_k5, color=colors_H[2], linestyle=line_styles[2], label=r'$[g^{(2)}(\infty, 0)]_{O}^{\kappa = 5}$', linewidth=3)
# ax_H.plot(instrument_responses, data_H_nfd_fast_jigg_k5, color=colors_H[3], linestyle=line_styles[3], label=r'$[g^{(2)}(\infty, 0)]_{(O,\textbf{q})}^{\kappa = 5}$', linewidth=3)

# # Set title and legend for H dimer
# ax_H.legend(
#     loc='lower center', 
#     bbox_to_anchor=(0.6, 0.45),  # Centered horizontally, slightly above the plot
#     fontsize=20, 
#     frameon=False, 
#     handletextpad=1, 
#     borderaxespad=0.5, 
#     labelspacing=1.5, 
#     ncol=1 
# )
# ax_H.tick_params(axis='x', labelsize=20, direction='in', length=6)
# ax_H.tick_params(axis='y', labelsize=20, direction='in', length=6)
# ax_H.set_ylim([0.3, 1.4])

# # Add insets for H dimer
# for inset in insets_H:
#     add_inset(ax_H, g2_corr_H, inset['ir'], inset['pos'], [-2, 2], [0.43, 0.9], color=colors_H[3])


# # Plot for J dimer
# ax_J.plot(instrument_responses, data_J_nfd_slow, color=colors_J[0], linestyle=line_styles[0], label=r'$[g^{(2)}(\infty, 0)]$', linewidth=3)
# ax_J.plot(instrument_responses, data_J_nfd_slow_jigg_k10, color=colors_J[1], linestyle=line_styles[1], label=r'$[g^{(2)}(\infty, 0)]_{O}^{\kappa = 10}$', linewidth=3)
# ax_J.plot(instrument_responses, data_J_nfd_slow_jigg_k5, color=colors_J[2], linestyle=line_styles[2], label=r'$[g^{(2)}(\infty, 0)]_{O}^{\kappa = 5}$', linewidth=3)
# ax_J.plot(instrument_responses, data_J_nfd_fast_jigg_k5, color=colors_J[3], linestyle=line_styles[3], label=r'$[g^{(2)}(\infty, 0)]_{(O,\textbf{q})}^{\kappa = 5}$', linewidth=3)

# # Set title and legend for J dimer
# ax_J.legend(
#     loc='lower center', 
#     bbox_to_anchor=(0.6, 0.45),  # Centered horizontally, slightly above the plot
#     fontsize=20, 
#     frameon=False, 
#     handletextpad=1, 
#     borderaxespad=0.5, 
#     labelspacing=1.5, 
#     ncol=1
# )
# fig.text(0.5, 0.04, r'Instrument response $\Delta \tau$ (ps)', ha='center', fontsize=25)
# ax_J.tick_params(axis='x', labelsize=20, direction='in', length=6)
# ax_J.tick_params(axis='y', labelsize=20, direction='in', length=6)
# ax_J.set_ylim([0.48, 0.9])

# # Add insets for J dimer
# for inset in insets_J:
#     add_inset(ax_J, g2_corr_J, inset['ir'], inset['pos'], [-2, 2], [0.6, 0.72], color=colors_J[3])

# # Add a shared y-axis label between the plots
# fig.text(0.03, 0.5, r'Zero delay coincidence', va='center', rotation='vertical', fontsize=25)
# plt.subplots_adjust(left=0.12, right=0.96, top=0.99, bottom=0.1)

# ax_H.set_yticks([0.4, 0.7, 1.0, 1.3])  # For H dimer plot
# ax_J.set_yticks([0.55, 0.65, 0.75, 0.85])  # For J dimer plot

# # Save the combined plot
# plt.savefig('zero_time_g2_H_J', dpi=300, bbox_inches='tight')
# plt.show()






# # Define instrument responses
# instrument_responses = np.linspace(0, 300, 10000)

# # Define folder paths
# folder_0 = "disorder kappa_k1 = kappa_k2 = 0"
# folder_5 = "disorder kappa_k1 = 0 and kappa_k2 = 5"

# # Define file names
# files = [
#     "H_strongC_2nm_symm_fast_tumb_jigg_convolved_smart_k5.npy",
#     "J_strongC_2nm_symm_fast_tumb_jigg_convolved_smart_k5.npy"
# ]

# # Plotting
# plt.figure(figsize=(12, 5))

# for i, file in enumerate(files):
#     # Load data from each folder
#     data_0 = np.load(os.path.join(folder_0, file))
#     data_5 = np.load(os.path.join(folder_5, file))
    
#     # Create subplot
#     plt.subplot(1, 2, i + 1)
#     plt.plot(instrument_responses, data_0, label=r'Second photon detection: concentration parameter kappa = 0')
#     plt.plot(instrument_responses, data_5, label=r'Second photon detection: concentration parameter kappa = 5')
#     plt.title(file.split('_')[0] + r'-dimer vs Instrument Response')
#     plt.xlabel(r'Instrument Response (ps)')
#     plt.ylabel(r'$g^{(2)}(\infty, 0)$')
#     plt.legend()
#     plt.grid(True)

# plt.tight_layout()
# plt.savefig('H_J_second_photon_sampling_disorder.png', dpi=300, bbox_inches='tight')
# plt.show()

# # Load data
# g2_corr_H = np.load('files/g2/H_strongC_2nm_symm_fast_tumb_jigg_normalized_smart_k5.npy')
# g2_corr_J = np.load('files/g2/J_strongC_2nm_symm_fast_tumb_jigg_normalized_smart_k5.npy')


# insets_H = [50, 100, 200]
# insets_J = [50, 100, 200]

# colors_H = ['midnightblue']
# colors_J = ['maroon']

# def plot_inset(g2_corr, instrument_response, title, color):
#     step_con = 1.0e4
#     dt = (5 * system_params.tau_L - (-5 * system_params.tau_L)) / step_con
#     sigma = instrument_response
    
#     g2_corr_appended = gaussian_filter1d(
#         np.append(np.flip(g2_corr[1:]), g2_corr), 
#         sigma=sigma / dt, 
#         mode='reflect'
#     )
#     times_appended = np.append(-np.flip((times_ns[1:] - times_ns[0])), times_ns)
    
#     plt.figure(figsize=(6, 4))
#     plt.plot(times_appended, g2_corr_appended, color=color)
#     plt.xlim([-2, 2])
#     plt.ylim([0.4, 0.9])
#     plt.xlabel('Time delay (ns)', fontsize=25)
#     plt.ylabel(r'$[g^{(2)}(\infty, 0)]_{(O,\textbf{q})}^{\kappa = 5}$', fontsize=25)
#     plt.xticks([-2, -1, 0, 1, 2], fontsize=20)  # Reduce number of x-ticks
#     plt.yticks(fontsize=20)
#     plt.grid(True, linestyle='--', alpha=0.6)
#     plt.tight_layout()  # Ensure the plot fits properly in the exported file
#     plt.savefig(f'{title.replace(" ", "_")}.svg', dpi=600)
#     plt.show()

# # Generate inset plots
# for ir in insets_H:
#     plot_inset(g2_corr_H, ir, f'H dimer inset - IR {ir} ps', colors_H[0])
    
# for ir in insets_J:
#     plot_inset(g2_corr_J, ir, f'J dimer inset - IR {ir} ps', colors_J[0])


# # np.save('files/g2/ortho_dimer/ortho_weakC_mic_T0.npy', g2_corr)

# # print(g2_corr)
# # np.save('g2_par', g2_corr)

'''g2 distribution'''

# # Create a Bloch-Redfield calculator instance
# BR_calculator = BlochRedfieldCalculator(H, a_ops_opt, Gamma_opt, rates, a_ops_vib = a_ops_vib, Gamma_vib = Gamma_vib, a_ops_coup = a_ops_coup, Gamma_coup = Gamma_coup, P_vib = P_vib, mode=mode, pump = pump)

# # Initialize the DistributionAnalysis class for calculating intensity
# quantum_analysis_g2_dist = DistributionAnalysis(BR_calculator, rates,quantity='g2 distribution', rho=rho0, final_time=tf, steps=steps)

# # Case 1: G^{(2)}(π/2, φ_p)
# n_phi_p = 100
# phi_vals_p = np.linspace(0, 2 * np.pi, n_phi_p)
# theta_fixed = np.pi / 2  # Fix theta at π/2
# g2_phi = np.zeros_like(phi_vals_p)

# for l, phi_p in enumerate(phi_vals_p):
#     g2_phi[l] = quantum_analysis_g2_dist.compute(theta_fixed, phi_p, theta_fixed, phi_p)

# # Normalize data for phi
# g2_phi_normalized = g2_phi #/ np.max(g2_phi)

# # Case 2: G^{(2)}(θ_p, 0)
# n_theta_p = 100
# theta_vals_p = np.linspace(-np.pi, np.pi, n_theta_p)
# phi_fixed = 0  # Fix phi at 0
# g2_theta = np.zeros_like(theta_vals_p)

# for k, theta_p in enumerate(theta_vals_p):
#     g2_theta[k] = quantum_analysis_g2_dist.compute(theta_p, phi_fixed, theta_p, phi_fixed)

# # Normalize data for theta
# g2_theta_normalized = g2_theta # / np.max(g2_theta)

# # Common legend configuration
# legend_fontsize = 50  # Font size for the legend
# tick_fontsize = 55    # Font size for tick labels

# # # Plot 1: Unnormalized data
# # fig1, ax1 = plt.subplots(figsize=(12, 9), subplot_kw={'projection': 'polar'})
# # ax1.plot(phi_vals_p, g2_phi, color='blue', linewidth=3, label=r'$G^{(2)}(\pi/2, \phi^{\prime})$')
# # ax1.plot(theta_vals_p, g2_theta, color='red', linewidth=3, linestyle='--', label=r'$G^{(2)}(\theta^{\prime}, 0)$')

# # # Customize tick labels
# # ticks = np.linspace(0, 2 * np.pi, 8, endpoint=False)
# # tick_labels = [r'$0$', r'$\frac{\pi}{4}$', r'$\frac{\pi}{2}$', r'$\frac{3\pi}{4}$',
# #                r'$\pi$', r'$\frac{5\pi}{4}$', r'$\frac{3\pi}{2}$', r'$\frac{7\pi}{4}$']
# # ax1.set_thetagrids(np.degrees(ticks), labels=tick_labels, fontsize=tick_fontsize)
# # ax1.set_rticks([0.5, 1])
# # ax1.tick_params(axis='both', labelsize=65, pad = 40) 
# # # Combine handles and labels for the common legend
# # handles1, labels1 = ax1.get_legend_handles_labels()

# # fig1.legend(
# #     handles1, labels1,
# #     loc='upper center',             # Align legends horizontally at the top
# #     bbox_to_anchor=(0.5, 1.2),     # Adjust position slightly above the plots
# #     fontsize=legend_fontsize,       # Set font size
# #     frameon=False,                  # Remove frame
# #     ncol=2,                          # Place all legend entries in one row
# #     handlelength=1,                 # Length of the legend line
# #     labelspacing=1,                 # Adjust the spacing between the label and the line
# #     markerscale=1.5,                # Control the size of the marker relative to the line
# # )

# # # Save each figure and show plots
# # fig1.savefig('unnorm_ortho_polar_radiation.png', dpi=300, bbox_inches='tight')



# # Plot 2: Normalized data
# fig2, ax2 = plt.subplots(figsize=(12, 9), subplot_kw={'projection': 'polar'})
# ax2.plot(phi_vals_p, g2_phi_normalized, color='blue', linewidth=3, label=r'$g^{(2)}(\pi/2, \phi^{\prime})$')
# ax2.plot(theta_vals_p, g2_theta_normalized, color='red', linewidth=3, linestyle='--', label=r'$g^{(2)}(\theta^{\prime}, 0)$')

# # Customize tick labels
# ticks = np.linspace(0, 2 * np.pi, 8, endpoint=False)
# tick_labels = [r'$0$', r'$\frac{\pi}{4}$', r'$\frac{\pi}{2}$', r'$\frac{3\pi}{4}$',
#                r'$\pi$', r'$\frac{5\pi}{4}$', r'$\frac{3\pi}{2}$', r'$\frac{7\pi}{4}$']
# ax2.set_thetagrids(np.degrees(ticks), labels=tick_labels, fontsize=tick_fontsize)
# ax2.set_rticks([0.5, 1.0])
# ax2.tick_params(axis='both', labelsize=65, pad = 40) 
# handles2, labels2 = ax2.get_legend_handles_labels()

# fig2.legend(
#     handles2, labels2,
#     loc='upper center',             # Align legends horizontally at the top
#     bbox_to_anchor=(0.5, 1.15),     # Adjust position slightly above the plots
#     fontsize=legend_fontsize,       # Set font size
#     frameon=False,                  # Remove frame
#     ncol=2,                          # Place all legend entries in one row
#     handlelength=1,                 # Length of the legend line
#     labelspacing=1,                 # Adjust the spacing between the label and the line
#     markerscale=1.5,                # Control the size of the marker relative to the line
# )

# fig2.savefig('norm_ortho_polar_radiation.png', dpi=300, bbox_inches='tight')

# plt.show()



# # Create a Bloch-Redfield calculator instance
# BR_calculator = BlochRedfieldCalculator(H, a_ops_opt, Gamma_opt, rates, a_ops_vib = a_ops_vib, Gamma_vib = Gamma_vib, a_ops_coup = a_ops_coup, Gamma_coup = Gamma_coup, P_vib = P_vib, mode=mode, pump = pump)

# # Initialize the DistributionAnalysis class for calculating intensity
# quantum_analysis_g2_dist = DistributionAnalysis(BR_calculator, rates,quantity='g2 distribution', rho=rho0, final_time=tf, steps=steps)


# n_phi = 50
# phi_vals = np.linspace(0, 2 * np.pi, n_phi)
# n_theta = 50
# theta_vals = np.linspace(-np.pi, np.pi, n_theta)

# phi, theta = np.meshgrid(phi_vals, theta_vals)

# # Initialize an array to store G^{(2)} values
# g2_values = np.zeros_like(phi)

# # Loop through the (theta, phi) grid and compute G^{(2)}(theta, phi)
# for i in range(phi.shape[0]):
#     for j in range(phi.shape[1]):
#         # Compute G^{(2)} for each pair of (theta, phi)
#         g2_values[i, j] = quantum_analysis_g2_dist.compute(theta[i, j], phi[i, j], theta[i, j], phi[i, j])

# # Now g2_values will contain the computed G^{(2)} values for each (theta, phi)

# # Convert spherical coordinates (r, theta, phi) to Cartesian coordinates (x, y, z)
# r = g2_values  # Radial distance is given by the G^{(2)} values
# x = r * np.sin(theta) * np.cos(phi)
# y = r * np.sin(theta) * np.sin(phi)
# z = r * np.cos(theta)

# surface = go.Surface(
#     x=x, 
#     y=y, 
#     z=z, 
#     surfacecolor=g2_values,  # Set surface color based on G^{(2)} values
#     colorscale='RdBu_r',      # Use the RdYlBu colormap
#     colorbar=dict(title='G^{(2)}(θ, φ)')  # Add a color bar with the title
# )

# # Layout for the plot
# layout = go.Layout(
#     title=r'3D Plot of $G^{(2)}(\theta, \phi)$',
#     scene=dict(
#         xaxis=dict(
#             gridcolor='rgb(255, 255, 255)',
#             zerolinecolor='rgb(255, 255, 255)',
#             showbackground=True,
#             backgroundcolor='rgb(230, 230,230)'
#         ),
#         yaxis=dict(
#             gridcolor='rgb(255, 255, 255)',
#             zerolinecolor='rgb(255, 255, 255)',
#             showbackground=True,
#             backgroundcolor='rgb(230, 230,230)'
#         ),
#         zaxis=dict(
#             gridcolor='rgb(255, 255, 255)',
#             zerolinecolor='rgb(255, 255, 255)',
#             showbackground=True,
#             backgroundcolor='rgb(230, 230,230)'
#         )
#     ),
#     height=800
# )

# # Create the figure and plot it
# fig = go.Figure(data=[surface], layout=layout)
# # fig.write_html("g2_J.html")
# fig.show()







# g2_corr_weak = np.load('files/J_weakC_2nm_symm_nfd_slow_tumb_nojigg.npy')
# g2_corr_strong = np.load('files/J_strongC_2nm_symm_nfd_slow_tumb_nojigg.npy')

# g2_corr_weak_unnorm = np.load('files/H_weakC_2nm_symm_nfd_slow_tumb_nojigg_unnorm.npy')
# g2_corr_strong_unnorm = np.load('files/H_strongC_2nm_symm_nfd_slow_tumb_nojigg_unnorm.npy')

# T_vib = np.linspace(0, 300, 100)

# # plt.plot(T_vib, g2_corr_weak)
# # plt.plot(T_vib, g2_corr_strong)

# plt.plot(T_vib, g2_corr_weak_unnorm)
# plt.plot(T_vib, g2_corr_strong_unnorm)

# plt.show()


# # # # Load data
# # # T_H = 90.41190640534819
# # # T_J = 180.82381281069638

# g2_corr_H0 = np.load('files/g2_temp_H.npy')
# g2_corr_J0 = np.load('files/g2_temp_J.npy')


# # Simulated data for demonstration (since actual files are not available)
# T_vib = np.linspace(0, 300, 500)
# # Create the figure and twin y-axis
# fig, ax_H = plt.subplots(figsize=(16, 10))
# ax_J = ax_H.twinx()

# # # Define a colormap gradient from blue to red
# # gradient = np.linspace(0, 1, 500).reshape(1, -1)
# # cmap = mcolors.LinearSegmentedColormap.from_list("gradient", ["blue", "red"])

# # # Add gradient background
# # ax_H.imshow(gradient, aspect="auto", extent=[0, 300, min(g2_corr_H0) - 5, max(g2_corr_H0) + 10],
# #             cmap=cmap, alpha=0.3, origin="lower")

# # Plot data for H and J dimers
# H_line, = ax_H.plot(T_vib, g2_corr_H0, color='blue', label='H dimer', linewidth=3)
# J_line, = ax_J.plot(T_vib, g2_corr_J0, color='red', label='J dimer', linewidth=3)

# # # Add vertical dashed lines
# # ax_H.axvline(T_H, color='black', linestyle='dashed', alpha=0.3)
# # ax_H.axvline(T_J, color='black', linestyle='dashed', alpha=0.3)

# # Set axis labels
# ax_H.set_xlabel(r'Temperature (K)', fontsize=40)
# ax_H.set_ylabel(r'$g^{(2)}(\infty, 0)$', fontsize=40, color='blue')
# ax_J.set_ylabel(r'$g^{(2)}(\infty, 0)$', fontsize=40, color='red')

# # Color the tick labels on right y-axis
# ax_J.yaxis.label.set_color('red')
# ax_H.yaxis.label.set_color('blue')
# ax_J.tick_params(axis='y', colors='red', labelsize=35, direction='in', length=6)
# ax_H.tick_params(axis='y', colors='blue', labelsize=35, direction='in', length=6)
# ax_H.tick_params(axis='x', labelsize=35, direction='in', length=6, pad=10)

# # # Add extra label for |T_vib / J_12| as a legend
# # extra_legend = plt.Line2D([0], [0], color='black', linestyle='dashed', label=r'$T_{\textrm{vib}} / \mathbf{J}_{12}$', alpha = 0.3)

# # # Combine all legends in lower left corner
# # ax_H.legend(handles=[H_line, J_line, extra_legend], loc='lower left', fontsize=40, frameon=False)
# ax_H.legend(handles=[H_line, J_line], loc='lower left', fontsize=40, frameon=False)

# ax_H.set_xlim(0, 300)

# # Set y-axis limits (optional, adjust as needed)
# ax_H.set_ylim(min(g2_corr_H0) - 5, max(g2_corr_H0) + 10)
# ax_J.set_ylim(min(g2_corr_J0) - 0.0005, max(g2_corr_J0) + 0.0005)

# # Set x-axis to log scale
# ax_H.set_xscale('log')
# ax_J.set_xscale('log')

# ax_H.set_xlim(1, 300)  

# # Save and display
# plt.savefig('g2_TempDep.png', dpi=300, bbox_inches='tight')
# plt.show()














# fig, ax = plt.subplots()
# ax.set_xlabel(r'$\phi$')
# ax.set_ylabel(r'$\theta$')
# img = ax.pcolormesh(phi_vals_p, theta_vals_p, g2_corr.T, cmap='hot', shading='auto')
# plt.colorbar(img)
# plt.show()

# fig, ax = plt.subplots()
# ax.set_xlabel(r'$\phi$')
# ax.set_ylabel(r'$\theta$')
# vimg = ax.imshow(g2_corr, cmap='hot', interpolation='nearest', 
#                 extent=[min(phi_vals_p), max(phi_vals_p), min(theta_vals_p), max(theta_vals_p)], 
#                 origin='lower')  # Set origin to 'lower'
# ax.set_aspect("auto")
# plt.colorbar(vimg)

# plt.show()

# # # # print(g2_corr)
# # Generate R_vals and Theta_vals for the polar plot
# R_vals, Theta_vals = np.meshgrid(np.sin(theta_vals_p), phi_vals_p)

# # Create the polar heatmap
# fig = plt.figure(figsize=(12, 8))  # Adjusted figure size
# ax1 = fig.add_subplot(projection='polar')

# # Create the colormap on the polar plot
# c = ax1.pcolormesh(Theta_vals, R_vals, g2_corr, cmap='coolwarm', shading='auto', vmin=0, vmax=1)

# # Set labels and formatting
# ax1.set_xlabel(r'$\theta$ (degrees)', fontsize=35, labelpad=0) 
# ax1.set_ylabel(r'$g^{(2)}(\theta)$', fontsize=35, labelpad=50)
# ax1.tick_params(axis='x', labelsize=30, pad=15, direction='in', length=6)  # pad=30 moves angle labels further outside
# ax1.tick_params(axis='y', labelsize=30, direction='in', length=6)

# # Move the text label to the upper-left corner
# ax1.text(-0.3, 1.05, r'(a)', size=30, transform=ax1.transAxes)  # Adjusted position to upper-left

# # Remove y-tick labels (optional based on your preference)
# ax1.set_yticklabels([])

# # Adjust layout to prevent cutting off theta labels
# plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.2)

# # Add colorbar with specific size and position adjustments
# cbar = plt.colorbar(c, ax=ax1, shrink=1.0, aspect=10, pad=0.1)  # Shrink and aspect to control size, pad to move it away

# # Adjust text size on colorbar
# cbar.ax.tick_params(labelsize=30)  # Change the number to adjust the colorbar tick label size

# # # Save the figure with the desired format and settings
# # plt.savefig('dist_g2_ortho.png', dpi=500, bbox_inches='tight')

# # Display the plot
# plt.show()

# '''Averaging over quantities'''

# #######################################################
# '''Averaging over k, k_prime (same) == The tumbling effect is slower that detection rate'''
# #######################################################

# # collapse operators (for initial incoherent pumping)
# c_ops_S = [np.sqrt(system_params.gamma_p_1) * sp1, np.sqrt(system_params.gamma_p_2) * sp2]
# # c_ops_S = [np.sqrt(system_params.gamma_p_1) * sp1, np.sqrt(system_params.gamma_p_2) * sp2, np.sqrt(system_params.gamma_q1) * sp1@sm1, np.sqrt(system_params.gamma_q1) * sp2@sm2]

# # Create a Bloch-Redfield calculator instance
# BR_calculator = BlochRedfieldCalculator(H, a_ops_opt, Gamma_opt, rates, a_ops_vib = a_ops_vib, Gamma_vib = Gamma_vib, a_ops_coup = a_ops_coup, Gamma_coup = Gamma_coup, P_vib = P_vib, mode=mode, pump = pump)

# # Initialize the DistributionAnalysis class for calculating g2
# quantum_analysis_g2_dist = DistributionAnalysis(BR_calculator, rates,quantity='g2 distribution', rho=rho0, final_time=tf, steps=steps)

# # Compute intensity distribution
# g2_corr = np.zeros((n_phi_p, n_theta_p, steps))
# for l, phi_p in enumerate(phi_vals_p):
#     for k, theta_p in enumerate(theta_vals_p):
#         theta, phi = theta_p, phi_p
#         g2_corr[l, k, :] = quantum_analysis_g2_dist.compute(theta, phi, theta_p, phi_p)

# g2_corr = g2_corr / np.max(g2_corr)
# g2_corr = np.mean(g2_corr, axis=(0, 1))


# # # print(norm_int_dir)
# # print(g2_corr[0])
# plt.plot(np.append(- np.flip((times_ns[1:] - times_ns[0])), (times_ns)), np.append(np.flip(g2_corr[1:]), g2_corr[:]))
# plt.grid(True)
# # plt.ylim([0.0,1.1])
# plt.show()

# plt.plot(np.append(- np.flip((times_ns[1:] - times_ns[0])), (times_ns)), np.append(np.flip(g2_corr[0, 0, 1:]), g2_corr[0, 0, :]))
# plt.show()






















# # File paths for H and J dimers
# files_H = [
#     'files/g2/H_2nm_symm_0p05meV.npy',
#     'files/g2/H_2nm_symm_0p5meV.npy',
#     'files/g2/H_2nm_symm_5meV.npy',
#     'files/g2/H_2nm_symm_50meV.npy'
# ]

# files_J = [
#     'files/g2/J_2nm_symm_0p05meV.npy',
#     'files/g2/J_2nm_symm_0p5meV.npy',
#     'files/g2/J_2nm_symm_5meV.npy',
#     'files/g2/J_2nm_symm_50meV.npy'
# ]

# # Reorganized energies (spaced along the axis going into the page)
# energies = [0.05, 0.5, 5, 50]  # Original values in meV
# log_energies = np.log10(energies)  # Log-transform energies for plotting

# # Plot H dimer curves (blue)
# fig = plt.figure(figsize=(20, 16))
# ax = fig.add_subplot(111, projection='3d')

# for i, file in enumerate(files_H):
#     g2 = np.load(file)
#     flipped_g2 = np.append(np.flip(g2[1:]), g2)  # Flip symmetrically
#     time_flipped = np.append(-np.flip(times_ns[1:]), times_ns)
#     log_energy = np.full_like(time_flipped, log_energies[i])  # Match length of time axis

#     ax.plot(log_energy, time_flipped, flipped_g2, color='blue', linewidth=2) #, label=None if i > 0 else "H dimer")

# # Set custom ticks for the x-axis (logarithmic scale)
# ax.set_xticks(log_energies)
# ax.set_xticklabels([f"{e}" for e in energies])  # Original energy values

# ax.set_ylim([-40, 40])
# # Set axis labels with additional padding
# ax.set_xlabel(r'$\lambda$ (meV)', fontsize=60, labelpad=55)
# ax.set_ylabel(r'$\tau$ (ns)', fontsize=60, labelpad=55)
# ax.set_zlabel(r'$g^{(2)}(\infty, \tau)$', fontsize=60, labelpad=55)

# # Customize tick sizes and distances
# ax.tick_params(axis='x', labelsize=45, direction='in', length=6, pad=20)
# ax.tick_params(axis='y', labelsize=45, direction='in', length=6, pad=20)
# ax.tick_params(axis='z', labelsize=45, direction='in', length=6, pad=30)

# # Adjust plot margins using tight_layout to avoid clipping
# plt.tight_layout()

# # Expand the canvas to avoid cropping when saving
# fig.subplots_adjust(left=0.2, right=0.8, top=0.9, bottom=0.3)

# # Add legend
# ax.legend(fontsize=40, frameon=False)

# # Set viewing angle
# ax.view_init(elev=30, azim=-60)

# # Save and show plot
# plt.savefig('H3d_g2_plot.png', dpi=600)
# plt.show()

# # Repeat for J dimers
# fig = plt.figure(figsize=(20, 16))
# ax = fig.add_subplot(111, projection='3d')

# for i, file in enumerate(files_J):
#     g2 = np.load(file)
#     flipped_g2 = np.append(np.flip(g2[1:]), g2)  # Flip symmetrically
#     time_flipped = np.append(-np.flip(times_ns[1:]), times_ns)
#     log_energy = np.full_like(time_flipped, log_energies[i])  # Match length of time axis

#     ax.plot(log_energy, time_flipped, flipped_g2, color='red', linewidth=2) #, label=None if i > 0 else "J dimer")

# # Set custom ticks for the x-axis (logarithmic scale)
# ax.set_xticks(log_energies)
# ax.set_xticklabels([f"{e}" for e in energies])  # Original energy values

# ax.set_ylim([-40, 40])
# # Set axis labels with additional padding
# ax.set_xlabel(r'$\lambda$ (meV)', fontsize=60, labelpad=55)
# ax.set_ylabel(r'$\tau$ (ns)', fontsize=60, labelpad=55)
# ax.set_zlabel(r'$g^{(2)}(\infty, \tau)$', fontsize=60, labelpad=55)

# # Customize tick sizes and distances
# ax.tick_params(axis='x', labelsize=45, direction='in', length=6, pad=20)
# ax.tick_params(axis='y', labelsize=45, direction='in', length=6, pad=20)
# ax.tick_params(axis='z', labelsize=45, direction='in', length=6, pad=30)

# # Adjust plot margins using tight_layout to avoid clipping
# plt.tight_layout()

# # Expand the canvas to avoid cropping when saving
# fig.subplots_adjust(left=0.2, right=0.8, top=0.9, bottom=0.3)

# # Add legend
# ax.legend(fontsize=40, frameon=False)

# # Set viewing angle
# ax.view_init(elev=30, azim=-60)

# # Save and show plot
# plt.savefig('J3d_g2_plot.png', dpi=600)
# plt.show()








t1 = time.time()

total = t1-t0
print(total)








# instrument_responses = np.linspace(1, 250, 100)
# H_fast_jigg = np.load('files/H_strongC_2nm_symm_fast_tumb_jigg_convolved_smart_5.npy')
# J_fast_jigg = np.load('files/J_strongC_2nm_symm_fast_tumb_jigg_convolved_smart_5.npy')

# plt.plot(instrument_responses, H_fast_jigg)
# plt.plot(instrument_responses, J_fast_jigg)
# plt.show()

































# indistinguishable_dimer = np.load('files/indistinguishable_strongC.npy')
# distinguishable_dimer = np.load('files/distinguishable_strongC.npy')
# # H_dimer_weak = np.load('files/intensity/H_dimer/H_polaronC.npy')
# # J_dimer_weak = np.load('files/intensity/J_dimer/J_polaronC.npy')
# # ortho_dimer_weak = np.load('files/intensity/ortho_dimer/ortho_polaronC.npy')
# H_dimer_strong = np.load('files/H_strongC.npy')
# J_dimer_strong = np.load('files/J_strongC.npy')
# ortho_dimer_strong = np.load('files/ortho_strongC.npy')
# # H_dimer_stronger = np.load('files/intensity/oppo_dimer/oppo_strongerC.npy')
# # J_dimer_stronger = np.load('files/intensity/J_dimer/J_strongerC.npy')
# # ortho_dimer = np.load('files/intensity/ortho_dimer/ortho_strongerC.npy')
# # anti_parallel = np.load('files/intensity/oppo_dimer/oppo_strongerC.npy')

# # Create main plot
# fig, ax = plt.subplots(figsize=(12, 8))

# # # Plot the data with appropriate labels
# # ax.plot(system_params.gamma_q1 * times, intensity_results, '-y', label='Anti-parallel', marker='s', markersize=5, markevery=4300)

# # line0, = ax.plot(system_params.gamma_q1 * times, indistinguishable_dimer, '--k', label='Orthogonal', linewidth = 3, alpha = 0.3)
# # line1, = ax.plot(system_params.gamma_q1 * times, distinguishable_dimer, '--k', label='Orthogonal', linewidth = 3, alpha = 0.3)
# # # line1, = ax.plot(system_params.gamma_q1 * times, anti_parallel, '-k', label='Anti-parallel', marker='s', markersize=6, markevery=25000)
# # line2, = ax.plot(system_params.gamma_q1 * times, H_dimer_strong, color = 'blue', label='H dimer', marker='s', markersize=0, markevery=10, linewidth = 3)
# # line3, = ax.plot(system_params.gamma_q1 * times, J_dimer_strong, color = 'tomato', label='J dimer', marker='^', markersize=0, markevery=13, linewidth = 3)
# # # line4, = ax.plot(system_params.gamma_q1 * times, H_dimer_strong, color = 'mediumblue', label='H dimer', marker='s', markersize=0, markevery=10)
# # # line5, = ax.plot(system_params.gamma_q1 * times, J_dimer_strong, color = 'red', label='J dimer', marker='^', markersize=0, markevery=13)
# # # line6, = ax.plot(system_params.gamma_q1 * times, H_dimer_stronger, color = 'darkblue', label='H dimer', marker='s', markersize=0, markevery=10)
# # # line7, = ax.plot(system_params.gamma_q1 * times, J_dimer_stronger, color = 'darkred', label='J dimer', marker='^', markersize=0, markevery=13)
# # line8, = ax.plot(system_params.gamma_q1 * times, ortho_dimer_strong, '--g', label='Orthogonal', linewidth = 3)
# end_y = indistinguishable_dimer[-1]  # Get the y-value at the end of the H_dimer_strong curve
# ax.axhline(y=end_y, color='gray', linestyle='--', linewidth=3, alpha=0.6)  # Customize appearance

# line0, = ax.plot(system_params.gamma_q1 * times, indistinguishable_dimer, 
#                  color='purple', label='Indistinguishable', 
#                  linestyle=':', linewidth = 3, alpha=0.8)

# # line1, = ax.plot(system_params.gamma_q1 * times, distinguishable_dimer, 
# #                  color='k', label='Distinguishable', 
# #                  linestyle='None', marker='o', markevery=100, markersize=6, alpha=1)

# # Keep H, J, and orthogonal with solid or dotted lines for contrast
# line2, = ax.plot(system_params.gamma_q1 * times, H_dimer_strong, 
#                  color='blue', label='H dimer', marker='s', markersize=0, linewidth=3)

# line3, = ax.plot(system_params.gamma_q1 * times, J_dimer_strong, 
#                  color='red', label='J dimer', marker='^', markersize=0, linewidth=3)

# line8, = ax.plot(system_params.gamma_q1 * times, ortho_dimer_strong, 
#                  linestyle='--', color='green', label='Orthogonal', linewidth=3)

# # Set labels and title
# ax.set_xlabel(r'$\gamma t$', fontsize=35)
# ax.set_ylabel(r'$I(t) / I_0$', fontsize=35)
# # ax.set_ylim([0, 2.1])
# # ax.set_xlim([system_params.gamma_q1 * times[0], system_params.gamma_q1 * times[-1]])
# # Set tick parameters
# ax.tick_params(axis='x', labelsize=30, direction='in', length=6)
# ax.tick_params(axis='y', labelsize=30, direction='in', length=6)
# # ax.text(0.05, 1.97, r'$T_{opt} = 5800$ K and $T_{vib} = 300$ K', size=30)

# # ax.set_xscale('log')
# # ax.set_title('Intensity vs time', fontsize=18)
# # Change number of ticks on y-axis
# ax.yaxis.set_major_locator(MaxNLocator(nbins=5))  # Change nbins to desired number of ticks

# # # Add legend and grid
# # ax.legend(handles=[line2, line3, line4, line5, line6, line7, line8], 
# #           labels=['H dimer-weak', 'J dimer-weak', 'H dimer-strong', 'J dimer-strong', 'H dimer-stronger', 'J dimer-stonger', 'Orthogonal'],
# #           loc='best', 
# #           fontsize=18, 
# #           frameon=False)

# # Add legend and grid
# # Dummy line for spacing in legend
# empty = Line2D([], [], linestyle='None')
# # First row: Indistinguishable, Distinguishable
# legend1 = ax.legend(handles=[line0],
#                     labels=['Indistinguishable'],
#                     loc='upper center',
#                     fontsize=25,
#                     frameon=False,
#                     ncol=2,
#                     bbox_to_anchor=(0.3, 0.93))  # Adjust horizontal position

# # Second row: H dimer, J dimer, Orthogonal
# legend2 = ax.legend(handles=[line2, line3, line8],
#                     labels=['H dimer', 'J dimer', 'Orthogonal'],
#                     loc='upper center',
#                     fontsize=25,
#                     frameon=False,
#                     ncol=3,
#                     bbox_to_anchor=(0.56, 1.02))  # Slightly lower than first row

# # Add first legend manually to the axes
# ax.add_artist(legend1)
# # ax.grid(True)

# # # # Inset plot
# # ax_inset = plt.axes([0.42, 0.3, 0.45, 0.25])  # [left, bottom, width, height]
# # diff_data_1 = H_dimer_weak - J_dimer_weak
# # # diff_data_2 = H_dimer_strong - J_dimer_strong
# # # diff_data_3 = H_dimer_stronger - J_dimer_stronger
# # ax_inset.plot(system_params.gamma_q1 * times, diff_data_1)
# # # ax_inset.plot(system_params.gamma_q1 * times, diff_data_2, label = r'$\lambda = 5 meV$, $\omega_{c} = 90 meV$')
# # # ax_inset.plot(system_params.gamma_q1 * times, diff_data_3, label = r'$\lambda = 20 meV$, $\omega_{c} = 90 meV$')
# # # diff_data = H_dimer - J_dimer
# # # ax_inset.plot(system_params.gamma_q1 * times, diff_data)
# # ax_inset.set_xlabel(r'$\gamma t$', fontsize=30)
# # ax_inset.set_ylabel(r'$\Delta I(t) / I_0$', fontsize=30)
# # ax_inset.tick_params(axis='x', labelsize=25)
# # ax_inset.tick_params(axis='y', labelsize=25)
# # ax_inset.set_xlim([0, 4])
# # ax_inset.legend(fontsize=22, frameon=False)
# # # ax_inset.set_ylim([-0.04, 0.0])
# # # ax_inset.grid(True)

# # # Save the plot
# plt.savefig("orient_intensity.png", dpi = 300, bbox_inches='tight')

# # Show the plot
# plt.show()






# times_sym = np.append(-np.flip(times_ns[1:] - times_ns[0]), times_ns)
# g2_perp = np.append(np.flip(np.load('files/g2/ortho_dimer/ortho_weakC_perp_T0.npy')[1:]), np.load('files/g2/ortho_dimer/ortho_weakC_perp_T0.npy'))
# g2_mic_deph = np.append(np.flip(np.load('files/g2/ortho_dimer/ortho_weakC_mic_T0.npy')[1:]), np.load('files/g2/ortho_dimer/ortho_weakC_mic_T0.npy'))
# g2_par = np.append(np.flip(np.load('files/g2/ortho_dimer/ortho_weakC_par_T0.npy')[1:]), np.load('files/g2/ortho_dimer/ortho_weakC_par_T0.npy'))

# # Create figure and axis
# fig, ax = plt.subplots(figsize=(14, 14))

# # Plot the data with appropriate labels
# ax.plot(times_sym, g2_perp, '-b', label='Perpendicular', linewidth = 1)
# ax.plot(times_sym, g2_mic_deph, 'r', label= r'$\theta = \pi/2, \phi = \pi/4$', linewidth = 1)
# ax.plot(times_sym, g2_par, '--g', label= 'Parallel', linewidth = 1)

# # Set labels and title
# ax.set_xlabel(r'$\tau$ (ns)', fontsize=40)
# ax.set_ylabel(r'$g^{(2)}(\infty, \tau)$', fontsize=40)
# ax.set_xlim([-7, 7])

# # Set tick parameters
# ax.tick_params(axis='x', labelsize=40, direction='in', length=6)
# ax.tick_params(axis='y', labelsize=40, direction='in', length=6)
# ax.text(-6.8, 0.91, '(b)', size=40)

# # Add legend and grid
# ax.legend(loc='best', fontsize=30, frameon=False)


# # Change number of ticks on y-axis
# ax.yaxis.set_major_locator(MaxNLocator(nbins=5))


# # # Inset plot
# # ax_inset = plt.axes([0.625, 0.3, 0.25, 0.25])  # [left, bottom, width, height]
# # ax_inset.plot(1.0e3 * times_sym, g2_mic_deph, '-r', markersize=0)
# # ax_inset.set_xlabel(r'$t$ (ps)', fontsize=14)
# # ax_inset.set_ylabel(r'$g^{(2)}(\infty, \tau)$', fontsize=14)
# # ax_inset.tick_params(axis='x', labelsize=14)
# # ax_inset.tick_params(axis='y', labelsize=14)
# # ax_inset.set_xlim([-1, 1])
# # # ax_inset.set_ylim([-0.04, 0.0])
# # # ax_inset.grid(True)


# # Save the plot
# plt.savefig("orient_g2_weakC_ortho_T0.png", dpi=500)

# # Show the plot
# plt.show()





















# # Plot the downsampled data
# fig, ax = plt.subplots(figsize=(14, 10))

# # Load data
# h_dimer = np.append(
#     np.flip(np.load('files/g2/H_strongC_2nm_symm.npy')[1:]),
#     np.load('files/g2/H_strongC_2nm_symm.npy')
# )
# j_dimer = np.append(
#     np.flip(np.load('files/g2/J_strongC_2nm_symm.npy')[1:]),
#     np.load('files/g2/J_strongC_2nm_symm.npy')
# )
# dist_dimer = np.append(
#     np.flip(np.load('files/ortho_strongC_2nm_symm_perp.npy')[1:]),
#     np.load('files/ortho_strongC_2nm_symm_perp.npy')
# )
# times = np.append(
#     -np.flip(np.linspace(0, 1.0e-3 * tf, 10000)[1:] - np.linspace(0, 1.0e-3 * tf, 10000)[0]),
#     np.linspace(0, 1.0e-3 * tf, 10000)
# )

# # sigma = 40
# # zero_delay_H_convolved_array = []
# # zero_delay_J_convolved_array = []
# # step_con = 5.0e5

# # dt = (5 * system_params.tau_L - (-5 * system_params.tau_L)) / step_con
# # g2_corr_H_appended = gaussian_filter1d(h_dimer, sigma=sigma / dt, mode='reflect')
# # g2_corr_J_appended = gaussian_filter1d(j_dimer, sigma=sigma / dt, mode='reflect')
# # times_appended = np.append(-np.flip((times_ns[1:] - times_ns[0])), times_ns)


# # # Main plot
# # line0, = ax.plot(times_appended, g2_corr_H_appended, '-b', linewidth=3)
# # line1, = ax.plot(times_appended, g2_corr_J_appended, '-r', linewidth=3)

# # # Set axis labels and limits
# # ax.set_xlabel(r'$\tau$ (ns)', fontsize=35)
# # ax.set_ylabel(r'$g^{(2)}(\infty, \tau)$', fontsize=35)
# # ax.set_xlim([-10, 10])
# # ax.set_ylim([0.4, 1.3])
# # ax.tick_params(axis='x', labelsize=35, direction='in', length=6)
# # ax.tick_params(axis='y', labelsize=35, direction='in', length=6)
# # ax.legend(handles=[line0, line1], labels=[r"H dimer", r"J dimer"], fontsize=30, frameon=False)

# # # plt.plot(times_appended, g2_corr_H_appended)
# # # plt.plot(times_appended, g2_corr_J_appended)
# # plt.savefig('orient_g2_symm_convolved', dpi=300)
# # plt.show()

# # Main plot
# line0, = ax.plot(times, h_dimer, '-b', linewidth=3, label="H dimer")
# line1, = ax.plot(times, j_dimer, '-r', linewidth=3, label="J dimer")
# line2, = ax.plot(times, dist_dimer, '--', color='black', alpha=0.3, linewidth=3, label="Independent Emitters")

# # Set axis labels and limits
# ax.set_xlabel(r'$\tau$ (ns)', fontsize=35)
# ax.set_ylabel(r'$g^{(2)}(\infty, \tau)$', fontsize=35)
# ax.set_xlim([-10, 10])
# ax.set_ylim([0.48, 1.3])
# ax.tick_params(axis='x', labelsize=35, direction='in', length=6)
# ax.tick_params(axis='y', labelsize=35, direction='in', length=6)
# ax.legend(fontsize=30, frameon=False)

# # Inset plot
# inset_xlim = [-0.5, 0.5]
# inset_ylim = [0.5, 1.25]
# ax_inset = plt.axes([0.22, 0.6, 0.25, 0.25])  # [left, bottom, width, height]
# ax_inset.plot(times, h_dimer, '-b', linewidth=3)
# ax_inset.plot(times, j_dimer, '-r', linewidth=3)
# # ax_inset.plot(times, dist_dimer, '--', color='black', alpha=0.3, linewidth=3)
# ax_inset.set_xlim(inset_xlim)
# ax_inset.set_ylim(inset_ylim)
# ax_inset.set_xlabel(r'$\tau$ (ns)', fontsize=35)
# ax_inset.set_ylabel(r'$g^{(2)}(\infty, \tau)$', fontsize=35)

# # Set a smaller number of ticks for x and y axes
# ax_inset.xaxis.set_major_locator(MaxNLocator(nbins=3))  # Fewer x ticks
# ax_inset.yaxis.set_major_locator(MaxNLocator(nbins=3))  # Fewer y ticks

# # Increase font size of tick labels
# ax_inset.tick_params(axis='x', labelsize=25, direction='in', length=6)  # bigger than before
# ax_inset.tick_params(axis='y', labelsize=25, direction='in', length=6)
# ax_inset.yaxis.set_major_locator(MaxNLocator(nbins=4))

# # Dashed lines connecting main plot to inset
# from matplotlib.patches import ConnectionPatch

# # Create connection lines
# # Arguments: (x1, y1) in inset coordinates, (x2, y2) in main axes coordinates

# # Bottom-left corner
# con1 = ConnectionPatch(
#     xyA=(inset_xlim[0], inset_ylim[0]), coordsA=ax_inset.transData,
#     xyB=(inset_xlim[0], inset_ylim[0] + 0.015), coordsB=ax.transData,
#     linestyle="--", color="black", linewidth=2
# )
# # Bottom-right corner
# con2 = ConnectionPatch(
#     xyA=(inset_xlim[1], inset_ylim[0]), coordsA=ax_inset.transData,
#     xyB=(inset_xlim[1], inset_ylim[0] + 0.015), coordsB=ax.transData,
#     linestyle="--", color="black", linewidth=2
# )

# ax.add_artist(con1)
# ax.add_artist(con2)
# # Save and show the plot
# plt.savefig('orient_g2_symm', dpi=300, bbox_inches='tight')
# plt.show()






# # Plot the downsampled data
# fig, ax = plt.subplots(figsize=(14, 10))

# # Load data
# h_dimer = np.append(
#     np.flip(np.load('files/g2/H_strongC_2nm_symm_fast_tumb_jigg_normalized_smart_5.npy')[1:]),
#     np.load('files/g2/H_strongC_2nm_symm_fast_tumb_jigg_normalized_smart_5.npy')
# )
# j_dimer = np.append(
#     np.flip(np.load('files/g2/J_strongC_2nm_symm_fast_tumb_jigg_normalized_smart_5.npy')[1:]),
#     np.load('files/g2/J_strongC_2nm_symm_fast_tumb_jigg_normalized_smart_5.npy')
# )
# times = np.append(
#     -np.flip(np.linspace(0, 1.0e-3 * tf, 10000)[1:] - np.linspace(0, 1.0e-3 * tf, 10000)[0]),
#     np.linspace(0, 1.0e-3 * tf, 10000)
# )

# # Main plot
# line0, = ax.plot(times, h_dimer, '-b', linewidth=3)
# line1, = ax.plot(times, j_dimer, '-r', linewidth=3)

# # Set axis labels and limits
# ax.set_xlabel(r'$\tau$ (ns)', fontsize=35)
# ax.set_ylabel(r'$g^{(2)}(\infty, \tau)$', fontsize=35)
# # ax.set_xlim([-10, 10])
# # ax.set_ylim([0.4, 1.3])
# ax.tick_params(axis='x', labelsize=35, direction='in', length=6)
# ax.tick_params(axis='y', labelsize=35, direction='in', length=6)
# ax.legend(handles=[line0, line1], labels=[r"H dimer", r"J dimer"], fontsize=30, frameon=False)

# plt.show()





















# # Create a Bloch-Redfield calculator instance
# BR_calculator = BlochRedfieldCalculator(H, a_ops_opt, Gamma_opt, rates, 
#                                         a_ops_vib=a_ops_vib, Gamma_vib=Gamma_vib, 
#                                         a_ops_coup=a_ops_coup, Gamma_coup=Gamma_coup, 
#                                         P_vib=P_vib, mode=mode, pump=pump)

# # Initialize the DistributionAnalysis class for calculating intensity
# quantum_analysis_g2_dist = DistributionAnalysis(BR_calculator, rates,
#                                                 quantity='g2 distribution', rho=rho0, 
#                                                 final_time=tf, steps=steps)

# # Case 1: G^{(2)}(π/2, φ_p) for 45-degree dimer
# n_phi_p = 100
# phi_vals_p = np.linspace(0, 2 * np.pi, n_phi_p)
# theta_fixed = np.pi / 2  # Fix theta at π/2
# g2_phi = np.zeros_like(phi_vals_p)

# for l, phi_p in enumerate(phi_vals_p):
#     g2_phi[l] = quantum_analysis_g2_dist.compute(theta_fixed, phi_p, theta_fixed, phi_p)
#     print(f'Completed iteration {l + 1} of {n_phi_p}')

# # Normalize data for phi
# g2_phi_normalized = g2_phi  # / np.max(g2_phi)

# # Case 2: G^{(2)}(θ_p, 0) for 45-degree dimer
# n_theta_p = 100
# theta_vals_p = np.linspace(-np.pi, np.pi, n_theta_p)
# phi_fixed = 0  # Fix phi at 0
# g2_theta = np.zeros_like(theta_vals_p)

# for k, theta_p in enumerate(theta_vals_p):
#     g2_theta[k] = quantum_analysis_g2_dist.compute(theta_p, phi_fixed, theta_p, phi_fixed)
#     print(f'Completed iteration {k + 1} of {n_theta_p}')

# # Normalize data for theta
# g2_theta_normalized = g2_theta  # / np.max(g2_theta)

# # Common legend configuration
# legend_fontsize = 50  # Font size for the legend
# tick_fontsize = 55    # Font size for tick labels

# # # Plot 1: Unnormalized data
# # fig1, ax1 = plt.subplots(figsize=(12, 9), subplot_kw={'projection': 'polar'})
# # ax1.plot(phi_vals_p, g2_phi, color='blue', linewidth=3, label=r'$G^{(2)}(\pi/2, \phi^{\prime})$')
# # ax1.plot(theta_vals_p, g2_theta, color='red', linewidth=3, linestyle='--', label=r'$G^{(2)}(\theta^{\prime}, 0)$')

# # Customize tick labels
# ticks = np.linspace(0, 2 * np.pi, 8, endpoint=False)
# tick_labels = [r'$0$', r'$\frac{\pi}{4}$', r'$\frac{\pi}{2}$', r'$\frac{3\pi}{4}$',
#                r'$\pi$', r'$\frac{5\pi}{4}$', r'$\frac{3\pi}{2}$', r'$\frac{7\pi}{4}$']
# # ax1.set_thetagrids(np.degrees(ticks), labels=tick_labels, fontsize=tick_fontsize)
# # ax1.set_rticks([0.5, 1])
# # ax1.tick_params(axis='both', labelsize=65, pad=40)

# # # Combine handles and labels for the common legend
# # handles1, labels1 = ax1.get_legend_handles_labels()
# # fig1.legend(
# #     handles1, labels1,
# #     loc='upper center',
# #     bbox_to_anchor=(0.5, 1.05),
# #     fontsize=legend_fontsize,
# #     frameon=False,
# #     ncol=2,
# #     handlelength=1,
# #     labelspacing=1,
# #     markerscale=1.5,
# # )

# # # Save each figure and show plots
# # fig1.savefig('unnorm_45_polar_radiation.png', dpi=300, bbox_inches='tight')

# # Plot 2: Normalized data
# fig2, ax2 = plt.subplots(figsize=(12, 9), subplot_kw={'projection': 'polar'})
# ax2.plot(phi_vals_p, g2_phi_normalized, color='blue', linewidth=3, label=r'$g^{(2)}(\pi/2, \phi^{\prime})$')
# ax2.plot(theta_vals_p, g2_theta_normalized, color='red', linewidth=3, linestyle='--', label=r'$g^{(2)}(\theta^{\prime}, 0)$')

# # Customize tick labels
# ax2.set_thetagrids(np.degrees(ticks), labels=tick_labels, fontsize=tick_fontsize)
# ax2.set_rticks([0.5, 1.0])
# ax2.tick_params(axis='both', labelsize=65, pad=40)

# # Add legend
# handles2, labels2 = ax2.get_legend_handles_labels()
# fig2.legend(
#     handles2, labels2,
#     loc='upper center',
#     bbox_to_anchor=(0.5, 1.05),
#     fontsize=legend_fontsize,
#     frameon=False,
#     ncol=2,
#     handlelength=1,
#     labelspacing=1,
#     markerscale=1.5,
# )

# # Save normalized plot
# fig2.savefig('norm_45_polar_radiation.png', dpi=300, bbox_inches='tight')

# plt.show()




# omega = 11.5

# g2_ps = np.load('files/g2/45_strongC_2nm_symm_d1.npy')[0:100]
# tf = 5 * system_params.tau_L
# steps = len(np.load('files/g2/45_strongC_2nm_symm_d1.npy')) #5000000
# times_ps = np.linspace(0, tf, steps)[0:100]

# plt.plot(times_ps, g2_ps)
# plt.plot(times_ps, 0.82 * np.sin(omega * times_ps)**2)
# plt.plot(times_ps, 0.82 * np.sin(system_params.J * kappa**2 * times_ps)**2)
# plt.show()




# def plot_combined_with_inset(files, times_ns, colors, labels, inset_xlim_ns, system_params, kappa):
#     # Convert inset x-axis limits to ps
#     inset_xlim_ps = [x * 1000 for x in inset_xlim_ns]

#     # Create figure with preferred formatting
#     fig, ax = plt.subplots(figsize=(14, 10))

#     # Iterate through each file and plot
#     for file, color, label in zip(files, colors, labels):
#         data = np.append(np.flip(np.load(file)[1:]), np.load(file))
#         times_ns_flipped = np.append(-np.flip((times_ns[1:] - times_ns[0])), times_ns)

#         ax.plot(times_ns_flipped, data, color, linewidth=2, label=label)

#     # Set main plot labels and limits
#     ax.set_xlabel(r'$\tau$ (ns)', fontsize=35)
#     ax.set_ylabel(r'$g^{(2)}(\infty, \tau)$', fontsize=35)
#     ax.set_xlim([-2, 2])
#     ax.set_ylim([-0.6, 1.1])
#     ax.tick_params(axis='both', labelsize=30, direction='in', length=6)

#     # # Adjust legend position inside plot
#     # ax.legend(fontsize=30, loc='upper right', frameon=False, handlelength=1)

#     # Inset plot
#     inset_ylim = [-0.5, 1.25]
#     ax_inset = fig.add_axes([0.23, 0.25, 0.25, 0.25])  # [left, bottom, width, height]

#     for file, color in zip(files, colors):
#         data = np.append(np.flip(np.load(file)[1:]), np.load(file))
#         times_ns_flipped = np.append(-np.flip((times_ns[1:] - times_ns[0])), times_ns)

#         ax_inset.plot(times_ns_flipped * 1000, data, color, linewidth=2)

#     # Add periodic plot to the inset
#     times_ps = times_ns * 1000
#     periodic = 1.2 * np.sin(system_params.J * kappa**2 * times_ps / 2)**2
#     periodic_flipped = np.append(np.flip(periodic[1:]), periodic)

#     ax_inset.plot(times_ns_flipped * 1000, periodic_flipped, 'k', linewidth=3, label=r'$\propto \sin(J^{\prime}_{1,2} t / 2)$', alpha=0.2)

#     # Set inset limits and labels
#     ax_inset.set_xlim(inset_xlim_ps)
#     ax_inset.set_ylim(inset_ylim)
#     ax_inset.set_xlabel(r'$\tau$ (ps)', fontsize=30)
#     ax_inset.set_ylabel(r'$g^{(2)}(\infty, \tau)$', fontsize=30)
#     ax_inset.tick_params(axis='both', labelsize=25, direction='in', length=6)
#     ax_inset.yaxis.set_major_locator(MaxNLocator(nbins=4))

#     # Add legend to inset
#     ax_inset.legend(
#     fontsize=25, 
#     frameon=False, 
#     # bbox_to_anchor=(0.5, 0.22),  # (x, y) in axes coordinates
#     handlelength=1)

#     # Save and show the plot
#     plt.savefig('orient_g2_45.png', dpi=300, bbox_inches='tight', pad_inches=0.05)
#     # plt.show()

# # Example usage:
# # Define inset x-limits in ns
# inset_xlim = [-0.0005, 0.0005]  

# # File paths, colors, and labels
# files = [
#     'files/g2/45_strongC_2nm_symm_d1.npy',
#     'files/g2/45_strongC_2nm_symm_d2_perp.npy',
#     'files/g2/45_strongC_2nm_symm_z.npy'
# ]
# colors = ['b', 'r', 'g']
# labels = [
#     r'$\textbf{q} = (\pm 1 , 0, 0)$',
#     r'$\textbf{q} = (\pm 1 , \mp 1, 0)$',
#     r'$\textbf{q} = (0 , 0, \pm 1)$'
# ]

# # # Parameters for the periodic plot
# # class SystemParams:
# #     J = 1.5  # Example value for J

# # kappa = 0.8  # Example value for kappa
# # times_ns = np.linspace(-0.002, 0.002, 1000)  # Example time range

# # Generate the plot
# plot_combined_with_inset(files, times_ns, colors, labels, inset_xlim, system_params, kappa)




# # Load your data
# g2_corr_H0 = np.load('files/g2_temp_H.npy')
# g2_corr_J0 = np.load('files/g2_temp_J.npy')

# # Simulated x-axis (Temperature)
# T_vib = np.linspace(0, 300, 500)

# # Create the figure
# fig, ax = plt.subplots(figsize=(16, 10))

# # Plot both H and J dimer data on the same axis
# H_line, = ax.plot(T_vib, g2_corr_H0, color='blue', label='H dimer', linewidth=3)
# J_line, = ax.plot(T_vib, g2_corr_J0, color='red', label='J dimer', linewidth=3)

# # Set axis labels and styling
# ax.set_xlabel(r'Temperature (K)', fontsize=40)
# ax.set_ylabel(r'$g^{(2)}(\infty, 0)$', fontsize=40)
# ax.tick_params(axis='both', labelsize=35, direction='in', length=6, pad=10)

# # Set x and y axes to log scale
# ax.set_xscale('log')
# ax.set_yscale('log')

# # Set axis limits
# ax.set_xlim(1, 300)
# ax.set_ylim(
#     min(min(g2_corr_H0), min(g2_corr_J0)) * 0.9,
#     max(max(g2_corr_H0), max(g2_corr_J0)) * 1.1
# )

# # Add legend
# ax.legend(handles=[H_line, J_line], loc='lower left', fontsize=40, frameon=False)

# # Show plot
# plt.show()












