# # from params import *
# from sys_params import *




# def calculate_spectral_densities():
#     """
#     Calculates the spectral densities for optical and vibrational modes.

#     Returns:
#         A tuple containing two functions:
#             - J_opt(w): Spectral density for the optical mode.
#             - J_vib(w): Spectral density for the vibrational mode.
#     """
    
#     constants = Constants()
#     dipole_params = DipoleParameters(config_name)  
#     system_params = SystemParameters(dipole_params=dipole_params, E_q1=E_q1*constants.eV, E_q2=E_q2*constants.eV)
#     system_params.calculate_parameters()

#     # Define spectral density functions
#     def J_opt(w):
#         k_opt_1 = 1 / system_params.tau_L 
#         return k_opt_1 
#         # return (w**3 * ((system_params.d_1)**2)) / (3 * np.pi * system_params.eps_0 * system_params.hbar_ps * (system_params.c_ps**3)) 

#     # omega_c = 900.0e-3 * constants.eV / system_params.hbar_ps
#     # E_reorg = 5.0e-6 * constants.eV / system_params.hbar_ps
#     # omega_c = 90.0e-3 * constants.eV / system_params.hbar_ps
#     # E_reorg = 500.0e-4 * constants.eV / system_params.hbar_ps
    
#     omega_c = 90.0e-3 * constants.eV / system_params.hbar_ps
#     # E_reorg = 5.0e-3 * constants.eV / system_params.hbar_ps

#     E_reorg = 62.0e-3 * constants.eV / system_params.hbar_ps
#     def J_vib(w):
#         k_vib_2 = E_reorg / (2 * omega_c**3)
#         # if w==0:
#         #     return 1*k_vib_2 # '''Somewhat arbitrary'''
#         # else:
#         #     return (k_vib_2 * np.abs(w)**3) * np.exp(- (np.abs(w) / omega_c))


#         # # return 1.0e3 * gamma_q1 #
#         return (k_vib_2 * np.abs(w)**3) * np.exp(- (np.abs(w) / omega_c))
#         # # return (k_vib_2 * (w)**3) * np.exp(- (w / omega_c)**2)

#     # omega_c = mp.mpf(90.0e-3 * constants.eV / system_params.hbar_ps)
#     # E_reorg = mp.mpf(20.0e-3 * constants.eV / system_params.hbar_ps)

#     # def J_vib(w):
#     #     k_vib_2 = E_reorg / (2 * omega_c**3)
        
   
#     #     if isinstance(w, np.ndarray):  # Check if w is a NumPy array
#     #         return np.array([k_vib_2 * mp.fabs(w_i)**3 * mp.exp(- (mp.fabs(w_i) / omega_c)) for w_i in w])
#     #     else:
#     #         return k_vib_2 * mp.fabs(w)**3 * mp.exp(- (mp.fabs(w) / omega_c))

#     return J_opt, J_vib

# # J_opt, J_vib = calculate_spectral_densities()

# # print(J_opt)




# # array([[-2824.93047719+0.j,     0.        +0.j,     0.        +0.j,
# #             0.        +0.j],
# #        [    0.        +0.j, -1412.46523859+0.j,   -46.4094714 +0.j,
# #             0.        +0.j],
# #        [    0.        +0.j,   -46.4094714 +0.j, -1412.46523859+0.j,
# #             0.        +0.j],
# #        [    0.        +0.j,     0.        +0.j,     0.        +0.j,
# #             0.        +0.j]])











# # ############################################################
# # '''Temperatures'''
# # ############################################################

# # if T_opt == 0.0:
# #     beta_opt = np.inf
# # elif T_opt == np.inf:
# #     beta_opt = 0.0
# # else:
# #     beta_opt = 1/(k_B_Js*T_opt) 


# # if T_vib == 0.0:
# #     beta_vib = np.inf
# # elif T_vib == np.inf:
# #     beta_vib = 0.0
# # else:
# #     beta_vib = 1/(k_B_Js * T_vib)     
    

# # ############################################################
# # '''Population'''
# # ############################################################

# # n_inv = lambda w, beta: (np.exp((constants.hbar * np.abs(w)) * beta) - 1)     # Bose Einstein Distribution


# # ############################################################
# # '''Spectral Densities'''
# # ############################################################

# # # Spectral Density for the Optical Modes
# # def J_opt(w):
# #     k_opt_1 = gamma_q1  
# #     return k_opt_1 


# # # Spectral Density for the Vibrational Modes
# # omega_c = 90e-3 * eV / hbar
# # E_reorg = 5.0e-3 * eV / hbar
# # def J_vib(w):
# #     k_vib_2 = E_reorg / (2 * omega_c**3)
# #     return (k_vib_2 * np.abs(w)**3) * np.exp(- (np.abs(w) / omega_c))


# # def Phi(t):
# #     integrand = lambda w: J_vib(w) * (np.cos(w * t) * (1 / np.tanh(beta_vib * w / 2)) - 1j * np.sin(w * t)) / w**2
# #     omega =  np.linspace(1.0e-4, 10000, 10000)
# #     inte = [integrand(w) for w in omega]
# #     result = simpson(inte, omega)
# #     return result 





# # def Phi(t):
# #     result, _ = quad_vec(lambda w: J_vib(w) * (np.cos(w * t) * (1 / np.tanh(beta_vib * w / 2)) - 1j * np.sin(w * t)) / w**2, 1.0e-4, np.inf)
# #     return np.array(result)


