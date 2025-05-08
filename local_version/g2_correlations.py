from system_ops import *
from spectral_dens import *
from br_tensor import *
from scipy.linalg import expm, null_space




def calculate_rho_qr(H, rho0, Gamma_coup, Gamma_opt, P_vib, theta, phi, k_q1):

    sx, sz, sp, sm = pauli_matrices()
    sx1, sx2, sz1, sz2, sp1, sp2, sm1, sm2 = two_qubit_pauli_matrices()

    N = steps
    # Integrate over times
    rho_qr = np.empty_like(rho0, dtype = complex) # correlation array
     
    # Coupling operator
    a_ops_coup = [sp1 @ sm2, sm1 @ sp2, sp2 @ sm1, sm2 @ sp1]
    a_ops_opt = [sp1, sm1, sp2, sm2]

    # Bloch-Redfield tensor
    R_coup, ekets = BR_tensor_coup(H, a_ops_coup, Gamma_coup)

    R_opt, ekets = BR_tensor_opt(H, a_ops_opt, Gamma_opt, P_vib)
    
    # calculate propagator
    P = expm((R_opt + R_coup) * dt) 

    # initialize state
    rho_ss = rho0
    dim = len(rho_ss[0])

    # Find the steady-state (null space)
    null = null_space(R_opt + R_coup)
    rho_ss = np.reshape(null, (dim, dim) ) 
    rho_ss /= np.trace(rho_ss)
    
    # k = k_q1 * np.array([np.cos(phi) * np.cos(theta), np.sin(phi) * np.cos(theta), -np.sin(theta)])
    
    e_1_vec = np.array([np.cos(phi) * np.cos(theta), np.sin(phi) * np.cos(theta), -np.sin(theta)])
    e_2_vec = np.array([np.sin(phi), -np.cos(phi), 0])
    e_vec = np.array([e_1_vec, e_2_vec])

    psi_g1g2 = np.array([0, 0, 0, 1])
    psi_e1g2 = np.array([0, 0, 1, 0])
    psi_e2g1 = np.array([0, 1, 0, 0])
    psi_e1e2 = np.array([1, 0, 0, 0])

    def compute_rho_qr(i):
        d1_pol = d_1_vec @ e_vec[i]
        d2_pol = d_2_vec @ e_vec[i]
    
        Norm = np.sqrt(d1_pol**2 + d2_pol**2) + 1.0e-35
        psi_k = ((d1_pol * np.exp(-1j * k_q1 * r * np.cos(theta) / 2)) * psi_e1g2 + 
                 (d2_pol * np.exp(1j * k_q2 * r * np.cos(theta) / 2)) * psi_e2g1) / (Norm)
    
        sk_m = ekets.conj().T @ (np.outer(psi_g1g2, psi_k) + np.outer(psi_k, psi_e1e2)) @ ekets
        sk_p = sk_m.conj().T
    
        rho_qr_0 = sk_m @ rho_ss @ sk_p

        return rho_qr_0
    
    rho_qr += Parallel(n_jobs=-1)(delayed(compute_rho_qr)(i) for i in range(len(e_vec)))
    
    return rho_qr

# def calculate_steady_state_g2_correlations(H, rho0, Gamma_coup, Gamma_opt, P_vib, theta, phi, theta_p, phi_p, k_q1):
    
#     sx, sz, sp, sm = pauli_matrices()
#     sx1, sx2, sz1, sz2, sp1, sp2, sm1, sm2 = two_qubit_pauli_matrices()
#     N = steps
#     # Integrate over times
#     g2_corr = 0 #np.zeros(N, dtype = complex) # correlation array
     

#     # Coupling operator
#     a_ops_coup = [sp1 @ sm2, sm1 @ sp2, sp2 @ sm1, sm2 @ sp1]
#     a_ops_opt = [sp1, sm1, sp2, sm2]

#     # Bloch-Redfield tensor
#     R_coup, ekets = BR_tensor_coup(H, a_ops_coup, Gamma_coup)

#     R_opt, ekets = BR_tensor_opt(H, a_ops_opt, Gamma_opt, P_vib)
    
#     # calculate propagator
#     P = expm((R_opt + R_coup) * dt) 

#     # initialize state
#     rho_ss = rho0
#     dim = len(rho_ss[0])

#     # # propagate for a_ops to steady state
#     # for _ in times_ss:   
#     #     rho_vec_ss = np.reshape(rho_ss, (dim**2, 1))
#     #     rho_vec_ss = P @ rho_vec_ss
#     #     rho_ss = np.reshape(rho_vec_ss, (dim, dim))   
    
#     # Find the steady-state (null space)
#     null = null_space(R_opt + R_coup)
#     rho_ss = np.reshape(null, (dim, dim) ) 
#     rho_ss /= np.trace(rho_ss)

#     e_1_vec = np.array([np.cos(phi) * np.cos(theta), np.sin(phi) * np.cos(theta), -np.sin(theta)])
#     e_2_vec = np.array([np.sin(phi), -np.cos(phi), 0])
#     e_vec = np.array([e_1_vec, e_2_vec])

#     psi_g1g2 = np.array([0, 0, 0, 1])
#     psi_e1g2 = np.array([0, 0, 1, 0])
#     psi_e2g1 = np.array([0, 1, 0, 0])
#     psi_e1e2 = np.array([1, 0, 0, 0])

#     def compute_g2(i):
#         d1_pol = d_1_vec @ e_vec[i]
#         d2_pol = d_2_vec @ e_vec[i]
    
#         Norm = np.sqrt(d1_pol**2 + d2_pol**2) + 1.0e-35
#         psi_k = ((d1_pol * np.exp(-1j * k_q1 * r * np.cos(theta) / 2)) * psi_e1g2 + 
#                  (d2_pol * np.exp(1j * k_q2 * r * np.cos(theta) / 2)) * psi_e2g1) / (Norm)
    
#         sk_m = ekets.conj().T @ (np.outer(psi_g1g2, psi_k) + np.outer(psi_k, psi_e1e2)) @ ekets
#         sk_p = sk_m.conj().T
    
#         rho_qr = sk_m @ rho_ss @ sk_p
    
#         # g2 = []

#         # for t in times:
#         #     g2.append(np.real(np.trace(np.dot(sk_p @ sk_m, rho_qr))))
#         #     rho_vec = np.reshape(rho_qr, (dim**2, 1))
#         #     rho_vec = P @ rho_vec
#         #     rho_qr = np.reshape(rho_vec, (dim, dim))
            
#         # g2.append(np.real(np.trace(np.dot(sk_p @ sk_m, rho_qr))))

#         g2 = np.real(np.trace(np.dot(sk_p @ sk_m, rho_qr)))
#         return g2 #/ (np.real(np.trace(np.dot(sk_p @ sk_m, rho_ss))))**2
    
#     g2_corr += Parallel(n_jobs=-1)(delayed(compute_g2)(i) for i in range(len(e_vec)))

#     # # Parallelize the loop
#     # results = Parallel(n_jobs=-1)(delayed(compute_g2)(i) for i in range(len(e_vec)))

#     # # Sum up the results
#     # g2_corr = sum(results)

#     return g2_corr

