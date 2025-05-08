from sys_params import *
from system_ops import *
from rates import *


class BlochRedfieldCalculator:
    def __init__(self, H, constants, dipole_params, system_params, a_ops_opt, Gamma_opt, rates, a_ops_vib=None, Gamma_vib=None, a_ops_coup=None, Gamma_coup=None, P_vib=None, mode=None, pump = None):
        self.H = H
        self.a_ops_opt = a_ops_opt
        self.Gamma_opt = Gamma_opt
        self.a_ops_vib = a_ops_vib
        self.Gamma_vib = Gamma_vib
        self.a_ops_coup = a_ops_coup
        self.Gamma_coup = Gamma_coup
        self.P_vib = P_vib
        self.mode = mode
        self.pump = pump
        self.dim = len(self.H)
        self.evals, self.ekets = np.linalg.eigh(self.H)
        self.constants = constants
        self.dipole_params = dipole_params
        self.system_params = system_params
        self.r_hat, self.d1_hat, self.d2_hat = self.dipole_params.get_dipole_parameters()


        # Calculate spectral densities and rates
        self.rates = rates #CombinedRates(constants, dipole_params, system_params, T_opt, T_vib, mode)
        self.kappa = self.rates.kappa()

        sx1, sx2, sz1, sz2, sp1, sp2, sm1, sm2 = two_qubit_pauli_matrices()

        if self.pump == 'sym' or self.pump == 'site':
            # self.c_ops_S = [0 * ((sp1 + sp2)), np.sqrt(system_params.gamma_p_1) * ((sp1 - sp2))]
            # self.c_ops_S = [np.sqrt(system_params.gamma_p_1) * ((sp1 + sp2)), 0 * ((sp1 - sp2))]
            self.c_ops_S = [np.sqrt(system_params.gamma_p_1) * sp1, np.sqrt(system_params.gamma_p_2) * sp2]
        elif self.pump == 'none':
            self.c_ops_S = [0 * sp1, 0 * sp2]

        self.P_lv = np.full((2,2), self.kappa**4)

    def _calc_a_ops_freq(self, H, a_ops):
        evals, ekets = np.linalg.eigh(H)
        a_ops_S = [ekets.conj().T @ a_op @ ekets for a_op in a_ops]
        trans_freq = np.array([])
        # for i in range(len(evals)):
        #     for j in range(len(evals)):
        #         if evals[i] - evals[j] == 0:
        #             trans_freq = np.append(trans_freq, 1.0e0)
        #         else:
        #             trans_freq = np.append(trans_freq, evals[i] - evals[j])
        trans_freq = np.array([evals[i] - evals[j] for i in range(len(evals)) for j in range(len(evals))])
        return a_ops_S, trans_freq

    def _Gamma_gen(self, Gamma, trans_freq):
        G = np.zeros((2, 2, len(trans_freq)))
        for (i, j), G_ab in np.ndenumerate(Gamma):
            G[i, j, :] = np.array([G_ab(w) for w in trans_freq])
        return G

    def _Gamma_pn(self, Gamma_coup, trans_freq):
        G = np.zeros((2, 2, len(trans_freq)), dtype=complex)

        def compute_G_ab(G_ab):
            with ThreadPoolExecutor(max_workers=10) as executor:
                return np.array(list(executor.map(G_ab, trans_freq)))

        for (i, j), G_ab in np.ndenumerate(Gamma_coup):
            G[i, j, :] = compute_G_ab(G_ab)
        print(G)
        return G

    # def _Gamma_pn(self, Gamma_coup, trans_freq):
    #     G = np.zeros((2, 2, len(trans_freq)), dtype = complex)
    #     for (i, j), G_ab in np.ndenumerate(Gamma_coup):
    #         G[i, j, :] = np.array([G_ab(w) for w in trans_freq])
    #     return G
    
    # def _Gamma_pn(self, Gamma_coup, trans_freq):
    #     G = np.zeros((2, 2, len(trans_freq)), dtype = complex)
    #     for (i, j), G_ab in np.ndenumerate(Gamma_coup):
    #         G[i, j, :] = np.array([G_ab(trans_freq)])
    #     return G

    def Unitary_Dynamics(self):
        I = np.eye(self.dim)
        trans_freq = np.array([self.evals[i]-self.evals[j] for i in range(self.dim) for j in range(self.dim)])
        # superH = -1j * (np.kron(I, self.H) - np.kron(self.H.T, I))
        # unitary part
        R_unitary = -1j * np.diag(trans_freq)
        return R_unitary
        # return superH

    # def Liouvillian(self):
    #     if self.c_ops_S is None:
    #         raise ValueError("Provide collapse operators")

    #     c_ops = self.ekets.conj().T @ self.c_ops_S @ self.ekets
    #     I = np.eye(self.dim)
    #     superL = sum([np.kron(c_op.conj(),c_op) 
    #               - 1/2 * ( np.kron(I, np.dot(c_op.conj().T, c_op)) +
    #                         np.kron(np.dot(c_op.T, c_op.conj()), I) 
    #                       ) for c_op in c_ops])
    #     return superL


    def Liouvillian(self):
        if self.c_ops_S is None:
            raise ValueError("Provide collapse operators")

        c_ops = self.ekets.conj().T @ self.c_ops_S @ self.ekets
        I = np.eye(self.dim)
        
        if self.pump == 'sym':
            superL = np.zeros((self.dim**2, self.dim**2), dtype=complex)

            for x in range(len(c_ops)):
                c_a = c_ops[x]
                for y in range(len(c_ops)):
                    c_b = c_ops[y]

                    P_lv = self.P_lv[x, y]

                    # Build the Liouvillian superoperator
                    superL += (np.kron(c_a.conj(), c_b) - 1/2 * (np.kron(I, np.dot(c_a.conj().T, c_b)) + np.kron(np.dot(c_a.T, c_b.conj()), I))) * P_lv 

        if self.pump == 'site' or self.pump == 'none':
            superL = sum([np.kron(c_op.conj(),c_op) 
                  - 1/2 * ( np.kron(I, np.dot(c_op.conj().T, c_op)) +
                            np.kron(np.dot(c_op.T, c_op.conj()), I) 
                          ) for c_op in c_ops])
        
                
        return superL

    def mathcal_F(self, k_q_r, alpha, beta, threshold=1e-8):
        
        mathcal_F = (alpha * ((np.sin(k_q_r)) / (k_q_r)) + 
                    beta * (((np.cos(k_q_r)) / (k_q_r)**2) - ((np.sin(k_q_r)) / (k_q_r)**3))) * (3 / 2)

        # Handle the case where k_q_r is very small (to avoid division by zero or extremely small values).
        mathcal_F[np.abs(k_q_r) < threshold] = np.dot(self.system_params.d1_hat, self.system_params.d2_hat)

        # mathcal_F = np.dot(system_params.d1_hat, system_params.d2_hat)
               
        # Create mathcal_F_prime as a matrix of ones
        mathcal_F_prime = np.ones_like(mathcal_F)
        
        # Construct the F tensor as specified
        F = np.array([[mathcal_F_prime, mathcal_F_prime, mathcal_F, mathcal_F],
                    [mathcal_F_prime, mathcal_F_prime, mathcal_F, mathcal_F],
                    [mathcal_F, mathcal_F, mathcal_F_prime, mathcal_F_prime],
                    [mathcal_F, mathcal_F, mathcal_F_prime, mathcal_F_prime]])
        
        return F

    def BR_tensor_opt(self):
        a_ops_S_opt, trans_freq = self._calc_a_ops_freq(self.H, self.a_ops_opt)
        Gamma_opt_0 = self._Gamma_gen(self.Gamma_opt, trans_freq)

        R_opt = np.zeros((self.dim**2, self.dim**2), dtype=complex)
        
        for x in range(len(a_ops_S_opt)):
            A_a = a_ops_S_opt[x]
            for y in range(len(a_ops_S_opt)):
                A_b = a_ops_S_opt[y]

                Gamma_Opt = np.tile(Gamma_opt_0, (2, 2, 1))
                G_ab = Gamma_Opt[x, y].reshape(self.dim, self.dim)

                tf = trans_freq.reshape(self.dim, self.dim)

                k_q = tf / self.system_params.c_ps

                alpha = (np.dot(self.system_params.d1_hat, self.system_params.d2_hat) - 
                        (np.dot(self.system_params.d1_hat, self.system_params.r_hat)) * 
                        (np.dot(self.system_params.d2_hat, self.system_params.r_hat)))
                beta = (np.dot(self.system_params.d1_hat, self.system_params.d2_hat) - 
                        3 * (np.dot(self.system_params.d1_hat, self.system_params.r_hat)) * 
                        (np.dot(self.system_params.d2_hat, self.system_params.r_hat)))

                k_q_r = k_q * self.system_params.r + 1e-16

                # Calculate F using the mathcal_F function
                F = self.mathcal_F(k_q_r, alpha, beta)

                G_ab = G_ab * F[x, y]
                
                P_ab = self.P_vib[x, y] if self.mode == 'polaron' else 1.0

                R_opt += -1/2 * (np.kron(np.dot(A_a, np.multiply(A_b, np.transpose(G_ab))), np.transpose(np.eye(self.dim, self.dim)))
                        - np.kron(np.multiply(np.transpose(G_ab), A_a), np.transpose(A_b))
                        + np.kron(np.eye(self.dim, self.dim), np.transpose(np.dot(np.multiply(A_a, G_ab), A_b)))
                        - np.kron(A_a, np.transpose(np.multiply(A_b, G_ab)))) * np.dot(self.d1_hat, self.d2_hat) * P_ab 
                        
        return R_opt

    def BR_tensor_vibcoup(self):

        if self.mode == 'weak':
          a_ops_S, trans_freq = self._calc_a_ops_freq(self.H, self.a_ops_vib)
          Gamma = self._Gamma_gen(self.Gamma_vib, trans_freq)

        elif self.mode == 'polaron':
          a_ops_S, trans_freq = self._calc_a_ops_freq(self.H, self.a_ops_coup)
          trans_freq[(np.abs(trans_freq) > 1000)] = np.nan
          Gamma = self._Gamma_pn(self.Gamma_coup, np.array(trans_freq))


        R = np.zeros((self.dim**2, self.dim**2), dtype=complex)
        for x in range(len(a_ops_S)):
            A_a = a_ops_S[x]
            for y in range(len(a_ops_S)):
                A_b = a_ops_S[y]
                
                G_ab = Gamma[x, y].reshape(self.dim, self.dim)

                R += -1/2 * (np.kron(np.dot(A_a, np.multiply(A_b, np.transpose(G_ab))), np.transpose(np.eye(self.dim, self.dim)))
                        - np.kron(np.multiply(np.transpose(G_ab), A_a), np.transpose(A_b))
                        + np.kron(np.eye(self.dim, self.dim), np.transpose(np.dot(np.multiply(A_a, G_ab), A_b)))
                        - np.kron(A_a, np.transpose(np.multiply(A_b, G_ab))))
  
        return R

    def total_BR_tensor(self):
        if self.mode == 'weak' and self.a_ops_vib is None:
            return self.Unitary_Dynamics() + self.BR_tensor_opt()
        if self.c_ops_S is None:
            return self.Unitary_Dynamics() + self.BR_tensor_opt() + self.BR_tensor_vibcoup()
        else:
            return self.Unitary_Dynamics() + self.BR_tensor_opt() + self.BR_tensor_vibcoup() + self.Liouvillian()









































# ############################################################
# '''Liouvillian'''
# ############################################################

# def Liouvillian(H, c_ops_S, ekets):   # There's no point in transforming c_ops, they're in the eigenbasis of the diag H.
#     # print(c_ops)
#     # print(ekets.conj().T @ c_ops @ ekets)
#     c_ops = ekets.conj().T @ c_ops_S @ ekets
#     dim = len(H[0]) # dimension of the system
#     I = np.eye(dim) # Identity matrix of the order of dim

#     superH = -1j * ( np.kron(I,H)-np.kron(H.T,I) ) # Hamiltonian part

#     superL = sum([np.kron(c_op.conj(),c_op) 
#                   - 1/2 * ( np.kron(I, np.dot(c_op.conj().T, c_op)) +
#                             np.kron(np.dot(c_op.T, c_op.conj()), I) 
#                           ) for c_op in c_ops])
#     return superL + superH 



# ############################################################
# '''Bloch-Redfield tensor in the Hamiltonian's basis'''
# ############################################################

# def BR_tensor_coup(H, a_ops_coup, Gamma_coup):

#     secular_cut_off = 0.01

#     dim = len(H) # dimension
#     evals,ekets = np.linalg.eigh(H) # HS's basis

#     # coupling operators in H basis and transition frequencies (w_ab)
#     a_ops_S_coup, trans_freq = calc_a_ops_freq(H, a_ops_coup)
    
#     # construct empty R
#     R = np.zeros((dim**2,dim**2),dtype = complex) 

#     Gamma = Gamma_gen(Gamma_coup, trans_freq, dim)
    
#     # unitary part
#     R = -1j * np.diag(trans_freq)

#     for x, in np.ndindex(dim): # loop over uncorrelated a_ops
#         A_a = a_ops_S_coup[x]

#         for y, in np.ndindex(dim):
#             A_b = a_ops_S_coup[y]

#             Gamma_c = np.tile(Gamma, (2, 2, 1))

#             G_ab = Gamma_c[x, y].reshape(dim, dim)


#             R+= - 1/2 * (np.kron(np.dot(A_a, np.multiply(A_b, np.transpose(G_ab))), np.transpose(np.eye(dim, dim)))
#                         - np.kron(np.multiply(np.transpose(G_ab), A_a), np.transpose(A_b))
#                         + np.kron(np.eye(dim, dim), np.transpose(np.dot(np.multiply(A_a, G_ab), A_b)))
#                         - np.kron(A_a, np.transpose(np.multiply(A_b, G_ab)))) 
            
#     return R, ekets



# def BR_tensor_opt(H, a_ops_opt, Gamma_opt, P_vib):

#     secular_cut_off = 0.01

#     dim = len(H) # dimension
#     evals,ekets = np.linalg.eigh(H) # HS's basis

#     # coupling operators in H basis and transition frequencies (w_ab)
#     a_ops_S_opt, trans_freq = calc_a_ops_freq(H, a_ops_opt)

#     Gamma = Gamma_gen(Gamma_opt, trans_freq, dim)
    
#     # construct empty R
#     R = np.zeros((dim**2,dim**2),dtype = complex) 

#     for x, in np.ndindex(dim): # loop over uncorrelated a_ops
#         A_a = a_ops_S_opt[x]

#         for y, in np.ndindex(dim):
#             A_b = a_ops_S_opt[y]

#             Gamma_o = np.tile(Gamma, (2, 2, 1))

#             G_ab = Gamma_o[x, y].reshape(dim, dim)
#             P_ab = P_vib[x, y]

#             R+= - 1/2 * (np.kron(np.dot(A_a, np.multiply(A_b, np.transpose(G_ab))), np.transpose(np.eye(dim, dim)))
#                         - np.kron(np.multiply(np.transpose(G_ab), A_a), np.transpose(A_b))
#                         + np.kron(np.eye(dim, dim), np.transpose(np.dot(np.multiply(A_a, G_ab), A_b)))
#                         - np.kron(A_a, np.transpose(np.multiply(A_b, G_ab)))) * np.dot(d1_hat, d2_hat) * P_ab
    
#     return R, ekets





# # compute the Bloch-Redfield tensor in the Hamiltonian's basis
# def BR_tensor_opt(H, a_ops_opt, Gamma_opt):

#     secular_cut_off = 0.01

#     dim = len(H) # dimension
#     evals,ekets = np.linalg.eigh(H) # HS's basis

#     # coupling operators in H basis and transition frequencies (w_ab)
#     a_ops_S_opt, trans_freq = calc_a_ops_freq(H, a_ops_opt)

#     Gamma = Gamma_gen(Gamma_opt, trans_freq)
    
#     # construct empty R
#     R = np.zeros((dim**2,dim**2),dtype = complex) 

#     # unitary part
#     R = -1j * np.diag(trans_freq)

#     for x, in np.ndindex(dim): # loop over uncorrelated a_ops
#         A_a = a_ops_S_opt[x]

#         for y, in np.ndindex(dim):
#             A_b = a_ops_S_opt[y]

#             Gamma_o = np.tile(Gamma, (2, 2, 1))

#             G_ab = Gamma_o[x, y].reshape(dim, dim)

#             R+= - 1/2 * (np.kron(np.dot(A_a, np.multiply(A_b, np.transpose(G_ab))), np.transpose(np.eye(dim, dim)))
#                         - np.kron(np.multiply(np.transpose(G_ab), A_a), np.transpose(A_b))
#                         + np.kron(np.eye(dim, dim), np.transpose(np.dot(np.multiply(A_a, G_ab), A_b)))
#                         - np.kron(A_a, np.transpose(np.multiply(A_b, G_ab)))) * np.dot(d1_hat, d2_hat)
    
#     return R, ekets



# # compute the Bloch-Redfield tensor in the Hamiltonian's basis
# def BR_tensor_vib(H, a_ops_vib, Gamma_vib):

#     secular_cut_off = 0.01

#     dim = len(H) # dimension
#     evals,ekets = np.linalg.eigh(H) # HS's basis

#     # coupling operators in H basis and transition frequencies (w_ab)
#     a_ops_S_vib, trans_freq = calc_a_ops_freq(H, a_ops_vib)
    
#     # construct empty R
#     R = np.zeros((dim**2,dim**2),dtype = complex) 

#     Gamma = Gamma_gen(Gamma_vib, trans_freq)
    
#     for x, in np.ndindex(2): # loop over uncorrelated a_ops
#         A_a = a_ops_S_vib[x]

#         for y, in np.ndindex(2):
#             A_b = a_ops_S_vib[y]

#             G_ab = Gamma[x, y].reshape(dim, dim)


#             R+= - 1/2 * (np.kron(np.dot(A_a, np.multiply(A_b, np.transpose(G_ab))), np.transpose(np.eye(dim, dim)))
#                         - np.kron(np.multiply(np.transpose(G_ab), A_a), np.transpose(A_b))
#                         + np.kron(np.eye(dim, dim), np.transpose(np.dot(np.multiply(A_a, G_ab), A_b)))
#                         - np.kron(A_a, np.transpose(np.multiply(A_b, G_ab)))) 
            
#     return R, ekets




















