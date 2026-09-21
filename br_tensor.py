from sys_params import *
from system_ops import *
from spectral_dens import *
from rates import *


import numpy as np
from concurrent.futures import ThreadPoolExecutor

class BlochRedfieldCalculator:
    """
    Evaluate the Bloch-Redfield master equation tensor in the system eigenbasis.
    """

    def __init__(self, system_params, rates, H: np.ndarray, a_ops_opt: list, Gamma_opt: np.ndarray, a_ops_vib: list = None, Gamma_vib: np.ndarray = None, a_ops_coup: list = None, Gamma_coup: np.ndarray = None, P_vib: np.ndarray = None):
        self.system_params = system_params
        self.rates = rates
        self.H = H
        self.dim = len(self.H)
        self.evals, self.ekets = np.linalg.eigh(self.H)

        self.a_ops_opt = a_ops_opt
        self.Gamma_opt = Gamma_opt
        self.a_ops_vib = a_ops_vib
        self.Gamma_vib = Gamma_vib
        self.a_ops_coup = a_ops_coup
        self.Gamma_coup = Gamma_coup
        self.P_vib = P_vib

        self.mode = self.system_params.config_params.mode
        self.pump = self.system_params.config_params.pump
        self.kappa = self.rates.kappa()

        sp = np.array([[0, 1], [0, 0]], dtype=complex)
        sp1 = np.kron(sp, np.identity(2))
        sp2 = np.kron(np.identity(2), sp)

        if self.pump in ['site', 'sym']:
            self.c_ops_S = [np.sqrt(self.system_params.gamma_p_1) * sp1, np.sqrt(self.system_params.gamma_p_2) * sp2]
        elif self.pump in ['coh', 'none']:
            self.c_ops_S = []
        else:
            raise ValueError(f"Unknown pump configuration: {self.pump}")

        self.P_lv = np.full((2, 2), self.kappa**4)

    def _calc_a_ops_freq(self, a_ops: list) -> tuple:
        a_ops_S = [self.ekets.conj().T @ a_op @ self.ekets for a_op in a_ops]
        trans_freq = np.array([self.evals[i] - self.evals[j] for i in range(self.dim) for j in range(self.dim)])
        return a_ops_S, trans_freq

    def _Gamma_gen(self, Gamma: np.ndarray, trans_freq: np.ndarray, shift_w: float = 0.0) -> np.ndarray:
        G = np.zeros((2, 2, len(trans_freq)), dtype=complex)
        for (i, j), G_ab in np.ndenumerate(Gamma):
            G[i, j, :] = np.array([G_ab(w + shift_w) for w in trans_freq])
        return G

    def _Gamma_pn(self, Gamma_coup: np.ndarray, trans_freq: np.ndarray) -> np.ndarray:
        G = np.zeros((2, 2, len(trans_freq)), dtype=complex)
        
        def compute_G_ab(G_ab):
            with ThreadPoolExecutor(max_workers=10) as executor:
                return np.array(list(executor.map(G_ab, trans_freq)))
                
        for (i, j), G_ab in np.ndenumerate(Gamma_coup):
            G[i, j, :] = compute_G_ab(G_ab)
        return G

    def Unitary_Dynamics(self) -> np.ndarray:
        trans_freq = np.array([self.evals[i] - self.evals[j] for i in range(self.dim) for j in range(self.dim)])
        return -1j * np.diag(trans_freq)

    def Liouvillian(self) -> np.ndarray:
        superL = np.zeros((self.dim**2, self.dim**2), dtype=complex)
        if not self.c_ops_S:
            return superL

        c_ops = [self.ekets.conj().T @ c_op @ self.ekets for c_op in self.c_ops_S]
        I = np.eye(self.dim)

        if self.pump == 'sym':
            for x, c_a in enumerate(c_ops):
                for y, c_b in enumerate(c_ops):
                    superL += (np.kron(c_a.conj(), c_b) - 0.5 * (np.kron(I, c_a.conj().T @ c_b) + np.kron(c_a.T @ c_b.conj(), I))) * self.P_lv[x, y]
        else:
            for c_op in c_ops:
                superL += (np.kron(c_op.conj(), c_op) - 0.5 * (np.kron(I, c_op.conj().T @ c_op) + np.kron(c_op.T @ c_op.conj(), I)))
                
        return superL

    def mathcal_F(self, k_q_r: np.ndarray, threshold: float = 1e-16) -> np.ndarray:
        alpha = lambda d1_hat, d2_hat: (np.dot(d1_hat, d2_hat) - (np.dot(d1_hat, self.system_params.r_hat)) * (np.dot(d2_hat, self.system_params.r_hat)))
        beta = lambda d1_hat, d2_hat: (np.dot(d1_hat, d2_hat) - 3.0 * (np.dot(d1_hat, self.system_params.r_hat)) * (np.dot(d2_hat, self.system_params.r_hat)))

        cross_F = (alpha(self.system_params.d1_hat, self.system_params.d2_hat) * (np.sin(k_q_r) / k_q_r) + 
                beta(self.system_params.d1_hat, self.system_params.d2_hat) * (np.cos(k_q_r) / (k_q_r**2) - np.sin(k_q_r) / (k_q_r**3))) * 1.5

        cross_F[np.abs(k_q_r) < threshold] = np.dot(self.system_params.d1_hat, self.system_params.d2_hat)
        
        self_F_1 = np.ones_like(k_q_r) * np.dot(self.system_params.d1_hat, self.system_params.d1_hat)
        self_F_2 = np.ones_like(k_q_r) * np.dot(self.system_params.d2_hat, self.system_params.d2_hat)

        F = np.array([
            [self_F_1, self_F_1, cross_F,  cross_F],
            [self_F_1, self_F_1, cross_F,  cross_F],
            [cross_F,  cross_F,  self_F_2, self_F_2],
            [cross_F,  cross_F,  self_F_2, self_F_2]
        ])

        return F

    def BR_tensor_opt(self) -> np.ndarray:
        a_ops_S_opt, trans_freq = self._calc_a_ops_freq(self.a_ops_opt)
        
        Gamma_opt_0 = self._Gamma_gen(self.Gamma_opt, trans_freq, shift_w=self.system_params.omega_L)
        
        R_opt = np.zeros((self.dim**2, self.dim**2), dtype=complex)

        k_q = np.reshape(np.ones_like(trans_freq) * self.system_params.k_q1, (self.dim, self.dim))
        k_q_r = k_q * self.system_params.r 
        F = self.mathcal_F(k_q_r) 

        for x in range(len(a_ops_S_opt)):
            A_a = a_ops_S_opt[x]
            for y in range(len(a_ops_S_opt)):
                A_b = a_ops_S_opt[y]

                Gamma_Opt = np.tile(Gamma_opt_0, (2, 2, 1))
                G_ab = Gamma_Opt[x, y].reshape(self.dim, self.dim)
                G_ab = G_ab * F[x, y]
                
                P_ab = self.P_vib[x, y] if self.mode == 'polaron' else 1.0
                
                R_opt += -0.5 * (np.kron(np.dot(A_a, np.multiply(A_b, np.transpose(G_ab))), np.transpose(np.eye(self.dim, self.dim)))
                        - np.kron(np.multiply(np.transpose(G_ab), A_a), np.transpose(A_b))
                        + np.kron(np.eye(self.dim, self.dim), np.transpose(np.dot(np.multiply(A_a, G_ab), A_b)))
                        - np.kron(A_a, np.transpose(np.multiply(A_b, G_ab)))) * P_ab 
                
        return R_opt

    def BR_tensor_vibcoup(self) -> np.ndarray:
        if self.mode == 'weak':
            a_ops_S, trans_freq = self._calc_a_ops_freq(self.a_ops_vib)
            Gamma = self._Gamma_gen(self.Gamma_vib, trans_freq, shift_w=0.0)
        elif self.mode == 'polaron':
            a_ops_S, trans_freq = self._calc_a_ops_freq(self.a_ops_coup)
            trans_freq[(np.abs(trans_freq) > 1000)] = np.nan
            Gamma = self._Gamma_pn(self.Gamma_coup, np.array(trans_freq))
        else:
            return np.zeros((self.dim**2, self.dim**2), dtype=complex)

        R = np.zeros((self.dim**2, self.dim**2), dtype=complex)
        
        for x in range(len(a_ops_S)):
            A_a = a_ops_S[x]
            for y in range(len(a_ops_S)):
                A_b = a_ops_S[y]
                
                G_ab = Gamma[x, y].reshape(self.dim, self.dim)

                R += -0.5 * (np.kron(np.dot(A_a, np.multiply(A_b, np.transpose(G_ab))), np.transpose(np.eye(self.dim, self.dim)))
                        - np.kron(np.multiply(np.transpose(G_ab), A_a), np.transpose(A_b))
                        + np.kron(np.eye(self.dim, self.dim), np.transpose(np.dot(np.multiply(A_a, G_ab), A_b)))
                        - np.kron(A_a, np.transpose(np.multiply(A_b, G_ab))))

        return R

    def total_BR_tensor(self) -> np.ndarray:
        tensor = self.Unitary_Dynamics() + self.BR_tensor_opt()
        
        if self.mode == 'polaron' or (self.mode == 'weak' and self.a_ops_vib is not None):
            tensor += self.BR_tensor_vibcoup()
            
        if self.c_ops_S:
            tensor += self.Liouvillian()
            
        return tensor










