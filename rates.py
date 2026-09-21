import numpy as np
from concurrent.futures import ProcessPoolExecutor, as_completed

def zero_func(w: float) -> float:
    return 0.0

def _compute_phi_chunk(start_idx: int, tau_chunk: np.ndarray, omega: np.ndarray, coeff_real: np.ndarray, coeff_imag: np.ndarray):
    """
    Module-level function for pickling by ProcessPoolExecutor.
    """
    wt = np.outer(omega, tau_chunk)
    integrand = coeff_real[:, None] * np.cos(wt) - 1j * coeff_imag[:, None] * np.sin(wt)
    return start_idx, tau_chunk.size, np.trapezoid(integrand, omega, axis=0)


class CombinedRates:
    """
    Evaluate spectral densities, thermal occupation factors, and transition rates for the 
    Bloch-Redfield dissipators.

    Parameters
    ----------
    system_params : SystemParameters
        Object containing all evaluated physical parameters, constants, and bath temperatures.

    Attributes
    ----------
    J : float
        Dipole-dipole coupling strength.
    gamma_q1 : float
        Bare spontaneous emission rate.
    T_opt : float
        Temperature of the optical bath.
    T_vib : float
        Temperature of the vibrational bath.
    mode : str
        Coupling regime mapping ('weak' or 'polaron').
    """

    def __init__(self, system_params):
        self.system_params = system_params; self.constants = system_params.constants; self.mode = system_params.config_params.mode; self.T_opt = system_params.T_opt; self.T_vib = system_params.T_vib; self.J = system_params.J; self.gamma_q1 = system_params.gamma_q1
        self.constants = system_params.constants
        self.mode = system_params.config_params.mode
        self.integration_method = system_params.config_params.integration_method
        self.T_opt = system_params.T_opt
        self.T_vib = system_params.T_vib
        self.J = system_params.J
        self.gamma_q1 = system_params.gamma_q1

        if self.integration_method == 'vectorised':
            self.omega_integral = np.linspace(1.0e-4, 1500, 10000)
            self.Jw = self.J_vib(self.omega_integral)
            if self.T_vib == 0.0:
                self.tfac = np.ones_like(self.omega_integral)
            else:
                self.tfac = 1.0 / np.tanh(self._beta(self.T_vib) * self.omega_integral / 2.0)
                
                
    def J_opt(self, w: float) -> float:
        """
        Evaluate the optical spectral density.
        """
        return self.gamma_q1


    def J_vib(self, w: float) -> float:
        """
        Evaluate the super-Ohmic vibrational spectral density.
        """
        omega_c = 90.0e-3 * self.constants.eV / self.constants.hbar_ps
        E_reorg = 5.0e-3 * self.constants.eV / self.constants.hbar_ps
        k_vib_2 = E_reorg / (2.0 * omega_c**3)
        return (k_vib_2 * np.abs(w)**3 * np.exp(-(np.abs(w) / omega_c)))

        # # --- GaAs Quantum Dot Spectral Density (Eq. 17 / S51) ---
        # mu = 5370.0  # Mass density in kg/m^3
        # c_sound = 5110.0  # Speed of sound MUST remain in m/s for unit cancellation
        #
        # # Deformation potentials MUST remain as energies (Joules), not frequencies
        # D_e = 7.0 * self.constants.eV
        # D_h = -3.5 * self.constants.eV
        #
        # # Cut-off frequencies correctly calculated in ps^-1
        # omega_e = 2.9e-3 * self.constants.eV / self.constants.hbar_ps 
        # omega_h = 4.4e-3 * self.constants.eV / self.constants.hbar_ps
        #
        # # Convert input w from ps^-1 to s^-1 for the prefactor calculation
        # w_SI = np.abs(w) * 1.0e12
        #
        # # Use self.constants.hbar (J s) and w_SI (s^-1) to ensure the physics units cancel
        # prefactor_SI = (w_SI**3) / (4.0 * np.pi**2 * mu * self.constants.hbar * (c_sound**5))
        #
        # # Bracket term computed in Joules. The exponents use ps^-1 / ps^-1, which correctly cancels out.
        # bracket_term = D_e * np.exp(-(np.abs(w)**2) / (omega_e**2)) - D_h * np.exp(-(np.abs(w)**2) / (omega_h**2))
        #
        # # Compute spectral density in SI units (s^-1) and converted back to ps^-1
        # J_w = prefactor_SI * (bracket_term**2) * 1.0e-12
        #
        # return J_w
        
    def _beta(self, T: float) -> float:
        """
        Calculate the inverse thermal energy scale.
        """
        if T == 0.0:
            return np.inf
        elif T == np.inf:
            return 0.0
        else:
            return 1.0 / (self.constants.k_B_ps * T)
        
    def _n_inv(self, w: float, beta: float) -> float:
        """
        Evaluate the inverse Bose-Einstein thermal occupation factor.
        """
        return (np.exp(np.abs(w) * beta) - 1.0)
    

    def _thermal_occupation(self, w: float, T: float) -> float:
        """
        Evaluate thermal occupation factor, explicitly handling the T=0 limit.
        """
        if T == 0.0:
            return 0.0
        
        beta = self._beta(T)
        _n_inv = self._n_inv(w, beta)
        if np.abs(_n_inv) > 1.0e-8:
            return 1.0 / _n_inv
        return 0.0

    def Opt_rate(self) -> np.ndarray:
        def G_opt_11(w):
            n_pt = self._thermal_occupation(w, self.T_opt)
            if w > 0.0:
                return (1.0 + n_pt) * self.J_opt(w)
            elif w == 0.0:
                return (1.0 + 2.0 * n_pt) * self.J_opt(w)
            else:
                return n_pt * self.J_opt(w)

        def G_opt_12(w):
            n_pt = self._thermal_occupation(w, self.T_opt)
            if w > 0.0:
                return np.sqrt((1.0 + n_pt) * self.J_opt(w)) * np.sqrt((1.0 + n_pt) * self.J_opt(w))
            elif w == 0.0:
                return np.sqrt((1.0 + 2.0 * n_pt) * self.J_opt(w)) * np.sqrt((1.0 + 2.0 * n_pt) * self.J_opt(w))
            else:
                return np.sqrt(n_pt * self.J_opt(w)) * np.sqrt(n_pt * self.J_opt(w))

        return np.array([[G_opt_11, G_opt_12], 
                         [G_opt_12, G_opt_11]])

    def Vib_rate(self) -> np.ndarray:
        if self.mode == 'weak':
            def G_vib_11(w):
                n_pn = self._thermal_occupation(w, self.T_vib)
                if w > 0.0:
                    return (1.0 + n_pn) * self.J_vib(w)
                elif w == 0.0:
                    return (1.0 + 2.0 * n_pn) * self.J_vib(w)
                else:
                    return n_pn * self.J_vib(w)

            return np.array([[G_vib_11, zero_func], 
                             [zero_func, G_vib_11]])
        return None
    
    
    def Phi(self, tau_values: np.ndarray, batch_size: int = 100) -> np.ndarray:
        """
        Evaluate the polaron displacement correlation function.
        """
        if self.integration_method == 'sequential':
            beta_vib = self._beta(self.T_vib)
            w = np.linspace(1.0e-4, 1500, 10000)
            
            tau_values = np.atleast_1d(tau_values)
            result = np.zeros_like(tau_values, dtype=complex)
            
            for idx, t in enumerate(tau_values):
                if self.T_vib == 0.0:
                    coth_factor = 1.0
                else:
                    coth_factor = 1.0 / np.tanh(beta_vib * w / 2.0)
                
                integrand = self.J_vib(w) * (np.cos(w * t) * coth_factor - 1j * np.sin(w * t)) / (w**2)
                result[idx] = np.trapezoid(integrand, w)
            
            return result[0] if result.size == 1 else result

        elif self.integration_method == 'vectorised':
            tau_values = np.atleast_1d(tau_values)
            result = np.zeros_like(tau_values, dtype=complex)
            
            coeff_real = self.Jw * self.tfac / (self.omega_integral**2)
            coeff_imag = self.Jw / (self.omega_integral**2)

            batch_starts = range(0, len(tau_values), batch_size)
            
            with ProcessPoolExecutor() as executor:
                futures = []
                for start in batch_starts:
                    tau_chunk = tau_values[start:start + batch_size]
                    future = executor.submit(
                        _compute_phi_chunk, start, tau_chunk, self.omega_integral, coeff_real, coeff_imag
                    )
                    futures.append(future)

                for future in as_completed(futures):
                    start, size, batch_result = future.result()
                    result[start:start + size] = batch_result

            return result[0] if result.size == 1 else result

    def kappa(self) -> float:
        """
        Evaluate the equilibrium polaron renormalization factor.
        """
        if self.mode == 'polaron':
            if self.integration_method == 'vectorised':
                integrand = self.Jw * self.tfac / (self.omega_integral**2)
                integral = np.trapezoid(integrand, self.omega_integral)
                return float(np.exp(-0.5 * integral))
            else:
                return float(np.exp(-0.5 * np.real(self.Phi(0.0))))
        else:
            return 1.0

    def Coup_rate(self) -> np.ndarray:
        """
        Generate the non-Markovian phonon-mediated coupling rate matrix.
        """
        if self.mode == 'polaron':
            if self.integration_method == 'vectorised':
                kappa_val = self.kappa()
                pref = ((self.J / 2.0)**2) * (kappa_val**4)
                
                # Precompute the displacement array outside the frequency loop
                tau_values = np.linspace(0, 0.5, 20000)
                phi_vals = self.Phi(tau_values)
                exp_term_minus_one = pref * (np.exp(2.0 * phi_vals) - 1.0)

                def G_c_vec(w):
                    if np.isnan(w):
                        return 0.0
                    elif w == 0.0:
                        integrand = exp_term_minus_one
                    else:
                        integrand = np.exp(1j * w * tau_values) * exp_term_minus_one

                    return np.trapezoid(integrand, tau_values).real
                
                return np.array([[zero_func, G_c_vec], [G_c_vec, zero_func]])
            
            else:
                cache_w0 = None
                def G_c_seq(w):
                    nonlocal cache_w0 
                    if np.isnan(w):
                        return 0.0
                    elif w == 0.0:
                        if cache_w0 is None:
                            tau_values = np.linspace(0, 0.5, 20000)
                            pref = ((self.J / 2.0)**2) * (self.kappa()**4)
                            integrand = np.array([pref * (np.exp(2.0 * self.Phi(t)) - 1.0) for t in tau_values])
                            cache_w0 = np.trapezoid(integrand, tau_values).real
                        return cache_w0  
                    else:
                        tau_values = np.linspace(0, 0.5, 20000)
                        pref = ((self.J / 2.0)**2) * (self.kappa()**4)
                        integrand = np.array([np.exp(1j * w * t) * pref * (np.exp(2.0 * self.Phi(t)) - 1.0) for t in tau_values])
                        return np.trapezoid(integrand, tau_values).real

                return np.array([[zero_func, G_c_seq], [G_c_seq, zero_func]])
        return None
    
    def P_weight(self) -> np.ndarray:
        """
        Generate the polaron tensor weighting matrix.
        """
        if self.mode == 'polaron':
            return np.array([[self.kappa()**4, 1.0, self.kappa()**2, self.kappa()**2],
                             [1.0, self.kappa()**4, self.kappa()**2, self.kappa()**2],
                             [self.kappa()**2, self.kappa()**2, self.kappa()**4, 1.0],
                             [self.kappa()**2, self.kappa()**2, 1.0, self.kappa()**4]])
        return None

