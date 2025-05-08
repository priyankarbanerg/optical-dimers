# from spectral_dens import *
# mp.dps = 20
from sys_params import *

constants = Constants()
dipole_params = DipoleParameters(config_name)  
system_params = SystemParameters(dipole_params=dipole_params, E_q1=E_q1*constants.eV, E_q2=E_q2*constants.eV)
system_params.calculate_parameters()
# J_opt, J_vib = calculate_spectral_densities()


############################################################
'''Calculate Rates Γ(ω) for Bloch-Redfield'''
############################################################

class CombinedRates:
    def __init__(self, T_opt, T_vib, mode):
        self.T_opt = T_opt
        self.T_vib = T_vib
        self.mode = mode
        self.J = system_params.J
        self.gamma_q1 = system_params.gamma_q1

    # Define spectral density functions
    def J_opt(self, w):
        k_opt_1 = self.gamma_q1
        return k_opt_1


    def J_vib(self, w):
        omega_c = 90.0e-3 * constants.eV / system_params.hbar_ps
        # E_reorg = 5.0e-3 * constants.eV / system_params.hbar_ps
        E_reorg = 5.0e-3 * constants.eV / system_params.hbar_ps

        k_vib_2 = E_reorg / (2 * omega_c**3)
        # if w==0:
        #     return k_vib_2
        # else:
        #     return (k_vib_2 * np.abs(w)**3) * np.exp(- (np.abs(w) / omega_c))
        # print(1.0e3 * self.gamma_q1)
        # print(k_vib_2 * np.abs(w)**3 * np.exp(- (np.abs(w) / omega_c)))
        # return 1.0e3 * self.gamma_q1 #
            # return (k_vib_2 * np.abs(w)**3) * np.exp(- (np.abs(w) / omega_c))
        return (k_vib_2 * np.abs(w)**3 * np.exp(- (np.abs(w) / omega_c)))

    # @jit
    def _beta(self, T):
        if T == 0.0:
            return np.inf
        elif T == np.inf:
            return 0.0
        else:
            return 1 / (constants.k_B_ps * T)

    def _n_inv(self, w, beta):
        return (np.exp((np.abs(w)) * beta) - 1)

    def Opt_rate(self):
        beta_opt = self._beta(self.T_opt)

        def G_opt_11(w):
            _n_pt = self._n_inv(w, beta_opt)
            if np.abs(_n_pt) > 1.0e-8:
                n_pt = 1 / _n_pt
                if w > 0.0:
                    return (1 + n_pt) * self.J_opt(w)
                elif w == 0.0:
                    return (1 + 2 * n_pt) * self.J_opt(w)
                else:
                    return (n_pt) * self.J_opt(w)
            else:
                return self.J_opt(w)

        def G_opt_12(w):
            _n_pt = self._n_inv(w, beta_opt)
            if np.abs(_n_pt) > 1.0e-8:
                n_pt = 1 / _n_pt
                if w > 0.0:
                    return np.sqrt((1 + n_pt) * self.J_opt(w)) * np.sqrt((1 + n_pt) * self.J_opt(w))
                elif w == 0.0:
                    return np.sqrt((1 + 2 * n_pt) * self.J_opt(w)) * np.sqrt((1 + 2 * n_pt) * self.J_opt(w))
                else:
                    return np.sqrt(n_pt * self.J_opt(w)) * np.sqrt(n_pt * self.J_opt(w)) 
            else:
                return np.sqrt(self.J_opt(w)) * np.sqrt(self.J_opt(w))

        def G_opt_21(w):
            _n_pt = self._n_inv(w, beta_opt)
            if np.abs(_n_pt) > 1.0e-8:
                n_pt = 1 / _n_pt
                if w > 0.0:
                    return np.sqrt((1 + n_pt) * self.J_opt(w)) * np.sqrt((1 + n_pt) * self.J_opt(w))
                elif w == 0.0:
                    return np.sqrt((1 + 2 * n_pt) * self.J_opt(w)) * np.sqrt((1 + 2 * n_pt) * self.J_opt(w))
                else:
                    return np.sqrt(n_pt * self.J_opt(w)) * np.sqrt(n_pt * self.J_opt(w)) 
            else:
                return np.sqrt(self.J_opt(w)) * np.sqrt(self.J_opt(w))

        def G_opt_22(w):
            _n_pt = self._n_inv(w, beta_opt)
            if np.abs(_n_pt) > 1.0e-8:
                n_pt = 1 / _n_pt
                if w > 0.0:
                    return (1 + n_pt) * self.J_opt(w)
                elif w == 0.0:
                    return (1 + 2 * n_pt) * self.J_opt(w)
                else:
                    return n_pt * self.J_opt(w)
            else:
                return self.J_opt(w)

        return np.array([[G_opt_11, G_opt_12], 
                         [G_opt_21, G_opt_22]])
    def Vib_rate(self):
        if self.mode == 'weak':
            beta_vib = self._beta(self.T_vib)

            def G_vib_11(w):
                _n_pn = self._n_inv(w, beta_vib)
                if np.abs(_n_pn) > 1.0e-8:
                    n_pn = 1 / _n_pn
                    if w > 0.0:
                        return (1 + n_pn) * self.J_vib(w)
                    elif w == 0.0:
                        return (1 + 2 * n_pn) * self.J_vib(w)
                    else:
                        return n_pn * self.J_vib(w)
                else:
                    return self.J_vib(w)

            return np.array([[G_vib_11, lambda w: 0], 
                             [lambda w: 0, G_vib_11]])
        else:
            return None


    # def Phi(self, t):
    #     beta_vib = self._beta(self.T_vib)
    #     integrand = lambda w: self.J_vib(w) * (np.cos(w * t)* (1 / np.tanh(beta_vib * w / 2)) - 1j * np.sin(w * t)) / w**2 # * (1 / np.tanh(beta_vib * w / 2))
    #     result = quad_vec(integrand, 1.0e-8, 500)
    #     return result[0]

    # def kappa(self):
    #     if self.mode == 'polaron':
    #         return np.exp(-0.5 * self.Phi(0))
    #     else:
    #         return 1.0
 

    # def Coup_rate(self):
    #     if self.mode == 'polaron':
    #         # Parameters
    #         t_max = 1  # Maximum time
    #         n_points = 10000  # Number of points in time domain

    #         # Time array
    #         t = np.linspace(0, t_max, n_points)
    #         dt = t[1] - t[0]  # Time step

    #         # Frequency array generated by FFT
    #         omega_fft = fftfreq(n_points, d=dt) * 2 * np.pi
    #         omega_fft_shifted = fftshift(omega_fft)

    #         J_scaled = self.J * self.kappa()**2
    #         pref = (self.kappa()**4) * ((J_scaled / 2) ** 2)
            
    #         def G_c(w):
    #             integrand = 2 * pref * (np.exp(2 * self.Phi(t)) - 1)
    #             # Fourier Transform using DFT (FFT)
    #             F_dft = ifft(integrand) 
    #             F_dft_shifted = ifftshift(F_dft)

    #             # Interpolate FFT result onto custom frequency array
    #             pchip_interpolator = PchipInterpolator(omega_fft_shifted, F_dft_shifted)
    #             F_dft_custom = pchip_interpolator(w)

    #             return np.real(F_dft_custom)
     

    #         return np.array([[lambda w: 0, G_c], [G_c, lambda w: 0]])
    #     else:
    #         return None




    def Phi(self, t):
        beta_vib = self._beta(self.T_vib)
        w = np.linspace(1.0e-4, 2000, 10000)
        integrand = self.J_vib(w) * (np.cos(w * t)* (1 / np.tanh(beta_vib * w / 2)) - 1j * np.sin(w * t)) / w**2 # * (1 / np.tanh(beta_vib * w / 2))
        return np.trapz(integrand, w)
    
    
    def kappa(self):
        if self.mode == 'polaron':
            return np.exp(-0.5 * self.Phi(0))
        else:
            return 1.0
 
    def Coup_rate(self):
        if self.mode == 'polaron':
            
            J_scaled = self.J * self.kappa()**2
            cache_w0 = None

            def G_c(w):
                nonlocal cache_w0  # Allow access to the cache

                if np.isnan(w):
                    return 0.0
                elif w == 0:
                    # Check if we already computed the value for w = 0
                    if cache_w0 is None:
                        tau_values = np.linspace(0, 20, 20000)
                        pref =  ((J_scaled / 2)**2) * (self.kappa()**4)
                        integrand = np.array([pref * (np.exp(2 * self.Phi(t)) - 1) for t in tau_values]) #np.array([(pref * np.exp(2 * self.Phi_re(t)) * (np.cos(w * t) * np.cos(2 * self.Phi_im(t)) - np.sin(w * t) * np.sin(2 * self.Phi_im(t))) - 1) for t in tau_values]) # np.array([np.exp(1j * w * t) * pref * (np.exp(2 * self.Phi(t)) - 1) for t in tau_values])
                        # print(w, np.trapz(integrand, tau_values).real)
                        cache_w0 = np.trapz(integrand, tau_values).real
                    return cache_w0  
                else:
                    tau_values = np.linspace(0, 20, 20000)
                    pref =  ((J_scaled / 2)**2) * (self.kappa()**4)
                    integrand = np.array([np.exp(1j * w * t) * pref * (np.exp(2 * self.Phi(t)) - 1) for t in tau_values]) #np.array([(pref * np.exp(2 * self.Phi_re(t)) * (np.cos(w * t) * np.cos(2 * self.Phi_im(t)) - np.sin(w * t) * np.sin(2 * self.Phi_im(t))) - 1) for t in tau_values]) # np.array([np.exp(1j * w * t) * pref * (np.exp(2 * self.Phi(t)) - 1) for t in tau_values])
                    result = np.trapz(integrand, tau_values).real
                    # if np.isclose(np.abs(w), 2 * np.abs(self.J), atol=1):
                    #     print(w, result)
                    # print(w, result)
                    return result 
            return np.array([[lambda w: 0, G_c], [G_c, lambda w: 0]])
        else:
            return None

        
    def P_weight(self):
        if self.mode == 'polaron':
            return np.array([[self.kappa()**4, 1, self.kappa()**2, self.kappa()**2],
                             [1, self.kappa()**4, self.kappa()**2, self.kappa()**2],
                             [self.kappa()**2, self.kappa()**2, self.kappa()**4, 1],
                             [self.kappa()**2, self.kappa()**2, 1, self.kappa()**4]])
        else:
            return None

############################################################
'''Built the rate matrix'''
############################################################

def Gamma_gen(Gamma, trans_freq):
    G= np.zeros((2, 2, len(trans_freq)))
    for (i, j), G_ab in np.ndenumerate(Gamma):
        G[i,j,:] = np.array([G_ab(w) for w in trans_freq])
    return G








































    # # Define the Phi(t) function
    # def Phi_re(self, t):
    #     beta_vib = self._beta(self.T_vib)
    #     integrand = lambda w: J_vib(w) * (np.cos(w * t) * (1 / np.tanh((beta_vib * w) / 2))) / w**2
    #     # result = quad(integrand, [1.0e-4, 500])
    #     result = quad_vec(integrand, 1.0e-4, 500)
    #     return result[0]
    
    # def Phi_im(self, t):
    #     # beta_vib = self._beta(self.T_vib)
    #     integrand = lambda w: J_vib(w) * (np.sin(w * t)) / w**2
    #     # result = quad(integrand, [1.0e-4, 500])
    #     result = quad_vec(integrand, 1.0e-4, 500)
    #     return result[0]


    # def kappa(self):
    #     if self.mode == 'polaron':
    #         return np.exp(-0.5 * (self.Phi_re(0) + 1j * self.Phi_im(0)))
    #     else:
    #         return 1.0
        


    # def Coup_rate_1(self):
    #     if self.mode == 'polaron':
    #         def G_c_1(w):
    #             tau_values = np.linspace(1.0e-4, 0.2, 100)
    #             integrand = lambda t: np.exp(1j * w * t) * (self.kappa()**4) * ((self.J / 2)**2) * (np.exp(2 * self.Phi_re(t)) * (np.cos(w * t) * np.cos(2 * self.Phi_im(t)) - np.sin(w * t) * np.sin(2 * self.Phi_im(t))) - 1)
    #             inte = np.array([integrand(t) for t in tau_values], dtype = complex)
    #             return simpson(inte).real

    #         def G_c_2(w):
    #             tau_values = np.linspace(1.0e-4, 0.2, 100)
    #             integrand = lambda t: np.exp(1j * w * t) * (self.kappa()**4) * ((self.J / 2)**2) *  (np.exp(-2 * self.Phi_re(t)) * (np.cos(w * t) * np.cos(2 * self.Phi_im(t)) + np.sin(w * t) * np.sin(2 * self.Phi_im(t))) - 1) 
    #             inte = np.array([integrand(t) for t in tau_values], dtype = complex)
    #             return simpson(inte).real

    #         return np.array([[G_c_1, G_c_2], [G_c_2, G_c_1]])
    #     else:
    #         return None

    # def Coup_rate_2(self):
    #     if self.mode == 'polaron':
    #         def G_c_1(w):
    #             tau_values = np.linspace(1.0e-4, 0.2, 100)
    #             integrand = lambda t: np.exp(-1j * w * t) * (self.kappa()**4) * ((self.J / 2)**2) * (np.exp(2 * self.Phi_re(-t)) * (np.cos(w * (-t)) * np.cos(2 * self.Phi_im(-t)) - np.sin(w * (-t)) * np.sin(2 * self.Phi_im(-t))) - 1)
    #             inte = np.array([integrand(t) for t in tau_values], dtype = complex)
    #             return simpson(inte).real

    #         def G_c_2(w):
    #             tau_values = np.linspace(1.0e-4, 0.2, 100)
    #             integrand = lambda t: np.exp(-1j * w * t) * (self.kappa()**4) * ((self.J / 2)**2) * (np.exp(-2 * self.Phi_re(-t)) * (np.cos(w * (-t)) * np.cos(2 * self.Phi_im(-t)) + np.sin(w * (-t)) * np.sin(2 * self.Phi_im(-t))) - 1)
    #             inte = np.array([integrand(t) for t in tau_values], dtype = complex)
    #             return simpson(inte).real

    #         return np.array([[G_c_1, G_c_2], [G_c_2, G_c_1]])
    #     else:
    #         return None




    # def Coup_rate_1(self):
    #     if self.mode == 'polaron':
    #         # Parameters
    #         t_max = 0.1  # Maximum time
    #         n_points = 10000  # Number of points in time domain

    #         # Time array
    #         t = np.linspace(0, t_max, n_points)
    #         dt = t[1] - t[0]  # Time step

    #         # Frequency array generated by FFT
    #         omega_fft = fftfreq(n_points, d=dt) * 2 * np.pi
    #         omega_fft_shifted = fftshift(omega_fft)

    #         pref = (self.kappa()** 4) * ((self.J / 2) ** 2) 
            
            
    #         def G_c_1(w):
    #             integrand = pref * (np.exp(2 * self.Phi_re(t)) - 1)
    #             # Fourier Transform using DFT (FFT)
    #             F_dft = fft(integrand) * dt
    #             F_dft_shifted = fftshift(F_dft)

    #             # Interpolate FFT result onto custom frequency array
    #             pchip_interpolator = PchipInterpolator(omega_fft_shifted, F_dft_shifted)
    #             F_dft_custom = pchip_interpolator(w)

    #             return np.real(F_dft_custom)
            
    #         def G_c_2(w):
    #             integrand = pref * (np.exp(-2 * self.Phi_re(t)) - 1) 
    #             # Fourier Transform using DFT (FFT)
    #             F_dft = fft(integrand) * dt
    #             F_dft_shifted = fftshift(F_dft)

    #             # Interpolate FFT result onto custom frequency array
    #             pchip_interpolator = PchipInterpolator(omega_fft_shifted, F_dft_shifted)
    #             F_dft_custom = pchip_interpolator(w)

    #             return np.real(F_dft_custom)



    #         return np.array([[G_c_1, G_c_2], [G_c_2, G_c_1]])
    #     else:
    #         return None
        
    # def Coup_rate_2(self):
    #     if self.mode == 'polaron':
    #         # Parameters
    #         t_max = 0.1  # Maximum time
    #         n_points = 10000  # Number of points in time domain

    #         # Time array
    #         t = np.linspace(0, t_max, n_points)
    #         dt = t[1] - t[0]  # Time step

    #         # Frequency array generated by FFT
    #         omega_fft = fftfreq(n_points, d=dt) * 2 * np.pi
    #         omega_fft_shifted = fftshift(omega_fft)

    #         pref = (self.kappa()** 4) * ((self.J / 2) ** 2) 
            
            
    #         def G_c_1(w):
    #             integrand = pref * (np.exp(2 * self.Phi_re(-t)) - 1)
    #             # Fourier Transform using DFT (FFT)
    #             F_dft = fft(integrand) * dt
    #             F_dft_shifted = fftshift(F_dft)

    #             # Interpolate FFT result onto custom frequency array
    #             pchip_interpolator = PchipInterpolator(omega_fft_shifted, F_dft_shifted)
    #             F_dft_custom = pchip_interpolator(w)

    #             return np.real(F_dft_custom)
            
    #         def G_c_2(w):
    #             integrand = pref * (np.exp(-2 * self.Phi_re(-t)) - 1)
    #             # Fourier Transform using DFT (FFT)
    #             F_dft = fft(integrand) * dt
    #             F_dft_shifted = fftshift(F_dft)

    #             # Interpolate FFT result onto custom frequency array
    #             pchip_interpolator = PchipInterpolator(omega_fft_shifted, F_dft_shifted)
    #             F_dft_custom = pchip_interpolator(w)

    #             return np.real(F_dft_custom)

    #         return np.array([[G_c_1, G_c_2], [G_c_2, G_c_1]])
    #     else:
    #         return None

    # def Coup_rate_1(self):
    #     if self.mode == 'polaron':
    #         # Parameters
    #         t_max = 1.0  # Maximum time
    #         n_points = 10000  # Number of points in time domain

    #         # Time array
    #         t = np.linspace(0, t_max, n_points)
    #         dt = t[1] - t[0]  # Time step

    #         # Frequency array generated by FFT
    #         omega_fft = fftfreq(n_points, d=dt) * 2 * np.pi
    #         omega_fft_shifted = fftshift(omega_fft)

    #         pref = (self.kappa()** 4) * ((self.J / 2) ** 2) 
            
            
    #         def G_c_1(w):
    #             integrand = pref * (np.exp(-2 * self.Phi(t)) - 1)
    #             # Fourier Transform using DFT (FFT)
    #             F_dft = fft(integrand) * dt
    #             F_dft_shifted = fftshift(F_dft)

    #             # Interpolate FFT result onto custom frequency array
    #             pchip_interpolator = PchipInterpolator(omega_fft_shifted, F_dft_shifted)
    #             F_dft_custom = pchip_interpolator(w)

    #             return np.real(F_dft_custom)
            
    #         def G_c_2(w):
    #             integrand = pref * (np.exp(2 * self.Phi(t)) - 1)
    #             # Fourier Transform using DFT (FFT)
    #             F_dft = fft(integrand) * dt
    #             F_dft_shifted = fftshift(F_dft)

    #             # Interpolate FFT result onto custom frequency array
    #             pchip_interpolator = PchipInterpolator(omega_fft_shifted, F_dft_shifted)
    #             F_dft_custom = pchip_interpolator(w)

    #             return np.real(F_dft_custom)



    #         return np.array([[G_c_1, G_c_2], [G_c_2, G_c_1]])
    #         # return np.array([[lambda w: 0, G_c_2], [G_c_2, lambda w: 0]])
    #     else:
    #         return None
        
    # def Coup_rate_2(self):
    #     if self.mode == 'polaron':
    #         # Parameters
    #         t_max = 1.0  # Maximum time
    #         n_points = 10000  # Number of points in time domain

    #         # Time array
    #         t = np.linspace(0, t_max, n_points)
    #         dt = t[1] - t[0]  # Time step

    #         # Frequency array generated by FFT
    #         omega_fft = fftfreq(n_points, d=dt) * 2 * np.pi
    #         omega_fft_shifted = fftshift(omega_fft)

    #         pref = (self.kappa()** 4) * ((self.J / 2) ** 2) 
            
            
    #         def G_c_1(w):
    #             integrand = pref * (np.exp(-2 * self.Phi(-t)) - 1)
    #             # Fourier Transform using DFT (FFT)
    #             F_dft = fft(integrand) * dt
    #             F_dft_shifted = fftshift(F_dft)

    #             # Interpolate FFT result onto custom frequency array
    #             pchip_interpolator = PchipInterpolator(omega_fft_shifted, F_dft_shifted)
    #             F_dft_custom = pchip_interpolator(w)

    #             return np.real(F_dft_custom)
            
    #         def G_c_2(w):
    #             integrand = pref * (np.exp(2 * self.Phi(-t)) - 1)
    #             # Fourier Transform using DFT (FFT)
    #             F_dft = fft(integrand) * dt
    #             F_dft_shifted = fftshift(F_dft)

    #             # Interpolate FFT result onto custom frequency array
    #             pchip_interpolator = PchipInterpolator(omega_fft_shifted, F_dft_shifted)
    #             F_dft_custom = pchip_interpolator(w)

    #             return np.real(F_dft_custom)


    #         # return np.array([[lambda w: 0, G_c_2], [G_c_2, lambda w: 0]])
    #         return np.array([[G_c_1, G_c_2], [G_c_2, G_c_1]])
    #     else:
    #         return None














    # def Coup_rate(self):
    #     if self.mode == 'polaron':
    #         def G_coup_11(w):
    #             tau_values = np.linspace(1.0e-4, 0.5, 3000)
    #             # Phi_values = self.Phi(tau_values)
    #             integrand = lambda t: np.exp(1j * w * t) * (self.kappa()**4) * ((self.J)**2) * (np.exp(2 * self.Phi(t)) - 1)
    #             inte = [integrand(t) for t in tau_values]
    #             return simpson(inte, tau_values).real

    #         def G_coup_22(w):
    #             tau_values = np.linspace(1.0e-4, 0.5, 3000)
    #             # Phi_values = self.Phi(tau_values)
    #             integrand = lambda t: np.exp(1j * w * t) * (self.kappa()**4) * ((self.J)**2) * (np.exp(2 * self.Phi(t)) - 1)
    #             inte = [integrand(t) for t in tau_values]
    #             return simpson(inte, tau_values).real

    #         def G_coup_12(w):
    #             tau_values = np.linspace(1.0e-4, 0.5, 3000)
    #             # Phi_values = self.Phi(tau_values)
    #             integrand = lambda t: np.exp(1j * w * t) * (self.kappa()**4) * ((self.J)**2) * (np.exp(-2 * self.Phi(t)) - 1)
    #             inte = [integrand(t) for t in tau_values]
    #             return simpson(inte, tau_values).real

    #         def G_coup_21(w):
    #             tau_values = np.linspace(1.0e-4, 0.5, 3000)
    #             # Phi_values = self.Phi(tau_values)
    #             integrand = lambda t: np.exp(1j * w * t) * (self.kappa()**4) * ((self.J)**2) * (np.exp(-2 * self.Phi(t)) - 1)
    #             inte = [integrand(t) for t in tau_values]
    #             return simpson(inte, tau_values).real

    #         return np.array([[G_coup_11, G_coup_12], 
    #                          [G_coup_21, G_coup_22]])
    #     else:
    #         return None













    # def Coup_rate(self):
    #     if self.mode == 'polaron':
    #         def G_coup_11(w):
    #             tau_values = np.linspace(1.0e-4, 0.5, 1000)
    #             Phi_values = self.Phi(tau_values)
    #             integrand = np.exp(1j * w * tau_values) * (self.kappa()**4) * ((self.J)**2) * (np.exp(2 * Phi_values) - 1)
    #             return simpson(integrand, tau_values).real

    #         def G_coup_22(w):
    #             tau_values = np.linspace(1.0e-4, 0.5, 1000)
    #             Phi_values = self.Phi(tau_values)
    #             integrand = np.exp(1j * w * tau_values) * (self.kappa()**4) * ((self.J)**2) * (np.exp(2 * Phi_values) - 1)
    #             return simpson(integrand, tau_values).real

    #         def G_coup_12(w):
    #             tau_values = np.linspace(1.0e-4, 0.5, 1000)
    #             Phi_values = self.Phi(tau_values)
    #             integrand = np.exp(1j * w * tau_values) * (self.kappa()**4) * ((self.J)**2) * (np.exp(-2 * Phi_values) - 1)
    #             return simpson(integrand, tau_values).real

    #         def G_coup_21(w):
    #             tau_values = np.linspace(1.0e-4, 0.5, 1000)
    #             Phi_values = self.Phi(tau_values)
    #             integrand = np.exp(1j * w * tau_values) * (self.kappa()**4) * ((self.J)**2) * (np.exp(-2 * Phi_values) - 1)
    #             return simpson(integrand, tau_values).real

    #         return np.array([[G_coup_11, G_coup_12], 
    #                          [G_coup_21, G_coup_22]])
    #     else:
    #         return None



    # def Coup_rate(self):
    #     if self.mode == 'polaron':
    #         def G_coup_11(w):
    #             tau_values = np.linspace(1.0e-4, 1, 1000)
    #             Phi_values = self.Phi(tau_values)
    #             integrand = np.exp(1j * w * tau_values) * (self.kappa()**4) * ((self.J)**2) * (np.exp(2 * Phi_values) - 1)
    #             return simpson(integrand, tau_values).real

    #         def G_coup_22(w):
    #             tau_values = np.linspace(1.0e-4, 1, 1000)
    #             Phi_values = self.Phi(tau_values)
    #             integrand = np.exp(1j * w * tau_values) * (self.kappa()**4) * ((self.J)**2) * (np.exp(2 * Phi_values) - 1)
    #             return simpson(integrand, tau_values).real

    #         def G_coup_12(w):
    #             tau_values = np.linspace(1.0e-4, 1, 1000)
    #             Phi_values = self.Phi(tau_values)
    #             integrand = np.exp(1j * w * tau_values) * (self.kappa()**4) * ((self.J)**2) * (np.exp(-2 * Phi_values) - 1)
    #             return simpson(integrand, tau_values).real

    #         def G_coup_21(w):
    #             tau_values = np.linspace(1.0e-4, 1, 1000)
    #             Phi_values = self.Phi(tau_values)
    #             integrand = np.exp(1j * w * tau_values) * (self.kappa()**4) * ((self.J)**2) * (np.exp(-2 * Phi_values) - 1)
    #             return simpson(integrand, tau_values).real

    #         return np.array([[G_coup_11, G_coup_12], 
    #                          [G_coup_21, G_coup_22]])
    #     else:
    #         return None







    # def Coup_rate(self):
    #     if self.mode == 'polaron':
    #         def G_coup_11(w):
    #             tau_values = np.linspace(1.0e-4, 1, 1000)
    #             Phi_values = self.Phi(tau_values)
    #             integrand = lambda t: np.exp(1j * w * t) * (self.kappa()**4) * ((self.J)**2) * (np.exp(2 * self.Phi(t)) - 1)
    #             return quad_vec(integrand, 1.0e-4, 1)[0].real

    #         def G_coup_22(w):
    #             tau_values = np.linspace(1.0e-4, 1, 1000)
    #             Phi_values = self.Phi(tau_values)
    #             integrand = lambda t: np.exp(1j * w * t) * (self.kappa()**4) * ((self.J)**2) * (np.exp(2 * self.Phi(t)) - 1)
    #             return quad_vec(integrand, 1.0e-4, 1)[0].real

    #         def G_coup_12(w):
    #             tau_values = np.linspace(1.0e-4, 1, 1000)
    #             Phi_values = self.Phi(tau_values)
    #             integrand = lambda t: np.exp(1j * w * t) * (self.kappa()**4) * ((self.J)**2) * (np.exp(-2 * self.Phi(t)) - 1)
    #             return quad_vec(integrand, 1.0e-4, 1)[0].real

    #         def G_coup_21(w):
    #             tau_values = np.linspace(1.0e-4, 1, 1000)
    #             Phi_values = self.Phi(tau_values)
    #             integrand = lambda t: np.exp(1j * w * t) * (self.kappa()**4) * ((self.J)**2) * (np.exp(-2 * self.Phi(t)) - 1)
    #             return quad_vec(integrand, 1.0e-4, 1)[0].real

    #         return np.array([[G_coup_11, G_coup_12], 
    #                          [G_coup_21, G_coup_22]])
    #     else:
    #         return None






















# ############################################################
# '''Calculate Rates Γ(ω) for Polaron Bloch-Redfield'''
# ############################################################

# # Additional coupling dissipator
# def Coup_rate(kappa):

#     def G_coup_11(w):
#         tau_values = np.linspace(1.0e-4, 1000, 10000)
#         Phi_values = Phi(tau_values)
#         integrand = [np.exp(1j * w * tau_values) * (kappa**4) * ((J/2)**2) * (np.exp(2 * Phi_values) - 1)]
#         return simpson(integrand, tau_values).real[0]

#     def G_coup_22(w):
#         tau_values = np.linspace(1.0e-4, 100, 10000)
#         Phi_values = Phi(tau_values)
#         integrand = [np.exp(1j * w * tau_values) * (kappa**4) * ((J/2)**2) * (np.exp(2 * Phi_values) - 1)]
#         return simpson(integrand, tau_values).real[0]
    
#     def G_coup_12(w):
#         tau_values = np.linspace(1.0e-4, 100, 10000)
#         Phi_values = Phi(tau_values)
#         integrand = [np.exp(1j * w * tau_values) * (kappa**4) * ((J/2)**2) * (np.exp(-2 * Phi_values) - 1)]
#         return simpson(integrand, tau_values).real[0]

#     def G_coup_21(w):
#         tau_values = np.linspace(1.0e-4, 100, 10000)
#         Phi_values = Phi(tau_values)
#         integrand = [np.exp(1j * w * tau_values) * (kappa**4) * ((J/2)**2) * (np.exp(-2 * Phi_values) - 1)]
#         return simpson(integrand, tau_values).real[0]

#     return np.array([[G_coup_11, G_coup_12], 
#                      [G_coup_21, G_coup_22]])



# # Optical dissipator
# def Opt_rate(T_opt):
    
#     def G_opt_11(w):

#         _n_pt = n_inv(w, beta_opt)
        
#         if np.abs(_n_pt) > 1.0e-8:
#             n_pt = (_n_pt**(-1))

#             if w >= 0.0:
#                 return (1 + n_pt) * J_opt(w)  # Spontaneous decay
#             else: 
#                 return n_pt * J_opt(w)        # Absorption
#         else:
#             if w >= 0.0:
#                 return J_opt(w)               # Spontaneous decay
#             else: 
#                 return J_opt(w)               # Absorption
    
#     def G_opt_12(w):

#         _n_pt = n_inv(w, beta_opt)

#         if np.abs(_n_pt) > 1.0e-8:
#             n_pt = (_n_pt**(-1))

#             if w >= 0.0:
#                 return np.sqrt((1 + n_pt) *  J_opt(w)) * np.sqrt((1 + n_pt) * J_opt(w))  # Spontaneous decay
#             else: 
#                 return np.sqrt(n_pt * J_opt(w)) * np.sqrt(n_pt * J_opt(w))        # Absorption
#         else:
#             if w >= 0.0:
#                 return np.sqrt(J_opt(w)) * np.sqrt(J_opt(w))        # Spontaneous decay
#             else: 
#                 return np.sqrt(J_opt(w)) * np.sqrt(J_opt(w))       # Absorption
    
#     def G_opt_21(w):

#         _n_pt = n_inv(w, beta_opt)
        
#         if np.abs(_n_pt) > 1.0e-8:
#             n_pt = (_n_pt**(-1))

#             if w >= 0.0:
#                 return np.sqrt((1 + n_pt) *  J_opt(w)) * np.sqrt((1 + n_pt) * J_opt(w))  # Spontaneous decay
#             else: 
#                 return np.sqrt(n_pt * J_opt(w)) * np.sqrt(n_pt * J_opt(w))        # Absorption
#         else:
#             if w >= 0.0:
#                 return np.sqrt(J_opt(w)) * np.sqrt(J_opt(w))        # Spontaneous decay
#             else: 
#                 return np.sqrt(J_opt(w)) * np.sqrt(J_opt(w))       # Absorption
            
#     def G_opt_22(w):

#         _n_pt = n_inv(w, beta_opt)
        
#         if np.abs(_n_pt) > 1.0e-8:
#             n_pt = (_n_pt**(-1))

#             if w >= 0.0:                
#                 return (1 + n_pt) * J_opt(w)  # Spontaneous decay
#             else: 
#                 return n_pt * J_opt(w)        # Absorption
#         else:
#             if w >= 0.0:
#                 return J_opt(w)        # Spontaneous decay
#             else: 
#                 return J_opt(w)        # Absorption
    
#     return np.array([[G_opt_11, G_opt_12], 
#                      [G_opt_21, G_opt_22]])


# # Scaling factor in the optical dissipator due to phonons
# def P_weight(kappa):
#     return np.array([[kappa**4, 1, kappa**2, kappa**2],
#                     [1, kappa**4, kappa**2, kappa**2],
#                     [kappa**2, kappa**2, kappa**4, 1],
#                     [kappa**2, kappa**2, 1, kappa**4]])























    # def Coup_rate_1(self):
    #     if self.mode == 'polaron':
    #         t_values = np.linspace(1.0e-4, 0.1, 1000)  # Discretize the t interval

    #         def G_c_1(w):
    #             # integrand_values = (
    #             #     (self.kappa()**4) * ((self.J / 2)**2) *
    #             #     (np.exp(2 * self.Phi_re(t_values)) *
    #             #      (np.cos(w * t_values) * np.cos(2 * self.Phi_im(t_values)) -
    #             #       np.sin(w * t_values) * np.sin(2 * self.Phi_im(t_values))) - 1)
    #             # )
    #             # result = np.trapz(integrand_values, t_values)  # Apply trapezoidal integration
    #             integrand = lambda t: (
    #                 mp.exp(1j * w * t) * (self.kappa()**4) * ((self.J / 2)**2) * (mp.exp(2 * self.Phi(t)) - 1)
    #             )
    #             # inte = np.array([integrand(t) for t in t_values], dtype = complex)
    #             # result = simpson(inte, t_values)
    #             result = mp.quad(integrand, [1.0e-4, 0.1])
    #             return float(result.real)

    #         def G_c_2(w):
    #             # integrand_values = (
    #             #     (self.kappa()**4) * ((self.J / 2)**2) *
    #             #     (np.exp(-2 * self.Phi_re(t_values)) *
    #             #      (np.cos(w * t_values) * np.cos(2 * self.Phi_im(t_values)) +
    #             #       np.sin(w * t_values) * np.sin(2 * self.Phi_im(t_values))) - 1)
    #             # )
    #             # result = np.trapz(integrand_values, t_values)  # Apply trapezoidal integration
    #             integrand = lambda t: (
    #                 mp.exp(1j * w * t) * (self.kappa()**4) * ((self.J / 2)**2) * (mp.exp(-2 * self.Phi(t)) - 1)
    #             )
    #             # inte = np.array([integrand(t) for t in t_values], dtype = complex)
    #             # result = simpson(inte, t_values)
    #             result = mp.quad(integrand, [1.0e-4, 0.1])
    #             return float(result.real)

    #         return np.array([[G_c_1, G_c_2], [G_c_2, G_c_1]])
    #     else:
    #         return None

    # def Coup_rate_2(self):
    #     if self.mode == 'polaron':
    #         t_values = np.linspace(1.0e-4, 0.2, 1000)  # Discretize the t interval

    #         def G_c_1(w):
    #             # integrand_values = (
    #             #     (self.kappa()**4) * ((self.J / 2)**2) *
    #             #     (np.exp(2 * self.Phi_re(-t_values)) *
    #             #      (np.cos(w * (-t_values)) * np.cos(2 * self.Phi_im(-t_values)) -
    #             #       np.sin(w * (-t_values)) * np.sin(2 * self.Phi_im(-t_values))) - 1)
    #             # )
    #             # result = np.trapz(integrand_values, t_values)  # Apply trapezoidal integration
    #             integrand = lambda t: (
    #                 mp.exp(1j * w * (-t)) * (self.kappa()**4) * ((self.J / 2)**2) * (mp.exp(2 * self.Phi(-t)) - 1)
    #             )
    #             # inte = np.array([integrand(t) for t in t_values], dtype = complex)
    #             # result = simpson(inte, t_values)
    #             result = mp.quad(integrand, [1.0e-4, 0.1])
    #             return float(result.real)

    #         def G_c_2(w):
    #             # integrand_values = (
    #             #     (self.kappa()**4) * ((self.J / 2)**2) *
    #             #     (np.exp(-2 * self.Phi_re(-t_values)) *
    #             #      (np.cos(w * (-t_values)) * np.cos(2 * self.Phi_im(-t_values)) +
    #             #       np.sin(w * (-t_values)) * np.sin(2 * self.Phi_im(-t_values))) - 1)
    #             # )
    #             # result = np.trapz(integrand_values, t_values)  # Apply trapezoidal integration
    #             integrand = lambda t: (
    #                 mp.exp(1j * w * (-t)) * (self.kappa()**4) * ((self.J / 2)**2) * (mp.exp(-2 * self.Phi(-t)) - 1)
    #             )
    #             # inte = np.array([integrand(t) for t in t_values], dtype = complex)
    #             # result = simpson(inte, t_values)
    #             result = mp.quad(integrand, [1.0e-4, 0.1])
    #             return float(result.real)

    #         return np.array([[G_c_1, G_c_2], [G_c_2, G_c_1]])
    #     else:
    #         return None





















    # def Phi(self, t):
    #     beta_vib = self._beta(self.T_vib)
    #     w = np.linspace(1.0e-4, 500, 10000)
    #     integrand = self.J_vib(w) * (np.cos(w * t)* (1 / np.tanh(beta_vib * w / 2)) - 1j * np.sin(w * t)) / w**2 # * (1 / np.tanh(beta_vib * w / 2))
    #     return np.trapz(integrand, w)
    #     # beta_vib = mp.mpf(self._beta(self.T_vib))
    #     # # Define the integrand using mpmath functions
    #     # def integrand(w, t):
    #     #     term1 = mp.cos(w * t) * (1 / mp.tanh(beta_vib * w / 2))
    #     #     term2 = -1j * mp.sin(w * t)
    #     #     return self.J_vib(w) * (term1 + term2) / w**2
    #     # # print(t)
    #     # result = mp.quadsubdiv(lambda w: integrand(w, t), [0, 2000], error = True)
    #     # # print('time:', t, float(result[0].real), result[1])
    #     # # Perform the integration using mpmath.quad
    #     # return result[0].real

    # def kappa(self):
    #     if self.mode == 'polaron':
    #         return np.exp(-0.5 * self.Phi(0))
    #     else:
    #         return 1.0
    

    # def Coup_rate_1(self):
    #     if self.mode == 'polaron':

    #         cache_w0 = None

    #         def G_c_1(w):
    #             nonlocal cache_w0  # Allow access to the cache

    #             if np.isnan(w):
    #                 return 0.0
    #             elif w == 0:
    #                 # Check if we already computed the value for w = 0
    #                 if cache_w0 is None:
    #                     tau_values = np.linspace(0, 1, 2000)
    #                     pref =  (self.kappa()**4) * ((self.J / 2)**2)
    #                     integrand = np.array([pref * (np.exp(2 * self.Phi(t)) - 1) for t in tau_values]) #np.array([(pref * np.exp(2 * self.Phi_re(t)) * (np.cos(w * t) * np.cos(2 * self.Phi_im(t)) - np.sin(w * t) * np.sin(2 * self.Phi_im(t))) - 1) for t in tau_values]) # np.array([np.exp(1j * w * t) * pref * (np.exp(2 * self.Phi(t)) - 1) for t in tau_values])
    #                     print(w, np.trapz(integrand, tau_values).real)
    #                     cache_w0 = np.trapz(integrand, tau_values).real
    #                 return cache_w0  
    #             else:
    #                 tau_values = np.linspace(0, 1, 2000)
    #                 pref =  (self.kappa()**4) * ((self.J / 2)**2)
    #                 integrand = np.array([np.exp(1j * w * t) * pref * (np.exp(2 * self.Phi(t)) - 1) for t in tau_values]) #np.array([(pref * np.exp(2 * self.Phi_re(t)) * (np.cos(w * t) * np.cos(2 * self.Phi_im(t)) - np.sin(w * t) * np.sin(2 * self.Phi_im(t))) - 1) for t in tau_values]) # np.array([np.exp(1j * w * t) * pref * (np.exp(2 * self.Phi(t)) - 1) for t in tau_values])
    #                 result = np.trapz(integrand, tau_values).real
    #                 print(w, result)
    #                 return result 

    #         def G_c_2(w):
    #             nonlocal cache_w0  # Allow access to the cache

    #             if np.isnan(w):
    #                 return 0.0
    #             elif w == 0:
    #                 # Check if we already computed the value for w = 0
    #                 if cache_w0 is None:
    #                     tau_values = np.linspace(0, 1, 2000)
    #                     pref =  (self.kappa()**4) * ((self.J / 2)**2)
    #                     integrand = np.array([pref * (np.exp(-2 * self.Phi(t)) - 1) for t in tau_values]) #np.array([(pref * np.exp(2 * self.Phi_re(t)) * (np.cos(w * t) * np.cos(2 * self.Phi_im(t)) - np.sin(w * t) * np.sin(2 * self.Phi_im(t))) - 1) for t in tau_values]) # np.array([np.exp(1j * w * t) * pref * (np.exp(2 * self.Phi(t)) - 1) for t in tau_values])
    #                     print(w, np.trapz(integrand, tau_values).real)
    #                     cache_w0 = np.trapz(integrand, tau_values).real
    #                 return cache_w0  
    #             else:
    #                 tau_values = np.linspace(0, 1, 2000)
    #                 pref =  (self.kappa()**4) * ((self.J / 2)**2)
    #                 integrand = np.array([np.exp(1j * w * t) * pref * (np.exp(-2 * self.Phi(t)) - 1) for t in tau_values]) #np.array([(pref * np.exp(2 * self.Phi_re(t)) * (np.cos(w * t) * np.cos(2 * self.Phi_im(t)) - np.sin(w * t) * np.sin(2 * self.Phi_im(t))) - 1) for t in tau_values]) # np.array([np.exp(1j * w * t) * pref * (np.exp(2 * self.Phi(t)) - 1) for t in tau_values])
    #                 result = np.trapz(integrand, tau_values).real
    #                 print(w, result)
    #                 return result 

    #         return np.array([[lambda w: 0, G_c_1], [G_c_1, lambda w: 0]])
    #         # return np.array([[G_c_2, G_c_1], [G_c_1, G_c_2]])
    #     else:
    #         return None

    # def Coup_rate_2(self):
    #     if self.mode == 'polaron':
    #         cache_w0 = None
    #         def G_c_1(w):
    #             nonlocal cache_w0  # Allow access to the cache

    #             if np.isnan(w):
    #                 return 0.0
    #             elif w == 0:
    #                 # Check if we already computed the value for w = 0
    #                 if cache_w0 is None:
    #                     tau_values = np.linspace(0, 1, 2000)
    #                     pref =  (self.kappa()**4) * ((self.J / 2)**2)
    #                     integrand = np.array([pref * (np.exp(2 * self.Phi(t)) - 1) for t in tau_values]) #np.array([(pref * np.exp(2 * self.Phi_re(t)) * (np.cos(w * t) * np.cos(2 * self.Phi_im(t)) - np.sin(w * t) * np.sin(2 * self.Phi_im(t))) - 1) for t in tau_values]) # np.array([np.exp(1j * w * t) * pref * (np.exp(2 * self.Phi(t)) - 1) for t in tau_values])
    #                     print(w, np.trapz(integrand, tau_values).real)
    #                     cache_w0 = np.trapz(integrand, tau_values).real
    #                 return cache_w0  
    #             else:
    #                 tau_values = np.linspace(0, 1, 2000)
    #                 pref =  (self.kappa()**4) * ((self.J / 2)**2)
    #                 integrand = np.array([np.exp(1j * w * t) * pref * (np.exp(2 * self.Phi(t)) - 1) for t in tau_values]) #np.array([(pref * np.exp(2 * self.Phi_re(t)) * (np.cos(w * t) * np.cos(2 * self.Phi_im(t)) - np.sin(w * t) * np.sin(2 * self.Phi_im(t))) - 1) for t in tau_values]) # np.array([np.exp(1j * w * t) * pref * (np.exp(2 * self.Phi(t)) - 1) for t in tau_values])
    #                 result = np.trapz(integrand, tau_values).real
    #                 print(w, result)
    #                 return result 
    
    #         def G_c_2(w):
    #             nonlocal cache_w0  # Allow access to the cache

    #             if np.isnan(w):
    #                 return 0.0
    #             elif w == 0:
    #                 # Check if we already computed the value for w = 0
    #                 if cache_w0 is None:
    #                     tau_values = np.linspace(0, 1, 2000)
    #                     pref =  (self.kappa()**4) * ((self.J / 2)**2)
    #                     integrand = np.array([pref * (np.exp(-2 * self.Phi(-t)) - 1) for t in tau_values]) #np.array([(pref * np.exp(2 * self.Phi_re(t)) * (np.cos(w * t) * np.cos(2 * self.Phi_im(t)) - np.sin(w * t) * np.sin(2 * self.Phi_im(t))) - 1) for t in tau_values]) # np.array([np.exp(1j * w * t) * pref * (np.exp(2 * self.Phi(t)) - 1) for t in tau_values])
    #                     print(w, np.trapz(integrand, tau_values).real)
    #                     cache_w0 = np.trapz(integrand, tau_values).real
    #                 return cache_w0  
    #             else:
    #                 tau_values = np.linspace(0, 1, 2000)
    #                 pref =  (self.kappa()**4) * ((self.J / 2)**2)
    #                 integrand = np.array([np.exp(1j * w * t) * pref * (np.exp(-2 * self.Phi(-t)) - 1) for t in tau_values]) #np.array([(pref * np.exp(2 * self.Phi_re(t)) * (np.cos(w * t) * np.cos(2 * self.Phi_im(t)) - np.sin(w * t) * np.sin(2 * self.Phi_im(t))) - 1) for t in tau_values]) # np.array([np.exp(1j * w * t) * pref * (np.exp(2 * self.Phi(t)) - 1) for t in tau_values])
    #                 result = np.trapz(integrand, tau_values).real
    #                 print(w, result)
    #                 return result 

    #         return np.array([[lambda w: 0, G_c_1], [G_c_1, lambda w: 0]])
    #         # return np.array([[G_c_2, G_c_1], [G_c_1, G_c_2]])
    #     else:
    #         return None


    # def Phi(self, t):
    #     # beta_vib = self._beta(self.T_vib)
    #     # w = np.linspace(1.0e-4, 2000, 5000)
    #     # integrand = self.J_vib(w) * (np.cos(w * t) * (1 / np.tanh(beta_vib * w / 2)) - 1j * np.sin(w * t)) / w**2 # * (1 / np.tanh(beta_vib * w / 2))
    #     # print(np.trapz(integrand, w))
    #     beta_vib = mp.mpf(self._beta(self.T_vib))
    #     # Define the integrand using mpmath functions
    #     def integrand(w, t):
    #         term1 = mp.cos(w * t) * (1 / mp.tanh(beta_vib * w / 2))
    #         term2 = -1j * mp.sin(w * t)
    #         return self.J_vib(w) * (term1 + term2) / w**2
    #     # print(t)
    #     result = mp.quadsubdiv(lambda w: integrand(w, t), [0, 2000], error = True)
    #     # print('time:', t, float(result[0].real), result[1])
    #     # Perform the integration using mpmath.quad
    #     return result[0].real

    #     # beta_vib = mp.mpf(self._beta(self.T_vib))
    #     # omega_c = mp.mpf(9.0e-3 * constants.eV / system_params.hbar_ps)
    #     # E_reorg = mp.mpf(20.0e-3 * constants.eV / system_params.hbar_ps)
    #     # E_reorg_J = mp.mpf(20.0e-3)

    #     # epsilon = 1 / (beta_vib * omega_c)

    #     # result = - E_reorg_J * ((1 / (1 - 1j * omega_c * t)**2) - (epsilon**2) * (mp.psi(1, epsilon + 1j * t / beta_vib) + mp.psi(1, epsilon - 1j * t / beta_vib)))
    #     # return result

    # def kappa(self):
    #     # beta_vib = mp.mpf(self._beta(self.T_vib))
    #     # omega_c = mp.mpf(9.0e-3 * constants.eV / system_params.hbar_ps)
    #     # E_reorg = mp.mpf(20.0e-3 * constants.eV / system_params.hbar_ps)
    #     # E_reorg_J = mp.mpf(20.0e-3)
        

    #     if self.mode == 'polaron':
    #         # epsilon = 1 / (beta_vib * omega_c)
    #         # return float(mp.exp(-0.5 * E_reorg_J * (2 * (epsilon**2) * mp.psi(1, epsilon) - 1)))
    #         return float(mp.exp(-0.5 * self.Phi(0).real))
    #     else:
    #         return 1.0

    # def Coup_rate_1(self):
    #     if self.mode == 'polaron':
    #         pref =  (self.kappa()**4) * ((self.J / 2)**2)
    #         beta_vib = mp.mpf(self._beta(self.T_vib))
    #         cache_w0 = None

    #         def G_c_1(w):
    #             nonlocal cache_w0  # Allow access to the cache

    #             if np.isnan(w):
    #                 return 0.0
    #             elif w == 0:
    #                 # Check if we already computed the value for w = 0
    #                 if cache_w0 is None:
    #                     integrand = lambda t: (
    #                         pref * (mp.exp(2 * self.Phi(t)) - 1) # * (mp.exp(2 * self.Phi_re(t)) * (mp.cos(2 * self.Phi_im(t))) - 1) 
    #                     )
    #                     result = mp.quad(integrand, [0, 0.2], error=True)
    #                     cache_w0 = float(result[0].real)  # Store the computed value in the cache
    #                     print('freq:', w, cache_w0, result[1])
    #                 return cache_w0  # Return the cached value
    #             else:      
    #                 integrand = lambda t: (
    #                     pref * mp.exp(1j * w * t) * (mp.exp(2 * self.Phi(t)) - 1) #(mp.exp(2 * self.Phi_re(t)) * (mp.cos(w * t) * mp.cos(2 * self.Phi_im(t)) - mp.sin(w * t) * mp.sin(2 * self.Phi_im(t))) - 1)
    #                 )
    #                 # result = mp.quadosc(integrand, [0, mp.inf], omega = mp.mpf(w))
    #                 # print(w, float(result.real))
    #                 # return float(result.real)
    #                 result = mp.quad(integrand, [0, 0.2], error = True)
    #                 print('freq:', w, float(result[0].real), result[1])
    #                 return float(result[0].real)

    #         return np.array([[lambda w: 0, G_c_1], [G_c_1, lambda w: 0]])
    #     else:
    #         return None

    # def Coup_rate_2(self):
    #     if self.mode == 'polaron':
    #         pref =  (self.kappa()**4) * ((self.J / 2)**2)
    #         cache_w0 = None

    #         def G_c_1(w):
    #             nonlocal cache_w0  # Allow access to the cache

    #             if np.isnan(w):
    #                 return 0.0
    #             elif w == 0:
    #                 # Check if we already computed the value for w = 0
    #                 if cache_w0 is None:
    #                     integrand = lambda t: (
    #                         pref * (mp.exp(2 * self.Phi(-t)) - 1)  # * (mp.exp(2 * self.Phi_re(t)) * (mp.cos(2 * self.Phi_im(t))) - 1) 
    #                     )
    #                     result = mp.quad(integrand, [0, 0.2], error=True)
    #                     cache_w0 = float(result[0].real)  # Store the computed value in the cache
    #                     print('freq:', w, cache_w0, result[1])
    #                 return cache_w0  # Return the cached value
    #             else:      
    #                 integrand = lambda t: (
    #                     mp.exp(1j * w * t) * pref * (mp.exp(2 * self.Phi(-t)) - 1)
    #                 )
    #                 # result = mp.quadosc(integrand, [0, mp.inf], omega = mp.mpf(w))
    #                 # print(w, float(result.real))
    #                 # return float(result.real)
    #                 result = mp.quad(integrand, [0, 0.2], error = True)
    #                 print('freq:', w, float(result[0].real), result[1])
    #                 return float(result[0].real)
                
    #         return np.array([[lambda w: 0, G_c_1], [G_c_1, lambda w: 0]])
    #     else:
    #         return None