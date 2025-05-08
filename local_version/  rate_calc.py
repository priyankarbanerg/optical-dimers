import numpy as np
import matplotlib.pyplot as plt

eV = 1.6 * 1.0e-19 #J
k_B = 1.380649e-23 / eV  # J/K
ns = 1.0e-9 #s

T_vib = 1.0 * 300 #K

omega_c = 90.0e-3 #* eV 
E_reorg = 5.0e-3 #* eV 

J = 0.01 #* eV

def J_vib(w):
    k_vib_2 = E_reorg / (2 * omega_c**3)
    return (k_vib_2 * np.abs(w)**3) * np.exp(- (np.abs(w) / omega_c))

def beta(T):
    if T == 0.0:
        return np.inf
    elif T == np.inf:
        return 0.0
    else:
        return 1 / (k_B * T)
    
def Phi(t):
    beta_vib = beta(T_vib)
    w = np.linspace(1.0e-4 , 1 , 5000)
    integrand = J_vib(w) * (np.cos(w * t)* (1 / np.tanh(beta_vib * w / 2)) - 1j * np.sin(w * t)) / w**2
    return np.trapz(integrand, w)

def kappa():
    return np.exp(-0.5 * Phi(0))
    



def Coup_rate_1():
    def G_c_1(w):
        tau_values = np.linspace(0, 100, 10000)
        pref =  (kappa()**4) * ((J / 2)**2)
        integrand = np.array([np.exp(1j * w * t) * pref * (np.exp(2 * Phi(t)) - 1) for t in tau_values]) 
        result = np.trapz(integrand, tau_values).real
        print(w, result)
        return result  

    def G_c_2(w):
        tau_values = np.linspace(0, 100, 10000)
        pref =  (kappa()**4) * ((J / 2)**2)
        integrand = np.array([np.exp(1j * w * t) * pref * (np.exp(-2 * Phi(t)) - 1) for t in tau_values]) 
        result = np.trapz(integrand, tau_values).real
        print(w, result)
        return (result.real)

    # return np.array([[lambda w: 0, G_c_1], [lambda w: 0, G_c_1]])
    return np.array([[G_c_1, G_c_2], [G_c_1, G_c_2]])


def Coup_rate_2():
    def G_c_1(w):
        tau_values = np.linspace(0, 1000, 10000)
        pref =  (kappa()**4) * ((J / 2)**2)
        integrand = np.array([np.exp(1j * w * t) * pref * (np.exp(2 * Phi(-t)) - 1) for t in tau_values]) 
        result = np.trapz(integrand, tau_values).real
        print(w, result)
        return result  

    def G_c_2(w):
        tau_values = np.linspace(0, 1000, 1000)
        pref =  (kappa()**4) * ((J / 2)**2)
        integrand = np.array([np.exp(1j * w * t) * pref * (np.exp(-2 * Phi(-t)) - 1) for t in tau_values]) 
        result = np.trapz(integrand, tau_values).real
        print(w, result)
        return (result.real)

    # return np.array([[lambda w: 0, G_c_1], [lambda w: 0, G_c_1]])
    return np.array([[G_c_1, G_c_2], [G_c_1, G_c_2]])


# trans_freq = np.linspace(-0.5, 0.5 , 10)

# Gamma_coup_1 = Coup_rate_1()
# Gamma_coup_2 = Coup_rate_2()


# G1 = np.zeros((2, 2, len(trans_freq)), dtype = complex)
# G2 = np.zeros((2, 2, len(trans_freq)), dtype = complex)

# for (i, j), G_ab_1 in np.ndenumerate(Gamma_coup_1):
    # G1[i, j, :] = np.array([G_ab_1(w) for w in trans_freq])
# for (i, j), G_ab_2 in np.ndenumerate(Gamma_coup_2):
#     G2[i, j, :] = np.array([G_ab_2(w) for w in trans_freq])


# w = np.linspace(1.0e-4, 2, 10000)
# plt.plot(w, J_vib(w))
# plt.xlabel('t')
# plt.ylabel(r'$J_{vib}$')
# plt.legend()
# plt.show()    


tau_values = np.linspace(0, 100, 1000)
w = 10
# plt.plot(tau_values, np.array([Phi(t) for t in tau_values])) #, label = r'$phi(t)$')
plt.plot(tau_values, np.array([Phi(t) for t in tau_values]).real )
plt.plot(tau_values, np.array([Phi(t) for t in tau_values]).imag )
plt.xlabel('t')
plt.ylabel(r'$phi(t)$')
plt.legend()
plt.show()



# # print(np.sin(1000))
tau_values = np.linspace(0, 100, 1000)
w = 10
# plt.plot(tau_values, np.array([Phi(t) for t in tau_values])) #, label = r'$phi(t)$')
plt.plot(tau_values, np.array([(kappa()**4) * ((J / 2)**2) * (np.exp(2 * Phi(t)) - 1) for t in tau_values]).real )
plt.plot(tau_values, np.array([(kappa()**4) * ((J / 2)**2) * (np.exp(2 * Phi(t)) - 1) for t in tau_values]).imag )
plt.xlabel('t')
plt.ylabel(r'$corr(t)$')
plt.legend()
plt.show()


# plt.plot(trans_freq, G1[0, 0], label = r'Rate - $\sigma_{1}^{+}(t)\sigma_{2}^{-}(t)\sigma_{1}^{+}(t-\tau)\sigma_{2}^{-}(t-\tau)$')
# # plt.plot(trans_freq, G1[0, 1], label = r'Rate - $\sigma_{1}^{+}(t)\sigma_{2}^{-}(t)\sigma_{2}^{+}(t-\tau)\sigma_{1}^{-}(t-\tau)$')
# # plt.plot(trans_freq, G1[1, 0], label = r'Rate - $\sigma_{2}^{+}(t)\sigma_{1}^{-}(t)\sigma_{1}^{+}(t-\tau)\sigma_{2}^{-}(t-\tau)$')
# # plt.plot(trans_freq, G1[1, 1], label = r'Rate - $\sigma_{2}^{+}(t)\sigma_{1}^{-}(t)\sigma_{2}^{+}(t-\tau)\sigma_{1}^{-}(t-\tau)$')
# plt.xlabel('w')
# plt.ylabel(r'Rates$')
# plt.legend()
# plt.show()