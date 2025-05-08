from rates import *



# constants = Constants()
# dipole_params = DipoleParameters(config_name)  
# system_params = SystemParameters(E_q1=1.8*constants.eV, E_q2=1.8*constants.eV, dipole_params=dipole_params)
# system_params.calculate_parameters()
# J_opt, J_vib = calculate_spectral_densities()
# rates = CombinedRates(T_opt, T_vib, J_opt, J_vib, mode)
# opt_rates = rates.Opt_rate()
# coup_rates = rates.Coup_rate()
# p_vib = rates.P_weight()

rates = CombinedRates(T_opt, T_vib, mode)

############################################################
'''Pauli matrices'''
############################################################

def pauli_matrices():
    sx = np.array([[0, 1], [1, 0]], dtype=complex)
    sz = np.array([[1, 0], [0, -1]], dtype=complex)
    sp = np.array([[0, 1], [0, 0]], dtype=complex)
    sm = np.array([[0, 0], [1, 0]], dtype=complex)
    
    return sx, sz, sp, sm

def two_qubit_pauli_matrices():
    sx, sz, sp, sm = pauli_matrices()
    sx1 = np.kron(sx, np.identity(2))
    sx2 = np.kron(np.identity(2), sx)
    sz1 = np.kron(sz, np.identity(2))
    sz2 = np.kron(np.identity(2), sz)
    sp1 = np.kron(sp, np.identity(2))
    sp2 = np.kron(np.identity(2), sp)
    sm1 = np.kron(sm, np.identity(2))
    sm2 = np.kron(np.identity(2), sm)
    
    return sx1, sx2, sz1, sz2, sp1, sp2, sm1, sm2


############################################################
'''Hamiltonian'''
############################################################

def hamiltonian(mode, kappa):
    sx1, sx2, sz1, sz2, sp1, sp2, sm1, sm2 = two_qubit_pauli_matrices()
    
    if mode == 'polaron':
        # shift = float(mp.quad(lambda w: J_vib(w) / w, [0, np.inf]))
        omega = np.linspace(1.0e-4, 2000, 10000)
        integrand = np.array([rates.J_vib(w) / w for w in omega]) # * (1 / np.tanh(beta_vib * w / 2))
        shift = float(np.trapz(integrand, omega))
    
    elif mode == 'weak':
        shift = 0.0

    # Considering both dipoles energies get polaron shifted by the same amount
    omega_1 = system_params.omega_q1 - shift
    omega_2 = system_params.omega_q2 - shift

    # # kappa_1 and kappa_2 scaling the dipole dipole coupling are also considered the same
    J_scaled = system_params.J * (kappa**2)

    H = omega_1 * sp1 @ sm1 + omega_2 * sp2 @ sm2 + (J_scaled / 2) * (sp1 @ sm2 + sm1 @ sp2)
    return H


# ############################################################
# '''Project system operators onto Hamiltonian's basis'''
# ############################################################

# def calc_a_ops_freq(H, a_op):
#     evals, ekets = eigvalvec(H)
#     a_ops = [np.dot(ekets.conj().T, np.dot(a, ekets)) for a in a_op]

#     dim = H.shape[0]
#     trans_freq = np.array([evals[i]-evals[j] for i in range(dim) for j in range(dim)])
#     return a_ops, trans_freq


############################################################
'''Other utilities'''
############################################################

def eigvalvec(H):
    evals, ekets = np.linalg.eigh(H)
    return evals, ekets

def psi2rho(psi0):
    # initial state as density operator
    rho0 = np.array(np.outer(psi0, psi0.conj().T), dtype=complex)
    return rho0

def initial_state(H, psi0):
    evals, ekets = eigvalvec(H)
    # initial state as density operator
    rho_sz = np.array(np.outer(psi0, psi0.conj().T), dtype=complex)
    rho0 = np.dot(ekets.conj().T, np.dot(rho_sz, ekets))
    H_diag = np.dot(ekets.conj().T, np.dot(H, ekets))
    return H_diag, rho0



def cart2sph(vector):
    x,y,z = vector[0], vector[1], vector[2]
    XsqPlusYsq = x**2 + y**2
    r = np.sqrt(XsqPlusYsq + z**2)               # r
    theta = np.arccos(z / r)     # theta
    phi = np.arctan2(y,x)                           # phi
    return r, theta, phi

def unit_vector(vector):
    """ Returns the unit vector of the vector.  """
    return vector / np.linalg.norm(vector)

def angle_between_pol(v1, v2): # Plane polar
    """ Returns the angle in radians between vectors 'v1' and 'v2'::
    """
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    return np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))


def angles_between_sphpol(u, v): # Spherical polar
    # Calculate dot product
    dot_product = np.dot(u, v)
    
    # Calculate magnitudes
    norm_u = np.linalg.norm(u)
    norm_v = np.linalg.norm(v)
    
    # Calculate theta (angle between the vectors)
    cos_theta = dot_product / (norm_u * norm_v)
    theta = np.arccos(np.clip(cos_theta, -1.0, 1.0))  # Clip values to avoid numerical errors

    # Calculate cross product
    cross_product = np.cross(u, v)
    
    # Calculate phi (angle between the vectors using cross product)
    sin_phi = np.linalg.norm(cross_product) / (norm_u * norm_v)
    phi = np.arcsin(np.clip(sin_phi, -1.0, 1.0))  # Clip values to avoid numerical errors
    
    # Convert radians to degrees
    theta_deg = np.degrees(theta)
    phi_deg = np.degrees(phi)
    
    return theta, phi


def spherical_to_plane_polar(r, theta, phi):
    # Spherical to Cartesian
    x = r * np.sin(theta) * np.cos(phi)
    y = r * np.sin(theta) * np.sin(phi)
    z = r * np.cos(theta)
    
    # Cartesian to Plane Polar
    R = np.sqrt(x**2 + y**2)
    Theta = np.arctan2(y, x)  # arctan2 handles the correct quadrant
    
    return R, Theta



def orthogonal_vectors(theta, phi):
    v = np.array([np.cos(phi) * np.sin(theta), np.sin(phi) * np.sin(theta), np.cos(theta)])
    
    # Check if v is aligned with the z-axis. If so, use x-axis as reference
    if np.allclose(v, np.array([0, 0, 1])) or np.allclose(v, np.array([0, 0, -1])):
        reference = np.array([1, 0, 0])
    else:
        reference = np.array([0, 0, 1])
    
    # First orthogonal vector
    orthogonal_1 = np.cross(v, reference)
    orthogonal_1 /= np.linalg.norm(orthogonal_1)  # Normalize the vector
    
    # Second orthogonal vector
    orthogonal_2 = np.cross(v, orthogonal_1)
    orthogonal_2 /= np.linalg.norm(orthogonal_2)  # Normalize the vector

    # orthogonal_1 = np.array([np.cos(phi) * np.cos(theta), np.sin(phi) * np.cos(theta), -np.sin(theta)])
    # orthogonal_2 = np.array([-np.sin(phi), np.cos(phi), 0])
    
    # return np.round(orthogonal_1), np.round(orthogonal_2)
    return orthogonal_1, orthogonal_2




def rotate_vector(vector, theta, phi):
    # Convert angles from degrees to radians
    theta_rad = np.radians(theta)
    phi_rad = np.radians(phi)
    
    # Rotation matrix around z-axis by theta (yaw)
    R_z = np.array([
        [np.cos(theta_rad), -np.sin(theta_rad), 0],
        [np.sin(theta_rad), np.cos(theta_rad), 0],
        [0, 0, 1]
    ])
    
    # Rotation matrix around y-axis by phi (pitch)
    R_y = np.array([
        [np.cos(phi_rad), 0, np.sin(phi_rad)],
        [0, 1, 0],
        [-np.sin(phi_rad), 0, np.cos(phi_rad)]
    ])
    
    # Combined rotation matrix
    R = np.dot(R_y, R_z)
    
    # Rotate the vector
    rotated_vector = np.dot(R, vector)
    
    return rotated_vector






# Function to load data and plot
def load_and_plot(ax, filename, times_ns, label, marker=None):
    data = np.load(filename)
    times = np.append(-np.flip(times_ns[1:] - times_ns[0]), times_ns)
    g2_data = np.append(np.flip(data[1:]), data)
    ax.plot(times, g2_data, label=label, marker = marker)

