import numpy as np

def pauli_matrices() -> tuple:
    """
    Construct the standard single-qubit Pauli matrices.

    Returns
    -------
    tuple
        A tuple containing (sx, sz, sp, sm) matrices of type complex.
    """
    sx = np.array([[0, 1], [1, 0]], dtype=complex)
    sz = np.array([[1, 0], [0, -1]], dtype=complex)
    sp = np.array([[0, 1], [0, 0]], dtype=complex)
    sm = np.array([[0, 0], [1, 0]], dtype=complex)
    
    return sx, sz, sp, sm

def two_qubit_pauli_matrices() -> tuple:
    """
    Construct the two-qubit Pauli matrices in the tensor product space.

    Returns
    -------
    tuple
        A tuple containing (sx1, sx2, sz1, sz2, sp1, sp2, sm1, sm2) matrices.
    """
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

class SystemOperators:
    """
    Construct the Hamiltonian and evaluate states in the system eigenbasis.

    Parameters
    ----------
    system_params : SystemParameters
        Object containing the physical system parameters and coherent driving amplitudes.
    rates : CombinedRates
        Object containing the spectral density functions and evaluated polaron parameters.
    """

    def __init__(self, system_params, rates):
        self.system_params = system_params
        self.rates = rates
        self.mode = self.system_params.config_params.mode
        self.pump = self.system_params.config_params.pump
        
        self.sx1, self.sx2, self.sz1, self.sz2, self.sp1, self.sp2, self.sm1, self.sm2 = two_qubit_pauli_matrices()
        
        self.H = self.hamiltonian()
        self.evals, self.ekets = np.linalg.eigh(self.H)

    def hamiltonian(self) -> np.ndarray:
        """
        Construct the time-independent system Hamiltonian in the rotating frame of the laser.

        Returns
        -------
        ndarray
            The 4x4 system Hamiltonian matrix.
        """
        if self.mode == 'polaron':
            omega = np.linspace(1.0e-4, 1500, 10000)
            integrand = np.array([self.rates.J_vib(w) / w for w in omega])
            shift = float(np.trapezoid(integrand, omega))
            DB_factor = self.rates.kappa()**2
        elif self.mode == 'weak':
            shift = 0.0
            DB_factor = 1.0
        else:
            raise ValueError(f"Unknown mode: {self.mode}")

        omega_1_shifted = self.system_params.omega_q1 - shift
        omega_2_shifted = self.system_params.omega_q2 - shift
        
        detun_laser = ((omega_1_shifted + omega_2_shifted) / 2.0) - self.system_params.omega_L
        detun_qubits = (omega_2_shifted - omega_1_shifted) / 2.0
        
        J_scaled = self.system_params.J * DB_factor

        H = (detun_laser - detun_qubits) * self.sp1 @ self.sm1 + \
            (detun_laser + detun_qubits) * self.sp2 @ self.sm2 + \
            (J_scaled / 2) * (self.sp1 @ self.sm2 + self.sm1 @ self.sp2)

        if self.pump == 'coh':
            if self.mode != 'weak':
                raise ValueError("Coherent driving is strictly limited to the 'weak' coupling mode.")
            
            H_d = (self.system_params.Omega_1) * (self.sp1 + self.sm1) + \
                  (self.system_params.Omega_2) * (self.sp2 + self.sm2)
            H += H_d

        return H
    
    def psi2rho(self, psi0: np.ndarray) -> np.ndarray:
        """
        Convert a state vector to a density matrix.

        Parameters
        ----------
        psi0 : ndarray
            The initial state vector.

        Returns
        -------
        ndarray
            The corresponding density matrix.
        """
        return np.array(np.outer(psi0, psi0.conj().T), dtype=complex)

    def initial_state(self, psi0: np.ndarray) -> tuple:
        """
        Transform the initial state and Hamiltonian into the Hamiltonian eigenbasis.

        Parameters
        ----------
        psi0 : ndarray
            The initial state vector in the site basis.

        Returns
        -------
        tuple
            A tuple containing the diagonalized Hamiltonian and the initial density matrix in the eigenbasis.
        """
        rho_sz = self.psi2rho(psi0)
        rho0 = self.ekets.conj().T @ rho_sz @ self.ekets
        H_diag = self.ekets.conj().T @ self.H @ self.ekets
        return H_diag, rho0


# ---------------------------------------------------------
# Standalone Geometric and Coordinate Utilities
# ---------------------------------------------------------

def cart2sph(vector: np.ndarray) -> tuple:
    """
    Convert Cartesian coordinates to spherical coordinates.

    Parameters
    ----------
    vector : ndarray
        A 3D vector [x, y, z].

    Returns
    -------
    tuple
        (r, theta, phi) coordinates.
    """
    x, y, z = vector[0], vector[1], vector[2]
    XsqPlusYsq = x**2 + y**2
    r = np.sqrt(XsqPlusYsq + z**2)
    theta = np.arccos(z / r)
    phi = np.arctan2(y, x)
    return r, theta, phi

def unit_vector(vector: np.ndarray) -> np.ndarray:
    """
    Return the normalized unit vector.

    Parameters
    ----------
    vector : ndarray
        Input vector.

    Returns
    -------
    ndarray
        Unit vector.
    """
    return vector / np.linalg.norm(vector)

def angle_between_pol(v1: np.ndarray, v2: np.ndarray) -> float:
    """
    Return the angle in radians between two vectors.

    Parameters
    ----------
    v1 : ndarray
        First vector.
    v2 : ndarray
        Second vector.

    Returns
    -------
    float
        Angle in radians.
    """
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    return np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))

def angles_between_sphpol(u: np.ndarray, v: np.ndarray) -> tuple:
    """
    Evaluate the spherical polar angles between two vectors.

    Parameters
    ----------
    u : ndarray
        First vector.
    v : ndarray
        Second vector.

    Returns
    -------
    tuple
        (theta_deg, phi_deg) in degrees.
    """
    dot_product = np.dot(u, v)
    norm_u = np.linalg.norm(u)
    norm_v = np.linalg.norm(v)
    
    cos_theta = dot_product / (norm_u * norm_v)
    theta = np.arccos(np.clip(cos_theta, -1.0, 1.0))

    cross_product = np.cross(u, v)
    sin_phi = np.linalg.norm(cross_product) / (norm_u * norm_v)
    phi = np.arcsin(np.clip(sin_phi, -1.0, 1.0))
    
    theta_deg = np.degrees(theta)
    phi_deg = np.degrees(phi)
    
    return theta_deg, phi_deg

def spherical_to_plane_polar(r: float, theta: float, phi: float) -> tuple:
    """
    Convert spherical coordinates to plane polar coordinates.

    Parameters
    ----------
    r : float
        Radial distance.
    theta : float
        Polar angle.
    phi : float
        Azimuthal angle.

    Returns
    -------
    tuple
        (R, Theta) plane polar coordinates.
    """
    x = r * np.sin(theta) * np.cos(phi)
    y = r * np.sin(theta) * np.sin(phi)
    
    R = np.sqrt(x**2 + y**2)
    Theta = np.arctan2(y, x)
    
    return R, Theta

def orthogonal_vectors(theta: float, phi: float) -> tuple:
    """
    Compute two orthogonal polarization unit vectors for a given propagation direction.

    Parameters
    ----------
    theta : float
        Polar angle of the propagation vector.
    phi : float
        Azimuthal angle of the propagation vector.

    Returns
    -------
    tuple
        (orthogonal_1, orthogonal_2) unit vectors.
    """
    v = np.array([np.cos(phi) * np.sin(theta), np.sin(phi) * np.sin(theta), np.cos(theta)])
    
    if np.allclose(v, np.array([0, 0, 1])) or np.allclose(v, np.array([0, 0, -1])):
        reference = np.array([1, 0, 0])
    else:
        reference = np.array([0, 0, 1])
    
    orthogonal_1 = np.cross(v, reference)
    orthogonal_1 /= np.linalg.norm(orthogonal_1)
    
    orthogonal_2 = np.cross(v, orthogonal_1)
    orthogonal_2 /= np.linalg.norm(orthogonal_2)

    return orthogonal_1, orthogonal_2

def rotate_vector(vector: np.ndarray, theta: float, phi: float) -> np.ndarray:
    """
    Rotate a vector by the specified pitch and yaw angles.

    Parameters
    ----------
    vector : ndarray
        The 3D vector to rotate.
    theta : float
        Yaw angle in degrees.
    phi : float
        Pitch angle in degrees.

    Returns
    -------
    ndarray
        The rotated vector.
    """
    theta_rad = np.radians(theta)
    phi_rad = np.radians(phi)
    
    R_z = np.array([
        [np.cos(theta_rad), -np.sin(theta_rad), 0],
        [np.sin(theta_rad), np.cos(theta_rad), 0],
        [0, 0, 1]
    ])
    
    R_y = np.array([
        [np.cos(phi_rad), 0, np.sin(phi_rad)],
        [0, 1, 0],
        [-np.sin(phi_rad), 0, np.cos(phi_rad)]
    ])
    
    R = np.dot(R_y, R_z)
    rotated_vector = np.dot(R, vector)
    
    return rotated_vector
 
 
def concurrence(rho: np.ndarray) -> float:
        """
        Compute the concurrence of a bipartite two-photon density matrix.
        To be added to system_ops.py.
        """
        trace_rho = np.trace(rho)
        if not np.isclose(trace_rho, 0.0):
            rho = rho / trace_rho
            
        sigma_y = np.array([[0, -1j], [1j, 0]], dtype=complex)
        rho_tilde = np.kron(sigma_y, sigma_y) @ rho.conj() @ np.kron(sigma_y, sigma_y)
        
        R = rho @ rho_tilde
        eigvals = np.real(np.linalg.eigvals(R))
        eigvals = np.clip(eigvals, 0, None)
        
        lambdas = np.sort(np.sqrt(eigvals))[::-1]
        C = max(0.0, float(lambdas[0] - lambdas[1] - lambdas[2] - lambdas[3]))
        
        return float(np.real_if_close(C))