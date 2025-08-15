import numpy as np
from scipy.interpolate import InterpolatedUnivariateSpline
from scipy.io import netcdf_file
from scipy.optimize import root

__all__ = [
    "boozer_to_cylindrical",
    "cylindrical_to_boozer",
    "boozer_to_vmec",
    "vmec_to_boozer",
    "vmec_to_cylindrical",
    "cylindrical_to_vmec",
]


def boozer_to_cylindrical(field, s, theta, zeta):
    r"""
    Convert from Boozer coordinates to cylindrical coordinates.

    Args:
        field : The :class:`BoozerMagneticField` instance used for field evaluation.
        s : A scalar or a numpy array of shape (npoints,) containing the
            normalized toroidal flux.
        theta : A scalar or a numpy array of shape (npoints,) containing the
            Boozer poloidal angle.
        zeta : A scalar or a numpy array of shape (npoints,) containing the
            Boozer toroidal angle.

    Returns:
        R : A scalar or a numpy array of shape (npoints,) containing the
            radial coordinate.
        phi : A scalar or a numpy array of shape (npoints,) containing the
            azimuthal angle.
        Z : A scalar or a numpy array of shape (npoints,) containing the
            vertical coordinate.
    """
    if not isinstance(s, np.ndarray):
        s = np.asarray(s)
    if not isinstance(theta, np.ndarray):
        theta = np.asarray(theta)
    if not isinstance(zeta, np.ndarray):
        zeta = np.asarray(zeta)

    # Handle scalar inputs - return scalars if any input is a scalar
    input_scalar = np.isscalar(s) or np.isscalar(theta) or np.isscalar(zeta)

    # Ensure all arrays have the same shape
    if s.shape != theta.shape or s.shape != zeta.shape:
        raise ValueError("s, theta, and zeta must have the same shape")

    npoints = s.size

    # Validate that arrays are not empty
    if npoints == 0:
        raise ValueError("Input arrays cannot be empty")

    points = np.zeros((npoints, 3))
    points[:, 0] = s.flatten()
    points[:, 1] = theta.flatten()
    points[:, 2] = zeta.flatten()

    field.set_points(points)

    R = field.R()[:, 0]
    Z = field.Z()[:, 0]
    nu = field.nu()[:, 0]
    phi = zeta - nu

    # Return scalars for scalar inputs, arrays for array inputs
    if input_scalar:
        return R[0], phi[0], Z[0]
    else:
        return R, phi, Z


def cylindrical_to_boozer(
    field,
    R,
    phi,
    Z,
    s_guess=None,
    theta_guess=None,
    zeta_guess=None,
    n_guesses=4,
    ftol=1e-6,
):
    r"""
    Convert from cylindrical coordinates to Boozer coordinates using root finding.

    Args:
        field : The :class:`BoozerMagneticField` instance used for field evaluation.
        R : A scalar or a numpy array of shape (npoints,) containing the
            radial coordinate.
        phi : A scalar or a numpy array of shape (npoints,) containing the
            azimuthal angle.
        Z : A scalar or a numpy array of shape (npoints,) containing the
            vertical coordinate.
        s_guess : float, optional
            Initial guess for s (default: None, uses 0.5). Must be a scalar.
        theta_guess : float, optional
            Initial guess for theta (default: None, uses 0.0). Must be a scalar.
        zeta_guess : float, optional
            Initial guess for zeta (default: None, uses phi). Must be a scalar.
        n_guesses : int, optional
            Number of initial guesses to try for each point (default: 4).
            Must be a positive integer.
        ftol : float, optional
            Tolerance for root finding convergence (default: 1e-6).

    Returns:
        s : A scalar or a numpy array of shape (npoints,) containing the
            normalized toroidal flux.
        theta : A scalar or a numpy array of shape (npoints,) containing the
            Boozer poloidal angle.
        zeta : A scalar or a numpy array of shape (npoints,) containing the
            Boozer toroidal angle.
    """
    if not isinstance(R, np.ndarray):
        R = np.asarray(R)
    if not isinstance(phi, np.ndarray):
        phi = np.asarray(phi)
    if not isinstance(Z, np.ndarray):
        Z = np.asarray(Z)

    # Handle scalar inputs - return scalars if any input is a scalar
    input_scalar = np.isscalar(R) or np.isscalar(phi) or np.isscalar(Z)

    # Ensure all arrays have the same shape
    if R.shape != phi.shape or R.shape != Z.shape:
        raise ValueError("R, phi, and Z must have the same shape")

    npoints = R.size

    # Validate that arrays are not empty
    if npoints == 0:
        raise ValueError("Input arrays cannot be empty")

    s = np.zeros(npoints)
    theta = np.zeros(npoints)
    zeta = np.zeros(npoints)

    # Set default guesses if not provided
    if s_guess is None:
        s_guess = 0.5
    if theta_guess is None:
        theta_guess = 0.0
    use_phi_for_zeta = zeta_guess is None

    # Validate that guesses are scalars
    if hasattr(s_guess, "__len__") and len(s_guess) > 1:
        raise ValueError("s_guess must be a scalar, not an array")
    if hasattr(theta_guess, "__len__") and len(theta_guess) > 1:
        raise ValueError("theta_guess must be a scalar, not an array")
    if hasattr(zeta_guess, "__len__") and len(zeta_guess) > 1:
        raise ValueError("zeta_guess must be a scalar, not an array")

    # Validate n_guesses parameter
    if not isinstance(n_guesses, int) or n_guesses <= 0:
        raise ValueError("n_guesses must be a positive integer")

    def objective_function(x, R_target, phi_target, Z_target):
        """Objective function for root finding."""
        s_val, theta_val, zeta_val = x

        # Ensure s is within bounds
        s_val = np.clip(s_val, 0.0, 1.0)

        points = np.zeros((1, 3))
        points[0, 0] = s_val
        points[0, 1] = theta_val
        points[0, 2] = zeta_val

        field.set_points(points)

        R_computed = field.R()[0, 0]
        Z_computed = field.Z()[0, 0]
        nu_computed = field.nu()[0, 0]
        phi_computed = zeta_val - nu_computed

        return [
            R_computed - R_target,
            np.arctan2(
                np.sin(phi_computed - phi_target), np.cos(phi_computed - phi_target)
            ),
            Z_computed - Z_target,
        ]

    for i in range(npoints):
        if use_phi_for_zeta:
            zeta_guess = phi[i]

        # Try multiple initial guesses for better convergence
        success = False
        initial_guesses = []

        # Add the user-provided guess as the first attempt
        initial_guesses.append([s_guess, theta_guess, zeta_guess])

        # Generate additional random guesses
        for _j in range(1, n_guesses):
            initial_guesses.append(
                [np.random.uniform(0, 1), np.random.uniform(0, 2 * np.pi), phi[i]]
            )

        for _attempt, x0 in enumerate(initial_guesses):
            sol = root(
                objective_function,
                x0,
                args=(R[i], phi[i], Z[i]),
                method="lm",
                options={"ftol": ftol},
            )
            if sol.success and np.all(np.abs(sol.fun) < ftol):
                s[i] = np.clip(sol.x[0], 0.0, 1.0)
                theta[i] = sol.x[1]
                zeta[i] = sol.x[2]
                success = True
                break

        if not success:
            raise RuntimeError(
                f"Root finding failed for point {i} with coordinates "
                f"R={R[i]}, phi={phi[i]}, Z={Z[i]}. "
                f"Tried {n_guesses} different initial guesses."
            )

    # Return scalars for scalar inputs, arrays for array inputs
    if input_scalar:
        return s[0], theta[0], zeta[0]
    else:
        return s, theta, zeta


def vmec_to_boozer(wout_filename, field, s_vmec, theta_vmec, phi_vmec, ftol=1e-6):
    r"""
    Convert from VMEC coordinates to Boozer coordinates.

    Args:
        wout_filename : str
            The name of the VMEC wout file.
        field : The :class:`BoozerMagneticField` instance used for field evaluation.
        s_vmec : A scalar or a numpy array of shape (npoints,) containing the
            normalized toroidal flux.
        theta_vmec : A scalar or a numpy array of shape (npoints,) containing
            the VMEC poloidal angle.
        phi_vmec : A scalar or a numpy array of shape (npoints,) containing
            the VMEC cylindrical angle.
        ftol : float, optional
            Tolerance for root finding convergence (default: 1e-6).

    Returns:
        theta_b : A numpy array of shape (npoints,) containing the Boozer
            poloidal angle.
        zeta_b : A numpy array of shape (npoints,) containing the Boozer toroidal angle.
    """
    # Validate that arrays are not empty
    if len(s_vmec) == 0:
        raise ValueError("Input arrays cannot be empty")

    # Load VMEC and booz_xform data
    f = netcdf_file(wout_filename, mmap=False)
    lmns = f.variables["lmns"][()]
    mnmax = f.variables["mnmax"][()]
    ns = f.variables["ns"][()]
    xm = f.variables["xm"][()]
    xn = f.variables["xn"][()]
    f.close()

    s_full_grid = np.linspace(0, 1, ns)
    s_half_grid = (s_full_grid[0:-1] + s_full_grid[1::]) / 2.0

    # Create splines for lmns
    lmns_splines = []
    for jmn in range(mnmax):
        lmns_splines.append(InterpolatedUnivariateSpline(s_half_grid, lmns[1::, jmn]))

    def vartheta_vmec(s, theta_vmec, phi_vmec):
        """Compute vartheta from VMEC data."""
        lmns = np.zeros((1, mnmax))
        for jmn in range(mnmax):
            lmns[:, jmn] = lmns_splines[jmn](s)

        angle = xm * theta_vmec - xn * phi_vmec
        sinangle = np.sin(angle)

        lambd = np.sum(lmns * sinangle)
        vartheta = theta_vmec + lambd
        return vartheta

    def vartheta_phi_vmec(s, theta_b, zeta_b):
        """Compute PEST angles from Boozer coordinates."""
        points = np.zeros((1, 3))
        points[:, 0] = s
        points[:, 1] = theta_b
        points[:, 2] = zeta_b
        field.set_points(points)
        nu = field.nu()[0, 0]
        iota = field.iota()[0, 0]
        vartheta = theta_b - iota * nu
        phi = zeta_b - nu
        return vartheta, phi

    def func_root(x, s, vartheta_target, phi_target):
        """Root finding function."""
        theta_b = x[0]
        zeta_b = x[1]
        vartheta, phi = vartheta_phi_vmec(s, theta_b, zeta_b)
        vartheta_diff = np.arctan2(
            np.sin(vartheta - vartheta_target), np.cos(vartheta - vartheta_target)
        )
        phi_diff = np.arctan2(np.sin(phi - phi_target), np.cos(phi - phi_target))
        return [vartheta_diff, phi_diff]

    theta_b = []
    zeta_b = []
    for i in range(len(s_vmec)):
        vartheta = vartheta_vmec(s_vmec[i], theta_vmec[i], phi_vmec[i])
        sol = root(
            func_root,
            [vartheta, phi_vmec[i]],
            args=(s_vmec[i], vartheta, phi_vmec[i]),
            method="lm",
            options={"ftol": ftol},
        )
        if sol.success and np.all(np.abs(sol.fun) < ftol):
            theta_b.append(sol.x[0])
            zeta_b.append(sol.x[1])
        else:
            raise RuntimeError(
                f"Root finding failed for point {i} with coordinates "
                f"s={s_vmec[i]}, theta_vmec={theta_vmec[i]}, phi_vmec={phi_vmec[i]}. "
            )

    return np.array(theta_b), np.array(zeta_b)


def boozer_to_vmec(wout_filename, field, s, theta_b, zeta_b, ftol=1e-6):
    r"""
    Convert from Boozer coordinates to VMEC coordinates.

    Args:
        wout_filename : str
            The name of the VMEC wout file.
        field : The :class:`BoozerMagneticField` instance used for field evaluation.
        s : A scalar or a numpy array of shape (npoints,) containing the
            normalized toroidal flux.
        theta_b : A scalar or a numpy array of shape (npoints,) containing
            the Boozer poloidal angle.
        zeta_b : A scalar or a numpy array of shape (npoints,) containing
            the Boozer toroidal angle.
        ftol : float, optional
            Tolerance for root finding convergence (default: 1e-6).

    Returns:
        theta_vmec : A scalar or a numpy array of shape (npoints,) containing
            the VMEC poloidal angle.
        phi_vmec : A scalar or a numpy array of shape (npoints,) containing
            the VMEC cylindrical angle.
    """
    # Handle scalar inputs - return scalars if any input is a scalar
    input_scalar = np.isscalar(s) or np.isscalar(theta_b) or np.isscalar(zeta_b)

    # Convert to arrays if needed
    s = np.asarray(s)
    theta_b = np.asarray(theta_b)
    zeta_b = np.asarray(zeta_b)

    # Ensure all inputs have the same shape
    if s.shape != theta_b.shape or s.shape != zeta_b.shape:
        raise ValueError("s, theta_b, and zeta_b must have the same shape")

    # Validate that arrays are not empty
    if s.size == 0:
        raise ValueError("Input arrays cannot be empty")

    # Load VMEC and booz_xform data
    f = netcdf_file(wout_filename, mmap=False)
    lmns = f.variables["lmns"][()]
    mnmax = f.variables["mnmax"][()]
    ns = f.variables["ns"][()]
    xm = f.variables["xm"][()]
    xn = f.variables["xn"][()]
    f.close()

    s_full_grid = np.linspace(0, 1, ns)
    s_half_grid = (s_full_grid[0:-1] + s_full_grid[1::]) / 2.0

    # Create splines for lmns
    lmns_splines = []
    for jmn in range(mnmax):
        lmns_splines.append(InterpolatedUnivariateSpline(s_half_grid, lmns[1::, jmn]))

    def vartheta_vmec(s, theta_vmec, phi_vmec):
        """Compute vartheta from VMEC data."""
        lmns = np.zeros((1, mnmax))
        for jmn in range(mnmax):
            lmns[:, jmn] = lmns_splines[jmn](s)

        angle = xm * theta_vmec - xn * phi_vmec
        sinangle = np.sin(angle)

        lambd = np.sum(lmns * sinangle)
        vartheta = theta_vmec + lambd
        return vartheta

    def vartheta_phi_vmec(s, theta_b, zeta_b):
        """Compute PEST angles from Boozer coordinates."""
        points = np.zeros((1, 3))
        points[:, 0] = s
        points[:, 1] = theta_b
        points[:, 2] = zeta_b
        field.set_points(points)
        nu = field.nu()[0, 0]
        iota = field.iota()[0, 0]
        vartheta = theta_b - iota * nu
        phi = zeta_b - nu
        return vartheta, phi

    def func_root(x, s, vartheta_boozer, zeta_boozer):
        """Root finding function."""
        theta_vmec = x[0]
        # Compute PEST angles from desired Boozer coordinates
        vartheta_target, phi_target = vartheta_phi_vmec(s, vartheta_boozer, zeta_boozer)
        # Compute PEST angles from VMEC coordinates
        vartheta = vartheta_vmec(s, theta_vmec, phi_target)
        vartheta_diff = np.arctan2(
            np.sin(vartheta - vartheta_target), np.cos(vartheta - vartheta_target)
        )
        return [vartheta_diff]

    theta_vmec = np.zeros_like(s)
    phi_vmec = np.zeros_like(s)
    for i in range(s.size):  # Use s.size instead of len(s)
        sol = root(
            func_root,
            [theta_b[i]],
            args=(s[i], theta_b[i], zeta_b[i]),
            method="lm",
            options={"ftol": ftol},
        )
        if sol.success and np.all(np.abs(sol.fun) < ftol):
            theta_vmec[i] = sol.x[0]
            vartheta, phi_vmec[i] = vartheta_phi_vmec(s[i], theta_b[i], zeta_b[i])
        else:
            raise RuntimeError(
                f"Root finding failed for point {i} with coordinates "
                f"s={s[i]}, theta_b={theta_b[i]}, zeta_b={zeta_b[i]}. "
            )

    # Return scalars for scalar inputs, arrays for array inputs
    if input_scalar:
        return theta_vmec[0], phi_vmec[0]
    else:
        return theta_vmec, phi_vmec


def vmec_to_cylindrical(wout_filename, s_vmec, theta_vmec, phi_vmec):
    r"""
    Convert from VMEC coordinates to cylindrical coordinates.

    Args:
        wout_filename : str
            The name of the VMEC wout file.
        s_vmec : A scalar or a numpy array of shape (npoints,) containing the
            normalized toroidal flux.
        theta_vmec : A scalar or a numpy array of shape (npoints,) containing
            the VMEC poloidal angle.
        phi_vmec : A scalar or a numpy array of shape (npoints,) containing
            the VMEC cylindrical angle.

    Returns:
        R : A scalar or a numpy array of shape (npoints,) containing the
            radial coordinate.
        phi_cyl : A scalar or a numpy array of shape (npoints,) containing
            the azimuthal angle.
        Z : A scalar or a numpy array of shape (npoints,) containing the
            vertical coordinate.
    """
    # Handle scalar inputs - return scalars if any input is a scalar
    input_scalar = (
        np.isscalar(s_vmec) or np.isscalar(theta_vmec) or np.isscalar(phi_vmec)
    )

    # Convert to arrays for processing
    s_vmec = np.asarray(s_vmec)
    theta_vmec = np.asarray(theta_vmec)
    phi_vmec = np.asarray(phi_vmec)

    # Validate that arrays are not empty
    if s_vmec.size == 0:
        raise ValueError("Input arrays cannot be empty")

    # Load VMEC data
    with netcdf_file(wout_filename, "r") as f:
        rmnc = f.variables["rmnc"][:]  # R harmonics (cos)
        zmns = f.variables["zmns"][:]  # Z harmonics (sin)
        xm = f.variables["xm"][:]  # poloidal mode numbers
        xn = f.variables["xn"][:]  # toroidal mode numbers
        ns = int(f.variables["ns"][()])  # number of radial surfaces (scalar)
        s_full = np.linspace(0, 1, ns)  # full radial grid

    # Initialize R and Z arrays
    R = np.zeros_like(s_vmec)
    Z = np.zeros_like(s_vmec)

    # For each point, compute R and Z using VMEC Fourier harmonics
    for i in range(s_vmec.size):
        s_i = s_vmec[i]
        theta_i = theta_vmec[i]
        phi_i = phi_vmec[i]

        # Interpolate harmonics to the desired s value
        rmnc_s = np.zeros_like(rmnc[0, :])
        zmns_s = np.zeros_like(zmns[0, :])

        for j in range(rmnc.shape[1]):  # Loop over mode numbers
            # Interpolate rmnc and zmns to s_i
            rmnc_s[j] = np.interp(s_i, s_full, rmnc[:, j])
            zmns_s[j] = np.interp(s_i, s_full, zmns[:, j])

        # Compute R and Z using Fourier series
        R[i] = 0.0
        Z[i] = 0.0

        for j in range(len(xm)):
            angle = xm[j] * theta_i - xn[j] * phi_i
            R[i] += rmnc_s[j] * np.cos(angle)
            Z[i] += zmns_s[j] * np.sin(angle)

    # phi_cyl is the same as phi_vmec
    phi_cyl = phi_vmec

    # Return scalars for scalar inputs, arrays for array inputs
    if input_scalar:
        return (
            R[0] if hasattr(R, "__len__") else R,
            phi_cyl[0] if hasattr(phi_cyl, "__len__") else phi_cyl,
            Z[0] if hasattr(Z, "__len__") else Z,
        )
    else:
        return R, phi_cyl, Z


def cylindrical_to_vmec(
    wout_filename,
    R,
    phi,
    Z,
    s_guess=None,
    theta_guess=None,
    zeta_guess=None,
    n_guesses=4,
    ftol=1e-6,
):
    r"""
    Convert from cylindrical coordinates to VMEC coordinates.

    Args:
        wout_filename : str
            The name of the VMEC wout file.
        R : A scalar or a numpy array of shape (npoints,) containing the
            radial coordinate.
        phi : A scalar or a numpy array of shape (npoints,) containing the
            azimuthal angle.
        Z : A scalar or a numpy array of shape (npoints,) containing the
            vertical coordinate.
        s_guess : float, optional
            Initial guess for s (default: None, uses 0.5). Must be a scalar.
        theta_guess : float, optional
            Initial guess for theta_vmec (default: None, uses arctan2(Z, R)).
            Must be a scalar.
        zeta_guess : float, optional
            Initial guess for phi_vmec (default: None, uses phi). Must be a scalar.
        n_guesses : int, optional
            Number of initial guesses to try for each point (default: 4).
            Must be a positive integer.
        ftol : float, optional
            Tolerance for root finding convergence (default: 1e-6).

    Returns:
        s_vmec : A scalar or a numpy array of shape (npoints,) containing the
            normalized toroidal flux.
        theta_vmec : A scalar or a numpy array of shape (npoints,) containing
            the VMEC poloidal angle.
        phi_vmec : A scalar or a numpy array of shape (npoints,) containing
            the VMEC cylindrical angle.
    """
    # Handle scalar inputs - return scalars if any input is a scalar
    input_scalar = np.isscalar(R) or np.isscalar(phi) or np.isscalar(Z)

    # Convert to arrays for processing
    R = np.asarray(R)
    phi = np.asarray(phi)
    Z = np.asarray(Z)

    # Validate that arrays are not empty
    if R.size == 0:
        raise ValueError("Input arrays cannot be empty")

    # Ensure all arrays have the same shape
    if R.shape != phi.shape or R.shape != Z.shape:
        raise ValueError("R, phi, and Z must have the same shape")

    # Load VMEC data
    with netcdf_file(wout_filename, "r") as f:
        rmnc = f.variables["rmnc"][:]  # R harmonics (cos)
        zmns = f.variables["zmns"][:]  # Z harmonics (sin)
        xm = f.variables["xm"][:]  # poloidal mode numbers
        xn = f.variables["xn"][:]  # toroidal mode numbers
        ns = int(f.variables["ns"][()])  # number of radial surfaces (scalar)
        s_full = np.linspace(0, 1, ns)  # full radial grid

    # Initialize output arrays
    s_vmec = np.zeros_like(R)
    theta_vmec = np.zeros_like(R)
    phi_vmec = np.zeros_like(R)

    # Set default guesses if not provided
    if s_guess is None:
        s_guess = 0.5
    if theta_guess is None:
        # Use arctan2(Z, R) as initial guess for theta
        theta_guess = 0.0  # Will be computed per point using arctan2
    use_phi_for_zeta = zeta_guess is None

    # Validate that guesses are scalars
    if hasattr(s_guess, "__len__") and len(s_guess) > 1:
        raise ValueError("s_guess must be a scalar, not an array")
    if hasattr(theta_guess, "__len__") and len(theta_guess) > 1:
        raise ValueError("theta_guess must be a scalar, not an array")
    if hasattr(zeta_guess, "__len__") and len(zeta_guess) > 1:
        raise ValueError("zeta_guess must be a scalar, not an array")

    # Validate n_guesses parameter
    if not isinstance(n_guesses, int) or n_guesses <= 0:
        raise ValueError("n_guesses must be a positive integer")

    def objective_function(x, R_target, phi_target, Z_target):
        """Objective function for root finding: R(x) - R_target, Z(x) - Z_target."""
        s_i, theta_i, phi_i = x

        # Interpolate harmonics to the desired s value
        rmnc_s = np.zeros_like(rmnc[0, :])
        zmns_s = np.zeros_like(zmns[0, :])

        for j in range(rmnc.shape[1]):  # Loop over mode numbers
            # Interpolate rmnc and zmns to s_i
            rmnc_s[j] = np.interp(np.clip(s_i, 0, 1), s_full, rmnc[:, j])
            zmns_s[j] = np.interp(np.clip(s_i, 0, 1), s_full, zmns[:, j])

        # Compute R and Z using Fourier series
        R_computed = 0.0
        Z_computed = 0.0

        for j in range(len(xm)):
            angle = xm[j] * theta_i - xn[j] * phi_i
            R_computed += rmnc_s[j] * np.cos(angle)
            Z_computed += zmns_s[j] * np.sin(angle)

        # Return residuals
        phi_diff = np.arctan2(np.sin(phi_i - phi_target), np.cos(phi_i - phi_target))
        return [R_computed - R_target, phi_diff, Z_computed - Z_target]

    # For each point, solve for VMEC coordinates using root finding
    for i in range(R.size):
        R_i = R[i]
        phi_i = phi[i]
        Z_i = Z[i]

        # Determine zeta_guess for this point
        zeta_guess_i = phi_i if use_phi_for_zeta else zeta_guess

        # Try multiple initial guesses for better convergence
        success = False
        initial_guesses = []

        # Add the user-provided guess as the first attempt
        initial_guesses.append([s_guess, theta_guess, zeta_guess_i])

        # Generate additional random guesses
        for _j in range(1, n_guesses):
            initial_guesses.append(
                [np.random.uniform(0, 1), np.random.uniform(0, 2 * np.pi), phi_i]
            )

        for _attempt, x0 in enumerate(initial_guesses):
            try:
                sol = root(
                    objective_function,
                    x0,
                    args=(R_i, phi_i, Z_i),
                    method="lm",
                    options={"ftol": ftol},
                )

                if sol.success and np.all(np.abs(sol.fun) < ftol):
                    s_vmec[i] = np.clip(sol.x[0], 0.0, 1.0)
                    theta_vmec[i] = sol.x[1]
                    phi_vmec[i] = sol.x[2]
                    success = True
                    break
            except Exception:
                continue

        if not success:
            raise RuntimeError(
                f"Root finding failed for point {i} with coordinates "
                f"R={R_i}, phi={phi_i}, Z={Z_i}. "
                f"Tried {n_guesses} different initial guesses."
            )

    # Return scalars for scalar inputs, arrays for array inputs
    if input_scalar:
        return s_vmec[0], theta_vmec[0], phi_vmec[0]
    else:
        return s_vmec, theta_vmec, phi_vmec
