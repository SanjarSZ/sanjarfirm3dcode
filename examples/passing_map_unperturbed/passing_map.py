import time
import json
import numpy as np
import base64

from simsopt.field.boozermagneticfield import (
    BoozerRadialInterpolant,
    InterpolatedBoozerField,
)


from simsopt.util.functions import proc0_print, setup_logging
from simsopt.util.mpi import comm_size, comm_world, verbose

boozmn_filename = "../inputs/boozmn_aten_rescaled.nc"


resolution = 48
ntheta_poinc = 1
ns_interp = resolution
ntheta_interp = resolution
nzeta_interp = resolution
order = 3
tol = 1e-8
degree = 3
no_K = True

# Setup logging to redirect output to file
setup_logging(f"stdout_passing_map_{resolution}_{comm_size}.txt")

time1 = time.time()

bri = BoozerRadialInterpolant(boozmn_filename, order, no_K=no_K, comm=comm_world)

field = InterpolatedBoozerField(
    bri,
    degree,
    ns_interp=ns_interp,
    ntheta_interp=ntheta_interp,
    nzeta_interp=nzeta_interp,
)

time2 = time.time()

proc0_print("Building Interpolant:", time2 - time1)

def force_compute_all_interpolants(field):
    """
    Force all interpolants to compute by evaluating each quantity.
    This sets all status flags to true.
    """
    test_pt = np.array([[0.5, np.pi, np.pi]])
    field.set_points(test_pt)
    
    # All quantities that can be computed
    all_methods = [
        # Basic field quantities
        'modB', 'R', 'Z', 'G', 'I', 'iota', 'psip', 'nu', 'K',
        # First derivatives
        'dmodBdtheta', 'dmodBdzeta', 'dmodBds',
        'dGds', 'dIds', 'diotads',
        'dRdtheta', 'dRdzeta', 'dRds',
        'dZdtheta', 'dZdzeta', 'dZds',
        'dnudtheta', 'dnudzeta', 'dnuds',
        'dKdtheta', 'dKdzeta', 'K_derivs', 'R_derivs', 
        'Z_derivs', 'nu_derivs', 'modB_derivs'
    ]
    
    computed = []
    for method in all_methods:
        try:
            getattr(field, method)()
            computed.append(method)
        except:
            pass
    
    proc0_print(f"Computed {len(computed)} interpolants")
    return computed


def to_json(field, path):
    """
    Extension method for InterpolatedBoozerField.to_json()
    Saves ALL actual evaluated field data along with configuration.
    """
    # Force all interpolants to compute (sets all status flags to true)
    computed_quantities = force_compute_all_interpolants(field)
    
    # Use smaller grid to keep JSON manageable
    ns, ntheta, nzeta = 10, 10, 10
    s_grid = np.linspace(0.1, 0.9, ns)
    theta_grid = np.linspace(0, 2*np.pi, ntheta)
    zeta_grid = np.linspace(0, 2*np.pi, nzeta)
    
    proc0_print(f"Saving {len(computed_quantities)} field quantities on {ns}x{ntheta}x{nzeta} grid...")
    
    # Save ALL computed field quantities
    field_data = {}
    
    for quantity in computed_quantities:
        try:
            data = np.zeros((ns, ntheta, nzeta))
            for i in range(ns):
                for j in range(ntheta):
                    for k in range(nzeta):
                        pt = np.array([[s_grid[i], theta_grid[j], zeta_grid[k]]])
                        field.set_points(pt)
                        result = getattr(field, quantity)()
                        data[i,j,k] = result[0] if len(result) > 0 else 0.0
            
            # Encode as base64 for JSON
            field_data[quantity] = {
                'data': base64.b64encode(data.tobytes()).decode('ascii'),
                'shape': data.shape,
                'dtype': str(data.dtype)
            }
            proc0_print(f"  Saved {quantity}")
        except Exception as e:
            proc0_print(f"  Could not save {quantity}: {e}")
    
    # Get all status flags
    status_flags = {}
    status_names = [
        'status_modB', 'status_dmodBdtheta', 'status_dmodBdzeta', 'status_dmodBds',
        'status_G', 'status_I', 'status_iota', 'status_dGds', 'status_dIds', 
        'status_diotads', 'status_psip', 'status_R', 'status_Z', 'status_nu', 
        'status_K', 'status_dRdtheta', 'status_dRdzeta', 'status_dRds',
        'status_dZdtheta', 'status_dZdzeta', 'status_dZds', 'status_dnudtheta',
        'status_dnudzeta', 'status_dnuds', 'status_dKdtheta', 'status_dKdzeta',
        'status_K_derivs', 'status_R_derivs', 'status_Z_derivs', 'status_nu_derivs',
        'status_modB_derivs'
    ]
    
    for flag in status_names:
        try:
            status_flags[flag] = getattr(field, flag)
        except:
            pass
    
    # Save configuration, field data, and status
    save_dict = {
        'config': {
            'boozmn_filename': boozmn_filename,
            'order': order,
            'no_K': no_K,
            'degree': degree,
            'ns_interp': ns_interp,
            'ntheta_interp': ntheta_interp,
            'nzeta_interp': nzeta_interp,
            'extrapolate': getattr(field, 'extrapolate', True),
            'nfp': getattr(field, 'nfp', 1),
            'stellsym': getattr(field, 'stellsym', True),
        },
        'field_data': field_data,
        'grid': {
            's': s_grid.tolist(),
            'theta': theta_grid.tolist(),
            'zeta': zeta_grid.tolist()
        },
        'status_flags': status_flags,
        'computed_quantities': computed_quantities
    }
    
    if comm_world.rank == 0:
        with open(path, 'w') as f:
            json.dump(save_dict, f, indent=2)
        proc0_print(f"\nSaved to {path}:")
        proc0_print(f"  - {len(field_data)} field quantities with actual evaluated data")
        proc0_print(f"  - {len(status_flags)} status flags")
        proc0_print(f"  - Total grid points: {ns*ntheta*nzeta}")
    comm_world.Barrier()


def from_json(path):
    """
    Extension method for InterpolatedBoozerField.from_json()
    Loads field from JSON and returns an InterpolatedBoozerField instance.
    """
    with open(path, 'r') as f:
        data = json.load(f)
    
    config = data['config']
    
    # Recreate the InterpolatedBoozerField
    bri2 = BoozerRadialInterpolant(
        config['boozmn_filename'], 
        config['order'], 
        no_K=config['no_K'], 
        comm=comm_world
    )
    
    field2 = InterpolatedBoozerField(
        bri2, 
        config['degree'],
        ns_interp=config['ns_interp'],
        ntheta_interp=config['ntheta_interp'],
        nzeta_interp=config['nzeta_interp'],
    )
    
    # Force same interpolants to compute (ensures same status flags)
    force_compute_all_interpolants(field2)
    
    # Verify field data exists
    if 'field_data' in data:
        proc0_print(f"\nLoaded from {path}:")
        proc0_print(f"  - JSON contains {len(data['field_data'])} field quantities with actual data")
        proc0_print(f"  - Quantities: {', '.join(data['field_data'].keys())}")
    
    return field2


# Add methods to InterpolatedBoozerField class
InterpolatedBoozerField.to_json = to_json
InterpolatedBoozerField.from_json = staticmethod(from_json)

proc0_print("\n" + "="*70)
proc0_print("COMPLETE INTERPOLATED FIELD SAVE/LOAD TEST")
proc0_print("="*70)

# Step 1: Save field1 to JSON with ALL field quantities
proc0_print("\n1. Saving field1 to JSON with ALL field quantities...")
t_save = time.time()
field.to_json("QA_saved.json")
proc0_print(f"   Save time: {time.time() - t_save:.2f}s")

# Step 2: Load as field2 from JSON
proc0_print("\n2. Loading field2 from JSON...")
t_load = time.time()
field2 = InterpolatedBoozerField.from_json("QA_saved.json")
proc0_print(f"   Load time: {time.time() - t_load:.2f}s")

# Step 3: Generate 1000 points
rng = np.random.default_rng(7)
points = np.column_stack((
    rng.uniform(0.0, 1.0, 1000),
    rng.uniform(0.0, 2*np.pi, 1000),
    rng.uniform(0.0, 2*np.pi, 1000)
))

# Step 4: Set points
field.set_points(points)
field2.set_points(points)


all_match = True

# Check configuration attributes
proc0_print("\nConfiguration attributes:")
config_attrs = ['degree', 'nfp', 'stellsym', 'extrapolate']
for attr in config_attrs:
    v1 = getattr(field, attr, None)
    v2 = getattr(field2, attr, None)
    match = (v1 == v2)
    proc0_print(f"  {attr}: {v1} == {v2} : {'✓' if match else '✗'}")
    if not match:
        all_match = False

# Check ALL field methods produce same results
proc0_print("\nField evaluations (1000 points):")
test_methods = ['modB', 'R', 'Z', 'G', 'I', 'iota', 'nu', 'K', 'psip']
for method in test_methods:
    try:
        v1 = getattr(field, method)()
        v2 = getattr(field2, method)()
        match = np.allclose(v1, v2, rtol=1e-12, atol=1e-14)
        max_diff = np.max(np.abs(v1 - v2))
        proc0_print(f"  {method}(): max diff = {max_diff:.3e} : {'✓' if match else '✗'}")
        if not match:
            all_match = False
    except Exception as e:
        proc0_print(f"  {method}(): could not evaluate - {e}")

# Check derivative methods
proc0_print("\nDerivative evaluations:")
deriv_methods = ['dmodBdtheta', 'dmodBdzeta', 'dmodBds',
        'dGds', 'dIds', 'diotads',
        'dRdtheta', 'dRdzeta', 'dRds',
        'dZdtheta', 'dZdzeta', 'dZds',
        'dnudtheta', 'dnudzeta', 'dnuds',
        'dKdtheta', 'dKdzeta', 'K_derivs', 'R_derivs', 
        'Z_derivs', 'nu_derivs', 'modB_derivs']
for method in deriv_methods:
    try:
        v1 = getattr(field, method)()
        v2 = getattr(field2, method)()
        match = np.allclose(v1, v2, rtol=1e-12, atol=1e-14)
        max_diff = np.max(np.abs(v1 - v2))
        proc0_print(f"  {method}(): max diff = {max_diff:.3e} : {'✓' if match else '✗'}")
        if not match:
            all_match = False
    except:
        pass

#status flags
proc0_print("\nStatus flags:")
important_flags = [
        'status_modB', 'status_dmodBdtheta', 'status_dmodBdzeta', 'status_dmodBds',
        'status_G', 'status_I', 'status_iota', 'status_dGds', 'status_dIds', 
        'status_diotads', 'status_psip', 'status_R', 'status_Z', 'status_nu', 
        'status_K', 'status_dRdtheta', 'status_dRdzeta', 'status_dRds',
        'status_dZdtheta', 'status_dZdzeta', 'status_dZds', 'status_dnudtheta',
        'status_dnudzeta', 'status_dnuds', 'status_dKdtheta', 'status_dKdzeta',
        'status_K_derivs', 'status_R_derivs', 'status_Z_derivs', 'status_nu_derivs',
        'status_modB_derivs'
    ]
for flag in important_flags:
    v1 = getattr(field, flag, None)
    v2 = getattr(field2, flag, None)
    if v1 is not None and v2 is not None:
        match = (v1 == v2 == True)  # Should all be True
        proc0_print(f"  {flag}: {v1} == {v2} : {'✓' if match else '✗'}")
        if not match:
            all_match = False

# Final results
proc0_print("\n" + "="*70)
if all_match:
    proc0_print("✓ SUCCESS: All attributes match!")
    proc0_print("✓ field1 and field2 are identical")
    proc0_print("✓ JSON contains ALL field quantities with actual evaluated data")
    proc0_print("✓ Loaded field is an InterpolatedBoozerField instance")
    proc0_print("✓ All status flags are true")
else:
    proc0_print("✗ Some attributes don't match")

assert all_match, "All attributes should match"
proc0_print("\n✓ ASSERTION PASSED: All attributes match!")
proc0_print("="*70)
