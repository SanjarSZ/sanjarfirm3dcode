import time
import json
import numpy as np
import pathlib

from simsopt.field.boozermagneticfield import (
    BoozerRadialInterpolant,
    InterpolatedBoozerField,
)
from simsopt.field.trajectory_helpers import PassingPoincare
from simsopt.util.constants import (
    ALPHA_PARTICLE_CHARGE,
    ALPHA_PARTICLE_MASS,
    FUSION_ALPHA_PARTICLE_ENERGY,
)
from simsopt.util.functions import proc0_print, setup_logging
from simsopt.util.mpi import comm_size, comm_world, verbose

boozmn_filename = "../inputs/boozmn_aten_rescaled.nc"

charge = ALPHA_PARTICLE_CHARGE
mass = ALPHA_PARTICLE_MASS
Ekin = FUSION_ALPHA_PARTICLE_ENERGY

resolution = 48
sign_vpar = 1.0
lam = 0.0
ntheta_poinc = 1
ns_poinc = 120
Nmaps = 1000
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

poinc = PassingPoincare(
    field,
    lam,
    sign_vpar,
    mass,
    charge,
    Ekin,
    ns_poinc=ns_poinc,
    ntheta_poinc=ntheta_poinc,
    Nmaps=Nmaps,
    comm=comm_world,
    solver_options={"reltol": tol, "abstol": tol},
)

if verbose:
    poinc.plot_poincare()

time2 = time.time()

proc0_print("poincare time: ", time2 - time1)

def save_field_json(path):
    cfg = dict(
        boozmn_filename=boozmn_filename,
        order=order,
        no_K=no_K,
        degree=degree,
        ns_interp=ns_interp,
        ntheta_interp=ntheta_interp,
        nzeta_interp=nzeta_interp,
    )
    if comm_world.rank == 0:  # avoid MPI races
        p = pathlib.Path(path)
        p.write_text(json.dumps(cfg))
        proc0_print(f"Saved to: {p.resolve()}")
    comm_world.Barrier()

def load_field_json(path):
    cfg = json.loads(pathlib.Path(path).read_text())
    bri2 = BoozerRadialInterpolant(
        cfg["boozmn_filename"], cfg["order"], no_K=cfg["no_K"], comm=comm_world
    )
    return InterpolatedBoozerField(
        bri2, cfg["degree"],
        ns_interp=cfg["ns_interp"],
        ntheta_interp=cfg["ntheta_interp"],
        nzeta_interp=cfg["nzeta_interp"],
    )

def _eq(a, b, rtol=1e-12, atol=0.0):
    return np.allclose(np.asarray(a), np.asarray(b), rtol=rtol, atol=atol, equal_nan=True)

# 1) save
json_path = "QA_saved.json"
save_field_json(json_path)

# 2) load (rebuild) the field
field2 = load_field_json(json_path)

# 1000 random Boozer points
rng = np.random.default_rng(7)
s_pts     = rng.uniform(0.0, 1.0,    1000)
theta_pts = rng.uniform(0.0, 2*np.pi, 1000)
zeta_pts  = rng.uniform(0.0, 2*np.pi, 1000)


pts = np.column_stack((s_pts, theta_pts, zeta_pts))

field.set_points(pts)
field2.set_points(pts)

points_ok = np.allclose(field.get_points(), field2.get_points())


attrs_to_check = ("field","degree", "nfp", "stellsym", "extrapolate", "s_range", "theta_range", "zeta_range")

missing_on_field = []
missing_on_field2 = []
shared_attrs = []
attrs_ok = True

for a in attrs_to_check:
    v1 = getattr(field, a, None)
    v2 = getattr(field2, a, None)
    if v1 is None:
        missing_on_field.append(a)
    if v2 is None:
        missing_on_field2.append(a)
    if (v1 is not None) and (v2 is not None):
        shared_attrs.append(a)
        if not _eq(v1, v2):
            attrs_ok = False
            proc0_print(f"Attribute mismatch: {a}: {v1} vs {v2}")
            break

# If nothing could be compared, fail for safety
if not shared_attrs:
    attrs_ok = False
    proc0_print("No shared attributes to compare; failing attribute check for safety.")

proc0_print(
    f"Compared attributes: {shared_attrs} "
    f"(missing on field: {missing_on_field}; missing on field2: {missing_on_field2})"
)

for a in shared_attrs: 
    v1, v2 = getattr(field, a), getattr(field2, a) 
    same = (np.isclose(v1, v2) if isinstance(v1, float) else v1 == v2) 
    if not same: 
        attrs_ok = False 
        break

# final confirmation print
if attrs_ok and points_ok:
    proc0_print("All attributes (and set points) match: OK")
else:
    proc0_print(f"Mismatch detected — attrs_ok={attrs_ok}, points_ok={points_ok}")
