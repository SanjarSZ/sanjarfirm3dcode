API Reference
============

This page provides detailed API documentation for the main FIRM3D modules and classes.

Magnetic Field Classes
---------------------

.. automodule:: simsopt.field.boozermagneticfield
   :members:
   :undoc-members:
   :show-inheritance:

.. automodule:: simsopt.field.coordinates
   :members:
   :undoc-members:
   :show-inheritance:

Trajectory Integration
---------------------

.. automodule:: simsopt.field.trajectory_helpers
   :members:
   :undoc-members:
   :show-inheritance:

Stopping Criteria
----------------

# Note: This module may not exist yet in the current codebase
# .. automodule:: simsopt.field.stopping_criteria
#    :members:
#    :undoc-members:
#    :show-inheritance:

Shear Alfvén Wave Classes
------------------------

# Note: This module may not exist yet in the current codebase
# .. automodule:: simsopt.field.shear_alfven_waves
#    :members:
#    :undoc-members:
#    :show-inheritance:

Utility Functions
----------------

.. automodule:: simsopt.util.constants
   :members:
   :undoc-members:
   :show-inheritance:

.. automodule:: simsopt.util.functions
   :members:
   :undoc-members:
   :show-inheritance:

Plotting Utilities
-----------------

.. automodule:: simsopt.plotting.plotting_helpers
   :members:
   :undoc-members:
   :show-inheritance:

Core Types and Utilities
-----------------------

.. automodule:: simsopt._core.types
   :members:
   :undoc-members:
   :show-inheritance:

.. automodule:: simsopt._core.util
   :members:
   :undoc-members:
   :show-inheritance:

SAW (Shear Alfvén Wave) Module
-----------------------------

.. automodule:: simsopt.saw.ae3d
   :members:
   :undoc-members:
   :show-inheritance:

.. automodule:: simsopt.saw.stellgap
   :members:
   :undoc-members:
   :show-inheritance:

Class Hierarchy
--------------

.. inheritance-diagram:: simsopt.field.boozermagneticfield
   :parts: 1

.. inheritance-diagram:: simsopt.field.trajectory_helpers
   :parts: 1

# .. inheritance-diagram:: simsopt.field.stopping_criteria
#    :parts: 1

# Function Index
# -------------
#
# .. autosummary::
#    :toctree: _autosummary
#    :template: function.rst
#    :recursive:
#
#    simsopt.field.trajectory_helpers.trace_particles_boozer
#    simsopt.field.trajectory_helpers.trace_particles_boozer_perturbed
#    simsopt.field.boozermagneticfield.BoozerAnalytic
#    simsopt.field.boozermagneticfield.BoozerRadialInterpolant
#    simsopt.field.boozermagneticfield.InterpolatedBoozerField
#
# Class Index
# ----------
#
# .. autosummary::
#    :toctree: _autosummary
#    :template: class.rst
#    :recursive:
#
#    simsopt.field.boozermagneticfield.BoozerMagneticField
#    simsopt.field.boozermagneticfield.BoozerAnalytic
#    simsopt.field.boozermagneticfield.BoozerRadialInterpolant
#    simsopt.field.boozermagneticfield.InterpolatedBoozerField
#    # simsopt.field.shear_alfven_waves.ShearAlfvenHarmonic
#    # simsopt.field.shear_alfven_waves.ShearAlfvenWavesSuperposition
#    # simsopt.field.stopping_criteria.StoppingCriterion
#    # simsopt.field.stopping_criteria.MaxToroidalFluxStoppingCriterion
#    # simsopt.field.stopping_criteria.MinToroidalFluxStoppingCriterion
#    # simsopt.field.stopping_criteria.ZetaStoppingCriterion
#    # simsopt.field.stopping_criteria.VparStoppingCriterion
#    # simsopt.field.stopping_criteria.ToroidalTransitStoppingCriterion
#    # simsopt.field.stopping_criteria.IterationStoppingCriterion
#    # simsopt.field.stopping_criteria.StepSizeStoppingCriterion
#
# Module Index
# -----------
#
# .. autosummary::
#    :toctree: _autosummary
#    :template: module.rst
#    :recursive:
#
#    simsopt.field
#    simsopt.field.boozermagneticfield
#    simsopt.field.coordinates
#    simsopt.field.trajectory_helpers
#    # simsopt.field.stopping_criteria
#    # simsopt.field.shear_alfven_waves
#    simsopt.util
#    simsopt.util.constants
#    simsopt.util.functions
#    simsopt.plotting
#    simsopt.plotting.plotting_helpers
#    simsopt._core
#    simsopt._core.types
#    simsopt._core.util
#    simsopt.saw
#    simsopt.saw.ae3d
#    simsopt.saw.stellgap
