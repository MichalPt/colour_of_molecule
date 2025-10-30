from contextlib import contextmanager
from threading import local
import numpy as np

# Global context for energy units
_energy_context = local()
_energy_context.current_units = "kcal/mol"

@contextmanager
def energy_units(units):
    """Context manager for setting energy units temporarily."""
    previous_units = getattr(_energy_context, 'current_units', "kcal/mol")
    _energy_context.current_units = units
    try:
        yield
    finally:
        _energy_context.current_units = previous_units

def get_current_energy_units():
    """Get the current energy units from context."""
    return getattr(_energy_context, 'current_units', "kcal/mol")

@np.vectorize
def convert_nm_to_rcm(x):
    if x == 0:
        return 0
    else:
        return 1e7 / x

