import builtins
import sys

from .mpi import verbose


def print(*args, **kwargs):
    r"""
    Overloaded print function to force flushing of stdout.
    All procs prints to stdout.
    """
    builtins.print(*args, **kwargs, flush=True, file=sys.stdout)


def proc0_print(*args, **kwargs):
    r"""
    Overloaded print function to force flushing of stdout.
    Only proc0 prints to stdout.
    """
    if verbose:
        builtins.print(*args, **kwargs, flush=True, file=sys.stdout)
