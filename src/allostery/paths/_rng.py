"""Numerical Recipes ``ran3`` random number generator.

The original path-sampling programs read/write the ``ran3`` state from a file
called ``ran3.dat``.  Here the state is an ``int32`` array of length 57:
``state[0:55]`` is ``ma(1:55)``, ``state[55]`` is ``inext`` and
``state[56]`` is ``inextp``.  The array is shared (by reference) with the
Fortran kernel, which advances it in place.
"""

from __future__ import annotations

import os

import numpy as np

MBIG = 1000000000
MSEED = 161803398
MZ = 0
FAC = 1.0 / MBIG

STATE_SIZE = 57


def seed_state(seed: int) -> np.ndarray:
    """Return a freshly initialized ``ran3`` state for an integer ``seed``.

    This is the initialization branch of ``ran3`` in Numerical Recipes
    (``idum < 0``), which the original programs did not include.
    """
    ma = [0] * 56  # 1-based
    mj = abs(MSEED - abs(int(seed)))
    mj %= MBIG
    ma[55] = mj
    mk = 1
    for i in range(1, 55):
        ii = (21 * i) % 55
        ma[ii] = mk
        mk = mj - mk
        if mk < MZ:
            mk += MBIG
        mj = ma[ii]
    for _ in range(4):
        for i in range(1, 56):
            ma[i] -= ma[1 + (i + 30) % 55]
            if ma[i] < MZ:
                ma[i] += MBIG
    state = np.zeros(STATE_SIZE, dtype=np.int32)
    state[:55] = ma[1:]
    state[55] = 0
    state[56] = 31
    return state


def random_seed() -> int:
    """A random positive seed (used when the user does not provide one)."""
    return int.from_bytes(os.urandom(4), "little") % (MSEED - 1) + 1


def as_state(seed_or_state) -> np.ndarray:
    """Coerce a seed (int), ``None`` (random) or an existing state array."""
    if seed_or_state is None:
        return seed_state(random_seed())
    if isinstance(seed_or_state, (int, np.integer)):
        return seed_state(int(seed_or_state))
    state = np.ascontiguousarray(seed_or_state, dtype=np.int32).copy()
    if state.shape != (STATE_SIZE,):
        raise ValueError(f"ran3 state must have shape ({STATE_SIZE},)")
    return state


def ran3(state: np.ndarray) -> float:
    """Advance ``state`` in place and return a uniform deviate in (0, 1)."""
    inext = int(state[55]) + 1
    if inext == 56:
        inext = 1
    inextp = int(state[56]) + 1
    if inextp == 56:
        inextp = 1
    mj = int(state[inext - 1]) - int(state[inextp - 1])
    if mj <= MZ:
        mj += MBIG - 1
    state[inext - 1] = mj
    state[55] = inext
    state[56] = inextp
    return mj * FAC


def read_ran3_file(filename) -> np.ndarray:
    """Read a legacy ``ran3.dat`` file (``inext inextp`` then 55 integers)."""
    with open(filename) as f:
        values = [int(v) for v in f.read().split()]
    state = np.zeros(STATE_SIZE, dtype=np.int32)
    state[55], state[56] = values[0], values[1]
    state[:55] = values[2:57]
    return state


def write_ran3_file(filename, state: np.ndarray) -> None:
    """Write a ``ran3`` state in the legacy ``ran3.dat`` format."""
    with open(filename, "w") as f:
        f.write(f"{int(state[55])} {int(state[56])}\n")
        for v in state[:55]:
            f.write(f"{int(v)}\n")
