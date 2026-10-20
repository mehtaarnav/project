"""Optional numba acceleration. `maybe_njit` is numba.njit(cache=True) when
numba is importable, otherwise the identity decorator."""
try:
    import numba as _nb
    HAVE_NUMBA = True

    def maybe_njit(fn):
        return _nb.njit(cache=True)(fn)
except ImportError:  # pragma: no cover
    HAVE_NUMBA = False

    def maybe_njit(fn):
        return fn
