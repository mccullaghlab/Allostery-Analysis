"""Pure-Python implementation of the path Monte Carlo kernel.

This is a line-by-line port of ``_twisp.f90`` (itself a consolidation of the
original ``paths_TWISP*.f`` programs).  It is slow -- use it only for small
systems, for testing, or when no Fortran compiler is available.  For the same
random-number state it produces exactly the same paths as the Fortran kernel.

IEEE semantics (overflow to Inf, x/0, NaN comparisons) are reproduced
explicitly, since Python raises exceptions where Fortran silently returns
Inf/NaN.
"""

from __future__ import annotations

import math

import numpy as np

from ._rng import ran3

INF = math.inf
NAN = math.nan


def _exp(x: float) -> float:
    try:
        return math.exp(x)
    except OverflowError:
        return INF


def _div(a: float, b: float) -> float:
    if b == 0.0:
        if a == 0.0 or a != a:
            return NAN
        return math.copysign(INF, a) * math.copysign(1.0, b)
    return a / b


def _pow(a: float, e: float) -> float:
    if a == 0.0 and e < 0:
        return INF
    return a ** e


def _neglog_abs(x: float) -> float:
    x = abs(x)
    if x == 0.0:
        return INF
    if x != x:
        return NAN
    return -math.log(x)


def cost_matrix(adjacency: np.ndarray) -> np.ndarray:
    """Edge costs ``C[i, j] = -log(A[i, j])``, ``-1`` where ``A[i, j] < 1e-10``.

    Computed element by element with :func:`math.log` and handed to *both*
    backends, so that the Fortran kernel and this implementation use
    bit-identical costs (a compiler may otherwise replace ``log`` by a
    vectorized version that differs in the last bit).
    """
    a = np.ascontiguousarray(adjacency, dtype=np.float64)
    vals = [-1.0 if v < 1e-10 else -math.log(v) for v in a.ravel().tolist()]
    return np.array(vals, dtype=np.float64).reshape(a.shape)


def dijkstra(cost_c: np.ndarray, source: int, sink: int):
    """Port of the O(n^2) Dijkstra routine; returns (1-based path, length).

    ``cost_c`` is the matrix returned by :func:`cost_matrix`.
    """
    n = cost_c.shape[0]
    cost = cost_c.T.tolist()  # cost[a-1][b-1] == Fortran cost(a, b)
    isrc, isnk = source + 1, sink + 1
    q = list(range(1, n + 1))
    dist = [9e99] * (n + 1)
    prev = [0] * (n + 1)
    dist[isrc] = 0.0
    mpar = n
    mdis = 0.0
    while mpar != 0:
        k = 1
        m = q[0]
        mdis = dist[m]
        for i in range(2, mpar + 1):
            j = q[i - 1]
            if dist[j] < mdis:
                mdis = dist[j]
                k = i
                m = j
        if m == isnk:
            mpar = 0
        else:
            q[k - 1] = q[mpar - 1]
            mpar -= 1
            cm = cost[m - 1]
            for i in range(1, n + 1):
                if i != m and cm[i - 1] > 0.0:
                    alt = mdis + cm[i - 1]
                    if alt < dist[i]:
                        dist[i] = alt
                        prev[i] = m
    path = []
    m = isnk
    while m != 0:
        path.append(m)
        m = prev[m]
    path.reverse()
    return np.array(path, dtype=np.int32), mdis


def _accept(prob: float, state) -> int:
    if prob < 1.0:
        return 1 if ran3(state) < prob else 0
    return 1


def _sign(x: int) -> int:
    return 1 if x > 0 else -1


def sample(cost_c, positions, angular, source, sink, temperature, nsweeps,
           path, npath, state, nodes_buf, lengths_buf, weights_buf):
    """Same interface and semantics as :func:`._fortran.sample`."""
    n = cost_c.shape[0]
    cost_t = cost_c.T.tolist()  # cost_t[a-1][b-1] == Fortran cost(a, b)

    def cost(a, b):  # Fortran cost(a, b), 1-based
        return cost_t[a - 1][b - 1]

    T = float(temperature)
    bet = 1.0 / T
    is3d = bool(angular)

    if is3d:  # see `scaled` in _twisp.f90
        def scaled(v):
            return v * bet
    else:
        def scaled(v):
            return v / T
    isrc, isnk = source + 1, sink + 1
    cap = nodes_buf.size
    status = 0

    rhat = None
    if is3d:
        x = np.asarray(positions, dtype=np.float64).tolist()
        rhat = {}
        for i in range(n):
            for j in range(n):
                if i != j:
                    d = [x[j][k] - x[i][k] for k in range(3)]
                    nrm = (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]) ** 0.5
                    rhat[(j + 1, i + 1)] = (d[0] / nrm, d[1] / nrm, d[2] / nrm)
                else:
                    rhat[(j + 1, i + 1)] = (0.0, 0.0, 0.0)

    def bcos(a, b, c):
        u = rhat[(a, b)]
        v = rhat[(b, c)]
        return u[0] * v[0] + u[1] * v[1] + u[2] * v[2]

    p = [0] + [int(v) for v in path[:npath]] + [0] * (n - npath + 1)  # 1-based
    ipath = [0] * (n + 1)
    for i in range(1, npath + 1):
        j = abs(p[i])
        ipath[j] = i * _sign(p[i])

    ndone = 0
    nnodes = 0
    ntry = [0] * (2 * n + 3)
    wtry = [0.0] * (2 * n + 3)
    for _ in range(int(nsweeps)):
        if nnodes + n > cap:
            break
        for _mgo in range(n):
            nn = int(ran3(state) * float(n)) + 1
            while nn == isrc or nn == isnk:
                nn = int(ran3(state) * float(n)) + 1
            ncon = 0
            i1 = abs(p[1])
            j1 = _sign(p[1])
            jgo = 1
            wmax = -1e12
            for i in range(2, npath + 1):
                i2 = abs(p[i])
                j2 = _sign(p[i])
                if j1 == j2:
                    if i2 != nn:
                        if cost(nn, i1) > -0.1 and cost(nn, i2) > -0.1:
                            ncon += 1
                            ntry[ncon] = j1 * i
                            wtry[ncon] = cost(i1, i2) - cost(nn, i1) - cost(nn, i2)
                            if wtry[ncon] > wmax:
                                wmax = wtry[ncon]
                        i1 = i2
                        j1 = j2
                    else:
                        i2 = abs(p[i + 1])
                        j2 = _sign(p[i + 1])
                        if j1 == j2 and cost(i1, i2) < -0.1:
                            jgo = 0
                else:
                    if i2 != nn:
                        if cost(nn, i1) > -0.1:
                            ncon += 1
                            ntry[ncon] = j1 * i
                            wtry[ncon] = -cost(nn, i1)
                            if wtry[ncon] > wmax:
                                wmax = wtry[ncon]
                        if cost(nn, i2) > -0.1:
                            ncon += 1
                            ntry[ncon] = j2 * i
                            wtry[ncon] = -cost(nn, i2)
                            if wtry[ncon] > wmax:
                                wmax = wtry[ncon]
                        i1 = i2
                        j1 = j2
            wtot = 0.0
            for i in range(1, ncon + 1):
                wtot = wtot + _exp(scaled(wtry[i] - wmax))
                wtry[i] = wtot

            if ipath[nn] == 0:
                if ncon != 0:
                    prob = _exp(scaled(wmax)) * wtot
                    igo = _accept(prob, state)
                    if igo == 1:
                        tr = ran3(state) * wtot
                        nswp = 1
                        while tr > wtry[nswp]:
                            nswp += 1
                            if nswp > ncon:
                                break
                        if nswp > ncon:
                            status = 1
                            nswp = ncon
                        nswp = ntry[nswp]
                        mswp = _sign(nswp)
                        nswp = abs(nswp)
                        if is3d:
                            i2 = abs(p[nswp - 1]); j2 = _sign(p[nswp - 1])
                            i3 = abs(p[nswp]); j3 = _sign(p[nswp])
                            prob = bcos(i2, nn, i3) if j2 == j3 else 1.0
                            if nswp != 2:
                                i1 = abs(p[nswp - 2]); j1 = _sign(p[nswp - 2])
                                if j1 == mswp:
                                    prob = prob * bcos(i1, i2, nn)
                                if j1 == j3:
                                    prob = _div(prob, bcos(i1, i2, i3))
                            if nswp != npath:
                                i4 = abs(p[nswp + 1]); j4 = _sign(p[nswp + 1])
                                if mswp == j4:
                                    prob = prob * bcos(nn, i3, i4)
                                if j2 == j4:
                                    prob = _div(prob, bcos(i2, i3, i4))
                            prob = _pow(abs(prob), bet)
                            igo = _accept(prob, state)
                        if igo == 1:
                            for i in range(npath, nswp - 1, -1):
                                j = abs(p[i])
                                k = _sign(p[i])
                                ipath[j] = k * (i + 1)
                                p[i + 1] = k * j
                            p[nswp] = mswp * nn
                            ipath[nn] = mswp * nswp
                            npath += 1
            elif jgo == 1:
                prob = _div(_exp(scaled(-wmax)), wtot)
                igo = _accept(prob, state)
                if igo == 1:
                    nswp = abs(ipath[nn])
                    mswp = _sign(ipath[nn])
                    if is3d:
                        i2 = abs(p[nswp - 1]); j2 = _sign(p[nswp - 1])
                        i3 = abs(p[nswp + 1]); j3 = _sign(p[nswp + 1])
                        prob = bcos(i2, nn, i3) if j2 == j3 else 1.0
                        if nswp != 2:
                            i1 = abs(p[nswp - 2]); j1 = _sign(p[nswp - 2])
                            if j1 == mswp:
                                prob = prob * bcos(i1, i2, nn)
                            if j1 == j3:
                                prob = _div(prob, bcos(i1, i2, i3))
                        if nswp != npath - 1:
                            i4 = abs(p[nswp + 2]); j4 = _sign(p[nswp + 2])
                            if mswp == j4:
                                prob = prob * bcos(nn, i3, i4)
                            if j2 == j4:
                                prob = _div(prob, bcos(i2, i3, i4))
                        prob = _pow(abs(prob), -bet)
                        igo = _accept(prob, state)
                    if igo == 1:
                        ipath[nn] = 0
                        for i in range(nswp + 1, npath + 1):
                            j = abs(p[i])
                            k = _sign(p[i])
                            p[i - 1] = k * j
                            ipath[j] = k * (i - 1)
                        npath -= 1

        # path weight
        i1 = abs(p[1]); j1 = _sign(p[1])
        wgt = 0.0
        for i in range(2, npath + 1):
            i2 = abs(p[i]); j2 = _sign(p[i])
            if j1 == j2:
                wgt = wgt + cost(i1, i2)
            i1, j1 = i2, j2
        wgt2 = 0.0
        if is3d and npath >= 3:
            i1 = abs(p[1]); j1 = _sign(p[1])
            i2 = abs(p[2]); j2 = _sign(p[2])
            for i in range(3, npath + 1):
                i3 = abs(p[i]); j3 = _sign(p[i])
                if j1 == j3:
                    wgt2 = wgt2 + _neglog_abs(bcos(i1, i2, i3))
                i1, j1, i2, j2 = i2, j2, i3, j3
        weights_buf[ndone, 0] = wgt + wgt2
        weights_buf[ndone, 1] = wgt
        weights_buf[ndone, 2] = wgt2
        lengths_buf[ndone] = npath
        nodes_buf[nnodes:nnodes + npath] = p[1:npath + 1]
        nnodes += npath
        ndone += 1

    path[:npath] = p[1:npath + 1]
    return npath, ndone, nnodes, status
