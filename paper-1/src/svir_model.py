"""SVIR compliance model and its partial-protection extension.

All epidemic parameters retain the manuscript's units. Efficacies es and ev
modify infection rates only; the payoff functions and thresholds are fixed.
"""
from dataclasses import dataclass, replace

import numpy as np
from scipy.optimize import brentq
from scipy.special import expit


@dataclass(frozen=True)
class Parameters:
    lam: float = 10.0
    pv: float = 0.897
    pb: float = 0.7934
    gamma: float = 0.143
    d: float = 0.007
    beta: float = 0.003
    ts: float = 56.0
    tv: float = 76.0
    Ts: float = 1.5
    Ss: float = -2.2
    Tv: float = 1.3
    Sv: float = -2.0
    es: float = 1.0
    ev: float = 1.0

    def __post_init__(self):
        assert min(self.lam, self.gamma, self.d, self.beta, self.ts) > 0
        assert 0 < self.pv < 1 and 0 < self.pb < 1
        assert self.ts < self.tv and 0 < self.es <= 1 and 0 < self.ev <= 1
        assert self.Ts > 0 and self.Tv > 0 and self.Ss < 0 and self.Sv < 0


def rhs_logit(t, z, p):
    S, V, I, R, zx, zy = z
    x, y, nx, ny = expit(zx), expit(zy), expit(-zx), expit(-zy)
    ks = p.beta * ((1 - p.es) + p.es * nx)
    kv = (1 - p.pb) * p.beta * ((1 - p.ev) + p.ev * ny)
    return np.array([
        (1 - p.pv) * p.lam - ks * S * I - p.d * S,
        p.pv * p.lam - kv * V * I - p.d * V,
        (ks * S + kv * V - p.gamma - p.d) * I,
        p.gamma * I - p.d * R,
        (p.Ts * x - p.Ss * nx) * (I - p.ts),
        (p.Tv * y - p.Sv * ny) * (I - p.tv),
    ])


def physical_state(z):
    u = np.array(z, dtype=float).copy()
    u[4:] = expit(u[4:])
    return u


def rhs_physical(u, p):
    S, V, I, R, x, y = u
    ks = p.beta * (1 - p.es * x)
    kv = (1 - p.pb) * p.beta * (1 - p.ev * y)
    gx, gy = p.Ts * x - p.Ss * (1 - x), p.Tv * y - p.Sv * (1 - y)
    return np.array([
        (1 - p.pv) * p.lam - ks * S * I - p.d * S,
        p.pv * p.lam - kv * V * I - p.d * V,
        (ks * S + kv * V - p.gamma - p.d) * I,
        p.gamma * I - p.d * R,
        x * (1 - x) * gx * (I - p.ts),
        y * (1 - y) * gy * (I - p.tv),
    ])


def jacobian(u, p):
    S, V, I, R, x, y = u
    ks = p.beta * (1 - p.es * x)
    kv = (1 - p.pb) * p.beta * (1 - p.ev * y)
    C, D = p.es * p.beta * S * I, p.ev * (1 - p.pb) * p.beta * V * I
    gx, gy = p.Ts * x - p.Ss * (1 - x), p.Tv * y - p.Sv * (1 - y)
    hx, hy = x * (1 - x) * gx, y * (1 - y) * gy
    dx = ((1 - 2 * x) * gx + x * (1 - x) * (p.Ts + p.Ss)) * (I - p.ts)
    dy = ((1 - 2 * y) * gy + y * (1 - y) * (p.Tv + p.Sv)) * (I - p.tv)
    return np.array([
        [-ks * I - p.d, 0, -ks * S, 0, C, 0],
        [0, -kv * I - p.d, -kv * V, 0, 0, D],
        [ks * I, kv * I, ks * S + kv * V - p.gamma - p.d, 0, -C, -D],
        [0, 0, p.gamma, -p.d, 0, 0],
        [0, 0, hx, 0, dx, 0],
        [0, 0, hy, 0, 0, dy],
    ])


def frozen_equilibrium(p, x, y):
    ks, kv = p.beta * (1 - p.es * x), (1 - p.pb) * p.beta * (1 - p.ev * y)

    def balance(I):
        return (ks * (1 - p.pv) * p.lam / (p.d + ks * I)
                + kv * p.pv * p.lam / (p.d + kv * I) - p.gamma - p.d)

    I = brentq(balance, 0, p.lam / (p.gamma + p.d), xtol=1e-13) if balance(0) > 0 else 0.0
    return np.array([(1 - p.pv) * p.lam / (p.d + ks * I),
                     p.pv * p.lam / (p.d + kv * I), I, p.gamma * I / p.d, x, y])


def pinned_equilibrium(p, active):
    k = p.gamma + p.d
    if active == 's':
        I = p.ts
        V = p.pv * p.lam / (p.d + (1 - p.pb) * p.beta * I)
        S = (p.lam - k * I) / p.d - V
        x = (1 - (k - (1 - p.pb) * p.beta * V) / (p.beta * S)) / p.es
        y = 0.0
    else:
        assert active == 'v'
        I = p.tv
        S = (1 - p.pv) * p.lam / (p.d + p.beta * (1 - p.es) * I)
        V = (p.lam - k * I) / p.d - S
        y = (1 - (k - p.beta * (1 - p.es) * S) / ((1 - p.pb) * p.beta * V)) / p.ev
        x = 1.0
    return np.array([S, V, I, p.gamma * I / p.d, x, y])


def predicted_equilibrium(p):
    e00 = frozen_equilibrium(p, 0, 0)
    e10 = frozen_equilibrium(p, 1, 0)
    e11 = frozen_equilibrium(p, 1, 1)
    if e00[2] == 0:
        return 'disease_free', e00
    if e00[2] < p.ts:
        return 'no_compliance', e00
    if e10[2] < p.ts:
        return 'susceptible_pinning', pinned_equilibrium(p, 's')
    if e10[2] < p.tv:
        return 'susceptible_compliance', e10
    if e11[2] < p.tv:
        return 'vaccinated_pinning', pinned_equilibrium(p, 'v')
    return 'both_compliant', e11


def self_check():
    p = Parameters(es=0.8, ev=0.6)
    z = np.array([103., 897., 10., 0., -0.4, 0.7])
    u = physical_state(z)
    transformed = rhs_logit(0, z, p)
    transformed[4:] *= u[4:] * (1 - u[4:])
    assert np.allclose(transformed, rhs_physical(u, p), rtol=1e-13, atol=1e-12)
    numeric = np.column_stack([(rhs_physical(u + np.eye(6)[j] * 1e-5, p)
                               - rhs_physical(u - np.eye(6)[j] * 1e-5, p)) / 2e-5
                              for j in range(6)])
    assert np.allclose(jacobian(u, p), numeric, rtol=2e-7, atol=2e-7)
    assert abs(frozen_equilibrium(Parameters(), 0, 0)[2] - 56.41995587716349) < 1e-10
    weak = Parameters(beta=0.0004, ts=2.0, tv=5.0)
    assert frozen_equilibrium(weak, 1, 0)[2] == 0
    assert abs(frozen_equilibrium(weak, 0, 0)[2] - 3.714559354) < 1e-6
    for p in [weak, Parameters(Ts=2.0, Ss=-2.0, Tv=3.0, Sv=-3.0)]:
        _, eq = predicted_equilibrium(p)
        assert np.max(np.abs(rhs_physical(eq, p))) < 1e-9
        assert np.max(np.linalg.eigvals(jacobian(eq, p)).real) < 0
    for efficacy in [0.25, 0.5, 0.75, 1.0]:
        for ts, tv in [(56, 76), (30, 40)]:
            p = Parameters(es=efficacy, ev=efficacy, ts=ts, tv=tv)
            name, eq = predicted_equilibrium(p)
            assert np.min(eq) >= -1e-10 and np.all(eq[4:] <= 1)
            assert np.max(np.abs(rhs_physical(eq, p))) < 1e-9
            assert np.max(np.linalg.eigvals(jacobian(eq, p)).real) < 0, name


if __name__ == '__main__':
    self_check()
    print('Logit chain rule, analytic Jacobian and reference equilibria PASS')
