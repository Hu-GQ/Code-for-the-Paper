"""Exact characteristic-polynomial and Hurwitz identities, then numeric checks."""
from pathlib import Path
from dataclasses import replace
import json

import numpy as np
import sympy as s

from svir_model import Parameters, frozen_equilibrium, pinned_equilibrium, rhs_physical, jacobian, self_check

ROOT = Path(__file__).resolve().parents[1]


def main():
    self_check()
    L, B, F, G, C, E, d = s.symbols('L B F G C E d', positive=True)
    z = s.Symbol('z')
    Q = C * E
    a, b, U, V = L + d, B + d, L * F, B * G
    J = s.Matrix([[-a, 0, -F, C], [0, -b, -G, 0],
                  [L, B, 0, -C], [0, 0, E, 0]])
    c1 = a + b
    c2 = a * b + U + V + Q
    c3 = b * U + a * V + (b + d) * Q
    c4 = b * d * Q
    assert s.expand(J.charpoly(z).as_expr() - (z**4+c1*z**3+c2*z**2+c3*z+c4)) == 0
    H = a * b * (a + b) + a * U + b * V
    delta2 = H + L * Q
    delta3 = ((b * U + a * V) * H
              + Q * (b**2 * (a + b) * L + (a * (b + d) + b * L) * U
                     + (b * (b + d) + a * L) * V) + (b + d) * L * Q**2)
    assert s.expand(c1*c2-c3-delta2) == 0
    assert s.expand(c1*c2*c3-c1*c1*c4-c3*c3-delta3) == 0
    cubic = z**3+a*z**2+(U+Q)*z+d*Q
    assert s.expand(J.charpoly(z).as_expr().subs({B: 0, G: 0})-(z+d)*cubic) == 0
    assert s.expand(a*(U+Q)-d*Q-(a*U+L*Q)) == 0
    rng = np.random.default_rng(20261019)
    rows = []
    for i in range(400):
        lam = float(10**rng.uniform(-0.3, 1.7))
        d0, gamma = float(10**rng.uniform(-3, -1)), float(10**rng.uniform(-1.3, -0.2))
        pv, pb = rng.uniform(0.1, 0.95, 2)
        R0 = float(rng.uniform(1.3, 30))
        beta = R0 * d0 * (gamma + d0) / (lam * (1 - pb * pv))
        es, ev = rng.uniform(0.15, 1, 2)
        if i % 5 == 0:
            es = ev = 1.0
        p0 = Parameters(lam=lam, d=d0, gamma=gamma, pv=float(pv), pb=float(pb),
                        beta=beta, es=float(es), ev=float(ev))
        I00 = frozen_equilibrium(p0, 0, 0)[2]
        I10 = frozen_equilibrium(p0, 1, 0)[2]
        I11 = frozen_equilibrium(p0, 1, 1)[2]
        for active, lower, upper in [('s', I10, I00), ('v', I11, I10)]:
            if upper - lower < 1e-8:
                continue
            theta = lower + float(rng.uniform(0.1, 0.9)) * (upper - lower)
            if active == 's':
                p = replace(p0, ts=theta, tv=(theta + lam/d0)/2)
            else:
                p = replace(p0, ts=theta/2, tv=theta)
            eq = pinned_equilibrium(p, active)
            residual = float(np.max(np.abs(rhs_physical(eq, p))))
            eig = np.linalg.eigvals(jacobian(eq, p))
            assert 0 < eq[4 if active == 's' else 5] < 1
            assert residual < 1e-7
            assert eig.real.max() < 0
            rows.append({'sample': i, 'active': active, 'es': es, 'ev': ev,
                         'max_real_eigenvalue': float(eig.real.max()), 'residual': residual})
    (ROOT/'proofs'/'symbolic_identities.txt').write_text(
        'Exact SymPy identities PASS\nCharacteristic polynomial\n'+str(J.charpoly(z).as_expr())
        +'\nDelta2\n'+str(delta2)+'\nDelta3\n'+str(delta3)+'\n', encoding='utf-8')
    result = {'symbolic_characteristic_polynomial': 'PASS', 'symbolic_delta2': 'PASS',
              'symbolic_delta3': 'PASS', 'cubic_reduction': 'PASS',
              'random_pinned_jacobians_checked': len(rows),
              'maximum_real_eigenvalue': max(r['max_real_eigenvalue'] for r in rows),
              'maximum_equilibrium_residual': max(r['residual'] for r in rows),
              'samples': rows}
    (ROOT/'results'/'stability_verification.json').write_text(json.dumps(result, indent=2), encoding='utf-8')
    print({k:v for k,v in result.items() if k != 'samples'})


if __name__ == '__main__':
    main()
