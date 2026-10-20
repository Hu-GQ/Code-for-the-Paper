"""Partial-protection branch tests and exact-logit integrations for Figure 5.

Run from any directory. All outputs stay below this project root.
"""
from dataclasses import replace
from pathlib import Path
import csv
import json
import time

import numpy as np
from scipy.integrate import solve_ivp
from scipy.special import logit

from svir_model import (Parameters, rhs_logit, physical_state, rhs_physical,
                        jacobian, frozen_equilibrium, predicted_equilibrium, self_check)

ROOT = Path(__file__).resolve().parents[1]
INITS = [(80, 900, 10, .2, .2), (300, 680, 20, .5, .3),
         (500, 480, 20, .8, .5), (60, 920, 20, .1, .9),
         (200, 780, 20, .95, .05)]


def integrate(p, init=(103, 897, 10, .3, .2), method='LSODA', tol=1e-10, horizon=20000):
    S, V, I, x, y = init
    z0 = [S, V, I, 0., logit(x), logit(y)]
    name, eq = predicted_equilibrium(p)
    initial_horizon = horizon
    initial_error = None
    for attempt in range(3):
        sol = solve_ivp(rhs_logit, (0, horizon), z0, args=(p,), method=method,
                        rtol=tol, atol=tol*.01, dense_output=True)
        assert sol.success, sol.message
        u = physical_state(sol.y[:, -1])
        err = float(np.max(np.abs(u-eq)))
        if initial_error is None:
            initial_error = err
        if err < 1e-6 or attempt == 2:
            break
        horizon *= 10
    assert np.min(sol.y[:4]) > -1e-6
    N0 = S + V + I
    expected_N = p.lam/p.d + (N0-p.lam/p.d)*np.exp(-p.d*sol.t)
    assert np.max(np.abs(sol.y[:4].sum(axis=0)-expected_N)) < 1e-5
    residual = float(np.max(np.abs(rhs_physical(u, p))))
    return sol, dict(es=p.es, ev=p.ev, theta_s=p.ts, theta_v=p.tv,
                     region=name, method=method, rtol=tol, horizon=horizon,
                     initial_horizon=initial_horizon, initial_state_error=initial_error,
                     I_pred=float(eq[2]), x_pred=float(eq[4]), y_pred=float(eq[5]),
                     I_end=float(u[2]), x_end=float(u[4]), y_end=float(u[5]),
                     max_state_error=err, residual=residual,
                     max_real_eigenvalue=float(np.linalg.eigvals(jacobian(eq,p)).real.max()),
                     evaluations=sol.nfev)


def save_csv(name, rows):
    with (ROOT/'results'/name).open('w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)


def main():
    self_check()
    started = time.perf_counter()
    base = Parameters()
    rows, paths = [], {}
    for efficacy in [.25, .5, .6, .75, .9, 1.0]:
        for pair in [(56.,76.), (30.,40.)]:
            p = replace(base, es=efficacy, ev=efficacy, ts=pair[0], tv=pair[1])
            sol, row = integrate(p)
            rows.append(row)
            paths[(efficacy,pair)] = sol
            print(f"eps={efficacy:.2f}, theta={pair}: {row['region']}, I={row['I_end']:.8f}, "
                  f"error={row['max_state_error']:.2e}", flush=True)
    for es, ev in [(.6,.9),(.9,.6)]:
        p = replace(base, es=es, ev=ev, ts=30, tv=40)
        _, row = integrate(p)
        rows.append(row)
    save_csv('robustness_integrations.csv', rows)
    # Test one interior point of each available branch for every efficacy.
    branch_rows = []
    for efficacy in [.25,.5,.75,1.]:
        p0 = replace(base, es=efficacy, ev=efficacy)
        I00,I10,I11 = [frozen_equilibrium(p0,*c)[2] for c in [(0,0),(1,0),(1,1)]]
        settings = [(1.1*I00,1.3*I00), ((I00+I10)/2,1.2*I00),
                    (.5*I10,(I10+I00)/2)]
        if I10 > 0:
            tv = (I11+I10)/2
            settings.append((tv/2,tv))
        if I11 > 0:
            settings.append((I11/3,2*I11/3))
        for ts,tv in settings:
            _, row = integrate(replace(p0,ts=ts,tv=tv),horizon=30000)
            branch_rows.append(row)
    save_csv('branch_checkpoints.csv',branch_rows)
    multi=[]
    for ts,tv in [(56,76),(70,100)]:
        for init in INITS:
            _, row = integrate(replace(base,ts=ts,tv=tv),init=init,horizon=30000)
            row.update(dict(zip(['S0','V0','I0','x0','y0'],init)))
            multi.append(row)
    save_csv('multistart_corrected.csv',multi)
    fig7=[]
    for beta,ts in [(.002,20),(.002,45),(.002,70),(.005,30),(.005,54),(.005,80)]:
        _, row = integrate(replace(base,beta=beta,ts=ts,tv=100),horizon=60000)
        row['beta']=beta
        fig7.append(row)
    save_csv('fig7_checkpoints_corrected.csv',fig7)
    controls=[]
    for efficacy,ts,tv in [(.75,56,76),(.75,30,40),(1.,56,76)]:
        p=replace(base,es=efficacy,ev=efficacy,ts=ts,tv=tv)
        for method,tol in [('LSODA',1e-11),('DOP853',1e-10)]:
            _,row=integrate(p,method=method,tol=tol,horizon=20000)
            controls.append(row)
    save_csv('solver_controls.csv',controls)

    grid=[]
    for efficacy in np.linspace(.05,1,191):
        p=replace(base,es=float(efficacy),ev=float(efficacy))
        i00,i10,i11=[frozen_equilibrium(p,*c)[2] for c in [(0,0),(1,0),(1,1)]]
        grid.append(dict(epsilon=efficacy,I00=i00,I10=i10,I11=i11,width=i00-i10))
    save_csv('efficacy_branches.csv',grid)
    from plot_robustness import make_figure
    make_figure(paths)
    allrows=rows+branch_rows+multi+fig7+controls
    summary=dict(integrations=len(allrows),robustness_runs=len(rows),branch_runs=len(branch_rows),
                 multistart_runs=len(multi),fig7_runs=len(fig7),solver_control_runs=len(controls),
                 max_state_error=max(r['max_state_error'] for r in allrows),
                 max_infected_error=max(abs(r['I_end']-r['I_pred']) for r in allrows),
                 max_physical_residual=max(r['residual'] for r in allrows),
                 max_fig7_infected_error=max(abs(r['I_end']-r['I_pred']) for r in fig7),
                 max_multistart_infected_error=max(abs(r['I_end']-r['I_pred']) for r in multi),
                 extended_horizon_runs=sum(r['horizon']>r['initial_horizon'] for r in allrows),
                 elapsed_seconds=time.perf_counter()-started)
    # Near-boundary transients must not silently be called converged.
    assert summary['max_infected_error'] < 1e-3, summary
    assert summary['max_physical_residual'] < 1e-4, summary
    (ROOT/'results'/'robustness_summary.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
    print(json.dumps(summary,indent=2),flush=True)


if __name__=='__main__':
    main()
