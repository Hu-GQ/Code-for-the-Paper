"""Enumerate feasible equilibria and verify the ten-equilibrium example."""
from pathlib import Path
import csv
import json
import sys

import numpy as np

ROOT=Path(__file__).resolve().parents[1]
from svir_model import Parameters, frozen_equilibrium, rhs_physical


def candidates(p):
    k=p.gamma+p.d
    found=[]

    def keep(name,u):
        if min(u[:4]) < -1e-9 or min(u[4:]) < -1e-9 or max(u[4:]) > 1+1e-9:
            return
        residual=float(max(abs(rhs_physical(u,p))))
        assert residual<1e-9,(name,residual)
        if not any(np.allclose(u,v,rtol=1e-9,atol=1e-8) for _,v,_ in found):
            found.append((name,u,residual))

    for (x,y),dfe,endemic in [((0,0),'P110','P112'),((0,1),'P121','P122'),
                              ((1,0),'P211','P212'),((1,1),'P22',None)]:
        keep(dfe,np.array([(1-p.pv)*p.lam/p.d,p.pv*p.lam/p.d,0,0,x,y]))
        u=frozen_equilibrium(p,x,y)
        if u[2]>1e-9:
            assert endemic is not None
            keep(endemic,u)
    for active,fixed,name in [('v',0,'P41'),('v',1,'P42'),('s',0,'P51'),('s',1,'P52')]:
        if active=='s':
            I=p.ts;y=fixed;kv=(1-p.pb)*p.beta*(1-y)
            V=p.pv*p.lam/(p.d+kv*I)
            S=(p.lam-k*I)/p.d-V
            if S<=0:
                continue
            x=1-(k-kv*V)/(p.beta*S)
            if not 0<x<1:
                continue
        else:
            I=p.tv;x=fixed;ks=p.beta*(1-x)
            S=(1-p.pv)*p.lam/(p.d+ks*I)
            V=(p.lam-k*I)/p.d-S
            if V<=0:
                continue
            y=1-(k-ks*S)/((1-p.pb)*p.beta*V)
            if not 0<y<1:
                continue
        keep(name,np.array([S,V,I,p.gamma*I/p.d,x,y]))
    assert len(found)<=10
    return found


def main():
    (ROOT/'results').mkdir(parents=True,exist_ok=True)
    cases=[(2,40),(30,40),(52,100)]
    rows=[]
    counts=[]
    for ts,tv in cases:
        p=Parameters(ts=ts,tv=tv)
        items=candidates(p)
        levels={f'J{x}{y}':float(frozen_equilibrium(p,x,y)[2]) for x,y in [(0,0),(0,1),(1,0),(1,1)]}
        counts.append({'theta_s':ts,'theta_v':tv,'count':len(items),'names':[n for n,_,_ in items],'levels':levels})
        for name,u,residual in items:
            rows.append(dict(theta_s=ts,theta_v=tv,equilibrium=name,S=u[0],V=u[1],I=u[2],R=u[3],x=u[4],y=u[5],residual=residual))
    assert [x['count'] for x in counts]==[10,9,8]
    # The two stable interior types require theta_v<J10<theta_s to coexist.
    # This contradicts theta_s<theta_v, so four interior candidates cannot coexist.
    tested=0
    rng=np.random.default_rng(20261020)
    for _ in range(1000):
        levels=sorted(rng.uniform(0.05,100,2))
        p=Parameters(ts=levels[0],tv=levels[1],pv=float(rng.uniform(.05,.95)),pb=float(rng.uniform(.05,.95)),beta=float(rng.uniform(.0001,.01)))
        candidates(p)
        tested+=1
    with (ROOT/'results/equilibrium_count_checks.csv').open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    result={'cases':counts,'maximum_residual':max(r['residual'] for r in rows),'random_parameter_checks':tested,
            'general_upper_bound':10,'bound_is_attained':True,
            'proof':'Four disease-free corners, at most three endemic corners, and at most three single-active equilibria. P42 and P51 cannot coexist under theta_s<theta_v. Collision duplicates cannot increase the count.'}
    (ROOT/'results/equilibrium_count_checks.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    print(json.dumps(result,indent=2))


if __name__=='__main__':
    main()
