"""Create the submission figure set from archived data and two VIR checks."""
from pathlib import Path
import csv
import json
import shutil
import subprocess
import sys

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.special import expit,logit

ROOT=Path(__file__).resolve().parents[1]
REPRO=ROOT
DATA=REPRO/'data/original'
RESULTS=REPRO/'results'
OUT=REPRO/'figures'


def load(name):
    with (DATA/name).open() as f:
        rows=list(csv.reader(f))
    try:float(rows[0][0])
    except ValueError:rows=rows[1:]
    return np.array(rows,dtype=float)


def full_vaccination():
    fig,axs=plt.subplots(1,3,figsize=(9.0,2.9),layout='constrained')
    colors=['#0072B2','#D55E00','#009E73','#AA4499','#CCBB44']
    for i,theta in enumerate([20,50,80]):
        a=load(f'rev_fig5a_{i+1}.csv')
        axs[0].plot(a[:,0],a[:,1],color=colors[i],lw=1,label=rf'$\theta_v={theta}$')
    for i,theta in enumerate([19,38,57,76,95]):
        a=load(f'rev_fig5b_{i+1}.csv')
        axs[1].plot(a[:,0],a[:,1],color=colors[i],lw=.75,label=rf'$\theta_v={theta}$')
    rows=[]
    for color,theta in zip(colors,[19.,38.]):
        def rhs(t,u):
            V,I,R,z=u
            y,ny=expit(z),expit(-z)
            incidence=(1-.7934)*.003*ny*V*I
            return [10-incidence-.007*V,incidence-.15*I,.143*I-.007*R,(1.3*y+2*ny)*(I-theta)]
        sol=solve_ivp(rhs,(0,20000),[1000,10,0,logit(.2)],method='LSODA',rtol=1e-11,atol=1e-13,dense_output=True)
        assert sol.success and np.min(sol.y[:3])>-1e-8
        expected=10/.007+(1010-10/.007)*np.exp(-.007*sol.t)
        assert max(abs(sol.y[:3].sum(axis=0)-expected))<1e-6
        t=np.unique(np.r_[np.geomspace(10,20000,4000),2000,20000])
        values=sol.sol(t)
        error=np.abs(values[1]-theta)
        assert error[-1]<1e-6
        axs[2].plot(t,np.maximum(error,1e-12),color=color,lw=.8,label=rf'$\theta_v={theta:g}$')
        with (RESULTS/f'vir_long_time_theta_{int(theta)}.csv').open('w',newline='') as f:
            w=csv.writer(f);w.writerow(['t','V','I','R','logit_y','absolute_infected_error']);w.writerows(zip(t,*values,error))
        rows.append(dict(theta_v=theta,I_2000=float(sol.sol(2000)[1]),I_20000=float(values[1,-1]),predicted_I=theta,terminal_error=float(error[-1]),evaluations=sol.nfev))
    axs[0].set(title='(a) Disease-free',xlabel='Time (days)',ylabel='Infected individuals')
    axs[1].set(title='(b) Endemic, early window',xlabel='Time (days)',xlim=(0,2000))
    axs[2].set(title='(c) Long-time agreement',xlabel='Time (days)',ylabel=r'$|I-\theta_v|$',xscale='log',yscale='log',xlim=(10,20000))
    axs[2].axvline(2000,color='#999999',ls=':',lw=.8)
    for ax in axs:ax.legend(fontsize=7,frameon=False)
    fig.savefig(OUT/'Figure_2.pdf');fig.savefig(OUT/'Figure_2.png',dpi=170);plt.close(fig)
    with (RESULTS/'full_vaccination_late_checks.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    theta=38;V=(10-.15*theta)/.007;y=1-.15/((1-.7934)*.003*V)
    L=.15/V*theta;F=.15;C=(1-.7934)*.003*V*theta;E=y*(1-y)*(1.3*y+2*(1-y))
    spectrum=np.roots([1,L+.007,L*F+C*E,.007*C*E])
    summary={'late_checks':rows,'P64_spectrum':[[float(z.real),float(z.imag)] for z in spectrum],
             'V':V,'y':y,'maximum_terminal_error':max(r['terminal_error'] for r in rows)}
    (RESULTS/'full_vaccination_plot_checks.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
    print(json.dumps(summary,indent=2))


def collision():
    with (DATA/'bif_thetas.csv').open() as f:rows=list(csv.DictReader(f))
    # Keep the original column order; the final column is a regime label.
    values=[]
    for row in rows:
        cells=list(row.values());values.append((float(cells[0]),float(cells[2]),cells[5]))
    fig,ax=plt.subplots(figsize=(5.8,3.4),layout='constrained')
    for tag,color,label in [('endemic','#0072B2','Endemic plateau'),('pinned_s','#D55E00','Susceptible pinning')]:
        points=[(t,i) for t,i,s in values if s==tag]
        ax.plot(*zip(*points),'o',ms=3.5,color=color,label=label)
    ax.axhline(48.50603420458216,color='#AA4499',lw=1,label=r'$I_{212}$')
    ax.axhline(56.41995587716349,color='#009E73',lw=1,label=r'$I_{112}$')
    ax.plot([30,130],[30,130],'--',color='#D55E00',lw=.8,label=r'$I=\theta_s$')
    ax.axvspan(100,130,color='#999999',alpha=.14)
    ax.axvline(100,color='#777777',lw=.8,ls=':')
    ax.text(115,71,'Beyond\nthreshold ordering',ha='center',va='top',fontsize=8,color='#555555')
    ax.set(xlabel=r'$\theta_s$',ylabel=r'$I_{\rm eq}$',xlim=(30,130),ylim=(0,75))
    ax.legend(fontsize=7.5,ncol=2,loc='lower right');ax.grid(alpha=.2)
    fig.savefig(OUT/'Figure_3.pdf');fig.savefig(OUT/'Figure_3.png',dpi=170);plt.close(fig)


def efficacy():
    with (RESULTS/'efficacy_branches.csv').open() as f:
        grid=[{k:float(v) for k,v in row.items()} for row in csv.DictReader(f)]
    eps=np.array([x['epsilon'] for x in grid])
    fig,axs=plt.subplots(2,2,figsize=(7.2,6.0),layout='constrained')
    for field,color,style in [('I00','#0072B2','-'),('I10','#D55E00','-'),('I11','#777777','--')]:
        axs[0,0].plot(eps,[r[field] for r in grid],color=color,ls=style,label='$J_{'+field[1:]+'}$')
    axs[0,0].fill_between(eps,[r['I10'] for r in grid],[r['I00'] for r in grid],color='#D55E00',alpha=.14)
    axs[0,0].set(xlabel=r'Protection efficacy $\varepsilon$',ylabel='Endemic infected level',title='(a) Frozen-strategy levels')
    axs[0,0].legend(frameon=False,loc='lower left')
    axs[0,1].plot(eps,[r['width'] for r in grid],color='#D55E00')
    axs[0,1].set(xlabel=r'Protection efficacy $\varepsilon$',ylabel=r'$J_{00}-J_{10}$',title='(b) Susceptible-pinning width')
    for ax,pair,panel in zip(axs[1],[(56,76),(30,40)],['c','d']):
        with (RESULTS/f'fig8_timeseries_theta_{pair[0]}_{pair[1]}.csv').open() as f:
            samples=list(csv.DictReader(f))
        t=[float(r['t']) for r in samples]
        for color,eff in zip(['#0072B2','#D55E00','#009E73'],[.5,.75,1.]):
            ax.plot(t,[float(r[f'I_eps_{eff:g}']) for r in samples],color=color,lw=.75,label=rf'$\varepsilon={eff:g}$')
        ax.axhline(pair[0] if pair[0]==56 else pair[1],color='#888888',ls=':',lw=.8)
        ax.set(xlabel='Time (days)',ylabel='Infected individuals',title=rf'({panel}) $\theta_s={pair[0]},\ \theta_v={pair[1]}$')
        ax.legend(frameon=False,fontsize=8)
    fig.savefig(OUT/'Figure_5.pdf');fig.savefig(OUT/'Figure_5.png',dpi=170);plt.close(fig)


def main():
    OUT.mkdir(exist_ok=True)
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'pdf.fonttype':42,'ps.fonttype':42})
    subprocess.run([sys.executable,str(REPRO/'src/plot_original_figures.py')],check=True)
    for original,destination in [('fig3.pdf','Figure_1.pdf'),('fig1.pdf','Figure_S1.pdf'),('fig2.pdf','Figure_S2.pdf'),('fig4.pdf','Figure_S3.pdf')]:
        shutil.copy2(OUT/'original'/original,OUT/destination)
    subprocess.run([sys.executable,str(REPRO/'src/plot_regime.py')],check=True)
    full_vaccination();collision();efficacy()
    print('Five main figures and three supporting figures exported')


if __name__=='__main__':
    main()
