"""Export the data and four-panel robustness figure at readable print size."""
from pathlib import Path
import csv
from dataclasses import replace

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1]


def make_figure(paths=None):
    (ROOT/'figures').mkdir(parents=True,exist_ok=True)
    if paths is None:
        from robustness import integrate
        from svir_model import Parameters
        paths={}
        for eff in [.5,.75,1.]:
            for pair in [(56.,76.),(30.,40.)]:
                paths[(eff,pair)]=integrate(replace(Parameters(),es=eff,ev=eff,ts=pair[0],tv=pair[1]))[0]
    with (ROOT/'results'/'efficacy_branches.csv').open() as f:
        grid=[{k:float(v) for k,v in row.items()} for row in csv.DictReader(f)]
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9.5,
                         'pdf.fonttype':42,'ps.fonttype':42,'axes.spines.top':False,
                         'axes.spines.right':False,'xtick.labelsize':8.5,'ytick.labelsize':8.5})
    fig,axs=plt.subplots(2,2,figsize=(7.2,6.5),layout='constrained')
    ep=np.array([r['epsilon'] for r in grid])
    ax=axs[0,0]
    for field,color,style in [('I00','#274C77','-'),('I10','#B15A2B','-'),('I11','#777777','--')]:
        ax.plot(ep,[r[field] for r in grid],color=color,ls=style,label='$J_{'+field[1:]+'}$')
    ax.fill_between(ep,[r['I10'] for r in grid],[r['I00'] for r in grid],color='#B15A2B',alpha=.15)
    ax.set(xlabel=r'Protection efficacy $\varepsilon$',ylabel='Endemic infected level',title='(a) Frozen-strategy branches')
    ax.legend(frameon=False,loc='lower left',fontsize=8)
    ax=axs[0,1]
    ax.plot(ep,[r['width'] for r in grid],color='#B15A2B')
    ax.set(xlabel=r'Protection efficacy $\varepsilon$',ylabel=r'$J_{00}-J_{10}$',title='(b) Susceptible-pinning width')
    colors=['#274C77','#B15A2B','#53816D']
    t=np.linspace(0,2000,4001)
    for ax,pair,panel in zip(axs[1],[(56.,76.),(30.,40.)],['c','d']):
        data={'t':t}
        for color,eff in zip(colors,[.5,.75,1.]):
            values=paths[(eff,pair)].sol(t)
            ax.plot(t,values[2],color=color,lw=.8,label=rf'$\varepsilon={eff:g}$')
            data[f'I_eps_{eff:g}']=values[2]
        file=ROOT/'results'/f'fig8_timeseries_theta_{int(pair[0])}_{int(pair[1])}.csv'
        with file.open('w',newline='') as f:
            writer=csv.writer(f)
            writer.writerow(data)
            writer.writerows(zip(*data.values()))
        ax.axhline(pair[0] if pair[0]==56 else pair[1],color='#999999',ls=':',lw=.9)
        ax.set(xlabel='Time (days)',ylabel='Infected individuals',
               title=rf'({panel}) $\theta_s={pair[0]:g},\ \theta_v={pair[1]:g}$')
        ax.legend(frameon=False,fontsize=8)
    fig.savefig(ROOT/'figures'/'Figure_5.pdf')
    fig.savefig(ROOT/'figures'/'Figure_5.png',dpi=180)
    plt.close(fig)
    print('Figure 8 and its sampled time-series data exported',flush=True)


if __name__=='__main__':
    make_figure()
