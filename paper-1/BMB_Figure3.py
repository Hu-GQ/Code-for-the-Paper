"""Render BMB Figure 3 as finite-time samples, without changing source data.

Run from the extracted Online Resource 2 directory:
    python BMB_Figure3.py
The script reads the archived 48-point scan and writes figures/Figure_3.pdf.
The original src/plot_manuscript_figures.py is preserved for provenance.
"""
from pathlib import Path
import csv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def render(data_file, pdf_file, png_file=None):
    with Path(data_file).open(encoding='utf-8') as handle:
        rows=list(csv.DictReader(handle))
    samples=[(float(row['thetas']),float(row['I_end']),row['type']) for row in rows]
    assert len(samples)==48
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'pdf.fonttype':42,'ps.fonttype':42})
    fig,ax=plt.subplots(figsize=(5.8,3.4),layout='constrained')
    for tag,color,label in [('endemic','#0072B2','Other terminal samples'),
                            ('pinned_s','#D55E00','Near susceptible threshold')]:
        points=[(theta,infected) for theta,infected,kind in samples if kind==tag]
        ax.plot(*zip(*points),marker='s' if tag=='pinned_s' else 'o',linestyle='none',ms=3.5,color=color,label=label)
    ax.axhline(48.50603420458216,color='#AA4499',lw=1,label=r'$I_{212}$')
    ax.axhline(56.41995587716349,color='#009E73',lw=1,label=r'$I_{112}$')
    ax.plot([30,130],[30,130],'--',color='#D55E00',lw=.8,label=r'$I=\theta_s$')
    ax.axvspan(100,130,color='#999999',alpha=.14)
    ax.axvline(100,color='#777777',lw=.8,ls=':')
    ax.text(115,71,'Beyond\nthreshold ordering',ha='center',va='top',fontsize=8,color='#555555')
    ax.set(xlabel=r'$\theta_s$',ylabel=r'$I(2000)$',xlim=(30,130),ylim=(0,75))
    ax.legend(fontsize=7.5,ncol=2,loc='lower right');ax.grid(alpha=.2)
    Path(pdf_file).parent.mkdir(parents=True,exist_ok=True)
    fig.savefig(pdf_file)
    if png_file:fig.savefig(png_file,dpi=170)
    plt.close(fig)
    return samples

if __name__=='__main__':
    root=Path(__file__).resolve().parent
    render(root/'data/original/bif_thetas.csv',root/'figures/Figure_3.pdf',root/'figures/Figure_3.png')
