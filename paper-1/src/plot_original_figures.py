# Plot baseline SVIR trajectories for Figure 1 and Supplementary Figures S1-S3.
import csv, os, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams["font.family"] = ["DejaVu Sans", "sans-serif"]
plt.rcParams["axes.unicode_minus"] = False
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams.update({"font.size": 8, "axes.labelsize": 8, "axes.titlesize": 8,
                     "xtick.labelsize": 7.5, "ytick.labelsize": 7.5,
                     "legend.fontsize": 7, "lines.linewidth": 1.0})
HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, "..", "data", "original") + os.sep
OUT = os.path.join(HERE, "..", "figures", "original") + os.sep
os.makedirs(OUT, exist_ok=True)
C = {"S":"#009E73","V":"#0072B2","I":"#D55E00","R":"#CC79A7","x":"#0072B2","y":"#E69F00"}

def _is_float(s):
    try:
        float(s); return True
    except ValueError:
        return False

def load(fn, cols=None):
    rows = list(csv.reader(open(SRC+fn)))
    if rows and not _is_float(rows[0][0]):
        rows = rows[1:]  # skip header only if present
    data = [[float(x) for x in r] for r in rows]
    if cols is None:
        return data
    return [[r[c] for c in cols] for r in data]

# ---------- fig1: P110 (2 panels) ----------
d = load("rev_P110.csv")
t=[r[0] for r in d]; S=[r[1] for r in d]; V=[r[2] for r in d]
I=[r[3] for r in d]; R=[r[4] for r in d]; x=[r[5] for r in d]; y=[r[6] for r in d]
fig, axes = plt.subplots(1,2,figsize=(6.3,2.5))
axes[0].plot(t,S,color=C["S"],lw=1.2,label="S"); axes[0].plot(t,V,color=C["V"],lw=1.2,label="V")
axes[0].plot(t,I,color=C["I"],lw=1.2,label="I"); axes[0].plot(t,R,color=C["R"],lw=1.2,label="R")
axes[0].set_xlabel("t"); axes[0].set_ylabel("population"); axes[0].legend(fontsize=8); axes[0].set_title("(a) S,V,I,R"); axes[0].set_xlim(0,1000)
axes[1].plot(t,x,color=C["x"],lw=1.2,label=r"$x_A$"); axes[1].plot(t,y,color=C["y"],lw=1.2,label=r"$y_A$")
axes[1].set_ylim(-0.05,1.1); axes[1].set_xlabel("t"); axes[1].set_ylabel("compliance propensity"); axes[1].legend(fontsize=8); axes[1].set_title(r"(b) $x_A,y_A$"); axes[1].set_xlim(0,1000)
plt.tight_layout(); plt.savefig(OUT+"fig1.pdf"); plt.close()

# ---------- fig2: P112 (3 panels) ----------
d = load("rev_P112.csv")
t=[r[0] for r in d]; S=[r[1] for r in d]; V=[r[2] for r in d]
I=[r[3] for r in d]; R=[r[4] for r in d]; x=[r[5] for r in d]; y=[r[6] for r in d]
fig, axes = plt.subplots(1,3,figsize=(6.3,2.35))
axes[0].plot(t,S,color=C["S"],lw=1.1,label="S"); axes[0].plot(t,V,color=C["V"],lw=1.1,label="V")
axes[0].plot(t,I,color=C["I"],lw=1.1,label="I"); axes[0].plot(t,R,color=C["R"],lw=1.1,label="R")
axes[0].set_xlabel("t"); axes[0].set_ylabel("population"); axes[0].legend(fontsize=7); axes[0].set_title("(a) S,V,I,R")
axes[1].plot(t,x,color=C["x"],lw=1.1,label=r"$x_A$"); axes[1].plot(t,y,color=C["y"],lw=1.1,label=r"$y_A$")
axes[1].set_ylim(-0.05,1.1); axes[1].set_xlabel("t"); axes[1].set_ylabel("compliance propensity"); axes[1].legend(fontsize=7); axes[1].set_title(r"(b) $x_A,y_A$")
ax = axes[2]; ax.plot(t,I,color=C["I"],lw=1.2,label="I"); ax.axhline(70,ls="--",lw=0.8,color="#2e7d32"); ax.axhline(100,ls="--",lw=0.8,color="#e65100")
ax2 = ax.twinx(); ax2.plot(t,x,color=C["x"],lw=0.9,alpha=0.8,label=r"$x_A$"); ax2.plot(t,y,color=C["y"],lw=0.9,alpha=0.8,label=r"$y_A$")
ax2.set_ylim(-0.05,1.1); ax.set_xlabel("t"); ax.set_ylabel("infected level"); ax2.set_ylabel("compliance propensity")
h1,l1=ax.get_legend_handles_labels(); h2,l2=ax2.get_legend_handles_labels(); ax.legend(h1+h2,l1+l2,fontsize=7,loc="upper right")
ax.set_title(r"(c) $I,x_A,y_A$")
plt.tight_layout(); plt.savefig(OUT+"fig2.pdf"); plt.close()

# ---------- fig3: P51 (3 panels) ----------
d = load("rev_P51.csv")
t=[r[0] for r in d]; S=[r[1] for r in d]; V=[r[2] for r in d]
I=[r[3] for r in d]; R=[r[4] for r in d]; x=[r[5] for r in d]; y=[r[6] for r in d]
fig, axes = plt.subplots(1,3,figsize=(6.3,2.35))
axes[0].plot(t,S,color=C["S"],lw=1.1,label="S"); axes[0].plot(t,V,color=C["V"],lw=1.1,label="V")
axes[0].plot(t,I,color=C["I"],lw=1.1,label="I"); axes[0].plot(t,R,color=C["R"],lw=1.1,label="R")
axes[0].set_xlabel("t"); axes[0].set_ylabel("population"); axes[0].legend(fontsize=7); axes[0].set_title("(a) S,V,I,R"); axes[0].set_xlim(0,1600)
axes[1].plot(t,x,color=C["x"],lw=1.1,label=r"$x_A$"); axes[1].plot(t,y,color=C["y"],lw=1.1,label=r"$y_A$")
axes[1].set_ylim(-0.05,1.1); axes[1].set_xlabel("t"); axes[1].set_ylabel("compliance propensity"); axes[1].legend(fontsize=7); axes[1].set_title(r"(b) $x_A,y_A$"); axes[1].set_xlim(0,1600)
ax = axes[2]; ax.plot(t,I,color=C["I"],lw=1.2,label="I"); ax.axhline(56,ls="--",lw=0.8,color="#2e7d32"); ax.axhline(76,ls="--",lw=0.8,color="#e65100")
ax2 = ax.twinx(); ax2.plot(t,x,color=C["x"],lw=0.9,alpha=0.8,label=r"$x_A$"); ax2.plot(t,y,color=C["y"],lw=0.9,alpha=0.8,label=r"$y_A$")
ax2.set_ylim(-0.05,1.1); ax.set_xlabel("t"); ax.set_ylabel("infected level"); ax2.set_ylabel("compliance propensity")
h1,l1=ax.get_legend_handles_labels(); h2,l2=ax2.get_legend_handles_labels(); ax.legend(h1+h2,l1+l2,fontsize=7,loc="upper right")
ax.set_title(r"(c) $I,x_A,y_A$"); ax.set_xlim(0,1600)
plt.tight_layout(); plt.savefig(OUT+"fig3.pdf"); plt.close()

# ---------- fig4: theta_s/theta_v effect (2 panels) ----------
fig, axes = plt.subplots(1,2,figsize=(6.3,2.5))
for i,th in enumerate([10,30,50]):
    d=load(f"rev_fig4a_{i+1}.csv"); axes[0].plot([r[0] for r in d],[r[1] for r in d],lw=1.2,label=rf"$\theta_s={th},\,\theta_v={th+10}$")
axes[0].set_xlabel("t"); axes[0].set_ylabel("infected level"); axes[0].legend(fontsize=7); axes[0].set_title("(a) disease-free"); axes[0].set_xlim(0,400)
for i,th in enumerate([40,56,70]):
    d=load(f"rev_fig4b_{i+1}.csv"); axes[1].plot([r[0] for r in d],[r[1] for r in d],lw=1.2,label=rf"$\theta_s={th},\,\theta_v={th+20}$")
axes[1].set_xlabel("t"); axes[1].set_ylabel("infected level"); axes[1].legend(fontsize=7); axes[1].set_title("(b) endemic"); axes[1].set_xlim(0,900)
plt.tight_layout(); plt.savefig(OUT+"fig4.pdf"); plt.close()

print('Baseline figures exported to', OUT)
