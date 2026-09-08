"""Produce figures, witness, CSV tables, and a report from the exact Ising counts."""
from __future__ import annotations
from collections import Counter
import csv
import hashlib
import itertools
import json
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
import numpy as np
from scipy.optimize import minimize_scalar

from ising_exact import min_fill_order, thermodynamics

ROOT = Path(__file__).resolve().parent


def max_cut_witness(n, edges):
    """Max-plus elimination with backtracking, independent of the counting digits."""
    order, _ = min_fill_order(n, edges)
    factors = [((v,), np.zeros(2, dtype=int)) for v in range(n)]
    for u, v in edges:
        if u != v:
            factors.append((tuple(sorted((u, v))), np.array([[0, 1], [1, 0]], dtype=int)))
    decisions = []
    for v in order:
        selected = [(s, a) for s, a in factors if v in s]
        factors = [(s, a) for s, a in factors if v not in s]
        scope = tuple(sorted(set().union(*(s for s, _ in selected))))
        table = np.zeros((2,) * len(scope), dtype=int)
        for s, a in selected:
            table += a.reshape(tuple(2 if w in s else 1 for w in scope))
        axis = scope.index(v)
        rest = tuple(w for w in scope if w != v)
        decisions.append((v, rest, table.argmax(axis=axis)))
        factors.append((rest, table.max(axis=axis)))
    assignment = {}
    for v, rest, choice in reversed(decisions):
        assignment[v] = int(choice[tuple(assignment[w] for w in rest)])
    optimum = sum(int(a) for _, a in factors)
    assert sum(assignment[u] != assignment[v] for u, v in edges) == optimum
    return [assignment[v] for v in range(n)], optimum


def graph_layout(n, edges, steps=900):
    """Deterministic spring diagram; positions have no metric or physical meaning."""
    rng = np.random.default_rng(143)
    pos = rng.uniform(-1, 1, (n, 2))
    nonloops = np.array([(u, v) for u, v in edges if u != v], dtype=int)
    k = math.sqrt(4/n)
    for step in range(steps):
        delta = pos[:, None, :] - pos[None, :, :]
        d2 = np.maximum(np.sum(delta*delta, axis=2), 1e-8)
        disp = np.sum(delta*(k*k/d2)[:, :, None], axis=1)
        u, v = nonloops.T
        difference = pos[u] - pos[v]
        distance = np.maximum(np.linalg.norm(difference, axis=1), 1e-8)
        force = difference*(distance/k)[:, None]
        np.add.at(disp, u, -force)
        np.add.at(disp, v, force)
        disp -= 0.025*pos
        temp = 0.075*(1-step/steps) + 0.0003
        norm = np.maximum(np.linalg.norm(disp, axis=1), 1e-12)
        pos += disp*(np.minimum(norm, temp)/norm)[:, None]
        pos -= pos.mean(axis=0)
    pos /= np.max(np.abs(pos))
    return pos


def probability(row, beta):
    counts = row["density_of_states"]
    ks = np.array([k for k, v in enumerate(counts) if v])
    logs = np.log([v for v in counts if v]) - 2*beta*ks
    weights = np.exp(logs - max(logs))
    return ks, weights/weights.sum()


def main():
    report = json.loads((ROOT / "results.json").read_text())
    rows = {r["N"]: r for r in report["levels"]}
    r = rows[143]
    independent = json.loads((ROOT / "cycle_counts143.json").read_text())
    assert independent == r["even_subgraph_counts"]
    betas = np.linspace(0, 1.4, 561)
    peaks = {}
    for N, row in rows.items():
        fn = lambda b: thermodynamics(row["density_of_states"], row["spins"], row["edges_count"], b)["heat_capacity_per_spin"]
        grid = np.linspace(.01, 3, 1000)
        ix = int(np.argmax([fn(b) for b in grid]))
        opt = minimize_scalar(lambda b: -fn(b), bounds=(grid[max(0,ix-1)], grid[min(999,ix+1)]),
                              method="bounded", options={"xatol": 1e-13})
        peaks[str(N)] = {"beta": float(opt.x), "C_per_spin": float(-opt.fun),
                         "search_range": [.01, 3.0], "class": "DIAGNOSTIC numerical peak"}
    edges6 = Counter(tuple(sorted(e)) for e in rows[6]["edges"])
    edges11 = Counter(tuple(sorted(e)) for e in rows[11]["edges"])
    iso = next(p for p in itertools.permutations(range(4))
               if Counter(tuple(sorted((p[u],p[v]))) for u,v in rows[6]["edges"]) == edges11)
    assert rows[6]["density_of_states"] == rows[11]["density_of_states"]
    spins, cut = max_cut_witness(r["spins"], r["edges"])
    assert cut == r["max_cut"]
    frustrated = [i for i,(u,v) in enumerate(r["edges"]) if spins[u] == spins[v]]
    analysis = {"peaks": peaks, "isomorphism_6_to_11": list(iso),
                "full_143_cycle_count_matches": True,
                "cycle_enumeration_log": (ROOT / "cycle_enumeration.log").read_text().strip(),
                "max_cut_witness_143": {"spin_bits": spins, "cut_edges": cut,
                   "frustrated_edge_indices": frustrated,
                   "antiferromagnetic_energy": r["edges_count"] - 2*cut,
                   "ground_state_degeneracy": r["max_cut_degeneracy"]}}
    (ROOT / "analysis.json").write_text(json.dumps(analysis, indent=2) + "\n")

    with (ROOT / "density_of_states_143.csv").open("w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["disagreeing_edges", "ferromagnetic_energy_J1", "exact_spin_count", "exact_even_subgraph_count"])
        writer.writerows((k, 2*k-r["edges_count"], count, r["even_subgraph_counts"][k])
                         for k,count in enumerate(r["density_of_states"]))
    with (ROOT / "thermodynamics.csv").open("w", newline="") as fh:
        writer = None
        for N in [6, 11, 15, 35, 55, 77, 105, 143]:
            row = rows[N]
            for beta in betas:
                obs = {"N": N, **thermodynamics(row["density_of_states"], row["spins"], row["edges_count"], beta)}
                if writer is None:
                    writer = csv.DictWriter(fh, fieldnames=list(obs)); writer.writeheader()
                writer.writerow(obs)

    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10.5,
                         "axes.spines.top": False, "axes.spines.right": False,
                         "axes.titleweight": "medium", "axes.labelcolor": "#243140",
                         "text.color": "#243140", "xtick.color": "#43505d", "ytick.color": "#43505d",
                         "svg.fonttype": "none"})
    teal, gold, purple = "#087e8b", "#cf6b1d", "#7755a6"
    fig, axs = plt.subplots(1, 2, figsize=(12.8, 4.9))
    fig.subplots_adjust(left=.075, right=.98, top=.78, bottom=.24, wspace=.3)
    fig.suptitle("An exact thermometer for the Manin graph", x=.075, y=.97, ha="left", fontsize=20)
    fig.text(.075,.895,"MTFT v0.26.0  ·  all 72,057,594,037,927,936 spin states counted at N = 143", fontsize=11)
    for N, color, style in [(143,teal,"-"),(105,purple,"--")]:
        row=rows[N]
        vals=[thermodynamics(row["density_of_states"],row["spins"],row["edges_count"],b)["heat_capacity_per_spin"] for b in betas]
        axs[0].plot(betas,vals,color=color,linestyle=style,lw=2.3,label=f"N = {N}  ·  {row['spins']} spins")
    pk=peaks["143"]
    axs[0].plot(pk["beta"],pk["C_per_spin"],"o",color=teal,ms=5)
    axs[0].annotate(f"N = 143 peak\nβ ≈ {pk['beta']:.4f}", xy=(pk["beta"],pk["C_per_spin"]),
                    xytext=(.88,.86), arrowprops={"arrowstyle":"-","color":teal},fontsize=10)
    axs[0].set(xlabel="Inverse temperature β = J/(kBT)",ylabel="Heat capacity per spin C/(n kB)",xlim=(0,1.4),ylim=(0,1.23))
    axs[0].legend(frameon=False,fontsize=9.5,loc="upper left")
    axs[0].grid(axis="y",alpha=.17)
    axs[0].set_title("Same genus, different thermal curves",loc="left",pad=13)
    max_probability = 0.
    for beta,col in [(.2,purple),(.64,gold),(1.,teal)]:
        k,p=probability(r,beta)
        max_probability = max(max_probability, float(p.max()))
        axs[1].plot(k,p,color=col,lw=2,label=f"β = {beta:.2f}")
        axs[1].fill_between(k,p,color=col,alpha=.1)
    axs[1].set(xlabel="Number of disagreeing edges k",ylabel="Thermal probability Pβ(k)",xlim=(-1,r["max_cut"]+1),ylim=(0,1.08*max_probability))
    axs[1].legend(frameon=False,loc="upper right",fontsize=10)
    axs[1].grid(axis="y",alpha=.17)
    axs[1].set_title("Cooling moves weight toward aligned spins",loc="left",pad=13)
    fig.text(.075,.062,"Counts: exact integers, independently checked. Curves: numerical evaluation of finite sums.\nThese are finite-graph peaks; no continuum critical point is asserted.",fontsize=9.5,color="#53616e")
    fig.savefig(ROOT / "MTFT_Ising_Thermometer.png",dpi=190,facecolor="white")
    fig.savefig(ROOT / "MTFT_Ising_Thermometer.svg",facecolor="white")
    plt.close(fig)

    pos=graph_layout(r["spins"],r["edges"])
    fig,ax=plt.subplots(figsize=(9,8.5))
    fig.subplots_adjust(left=.05,right=.95,top=.84,bottom=.20)
    fig.suptitle("When neighbors prefer to disagree",x=.07,y=.965,ha="left",fontsize=20)
    fig.text(.07,.915,"N = 143  ·  antiferromagnetic ground state  ·  75 of 84 edges satisfied",fontsize=11)
    for i,(u,v) in enumerate(r["edges"]):
        fail=i in frustrated
        color="#c44445" if fail else "#b5c1ca"
        if u==v:
            direction=pos[u]/max(np.linalg.norm(pos[u]),1e-8)
            ax.add_patch(Circle(pos[u]+.054*direction,.054,fill=False,color=color,lw=2.3,zorder=3))
        else:
            ax.plot(*pos[[u,v]].T,color=color,lw=2.3 if fail else 1,alpha=1 if fail else .65,zorder=1,
                    linestyle="--" if fail else "-")
    for spin,col,label in [(0,teal,"Spin −1"),(1,gold,"Spin +1")]:
        mask=np.array(spins)==spin
        ax.scatter(*pos[mask].T,s=125,color=col,edgecolors="white",linewidths=.8,zorder=4,label=label)
    ax.plot([],[],"--",color="#c44445",lw=2.3,label="Unsatisfied edge")
    handles, labels = ax.get_legend_handles_labels()
    fig.legend(handles,labels,frameon=False,loc="lower center",bbox_to_anchor=(.5,.11),ncol=3)
    ax.set_aspect("equal")
    ax.set_xlim(pos[:,0].min()-.14,pos[:,0].max()+.14)
    ax.set_ylim(pos[:,1].min()-.14,pos[:,1].max()+.14)
    ax.axis("off")
    fig.text(.07,.045,"Exactly two ground states: this coloring and its global reversal.\nNine unsatisfied edges include one self-loop. Spring layout is schematic, not surface coordinates.",fontsize=9.5,color="#53616e")
    fig.savefig(ROOT / "MTFT_Ising_Antiferromagnet.png",dpi=180,facecolor="white")
    fig.savefig(ROOT / "MTFT_Ising_Antiferromagnet.svg",facecolor="white")
    plt.close(fig)
    print(json.dumps(analysis,indent=2))


if __name__ == "__main__":
    main()
