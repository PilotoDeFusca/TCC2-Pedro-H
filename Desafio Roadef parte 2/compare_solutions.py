#!/usr/bin/env python3
"""
compare_solutions.py

Read the AMNS, VND and (if present) trivial solution JSON files and produce
comparison figures. Each instance gets its own subfolder under --outdir.

Expected layout (matching the C++ exporters):

    amns_solution/<instance>-solution.json
    vnd_solution/<instance>-solution.json
    trivial_solution/<instance>-solution.json      (optional)

Output layout:

    figures/
    |-- summary-instances.png   every instance: worst-slot load, change vs
    |                           trivial, solve time (one row per instance)
    |-- summary-by-set.png      mean change vs trivial per instance set
    |-- summary.csv             the numbers behind the summary figures
    |-- setA-01/
    |   |-- overview.png        max / mean load per slot, load distribution,
    |   |                       waypoints, reduction by move, solve time
    |   `-- history.png         search progress, one panel per slot
    `-- ...

Usage:
    python compare_solutions.py
    python compare_solutions.py --root path/to/dir --outdir figures
    python compare_solutions.py --no-history      (summary + overviews only)
"""

import argparse
import json
import locale
import sys
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional

import matplotlib.pyplot as plt
import numpy as np


# ---------------------------------------------------------------------------
# Data structures
# ---------------------------------------------------------------------------

@dataclass
class DemandPath:
    demand_index: int
    time_slot: int
    waypoints: List[int]
    links: List[int]
    routed_demand: float
    feasible: bool
    total_metric: float


@dataclass
class LinkLoad:
    time_slot: int
    link_id: int
    capacity: float
    total_demand: float
    load: float


@dataclass
class MaxLoadPoint:
    iteration: int
    max_load: float
    op: str = "start"          # move that produced this drop


@dataclass
class NbhdSwitch:
    iteration: int             # this neighborhood is active from here on
    op: str


@dataclass
class Solution:
    instance_name: str
    max_link_load: float
    max_link_id: int
    max_link_slot: int
    total_metric: float
    solve_time_seconds: float = 0.0
    paths: List[DemandPath] = field(default_factory=list)
    link_loads: List[LinkLoad] = field(default_factory=list)
    # max_load_history[slot] = list of MaxLoadPoint
    max_load_history: List[List["MaxLoadPoint"]] = field(default_factory=list)
    # nbhd_switches[slot] = list of NbhdSwitch (empty for older JSON files)
    nbhd_switches: List[List["NbhdSwitch"]] = field(default_factory=list)

    def num_slots(self) -> int:
        slots = {p.time_slot for p in self.paths}
        slots |= {l.time_slot for l in self.link_loads}
        return (max(slots) + 1) if slots else 0

    def num_demands(self) -> int:
        return max((p.demand_index for p in self.paths), default=-1) + 1

    def path_at(self, slot: int, demand: int) -> Optional[DemandPath]:
        for p in self.paths:
            if p.time_slot == slot and p.demand_index == demand:
                return p
        return None

    def link_loads_at(self, slot: int) -> List[LinkLoad]:
        return [l for l in self.link_loads if l.time_slot == slot]

    def max_load_per_slot(self) -> List[float]:
        out = []
        for t in range(self.num_slots()):
            loads = [l.load for l in self.link_loads_at(t)
                     if np.isfinite(l.load)]
            out.append(max(loads) if loads else 0.0)
        return out

    def avg_load_per_slot(self) -> List[float]:
        out = []
        for t in range(self.num_slots()):
            loads = [l.load for l in self.link_loads_at(t)
                     if np.isfinite(l.load)]
            out.append(float(np.mean(loads)) if loads else 0.0)
        return out

    def overflow_per_slot(self) -> List[float]:
        out = []
        for t in range(self.num_slots()):
            ov = sum(max(0.0, l.load - 1.0)
                     for l in self.link_loads_at(t)
                     if np.isfinite(l.load))
            out.append(ov)
        return out

    def total_waypoints_per_slot(self) -> List[int]:
        out = []
        for t in range(self.num_slots()):
            out.append(sum(len(p.waypoints)
                           for p in self.paths
                           if p.time_slot == t and p.feasible))
        return out

    def all_loads(self) -> List[float]:
        return [l.load for l in self.link_loads if np.isfinite(l.load)]


# ---------------------------------------------------------------------------
# Loading
# ---------------------------------------------------------------------------

def load_solution(path: Path) -> Solution:
    with open(path, "r") as f:
        j = json.load(f)

    paths = [
        DemandPath(
            demand_index=e["d"],
            time_slot=e["t"],
            waypoints=list(e.get("w", [])),
            links=list(e.get("links", [])),
            routed_demand=float(e.get("routed_demand", 0.0)),
            feasible=bool(e.get("feasible", False)),
            total_metric=float(e.get("total_metric", 0.0)),
        )
        for e in j.get("paths", [])
    ]

    loads = [
        LinkLoad(
            time_slot=e["t"],
            link_id=e["link_id"],
            capacity=float(e["capacity"]),
            total_demand=float(e["total_demand"]),
            load=float(e["load"]),
        )
        for e in j.get("link_loads", [])
    ]

    raw_history = j.get("max_load_history", [])
    history: List[List[MaxLoadPoint]] = []
    for slot_arr in raw_history:
        slot_pts = [
            MaxLoadPoint(
                iteration=int(e["iter"]),
                max_load=float(e["max_load"]),
                op=str(e.get("op", "start")),
            )
            for e in slot_arr
        ]
        history.append(slot_pts)

    switches: List[List[NbhdSwitch]] = [
        [NbhdSwitch(iteration=int(e["iter"]), op=str(e["op"])) for e in slot_arr]
        for slot_arr in j.get("nbhd_switches", [])
    ]

    return Solution(
        instance_name=j.get("instance_name", path.stem),
        max_link_load=float(j.get("max_link_load", 0.0)),
        max_link_id=int(j.get("max_link_id", -1)),
        max_link_slot=int(j.get("max_link_slot", -1)),
        total_metric=float(j.get("total_metric", 0.0)),
        solve_time_seconds=float(j.get("solve_time_seconds", 0.0)),
        paths=paths,
        link_loads=loads,
        max_load_history=history,
        nbhd_switches=switches,
    )


def discover_instances(root: Path) -> List[str]:
    amns_dir = root / "amns_solution"
    vnd_dir  = root / "vnd_solution"

    if not amns_dir.is_dir() or not vnd_dir.is_dir():
        print(f"[warn] expected folders missing under {root}: "
              f"{amns_dir.name}/ and {vnd_dir.name}/", file=sys.stderr)
        return []

    def stems(d: Path) -> Dict[str, Path]:
        out = {}
        for p in d.glob("*-solution.json"):
            stem = p.name[:-len("-solution.json")]
            out[stem] = p
        return out

    a = stems(amns_dir)
    v = stems(vnd_dir)
    common = sorted(set(a) & set(v))
    if not common:
        print("[warn] no instance has both AMNS and VND solutions.",
              file=sys.stderr)
    return common


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _savefig(fig, out: Path, dpi: int, tries: int = 6):
    """fig.savefig with a few retries: in a OneDrive-synced folder a file can
    be briefly locked while it syncs (OSError: Invalid argument)."""
    for k in range(tries):
        try:
            fig.savefig(out, dpi=dpi, bbox_inches="tight")
            return
        except OSError:
            if k == tries - 1:
                raise
            time.sleep(0.5 * (k + 1))


# ---------------------------------------------------------------------------
# Style
# ---------------------------------------------------------------------------
#
# One encoding everywhere: AMNS blue / filled circles, VND orange / hollow
# squares (so the two stay visible when they coincide), trivial baseline gray /
# dotted. Colors are slots 1-2 of a CVD-checked categorical palette; the gray
# baseline is a reference, not a third hue.

C_AMNS = "#2a78d6"
C_VND = "#eb6834"
C_TRIV = "#8a8984"
INK = "#0b0b0b"
INK2 = "#52514e"
GRID = "#e4e3df"

METHOD_STYLE = {
    "Trivial": dict(color=C_TRIV, linestyle=":", marker="D", markersize=4,
                    markerfacecolor=C_TRIV, linewidth=1.6),
    "AMNS": dict(color=C_AMNS, linestyle="-", marker="o", markersize=5,
                 markerfacecolor=C_AMNS, linewidth=2.0),
    "VND": dict(color=C_VND, linestyle="-", marker="s", markersize=6,
                markerfacecolor="white", markeredgewidth=1.5, linewidth=2.0),
}

plt.rcParams.update({
    "font.size": 9,
    "axes.titlesize": 10,
    "axes.titleweight": "bold",
    "axes.titlelocation": "left",
    "axes.labelsize": 9,
    "axes.labelcolor": INK2,
    "axes.edgecolor": INK2,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.grid": True,
    "axes.axisbelow": True,
    "grid.color": GRID,
    "grid.linewidth": 0.8,
    "xtick.color": INK2,
    "ytick.color": INK2,
    "legend.frameon": False,
    "figure.facecolor": "white",
})

# Decimal comma on the axes (Portuguese); falls back to the default locale.
for _loc in ("pt_BR.UTF-8", "pt_BR.utf8", "Portuguese_Brazil.1252", "pt-BR"):
    try:
        locale.setlocale(locale.LC_NUMERIC, _loc)
        plt.rcParams["axes.formatter.use_locale"] = True
        break
    except locale.Error:
        pass

# Moves credited with max-load drops (same order in every figure), with the
# names shown in the figures.
OP_ORDER = ["insert", "change", "delete", "extend", "robust", "shake", "kick"]
OP_NAME = {"insert": "inserção", "change": "troca", "delete": "remoção",
           "extend": "extensão", "robust": "robusta", "shake": "perturbação",
           "kick": "kick"}


def br(x: float, d: int) -> str:
    """Number with d decimals and a decimal comma."""
    return f"{x:.{d}f}".replace(".", ",")


def _write(fig, path: Path, dpi: int):
    _savefig(fig, path, dpi)
    plt.close(fig)
    print(f"  wrote {path}")


def _capacity_line(ax, values: List[float]):
    """Dashed reference at load 1, only when the data comes near it."""
    finite = [v for v in values if np.isfinite(v)]
    if finite and max(finite) >= 0.8:
        ax.axhline(1.0, color=INK2, linestyle="--", linewidth=1.0, zorder=1)
        ax.text(1.0, 1.0, "capacidade", transform=ax.get_yaxis_transform(),
                ha="right", va="bottom", fontsize=8, color=INK2)


def _methods(amns: Solution, vnd: Solution, trivial: Optional[Solution]):
    out = []
    if trivial is not None:
        out.append(("Trivial", trivial))
    out += [("AMNS", amns), ("VND", vnd)]
    return out


def op_reduction(sol: Solution) -> Dict[str, float]:
    """Total max-load reduction credited to each move, over all slots."""
    red = {op: 0.0 for op in OP_ORDER}
    for pts in sol.max_load_history:
        for prev, cur in zip(pts, pts[1:]):
            if cur.op in red:
                red[cur.op] += prev.max_load - cur.max_load
    return red


# ---------------------------------------------------------------------------
# Per-instance figures
# ---------------------------------------------------------------------------

def save_overview(name: str, amns: Solution, vnd: Solution,
                  trivial: Optional[Solution], outdir: Path, dpi: int):
    """
    overview.png: six small panels comparing the methods on one instance.
      max link load per slot | mean link load per slot | load distribution
      waypoints per slot     | reduction by move       | solve time
    """
    methods = _methods(amns, vnd, trivial)
    n_slots = max(s.num_slots() for _, s in methods)
    slots = np.arange(n_slots)

    fig, axes = plt.subplots(2, 3, figsize=(15, 8))
    (a_max, a_mean, a_cdf), (a_wp, a_ops, a_time) = axes

    def per_slot(ax, getter, title, ylabel):
        allv = []
        for label, sol in methods:
            v = getter(sol)
            v = v + [np.nan] * (n_slots - len(v))
            allv += v
            ax.plot(slots, v, label=label, **METHOD_STYLE[label])
        ax.set_title(title)
        ax.set_xlabel("período")
        ax.set_ylabel(ylabel)
        ax.set_xticks(slots if n_slots <= 12 else slots[::2])
        return allv

    vals = per_slot(a_max, lambda s: s.max_load_per_slot(),
                    "Carga máxima por período", "carga máxima")
    _capacity_line(a_max, vals)
    a_max.set_ylim(bottom=0)
    per_slot(a_mean, lambda s: s.avg_load_per_slot(),
             "Carga média por período", "carga média")
    a_mean.set_ylim(bottom=0)
    per_slot(a_wp, lambda s: s.total_waypoints_per_slot(),
             "Waypoints em uso por período", "total de waypoints")
    a_wp.set_ylim(bottom=0)

    # The most loaded (slot, link) pairs, ranked: the official objective
    # minimizes these in order (highest first, then the second, ...).
    top_n = 30
    tops = []
    for label, sol in methods:
        loads = np.sort(np.asarray(sol.all_loads()))[::-1][:top_n]
        tops += list(loads)
        st = dict(METHOD_STYLE[label])
        st["markersize"] = st["markersize"] - 1
        a_cdf.plot(np.arange(1, len(loads) + 1), loads, label=label, **st)
    _capacity_line(a_cdf, tops)
    a_cdf.set_title(f"{top_n} enlaces mais carregados (todos os períodos)")
    a_cdf.set_xlabel("posição (1 = enlace mais carregado)")
    a_cdf.set_ylabel("carga do enlace")
    a_cdf.set_ylim(bottom=0)
    a_cdf.set_xticks([1, 5, 10, 15, 20, 25, 30])

    # Share of the max-load reduction credited to each move.
    y = np.arange(len(OP_ORDER))
    h = 0.38
    for k, (label, sol, color, hatch_face) in enumerate(
            [("AMNS", amns, C_AMNS, C_AMNS), ("VND", vnd, C_VND, "white")]):
        red = op_reduction(sol)
        total = sum(red.values())
        share = [100.0 * red[o] / total if total > 1e-12 else 0.0 for o in OP_ORDER]
        a_ops.barh(y + (k - 0.5) * h, share, height=h * 0.9, color=hatch_face,
                   edgecolor=color, linewidth=1.5, label=label)
        for yi, s in zip(y + (k - 0.5) * h, share):
            if s >= 1.0:
                a_ops.text(s, yi, f" {br(s, 0)}%", va="center", fontsize=7, color=INK2)
    a_ops.set_yticks(y)
    a_ops.set_yticklabels([OP_NAME[o] for o in OP_ORDER])
    a_ops.invert_yaxis()
    a_ops.set_xlim(0, 105)
    a_ops.set_title("Redução da carga máxima por movimento")
    a_ops.set_xlabel("parcela da redução total (%)")
    a_ops.grid(axis="y", visible=False)

    # Solve time.
    labels = [m for m, _ in methods]
    times = [s.solve_time_seconds for _, s in methods]
    yy = np.arange(len(labels))
    a_time.barh(yy, times, height=0.55,
                color=[METHOD_STYLE[m]["color"] if m != "VND" else "white" for m in labels],
                edgecolor=[METHOD_STYLE[m]["color"] for m in labels], linewidth=1.5)
    for yi, t in zip(yy, times):
        a_time.text(t, yi, f"  {br(t, 1)} s", va="center", fontsize=8, color=INK)
    a_time.set_yticks(yy)
    a_time.set_yticklabels(labels)
    a_time.invert_yaxis()
    a_time.set_title("Tempo de execução")
    a_time.set_xlabel("segundos")
    a_time.set_xlim(0, max(times + [1e-9]) * 1.25)
    a_time.grid(axis="y", visible=False)

    worst = "   ".join(f"{m} {br(s.max_link_load, 3)}" for m, s in methods)
    fig.suptitle(f"{name}   ·   carga máxima do pior período:   {worst}",
                 x=0.01, ha="left", fontsize=12, fontweight="bold", color=INK)
    handles, leg_labels = a_max.get_legend_handles_labels()
    fig.legend(handles, leg_labels, loc="upper right", ncol=len(methods),
               bbox_to_anchor=(0.99, 0.995), fontsize=9)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    _write(fig, outdir / "overview.png", dpi)


def save_history(name: str, amns: Solution, vnd: Solution,
                 outdir: Path, dpi: int):
    """
    history.png: one small panel per slot, best max link load against the
    search iteration, AMNS and VND. Long flat tails use a symlog x-axis.
    """
    total_slots = max(len(amns.max_load_history), len(vnd.max_load_history))

    # Only slots where either search lowered the max load get a panel.
    def improved(sol, t):
        return t < len(sol.max_load_history) and len(sol.max_load_history[t]) > 1
    active = [t for t in range(total_slots) if improved(amns, t) or improved(vnd, t)]
    flat = [t for t in range(total_slots) if t not in active]
    if not active:
        return
    n_slots = len(active)
    ncols = min(4, n_slots)
    nrows = (n_slots + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=(4.2 * max(ncols, 2), 2.7 * nrows + 0.8),
                             squeeze=False)

    for k in range(nrows * ncols):
        ax = axes[k // ncols][k % ncols]
        if k >= n_slots:
            ax.set_visible(False)
            continue
        slot = active[k]
        finals = []
        long_tail = False
        for label, sol in (("AMNS", amns), ("VND", vnd)):
            if slot >= len(sol.max_load_history) or not sol.max_load_history[slot]:
                continue
            pts = sol.max_load_history[slot]
            sw = sol.nbhd_switches[slot] if slot < len(sol.nbhd_switches) else []
            iters = [p.iteration for p in pts]
            vals = [p.max_load for p in pts]
            last = max([iters[-1]] + [s.iteration for s in sw]) + 1
            long_tail |= iters[-1] < 0.25 * last and last > 200
            st = METHOD_STYLE[label]
            ax.step(iters + [last], vals + [vals[-1]], where="post", label=label,
                    color=st["color"], linewidth=st["linewidth"],
                    linestyle="--" if label == "VND" else "-")
            ax.plot(iters[-1], vals[-1], marker=st["marker"], color=st["color"],
                    markerfacecolor=st["markerfacecolor"], markersize=st["markersize"],
                    markeredgewidth=st.get("markeredgewidth", 1.0), linestyle="none")
            finals.append(f"{label} {br(vals[-1], 3)}")
        if long_tail:
            ax.set_xscale("symlog", linthresh=10)
        ax.set_title(f"período {slot}   " + "  ·  ".join(finals), fontsize=9)
        if k % ncols == 0:
            ax.set_ylabel("melhor carga máxima")
        if k + ncols >= n_slots:
            ax.set_xlabel("iteração" + (" (escala symlog)" if long_tail else ""))

    height = 2.7 * nrows + 0.8
    fig.suptitle(f"{name}   ·   evolução da busca por período", x=0.01, y=0.995, ha="left",
                 va="top", fontsize=12, fontweight="bold", color=INK)
    if flat:
        fig.text(0.01, 1 - 0.42 / height,
                 f"Sem melhora no{'s' if len(flat) > 1 else ''} período{'s' if len(flat) > 1 else ''} "
                 + ", ".join(str(t) for t in flat) + " (não mostrado"
                 + ("s" if len(flat) > 1 else "") + ").",
                 ha="left", va="top", fontsize=9, color=INK2)
    handles = [plt.Line2D([], [], color=C_AMNS, linewidth=2, marker="o", label="AMNS"),
               plt.Line2D([], [], color=C_VND, linewidth=2, linestyle="--", marker="s",
                          markerfacecolor="white", label="VND")]
    fig.legend(handles=handles, loc="lower center", ncol=2, bbox_to_anchor=(0.5, 0.0))
    fig.tight_layout(rect=(0, 0.35 / height, 1, 1 - 0.7 / height))
    _write(fig, outdir / "history.png", dpi)


# ---------------------------------------------------------------------------
# Cross-instance summary
# ---------------------------------------------------------------------------

def save_summary(rows: List[dict], outdir: Path, dpi: int):
    """
    summary-instances.png: one row per instance, three aligned panels
      worst-slot max load (trivial -> AMNS / VND) | change vs trivial | time
    summary-by-set.png: mean change vs trivial per instance set.
    summary.csv: the numbers behind both.
    """
    n = len(rows)
    y = np.arange(n)
    has_triv = [np.isfinite(r["triv"]) for r in rows]

    fig, (a_load, a_change, a_time) = plt.subplots(
        1, 3, figsize=(15, 0.26 * n + 2.0), sharey=True,
        gridspec_kw={"width_ratios": [1.3, 1.0, 0.9]})

    # Panel 1: dumbbell from the trivial value to the solvers' values.
    for yi, r in zip(y, rows):
        lo = min(r["amns"], r["vnd"], r["triv"] if np.isfinite(r["triv"]) else np.inf)
        hi = max(r["amns"], r["vnd"], r["triv"] if np.isfinite(r["triv"]) else -np.inf)
        a_load.plot([lo, hi], [yi, yi], color=GRID, linewidth=3, zorder=1,
                    solid_capstyle="round")
    if any(has_triv):
        a_load.scatter([r["triv"] for r in rows], y, s=26, marker="D", color=C_TRIV,
                       zorder=2, label="Trivial")
    a_load.scatter([r["amns"] for r in rows], y, s=36, marker="o", color=C_AMNS,
                   zorder=3, label="AMNS")
    a_load.scatter([r["vnd"] for r in rows], y, s=46, marker="s", facecolor="white",
                   edgecolor=C_VND, linewidth=1.5, zorder=4, label="VND")
    a_load.axvline(1.0, color=INK2, linestyle="--", linewidth=1.0, zorder=0)
    a_load.text(1.0, 1.0, " capacidade", transform=a_load.get_xaxis_transform(),
                ha="left", va="bottom", fontsize=8, color=INK2)
    a_load.set_title("Carga máxima do pior período")
    a_load.set_xlabel("carga máxima (menor é melhor)")
    a_load.set_xlim(left=0)

    # Panel 2: change vs trivial (%).
    if any(has_triv):
        ch_a = [100 * (r["amns"] - r["triv"]) / r["triv"] if r["triv"] > 0 else np.nan
                for r in rows]
        ch_v = [100 * (r["vnd"] - r["triv"]) / r["triv"] if r["triv"] > 0 else np.nan
                for r in rows]
        for yi, a, v in zip(y, ch_a, ch_v):
            if np.isfinite(a):
                a_change.plot([0, min(a, v)], [yi, yi], color=GRID, linewidth=3, zorder=1,
                              solid_capstyle="round")
        a_change.scatter(ch_a, y, s=36, marker="o", color=C_AMNS, zorder=3)
        a_change.scatter(ch_v, y, s=46, marker="s", facecolor="white", edgecolor=C_VND,
                         linewidth=1.5, zorder=4)
        for yi, a, v in zip(y, ch_a, ch_v):
            if np.isfinite(a):
                a_change.text(min(a, v), yi, f"{br(min(a, v), 0)}%  ", ha="right",
                              va="center", fontsize=7, color=INK2)
        a_change.axvline(0, color=INK2, linewidth=1.0)
        a_change.set_xlim(-110, 10)
    a_change.set_title("Variação em relação à trivial")
    a_change.set_xlabel("variação da carga do pior período, em % (negativo é melhor)")

    # Panel 3: solve time (log scale).
    for key, marker, kw in (("t_triv", "D", dict(color=C_TRIV, s=22)),
                            ("t_amns", "o", dict(color=C_AMNS, s=36)),
                            ("t_vnd", "s", dict(facecolor="white", edgecolor=C_VND,
                                                linewidth=1.5, s=46))):
        xs = [r[key] for r in rows]
        if any(np.isfinite(x) and x > 0 for x in xs):
            a_time.scatter(xs, y, marker=marker, zorder=3, **kw)
    a_time.set_xscale("log")
    a_time.set_title("Tempo de execução")
    a_time.set_xlabel("segundos (escala log)")

    # Rows: instance names, set separators.
    a_load.set_yticks(y)
    a_load.set_yticklabels([r["name"] for r in rows], fontsize=8)
    a_load.invert_yaxis()
    sets = [r["name"].split("-")[0] for r in rows]
    for k in range(1, n):
        if sets[k] != sets[k - 1]:
            for ax in (a_load, a_change, a_time):
                ax.axhline(k - 0.5, color=INK2, linewidth=0.8)
    for ax in (a_load, a_change, a_time):
        ax.grid(axis="y", visible=False)
        ax.set_ylim(n - 0.5, -0.5)

    fig.suptitle("Todas as instâncias: solução trivial vs AMNS vs VND", x=0.01, ha="left",
                 fontsize=12, fontweight="bold", color=INK)
    fig.text(0.01, 1 - 0.55 / (0.26 * n + 2.0),
             "Onde só o quadrado laranja aparece, AMNS e VND atingem o mesmo valor.",
             ha="left", va="top", fontsize=9, color=INK2)
    handles, labels = a_load.get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper right", ncol=3, bbox_to_anchor=(0.99, 0.995))
    fig.tight_layout(rect=(0, 0, 1, 1 - 0.6 / (0.26 * n + 2.0)))
    _write(fig, outdir / "summary-instances.png", dpi)

    # Mean change vs trivial per set.
    if any(has_triv):
        groups, mean_a, mean_v = [], [], []
        for g in sorted(set(sets)) + ["todos"]:
            sel = [r for r, s in zip(rows, sets) if (g == "todos" or s == g)
                   and np.isfinite(r["triv"]) and r["triv"] > 0]
            if not sel:
                continue   # no trivial solution for this set
            groups.append(f"{g} ({len(sel)})")
            mean_a.append(np.mean([100 * (r["amns"] - r["triv"]) / r["triv"] for r in sel]))
            mean_v.append(np.mean([100 * (r["vnd"] - r["triv"]) / r["triv"] for r in sel]))
        fig, ax = plt.subplots(figsize=(7, 3.6))
        x = np.arange(len(groups))
        w = 0.36
        ax.bar(x - w / 2, mean_a, width=w * 0.92, color=C_AMNS, label="AMNS")
        ax.bar(x + w / 2, mean_v, width=w * 0.92, color="white", edgecolor=C_VND,
               linewidth=1.5, label="VND")
        for xi, v in list(zip(x - w / 2, mean_a)) + list(zip(x + w / 2, mean_v)):
            ax.text(xi, v, f"{br(v, 0)}%", ha="center", va="top", fontsize=8, color=INK,
                    bbox=dict(facecolor="white", edgecolor="none", pad=1))
        ax.axhline(0, color=INK2, linewidth=1.0)
        ax.set_xticks(x)
        ax.set_xticklabels(groups)
        ax.set_ylabel("variação média vs trivial (%)")
        ax.set_title("Variação média da carga do pior período vs trivial, por conjunto")
        ax.grid(axis="x", visible=False)
        ax.legend(loc="lower right")
        fig.tight_layout()
        _write(fig, outdir / "summary-by-set.png", dpi)

    # Table view.
    csv = outdir / "summary.csv"
    lines = ["instance,trivial_max,amns_max,vnd_max,amns_vs_trivial_pct,"
             "vnd_vs_trivial_pct,amns_overflow,vnd_overflow,trivial_s,amns_s,vnd_s"]
    for r in rows:
        def pct(v):
            return f"{100 * (v - r['triv']) / r['triv']:.2f}" if r["triv"] > 0 else ""
        lines.append(",".join([
            r["name"], f"{r['triv']:.6f}" if np.isfinite(r["triv"]) else "",
            f"{r['amns']:.6f}", f"{r['vnd']:.6f}", pct(r["amns"]), pct(r["vnd"]),
            f"{r['ov_amns']:.6f}", f"{r['ov_vnd']:.6f}",
            f"{r['t_triv']:.3f}" if np.isfinite(r["t_triv"]) else "",
            f"{r['t_amns']:.3f}", f"{r['t_vnd']:.3f}"]))
    for k in range(6):
        try:
            csv.write_text("\n".join(lines) + "\n")
            break
        except OSError:
            if k == 5:
                raise
            time.sleep(0.5 * (k + 1))
    print(f"  wrote {csv}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", type=Path, default=Path("."),
                    help="folder containing amns_solution/ and vnd_solution/ "
                         "(and optionally trivial_solution/)")
    ap.add_argument("--outdir", type=Path, default=Path("figures"),
                    help="where to save the figures (one subfolder per instance)")
    ap.add_argument("--no-history", action="store_true",
                    help="skip the per-slot search history figures")
    ap.add_argument("--show", action="store_true",
                    help="also call plt.show() at the end")
    ap.add_argument("--dpi", type=int, default=150,
                    help="DPI for saved figures")
    args = ap.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)

    instances = discover_instances(args.root)
    if not instances:
        print("Nothing to plot.", file=sys.stderr)
        return 1
    print(f"Found {len(instances)} instance(s) with both solutions.")

    rows = []
    for name in instances:
        amns = load_solution(args.root / "amns_solution" / f"{name}-solution.json")
        vnd = load_solution(args.root / "vnd_solution" / f"{name}-solution.json")
        t_path = args.root / "trivial_solution" / f"{name}-solution.json"
        trivial = load_solution(t_path) if t_path.is_file() else None

        inst_dir = args.outdir / name
        inst_dir.mkdir(parents=True, exist_ok=True)
        save_overview(name, amns, vnd, trivial, inst_dir, args.dpi)
        if not args.no_history:
            save_history(name, amns, vnd, inst_dir, args.dpi)

        rows.append(dict(
            name=name,
            triv=trivial.max_link_load if trivial else np.nan,
            amns=amns.max_link_load, vnd=vnd.max_link_load,
            ov_amns=sum(amns.overflow_per_slot()), ov_vnd=sum(vnd.overflow_per_slot()),
            t_triv=trivial.solve_time_seconds if trivial else np.nan,
            t_amns=amns.solve_time_seconds, t_vnd=vnd.solve_time_seconds,
        ))

    save_summary(rows, args.outdir, args.dpi)

    print(f"\n{'instance':<12} {'trivial':>9} {'AMNS':>9} {'VND':>9} "
          f"{'AMNS/Tr':>8} {'VND/Tr':>8}")
    for r in rows:
        def pct(v):
            return f"{100 * (v - r['triv']) / r['triv']:+.1f}%" if r["triv"] > 0 else "n/a"
        print(f"{r['name']:<12} {r['triv']:>9.4f} {r['amns']:>9.4f} {r['vnd']:>9.4f} "
              f"{pct(r['amns']):>8} {pct(r['vnd']):>8}")

    if args.show:
        plt.show()
    return 0


if __name__ == "__main__":
    sys.exit(main())
