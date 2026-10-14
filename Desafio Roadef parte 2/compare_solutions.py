#!/usr/bin/env python3
"""
compare_solutions.py

Read the AMNS and VND solution JSON files and produce comparison plots.
Each instance's plots go into its own subfolder under --outdir.

Expected layout (matching the C++ exporters):

    amns_solution/setA-01-solution.json
    amns_solution/setA-02-solution.json
    ...
    vnd_solution/setA-01-solution.json
    vnd_solution/setA-02-solution.json
    ...

Output layout:

    plots/
    ├── setA-01/
    │   ├── max_load_per_slot.png
    │   ├── max_load_difference.png
    │   ├── avg_load_per_slot.png
    │   ├── overflow_per_slot.png
    │   ├── waypoints_per_slot.png
    │   ├── load_distribution_amns.png
    │   ├── load_distribution_vnd.png
    │   ├── max_load_history_amns.png
    │   ├── max_load_history_vnd.png
    │   └── solve_time.png
    ├── setA-02/
    │   └── ...
    ├── setB-01/
    │   └── ...
    ├── summary-all-instances.png
    └── summary-solve-time.png

Usage:
    python compare_solutions.py
    python compare_solutions.py --root path/to/dir
    python compare_solutions.py --outdir plots
"""

import argparse
import json
import sys
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
            )
            for e in slot_arr
        ]
        history.append(slot_pts)

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

def _finalize(fig, ax, outdir: Path, fname: str, dpi: int,
              legend: bool = True):
    if legend:
        ax.legend()
    fig.tight_layout()
    out = outdir / fname
    fig.savefig(out, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {out}")


# ---------------------------------------------------------------------------
# Per-slot comparison plots
# ---------------------------------------------------------------------------

def save_max_load_history_slot(sol: Solution,
                               slot: int,
                               method_label: str,
                               color: str,
                               outdir: Path, dpi: int):
    """
    Plot the moments where the accepted solution's max link load dropped
    within a single slot, as a step plot vs. the iteration counter inside
    that slot.
    """
    if slot >= len(sol.max_load_history):
        return
    pts = sol.max_load_history[slot]
    if not pts:
        # Nothing to draw; skip silently.
        return

    iters = [p.iteration for p in pts]
    vals  = [p.max_load  for p in pts]

    fig, ax = plt.subplots(figsize=(9, 5))
    ax.step(iters, vals, where="post", color=color, linewidth=1.5,
            label=f"{method_label} max load")
    ax.scatter(iters, vals, s=20, color=color, zorder=3)

    ax.annotate(f"start = {vals[0]:.4f}",
                xy=(iters[0], vals[0]),
                xytext=(5, 5), textcoords="offset points",
                fontsize=9, color=color)
    ax.annotate(f"final = {vals[-1]:.4f}",
                xy=(iters[-1], vals[-1]),
                xytext=(5, -15), textcoords="offset points",
                fontsize=9, color=color)

    ax.axhline(1.0, color="red", linestyle=":", linewidth=1.0,
               label="capacity limit")

    ax.set_title(f"{sol.instance_name}: max link load history — "
                 f"{method_label}, slot {slot}")
    ax.set_xlabel("iteration within slot")
    ax.set_ylabel("max link load")
    ax.grid(True, alpha=0.3)
    ax.legend(loc="upper right")

    out = outdir / f"max_load_history_{method_label.lower()}.png"
    fig.tight_layout()
    fig.savefig(out, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {out}")

def save_max_load_per_slot(amns: Solution, vnd: Solution,
                           outdir: Path, dpi: int):
    """
    Paired vertical bars per time slot: AMNS (blue) next to VND (orange).
    The y-axis is the max link load in that slot.
    """
    n = max(amns.num_slots(), vnd.num_slots())
    a = amns.max_load_per_slot()
    v = vnd.max_load_per_slot()
    a = a + [np.nan] * (n - len(a))
    v = v + [np.nan] * (n - len(v))

    x = np.arange(n)
    w = 0.40

    fig, ax = plt.subplots(figsize=(max(8, 0.7 * n + 2), 5))
    ax.bar(x - w/2, a, width=w, label="AMNS", color="C0", alpha=0.9)
    ax.bar(x + w/2, v, width=w, label="VND",  color="C1", alpha=0.9)

    ax.axhline(1.0, color="red", linestyle=":", linewidth=1.0,
               label="capacity limit")

    if n <= 30:
        for xi, val in zip(x - w/2, a):
            if np.isfinite(val):
                ax.text(xi, val, f"{val:.3f}", ha="center", va="bottom",
                        fontsize=7, color="C0")
        for xi, val in zip(x + w/2, v):
            if np.isfinite(val):
                ax.text(xi, val, f"{val:.3f}", ha="center", va="bottom",
                        fontsize=7, color="C1")

    ax.set_title(f"{amns.instance_name}: max link load per slot "
                 "(AMNS vs VND)")
    ax.set_xlabel("time slot")
    ax.set_ylabel("max link load")
    ax.set_xticks(x)
    ax.grid(True, axis="y", alpha=0.3)
    ax.legend(loc="best")

    _finalize(fig, ax, outdir, "max_load_per_slot.png", dpi,
              legend=False)


def save_max_load_difference(amns: Solution, vnd: Solution,
                             outdir: Path, dpi: int):
    """
    Difference of max load per slot: AMNS - VND.
    Negative => AMNS is better (lower max load).
    Positive => VND  is better.
    """
    n = max(amns.num_slots(), vnd.num_slots())
    a = amns.max_load_per_slot() + [np.nan] * (n - amns.num_slots())
    v = vnd.max_load_per_slot()  + [np.nan] * (n - vnd.num_slots())
    diff = np.array(a, dtype=float) - np.array(v, dtype=float)

    fig, ax = plt.subplots(figsize=(max(8, 0.7 * n + 2), 5))
    colors = ["C0" if d < 0 else ("C1" if d > 0 else "gray")
              for d in np.nan_to_num(diff)]
    ax.bar(range(n), diff, color=colors, alpha=0.85)
    ax.axhline(0.0, color="black", linewidth=0.8)
    ax.set_xlabel("time slot")
    ax.set_ylabel("AMNS max load - VND max load")
    ax.set_xticks(range(n))

    from matplotlib.patches import Patch
    handles = [
        Patch(color="C0", label="AMNS better (lower max load)"),
        Patch(color="C1", label="VND better (lower max load)"),
        Patch(color="gray", label="tie"),
    ]
    ax.legend(handles=handles, loc="best", fontsize=9)

    ax.set_title(f"{amns.instance_name}: max load difference per slot "
                 "(AMNS - VND)")
    ax.grid(True, axis="y", alpha=0.3)

    _finalize(fig, ax, outdir, "max_load_difference.png", dpi,
              legend=False)


def save_avg_load_per_slot(amns: Solution, vnd: Solution,
                           outdir: Path, dpi: int):
    fig, ax = plt.subplots(figsize=(8, 5))
    t_a = list(range(amns.num_slots()))
    t_v = list(range(vnd.num_slots()))
    ax.plot(t_a, amns.avg_load_per_slot(), marker="o", label="AMNS")
    ax.plot(t_v, vnd.avg_load_per_slot(),  marker="s", ls="--", label="VND")
    ax.set_title(f"{amns.instance_name}: average link load per slot")
    ax.set_xlabel("time slot")
    ax.set_ylabel("average link load")
    ax.grid(True, alpha=0.3)
    _finalize(fig, ax, outdir, "avg_load_per_slot.png", dpi)


def save_overflow_per_slot(amns: Solution, vnd: Solution,
                           outdir: Path, dpi: int):
    fig, ax = plt.subplots(figsize=(8, 5))
    t_a = list(range(amns.num_slots()))
    t_v = list(range(vnd.num_slots()))
    ax.plot(t_a, amns.overflow_per_slot(), marker="o", label="AMNS")
    ax.plot(t_v, vnd.overflow_per_slot(),  marker="s", ls="--", label="VND")
    ax.set_title(f"{amns.instance_name}: overflow per slot "
                 "(sum of load-1 over congested links)")
    ax.set_xlabel("time slot")
    ax.set_ylabel("overflow")
    ax.grid(True, alpha=0.3)
    _finalize(fig, ax, outdir, "overflow_per_slot.png", dpi)


def save_waypoints_per_slot(amns: Solution, vnd: Solution,
                            outdir: Path, dpi: int):
    fig, ax = plt.subplots(figsize=(8, 5))
    t_a = list(range(amns.num_slots()))
    t_v = list(range(vnd.num_slots()))
    ax.plot(t_a, amns.total_waypoints_per_slot(), marker="o", label="AMNS")
    ax.plot(t_v, vnd.total_waypoints_per_slot(),  marker="s", ls="--",
            label="VND")
    ax.set_title(f"{amns.instance_name}: waypoints per slot "
                 "(sum over demands)")
    ax.set_xlabel("time slot")
    ax.set_ylabel("# waypoints")
    ax.grid(True, alpha=0.3)
    _finalize(fig, ax, outdir, "waypoints_per_slot.png", dpi)


# ---------------------------------------------------------------------------
# One load-distribution figure per method
# ---------------------------------------------------------------------------

def save_load_distribution_one(sol: Solution,
                               method_label: str,
                               color: str,
                               outdir: Path, dpi: int):
    loads = sol.all_loads()
    max_reported = sol.max_link_load
    max_emp = max(loads) if loads else 0.0

    hi = max(1.2, (max(loads + [0.0]) * 1.05) if loads else 1.2)
    bins = np.linspace(0.0, hi, 40)

    fig, ax = plt.subplots(figsize=(9, 5))
    ax.hist(loads, bins=bins, alpha=0.85, color=color, label=method_label)
    ax.axvline(1.0, color="red", linestyle=":", linewidth=1.0,
               label="capacity limit")
    ax.axvline(max_reported, color="black", linestyle="--", linewidth=1.0,
               label=f"max load = {max_reported:.4f}")

    ymax = ax.get_ylim()[1]
    ax.annotate(f"max = {max_reported:.4f}\n"
                f"(link {sol.max_link_id}, slot {sol.max_link_slot})",
                xy=(max_reported, ymax * 0.90),
                xytext=(5, 0), textcoords="offset points",
                color="black", fontsize=9, va="top")

    info = (f"{method_label} loads: n={len(loads)}\n"
            f"reported max = {max_reported:.4f}\n"
            f"empirical max = {max_emp:.4f}\n"
            f"mean load = {np.mean(loads) if loads else 0.0:.4f}")
    ax.text(0.99, 0.99, info, transform=ax.transAxes,
            ha="right", va="top", fontsize=8,
            bbox=dict(boxstyle="round", facecolor="white", alpha=0.85))

    ax.set_title(f"{sol.instance_name}: link load distribution — "
                 f"{method_label} (all slots)")
    ax.set_xlabel("link load")
    ax.set_ylabel("# (slot, link) samples")
    ax.grid(True, alpha=0.3)
    ax.legend(loc="upper right")

    out = outdir / f"load_distribution_{method_label.lower()}.png"
    fig.tight_layout()
    fig.savefig(out, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {out}")


# ---------------------------------------------------------------------------
# Max-load history (one figure per method)
# ---------------------------------------------------------------------------

def save_max_load_history(sol: Solution,
                          method_label: str,
                          color: str,
                          outdir: Path, dpi: int):
    """
    Plot the moments where the accepted solution's max link load dropped,
    as a step plot vs. the global iteration counter.
    """
    pts = sol.max_load_history
    if not pts:
        print(f"  [skip] {sol.instance_name}: no max_load_history for "
              f"{method_label}")
        return

    iters = [p.iteration for p in pts]
    vals  = [p.max_load  for p in pts]

    fig, ax = plt.subplots(figsize=(9, 5))
    # Step plot: each new value holds until the next drop.
    ax.step(iters, vals, where="post", color=color, linewidth=1.5,
            label=f"{method_label} max load")
    ax.scatter(iters, vals, s=18, color=color, zorder=3)

    # Annotate the final value.
    ax.annotate(f"final = {vals[-1]:.4f}",
                xy=(iters[-1], vals[-1]),
                xytext=(5, 5), textcoords="offset points",
                fontsize=9, color=color)

    ax.axhline(1.0, color="red", linestyle=":", linewidth=1.0,
               label="capacity limit")

    ax.set_title(f"{sol.instance_name}: max link load history — "
                 f"{method_label}")
    ax.set_xlabel("global iteration (across all slots)")
    ax.set_ylabel("max link load")
    ax.grid(True, alpha=0.3)
    ax.legend(loc="upper right")

    out = outdir / f"max_load_history_{method_label.lower()}.png"
    fig.tight_layout()
    fig.savefig(out, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {out}")


# ---------------------------------------------------------------------------
# Solve-time plots
# ---------------------------------------------------------------------------

def save_solve_time(amns: Solution, vnd: Solution,
                    outdir: Path, dpi: int):
    """
    Grouped vertical bars: AMNS vs VND solve time (seconds) for this instance.
    """
    labels = ["AMNS", "VND"]
    times  = [amns.solve_time_seconds, vnd.solve_time_seconds]

    fig, ax = plt.subplots(figsize=(6, 5))
    bars = ax.bar(labels, times, color=["C0", "C1"], alpha=0.9)
    ax.set_title(f"{amns.instance_name}: solve time (seconds)")
    ax.set_ylabel("seconds")
    ax.grid(True, axis="y", alpha=0.3)

    for b, t in zip(bars, times):
        ax.text(b.get_x() + b.get_width() / 2, t,
                f"{t:.3f}s", ha="center", va="bottom", fontsize=9)

    out = outdir / "solve_time.png"
    fig.tight_layout()
    fig.savefig(out, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {out}")


# ---------------------------------------------------------------------------
# Cross-instance summary
# ---------------------------------------------------------------------------

def save_instance_summary(instances: List[str],
                          amns_max: List[float], vnd_max: List[float],
                          amns_ov: List[float],  vnd_ov: List[float],
                          outdir: Path, dpi: int):
    fig, ax = plt.subplots(figsize=(max(6, 1.4 * len(instances)), 5))
    x = np.arange(len(instances))
    w = 0.35

    ax2 = ax.twinx()
    ax.bar(x - w/2, amns_max, width=w, label="AMNS max load",
           color="C0", alpha=0.85)
    ax.bar(x + w/2, vnd_max,  width=w, label="VND max load",
           color="C1", alpha=0.85)
    ax2.plot(x, amns_ov, marker="o", linestyle="--", color="C2",
             label="AMNS overflow")
    ax2.plot(x, vnd_ov,  marker="s", linestyle="--", color="C3",
             label="VND overflow")

    ax.set_xticks(x)
    ax.set_xticklabels(instances, rotation=30, ha="right")
    ax.set_ylabel("max link load")
    ax2.set_ylabel("overflow (sum of load-1)")
    ax.set_title("Per-instance summary: AMNS vs VND")
    ax.grid(True, alpha=0.3)

    h1, l1 = ax.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    ax.legend(h1 + h2, l1 + l2, loc="upper left", fontsize=8)

    out = outdir / "summary-all-instances.png"
    fig.tight_layout()
    fig.savefig(out, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {out}")


def save_time_summary(instances: List[str],
                      amns_t: List[float], vnd_t: List[float],
                      outdir: Path, dpi: int):
    """
    Grouped bars per instance: AMNS solve time vs VND solve time (seconds).
    """
    x = np.arange(len(instances))
    w = 0.35

    fig, ax = plt.subplots(figsize=(max(6, 1.4 * len(instances)), 5))
    ax.bar(x - w/2, amns_t, width=w, label="AMNS", color="C0", alpha=0.9)
    ax.bar(x + w/2, vnd_t,  width=w, label="VND",  color="C1", alpha=0.9)

    ax.set_xticks(x)
    ax.set_xticklabels(instances, rotation=30, ha="right")
    ax.set_ylabel("solve time (seconds)")
    ax.set_title("Solve time per instance: AMNS vs VND")
    ax.grid(True, axis="y", alpha=0.3)
    ax.legend(loc="best")

    for xi, t in zip(x - w/2, amns_t):
        ax.text(xi, t, f"{t:.2f}", ha="center", va="bottom",
                fontsize=7, color="C0")
    for xi, t in zip(x + w/2, vnd_t):
        ax.text(xi, t, f"{t:.2f}", ha="center", va="bottom",
                fontsize=7, color="C1")

    out = outdir / "summary-solve-time.png"
    fig.tight_layout()
    fig.savefig(out, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {out}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", type=Path, default=Path("."),
                    help="folder containing amns_solution/ and vnd_solution/")
    ap.add_argument("--outdir", type=Path, default=Path("plots"),
                    help="where to save the generated plots "
                         "(one subfolder per instance)")
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

    print(f"Found {len(instances)} instance(s) with both solutions:")
    for name in instances:
        print(f"  - {name}")

    summary_amns_max: List[float] = []
    summary_vnd_max:  List[float] = []
    summary_amns_ov:  List[float] = []
    summary_vnd_ov:   List[float] = []
    summary_amns_time: List[float] = []
    summary_vnd_time:  List[float] = []

    for name in instances:
        a_path = args.root / "amns_solution" / f"{name}-solution.json"
        v_path = args.root / "vnd_solution"  / f"{name}-solution.json"

        amns = load_solution(a_path)
        vnd  = load_solution(v_path)

        summary_amns_max.append(amns.max_link_load)
        summary_vnd_max.append(vnd.max_link_load)
        summary_amns_ov.append(sum(amns.overflow_per_slot()))
        summary_vnd_ov.append(sum(vnd.overflow_per_slot()))
        summary_amns_time.append(amns.solve_time_seconds)
        summary_vnd_time.append(vnd.solve_time_seconds)

        # Per-instance output folder.
        inst_dir = args.outdir / name
        inst_dir.mkdir(parents=True, exist_ok=True)
        

        # Per-slot comparison plots.
        save_max_load_per_slot(amns, vnd, inst_dir, args.dpi)
        save_max_load_difference(amns, vnd, inst_dir, args.dpi)
        save_avg_load_per_slot(amns, vnd, inst_dir, args.dpi)
        save_overflow_per_slot(amns, vnd, inst_dir, args.dpi)
        save_waypoints_per_slot(amns, vnd, inst_dir, args.dpi)
        save_solve_time(amns, vnd, inst_dir, args.dpi)

        # One load-distribution figure per method.
        save_load_distribution_one(amns, "AMNS", "C0",
                                   inst_dir, args.dpi)
        save_load_distribution_one(vnd,  "VND",  "C1",
                                   inst_dir, args.dpi)

        # Max-load history, one figure per method.
                # Max-load history, one figure per method PER SLOT, each in its own
        # Load_history_slotX/ folder under the instance directory.
        n_slots = max(amns.num_slots(), vnd.num_slots())
        for slot in range(n_slots):
            slot_dir = inst_dir / f"Load_history_slot{slot}"
            slot_dir.mkdir(parents=True, exist_ok=True)
            save_max_load_history_slot(amns, slot, "AMNS", "C0",
                                       slot_dir, args.dpi)
            save_max_load_history_slot(vnd,  slot, "VND",  "C1",
                                       slot_dir, args.dpi)

    save_instance_summary(
        instances,
        summary_amns_max, summary_vnd_max,
        summary_amns_ov,  summary_vnd_ov,
        args.outdir, args.dpi,
    )
    save_time_summary(
        instances,
        summary_amns_time, summary_vnd_time,
        args.outdir, args.dpi,
    )

    print("\nPer-instance totals")
    print(f"{'instance':<20} {'AMNS max':>10} {'VND max':>10} "
          f"{'AMNS ov':>10} {'VND ov':>10} "
          f"{'AMNS s':>8} {'VND s':>8}")
    for i, name in enumerate(instances):
        print(f"{name:<20} {summary_amns_max[i]:>10.6f} "
              f"{summary_vnd_max[i]:>10.6f} "
              f"{summary_amns_ov[i]:>10.6f} "
              f"{summary_vnd_ov[i]:>10.6f} "
              f"{summary_amns_time[i]:>8.3f} "
              f"{summary_vnd_time[i]:>8.3f}")

    if args.show:
        plt.show()

    return 0


if __name__ == "__main__":
    sys.exit(main())