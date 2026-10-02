#!/usr/bin/env python3
"""Standalone plotter for aoa_sweep_results.csv (does not touch the solver scripts).
Usage:  python plot_results.py [path/to/aoa_sweep_results.csv] [output.png]
Draws: Cl vs AoA, Cd vs AoA, drag polar (Cl vs Cd), Cl/Cd vs AoA."""
import csv, sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

csv_path = sys.argv[1] if len(sys.argv) > 1 else "aoa_sweep_results.csv"
out_path = sys.argv[2] if len(sys.argv) > 2 else "sweep_plots.png"

rows = []
with open(csv_path, newline="") as f:
    for r in csv.DictReader(f):
        if r["Cl"] == "" or r["Cd"] == "":
            continue
        rows.append((float(r["AoA_deg"]), float(r["Cl"]), float(r["Cd"]),
                     r.get("converged", "True").strip().lower() == "true"))
rows.sort()
a = [r[0] for r in rows]; cl = [r[1] for r in rows]; cd = [r[2] for r in rows]
ld = [c / d if d else float("nan") for c, d in zip(cl, cd)]
bad = [i for i, r in enumerate(rows) if not r[3]]

fig, ax = plt.subplots(2, 2, figsize=(11, 8))
def draw(axis, x, y, xl, yl, title, color):
    axis.plot(x, y, marker="o", color=color)
    if bad:
        axis.plot([x[i] for i in bad], [y[i] for i in bad], "o", mfc="none", mec="red", ms=11,
                  label="not converged")
        axis.legend()
    axis.set_xlabel(xl); axis.set_ylabel(yl); axis.set_title(title); axis.grid(True, alpha=0.3)

draw(ax[0][0], a, cl, "Angle of attack (deg)", "Cl", "Lift coefficient vs AoA", "tab:blue")
draw(ax[0][1], a, cd, "Angle of attack (deg)", "Cd", "Drag coefficient vs AoA", "tab:orange")
draw(ax[1][0], cd, cl, "Cd", "Cl", "Drag polar", "tab:green")
draw(ax[1][1], a, ld, "Angle of attack (deg)", "Cl / Cd", "Lift-to-drag ratio vs AoA", "tab:red")
fig.suptitle("NACA 0012 - OpenFOAM sweep", fontsize=13)
fig.tight_layout()
fig.savefig(out_path, dpi=150)
print("Saved:", out_path)
