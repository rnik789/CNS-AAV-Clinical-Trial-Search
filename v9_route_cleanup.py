#!/usr/bin/env python3
"""
Remove "Explicit CNS delivery" from cns_route_categories in v9 CNS CSV.
Saves cleaned file and regenerates route_breakdown.csv / route_breakdown.png
in cns_landscape_analysis/.
"""

import os
import csv
from collections import Counter

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

WORK_DIR  = r"D:\ClaudeCode\FDA_AAV_Search_061826"
OUT_DIR   = os.path.join(WORK_DIR, "cns_landscape_analysis")
IN_CSV    = os.path.join(WORK_DIR, "aav_cns_confirmed_v9_2026-06-18.csv")
OUT_CSV   = os.path.join(WORK_DIR, "aav_cns_confirmed_v9_2026-06-18_no_explicit_route.csv")
EXPECTED_N = 92
BAR_COLOR  = "#2C7BB6"
REMOVE_TERM = "Explicit CNS delivery"

os.makedirs(OUT_DIR, exist_ok=True)

# ── Load ───────────────────────────────────────────────────────────────────────
with open(IN_CSV, encoding="utf-8", newline="") as f:
    reader = csv.DictReader(f)
    fieldnames = reader.fieldnames
    rows = list(reader)

assert len(rows) == EXPECTED_N, f"Expected {EXPECTED_N} rows, got {len(rows)}"
print(f"Loaded {len(rows)} rows from {os.path.basename(IN_CSV)}")


def strip_term(route_str, term):
    """Remove `term` from a pipe-delimited route string; collapse empty segments."""
    parts = [p.strip() for p in route_str.split("|") if p.strip() != term]
    return "|".join(parts)


# ── Edit cns_route_categories ──────────────────────────────────────────────────
affected = []
for row in rows:
    old = row["cns_route_categories"]
    if REMOVE_TERM in old:
        new = strip_term(old, REMOVE_TERM)
        affected.append((row["nct_id"], old, new))
        row["cns_route_categories"] = new

# ── Before/after report ────────────────────────────────────────────────────────
print(f"\nAffected rows ({len(affected)}):")
print(f"{'nct_id':<15} {'old cns_route_categories':<65} {'new cns_route_categories'}")
print("-" * 130)
for nct, old, new in affected:
    print(f"{nct:<15} {old:<65} {new}")

# ── Save cleaned CSV ───────────────────────────────────────────────────────────
with open(OUT_CSV, "w", encoding="utf-8", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=fieldnames)
    writer.writeheader()
    writer.writerows(rows)
print(f"\nSaved: {os.path.basename(OUT_CSV)}  ({len(rows)} rows, columns unchanged)")

# ── Regenerate route breakdown ─────────────────────────────────────────────────
route_counter: Counter = Counter()
for row in rows:
    val = row["cns_route_categories"]
    if not val or not val.strip():
        route_counter["NOT REPORTED"] += 1
    else:
        for seg in val.split("|"):
            seg = seg.strip()
            if seg:
                route_counter[seg] += 1

route_items = sorted(route_counter.items(), key=lambda x: x[1], reverse=True)

# CSV
route_csv = os.path.join(OUT_DIR, "route_breakdown.csv")
with open(route_csv, "w", encoding="utf-8", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=["route", "trial_count", "percent_of_92_trials"])
    writer.writeheader()
    for route, cnt in route_items:
        writer.writerow({"route": route, "trial_count": cnt,
                         "percent_of_92_trials": round(cnt / EXPECTED_N * 100, 1)})

total = sum(c for _, c in route_items)
print(f"\nRoute breakdown (total occurrences: {total}, exceeds 92 due to multi-route trials):")
print(f"{'route':<35} {'trial_count':>12} {'pct_of_92':>10}")
print("-" * 60)
for route, cnt in route_items:
    print(f"{route:<35} {cnt:>12} {cnt/EXPECTED_N*100:>9.1f}%")

# Chart
route_chart = os.path.join(OUT_DIR, "route_breakdown.png")
labels = [r for r, _ in route_items]
vals   = [c for _, c in route_items]
labels_plot = labels[::-1]
vals_plot   = vals[::-1]

fig, ax = plt.subplots(figsize=(8, 5))
bars = ax.barh(labels_plot, vals_plot, color=BAR_COLOR)
x_max = max(vals_plot) if vals_plot else 1
for bar, val in zip(bars, vals_plot):
    ax.text(bar.get_width() + x_max * 0.01,
            bar.get_y() + bar.get_height() / 2,
            str(val), va="center", fontsize=9)
ax.set_xlabel("Number of Trials")
ax.set_title("CNS Delivery Routes Among CNS-Targeted AAV Trials\n(n=92, routes not mutually exclusive)")
ax.set_xlim(0, x_max * 1.15)
plt.tight_layout()
plt.savefig(route_chart, dpi=150)
plt.close()
print(f"\nSaved: {route_csv}")
print(f"Saved: {route_chart}")
