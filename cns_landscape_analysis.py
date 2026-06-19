#!/usr/bin/env python3
"""
CNS AAV Trial Landscape Analysis
Inputs:  aav_cns_confirmed_v9_2026-06-18.csv (92 rows)
         aav_serotype_frequency_v9_2026-06-18.csv
Outputs: cns_landscape_analysis/ (CSVs + PNGs)
"""

import os
import sys
import csv
import re
from collections import Counter, defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

# ── Paths ──────────────────────────────────────────────────────────────────────
WORK_DIR = r"D:\ClaudeCode\FDA_AAV_Search_061826"
OUT_DIR  = os.path.join(WORK_DIR, "cns_landscape_analysis")
CNS_CSV  = os.path.join(WORK_DIR, "aav_cns_confirmed_v9_2026-06-18.csv")
SER_CSV  = os.path.join(WORK_DIR, "aav_serotype_frequency_v9_2026-06-18.csv")

EXPECTED_N = 92
BAR_COLOR  = "#2C7BB6"

os.makedirs(OUT_DIR, exist_ok=True)

# ── Load & Validate ────────────────────────────────────────────────────────────
df = pd.read_csv(CNS_CSV, dtype=str, keep_default_na=False)
actual_n = len(df)
if actual_n != EXPECTED_N:
    print(f"STOP: expected {EXPECTED_N} rows, found {actual_n}.")
    sys.exit(1)
print(f"Loaded {actual_n} rows from {os.path.basename(CNS_CSV)}")

if "sponsor" not in df.columns:
    print("STOP: 'sponsor' column not found. Actual columns:")
    for col in df.columns:
        print(f"  {col}")
    sys.exit(1)
print("'sponsor' column confirmed.\n")

# ── Helpers ────────────────────────────────────────────────────────────────────
def save_csv(records, path):
    if not records:
        return
    with open(path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(records[0].keys()))
        writer.writeheader()
        writer.writerows(records)


def hbar(series, chart_path, title, xlabel="Number of Trials"):
    """Horizontal bar chart sorted descending with count labels at bar ends."""
    pairs = sorted(zip(series.values, series.index), reverse=True)
    vals   = [p[0] for p in pairs]
    labels = [p[1] for p in pairs]
    # Reverse so highest is at top (barh plots bottom=>top)
    vals_plot   = vals[::-1]
    labels_plot = labels[::-1]

    fig, ax = plt.subplots(figsize=(8, 5))
    bars = ax.barh(labels_plot, vals_plot, color=BAR_COLOR)
    x_max = max(vals_plot) if vals_plot else 1
    for bar, val in zip(bars, vals_plot):
        ax.text(
            bar.get_width() + x_max * 0.01,
            bar.get_y() + bar.get_height() / 2,
            str(val), va="center", fontsize=9
        )
    ax.set_xlabel(xlabel)
    ax.set_title(title)
    ax.set_xlim(0, x_max * 1.15)
    plt.tight_layout()
    plt.savefig(chart_path, dpi=150)
    plt.close()


summary_rows = []


# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 1 — STATUS
# ═══════════════════════════════════════════════════════════════════════════════
print("=" * 65)
print("SECTION 1: STATUS")
print("=" * 65)

status_counts = df["overall_status"].value_counts().sort_values(ascending=False)
status_records = [
    {"status": s, "count": int(c), "percent_of_total": round(c / EXPECTED_N * 100, 1)}
    for s, c in status_counts.items()
]
csv_path = os.path.join(OUT_DIR, "status_breakdown.csv")
save_csv(status_records, csv_path)

print(f"{'status':<30} {'count':>6} {'percent':>9}")
for r in status_records:
    print(f"{r['status']:<30} {r['count']:>6} {r['percent_of_total']:>8.1f}%")

chart_path = os.path.join(OUT_DIR, "status_breakdown.png")
hbar(status_counts, chart_path, "CNS-Targeted AAV Trials by Status (n=92)")
summary_rows.append({"section": "1. Status", "csv": "status_breakdown.csv",
                     "chart": "status_breakdown.png", "rows": len(status_records)})


# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 2 — PHASE
# ═══════════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 65)
print("SECTION 2: PHASE")
print("=" * 65)
print("Note: combined phases (e.g. 'PHASE1; PHASE2') treated as a single"
      " category; 'NA'/'N/A'/blank => 'NOT REPORTED'")


def normalize_phase(val):
    if not val or val.strip() in ("NA", "N/A", ""):
        return "NOT REPORTED"
    return val.strip()


df["phase_norm"] = df["phase"].apply(normalize_phase)
phase_counts = df["phase_norm"].value_counts().sort_values(ascending=False)
phase_records = [
    {"phase": p, "count": int(c), "percent_of_total": round(c / EXPECTED_N * 100, 1)}
    for p, c in phase_counts.items()
]
csv_path = os.path.join(OUT_DIR, "phase_breakdown.csv")
save_csv(phase_records, csv_path)

print(f"\n{'phase':<30} {'count':>6} {'percent':>9}")
for r in phase_records:
    print(f"{r['phase']:<30} {r['count']:>6} {r['percent_of_total']:>8.1f}%")

chart_path = os.path.join(OUT_DIR, "phase_breakdown.png")
hbar(phase_counts, chart_path, "CNS-Targeted AAV Trials by Phase (n=92)")
summary_rows.append({"section": "2. Phase", "csv": "phase_breakdown.csv",
                     "chart": "phase_breakdown.png", "rows": len(phase_records)})


# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 3 — START YEAR
# ═══════════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 65)
print("SECTION 3: START YEAR")
print("=" * 65)


def extract_year(val):
    if not val:
        return None
    m = re.match(r"(\d{4})", str(val).strip())
    return int(m.group(1)) if m else None


years = df["start_date"].apply(extract_year)
null_years = int(years.isna().sum())
if null_years:
    print(f"NOTE: {null_years} row(s) with no parseable start year — listed as NOT REPORTED")

year_counts = years.dropna().astype(int).value_counts().sort_index()
min_yr, max_yr = int(year_counts.index.min()), int(year_counts.index.max())

year_records = [{"year": yr, "count": int(year_counts.get(yr, 0))}
                for yr in range(min_yr, max_yr + 1)]
if null_years:
    year_records.append({"year": "NOT REPORTED", "count": null_years})

csv_path = os.path.join(OUT_DIR, "start_year_breakdown.csv")
save_csv(year_records, csv_path)

print(f"{'year':<14} {'count':>6}")
for r in year_records:
    print(f"{str(r['year']):<14} {r['count']:>6}")

# Line chart (numeric years only)
chart_path = os.path.join(OUT_DIR, "start_year_trend.png")
x_vals = list(range(min_yr, max_yr + 1))
y_vals = [int(year_counts.get(yr, 0)) for yr in x_vals]
fig, ax = plt.subplots(figsize=(8, 5))
ax.plot(x_vals, y_vals, marker="o", color=BAR_COLOR)
for x, y in zip(x_vals, y_vals):
    if y > 0:
        ax.annotate(str(y), (x, y), textcoords="offset points",
                    xytext=(0, 7), ha="center", fontsize=8)
ax.set_xlabel("Year")
ax.set_ylabel("Number of Trials Started")
ax.set_title("CNS-Targeted AAV Trial Starts by Year")
ax.set_xticks(x_vals)
ax.tick_params(axis="x", rotation=45)
plt.tight_layout()
plt.savefig(chart_path, dpi=150)
plt.close()

summary_rows.append({"section": "3. Start Year", "csv": "start_year_breakdown.csv",
                     "chart": "start_year_trend.png", "rows": len(year_records)})


# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 4 — CONDITION / DISEASE CATEGORY
# ═══════════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 65)
print("SECTION 4: CONDITION / DISEASE CATEGORY")
print("=" * 65)

# ── 4a  Raw counts ─────────────────────────────────────────────────────────────
raw_counts = df["conditions"].value_counts().sort_values(ascending=False)
raw_records = [{"raw_condition": c, "count": int(n)} for c, n in raw_counts.items()]
csv_path = os.path.join(OUT_DIR, "raw_condition_counts.csv")
save_csv(raw_records, csv_path)

print("\n4a. Raw condition value_counts (unmodified):")
for r in raw_records:
    print(f"  {r['count']:>3} | {r['raw_condition']}")

# ── 4b  Manual mapping dict ───────────────────────────────────────────────────
CONDITION_MAP = {
    # ── Parkinson's Disease ─────────────────────────────────────────────────
    "Parkinson's Disease":                                                     "Parkinson's Disease",
    "Parkinson Disease":                                                       "Parkinson's Disease",
    "Idiopathic Parkinson's Disease":                                          "Parkinson's Disease",
    "Parkinson Disease; Parkinson's Disease":                                  "Parkinson's Disease",
    "PD":                                                                      "Parkinson's Disease",
    "Parkinson Disease (PD)":                                                  "Parkinson's Disease",
    # ── Spinal Muscular Atrophy ─────────────────────────────────────────────
    "Spinal Muscular Atrophy":                                                 "Spinal Muscular Atrophy",
    "Spinal Muscular Atrophy (SMA)":                                           "Spinal Muscular Atrophy",
    "SMA - Spinal Muscular Atrophy":                                           "Spinal Muscular Atrophy",
    "Spinal Muscular Atrophy Type I":                                          "Spinal Muscular Atrophy",
    "Spinal Muscular Atrophy 1":                                               "Spinal Muscular Atrophy",
    "SMA - Spinal Muscular Atrophy; Gene Therapy":                             "Spinal Muscular Atrophy",
    "Spinal Muscular Atrophy Type I; Spinal Muscular Atrophy Type II; "
    "Spinal Muscular Atrophy Type III; SMA":                                   "Spinal Muscular Atrophy",
    "Spinal Muscular Atrophy Type 3":                                          "Spinal Muscular Atrophy",
    "Spinal Muscular Atrophy Type 2":                                          "Spinal Muscular Atrophy",
    "SMA II":                                                                  "Spinal Muscular Atrophy",
    # ── Rett Syndrome ───────────────────────────────────────────────────────
    "Rett Syndrome":                                                           "Rett Syndrome",
    "RETT Syndrome With Proven MECP2 Mutation":                                "Rett Syndrome",
    # ── Huntington's Disease ────────────────────────────────────────────────
    "Huntington Disease":                                                      "Huntington's Disease",
    "Huntington's Disease":                                                    "Huntington's Disease",
    # ── Amyotrophic Lateral Sclerosis ───────────────────────────────────────
    "Amyotrophic Lateral Sclerosis":                                           "Amyotrophic Lateral Sclerosis",
    "ALS - Amyotrophic Lateral Sclerosis; ALS (Amyotrophic Lateral Sclerosis)": "Amyotrophic Lateral Sclerosis",
    # ── AADC Deficiency ─────────────────────────────────────────────────────
    "AADC Deficiency":                                                         "AADC Deficiency",
    "Aromatic L-amino Acid Decarboxylase (AADC) Deficiency":                  "AADC Deficiency",
    "Aromatic Amino Acid Decarboxylase Deficiency":                            "AADC Deficiency",
    # ── GM1 Gangliosidosis ──────────────────────────────────────────────────
    "Lysosomal Diseases; Gangliosidosis; GM1":                                 "GM1 Gangliosidosis",
    "GM1 Gangliosidosis":                                                      "GM1 Gangliosidosis",
    "GM1 Gangliosidosis; GM1 Gangliosidosis, Type I; GM1 Gangliosidosis, "
    "Type 2; Beta-Galactosidase-1 (GLB1) Deficiency":                         "GM1 Gangliosidosis",
    # ── Alzheimer's Disease ─────────────────────────────────────────────────
    "Alzheimer's Disease":                                                     "Alzheimer's Disease",
    "Alzheimer's Disease; Mild Cognitive Impairment":                          "Alzheimer's Disease",
    # ── Giant Axonal Neuropathy ─────────────────────────────────────────────
    "Giant Axonal Neuropathy (GAN)":                                           "Giant Axonal Neuropathy",
    "Giant Axonal Neuropathy; Gene Transfer":                                  "Giant Axonal Neuropathy",
    # ── Gaucher Disease Type II ─────────────────────────────────────────────
    "Type II Gaucher Disease":                                                 "Gaucher Disease Type II",
    "Gaucher Disease, Type 2":                                                 "Gaucher Disease Type II",
    # ── Frontotemporal Dementia ─────────────────────────────────────────────
    "Frontotemporal Dementia; FTD; FTD-GRN; Dementia, Frontotemporal":        "Frontotemporal Dementia",
    "Frontotemporal Dementia":                                                 "Frontotemporal Dementia",
    # ── MPS II (Hunter Syndrome) ────────────────────────────────────────────
    "Mucopolysaccharidosis Type II (MPS II)":                                  "MPS II (Hunter Syndrome)",
    "MPS II; Hunter Syndrome (MPS II)":                                        "MPS II (Hunter Syndrome)",
    # ── Sanfilippo Syndrome A (MPS IIIA) ───────────────────────────────────
    "MPS IIIA; Sanfilippo Syndrome; Sanfilippo A; Mucopolysaccharidosis III":  "Sanfilippo Syndrome (MPS IIIA)",
    # ── Sanfilippo Syndrome B (MPS IIIB) ───────────────────────────────────
    "Sanfilippo Syndrome B":                                                   "Sanfilippo Syndrome (MPS IIIB)",
    "Mucopolysaccharidosis Type 3 B":                                          "Sanfilippo Syndrome (MPS IIIB)",
    # ── Krabbe Disease ──────────────────────────────────────────────────────
    "Krabbe Disease":                                                          "Krabbe Disease",
    # ── Neuronal Ceroid Lipofuscinosis / Batten Disease ─────────────────────
    "Batten Disease; Late-Infantile Neuronal Ceroid Lipofuscinosis":           "Neuronal Ceroid Lipofuscinosis (Batten Disease)",
    "Batten Disease; Late Infantile Neuronal Ceroid Lipofuscinosis":           "Neuronal Ceroid Lipofuscinosis (Batten Disease)",
    # ── MPS I (Hurler Syndrome) ─────────────────────────────────────────────
    "Mucopolysaccharidosis Type I (MPS I); Hurler Syndrome; Hurler-Scheie Syndrome": "MPS I (Hurler Syndrome)",
    # ── Neuronal Ceroid Lipofuscinosis / Batten Disease (all CLN subtypes) ──
    "Neuronal Ceroid Lipofuscinosis CLN5":                                     "Neuronal Ceroid Lipofuscinosis (Batten Disease)",
    "CLN7":                                                                    "Neuronal Ceroid Lipofuscinosis (Batten Disease)",
    "CLN3; Batten Disease":                                                    "Neuronal Ceroid Lipofuscinosis (Batten Disease)",
    "Variant Late-Infantile Neuronal Ceroid Lipofuscinosis":                   "Neuronal Ceroid Lipofuscinosis (Batten Disease)",
    # ── GM2 Gangliosidosis (Tay-Sachs / Sandhoff) ───────────────────────────
    "Tay-Sachs Disease; Sandhoff Disease":                                     "GM2 Gangliosidosis (Tay-Sachs/Sandhoff)",
    "Infantile GM2 Gangliosidosis (Disorder)":                                 "GM2 Gangliosidosis (Tay-Sachs/Sandhoff)",
    # ── Spastic Paraplegia Type 50 ───────────────────────────────────────────
    "Spastic Paraplegia Type 50":                                              "Spastic Paraplegia Type 50",
    "Spasticity, Muscle; Microcephaly; Intellectual Deficiency; Growth Retardation; SPG50; Spastic Paraplegia":
                                                                               "Spastic Paraplegia Type 50",
}

# Apply mapping; unmapped strings pass through unchanged
df["condition_mapped"] = df["conditions"].apply(lambda c: CONDITION_MAP.get(c, c))
unmapped_unique = sorted(set(c for c in df["conditions"] if c not in CONDITION_MAP))

print("\n4b. Strings NOT in mapping dict (pass-through — review for additional merges):")
if unmapped_unique:
    for s in unmapped_unique:
        print(f"  [{s}]")
else:
    print("  (all strings mapped)")

# ── 4c  Mapped breakdown ───────────────────────────────────────────────────────
mapped_counts = df["condition_mapped"].value_counts().sort_values(ascending=False)
mapped_records = [
    {"condition": c, "count": int(n), "percent_of_total": round(n / EXPECTED_N * 100, 1)}
    for c, n in mapped_counts.items()
]
csv_path = os.path.join(OUT_DIR, "condition_breakdown.csv")
save_csv(mapped_records, csv_path)

print("\n4c. Mapped condition breakdown:")
print(f"{'condition':<58} {'count':>6} {'percent':>9}")
for r in mapped_records:
    print(f"{r['condition']:<58} {r['count']:>6} {r['percent_of_total']:>8.1f}%")

chart_path = os.path.join(OUT_DIR, "condition_breakdown.png")
top15 = mapped_counts.head(15)
hbar(top15, chart_path, "Top 15 Conditions Among CNS-Targeted AAV Trials (n=92)")
summary_rows.append({"section": "4. Condition", "csv": "condition_breakdown.csv",
                     "chart": "condition_breakdown.png", "rows": len(mapped_records)})


# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 5 — DELIVERY ROUTE
# ═══════════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 65)
print("SECTION 5: DELIVERY ROUTE")
print("=" * 65)
print("NOTE: route counts sum to more than 92 due to multi-route trials.")

route_counter: Counter = Counter()
for val in df["cns_route_categories"]:
    if not val or not val.strip():
        route_counter["NOT REPORTED"] += 1
    else:
        for seg in val.split("|"):
            seg = seg.strip()
            if seg:
                route_counter[seg] += 1

route_items = sorted(route_counter.items(), key=lambda x: x[1], reverse=True)
route_records = [
    {"route": r, "trial_count": c, "percent_of_92_trials": round(c / EXPECTED_N * 100, 1)}
    for r, c in route_items
]
csv_path = os.path.join(OUT_DIR, "route_breakdown.csv")
save_csv(route_records, csv_path)

total_route_occs = sum(r["trial_count"] for r in route_records)
print(f"\n{'route':<35} {'trial_count':>12} {'pct_of_92':>10}")
for r in route_records:
    print(f"{r['route']:<35} {r['trial_count']:>12} {r['percent_of_92_trials']:>9.1f}%")
print(f"  Total route occurrences: {total_route_occs}  (exceeds 92 as expected)")

chart_path = os.path.join(OUT_DIR, "route_breakdown.png")
route_series = pd.Series({r["route"]: r["trial_count"] for r in route_records})
hbar(route_series, chart_path,
     "CNS Delivery Routes Among CNS-Targeted AAV Trials\n(n=92, routes not mutually exclusive)")
summary_rows.append({"section": "5. Route", "csv": "route_breakdown.csv",
                     "chart": "route_breakdown.png", "rows": len(route_records)})


# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 6 — SEROTYPE
# ═══════════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 65)
print("SECTION 6: SEROTYPE")
print("=" * 65)

ser_df = pd.read_csv(SER_CSV, dtype=str, keep_default_na=False)
ser_df["_cnt"] = pd.to_numeric(ser_df["count_cns_trials"], errors="coerce").fillna(0).astype(int)
ser_df = ser_df[ser_df["_cnt"] > 0].sort_values("_cnt", ascending=False)

ser_records = [
    {"serotype": row["serotype"],
     "count_cns_trials": row["count_cns_trials"],
     "percent_of_cns_trials": row["percent_of_cns_trials"]}
    for _, row in ser_df.iterrows()
]
csv_path = os.path.join(OUT_DIR, "serotype_breakdown.csv")
save_csv(ser_records, csv_path)

print(f"{'serotype':<15} {'count_cns_trials':>17} {'pct_of_cns_trials':>18}")
for r in ser_records:
    print(f"{r['serotype']:<15} {r['count_cns_trials']:>17} {r['percent_of_cns_trials']:>18}")

chart_path = os.path.join(OUT_DIR, "serotype_breakdown.png")
ser_series = pd.Series({r["serotype"]: int(r["count_cns_trials"]) for r in ser_records})
hbar(ser_series, chart_path,
     "Serotype Usage Among CNS-Targeted AAV Trials\n(n=92, not mutually exclusive)")
summary_rows.append({"section": "6. Serotype", "csv": "serotype_breakdown.csv",
                     "chart": "serotype_breakdown.png", "rows": len(ser_records)})


# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 7 — SPONSOR
# ═══════════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 65)
print("SECTION 7: SPONSOR")
print("=" * 65)


# Explicit alias merges (beyond whitespace/case), applied per user instruction
SPONSOR_MAP = {
    "Novartis Pharmaceuticals": "Novartis Gene Therapies",
}


def normalize_ws(s: str) -> str:
    """Collapse internal whitespace and strip ends."""
    return " ".join(s.split())


raw_sponsors = df["sponsor"].tolist()
# Apply explicit alias merges first, then whitespace-normalize
ws_normalized = [normalize_ws(SPONSOR_MAP.get(s.strip(), s)) for s in raw_sponsors]

# Report alias merges applied
alias_applied = [(raw.strip(), SPONSOR_MAP[raw.strip()])
                 for raw in raw_sponsors if raw.strip() in SPONSOR_MAP]
if alias_applied:
    print("Explicit sponsor alias merges (user-directed):")
    for raw, mapped in sorted(set(alias_applied)):
        print(f"  [{raw}] => [{mapped}]")

# Report whitespace changes
ws_changes = [(raw, norm) for raw, norm in zip(raw_sponsors, ws_normalized) if raw != norm]
if ws_changes:
    print("Whitespace normalization changes (raw => normalized):")
    for raw, norm in sorted(set(ws_changes)):
        print(f"  [{raw}] => [{norm}]")
else:
    print("Whitespace normalization: no changes found.")

# Case normalization: group by lowercase key, display canonical (most-common) form
case_groups: dict[str, list[str]] = defaultdict(list)
for s in ws_normalized:
    case_groups[s.lower()].append(s)

case_merges = {k: sorted(set(v)) for k, v in case_groups.items() if len(set(v)) > 1}
if case_merges:
    print("Case normalization merges:")
    for key, variants in case_merges.items():
        print(f"  key=[{key}]  variants={variants}")
else:
    print("Case normalization: no merges found.")


def canonical(s: str) -> str:
    key = s.lower()
    return Counter(case_groups[key]).most_common(1)[0][0]


df["sponsor_norm"] = [canonical(s) for s in ws_normalized]
sponsor_counts = df["sponsor_norm"].value_counts().sort_values(ascending=False)

sponsor_records = []
for spons, cnt in sponsor_counts.items():
    display = spons if spons else "NOT REPORTED"
    sponsor_records.append({
        "sponsor": display,
        "count": int(cnt),
        "percent_of_total": round(cnt / EXPECTED_N * 100, 1)
    })
csv_path = os.path.join(OUT_DIR, "sponsor_breakdown.csv")
save_csv(sponsor_records, csv_path)

print(f"\n{'sponsor':<52} {'count':>6} {'percent':>9}")
for r in sponsor_records:
    print(f"{r['sponsor']:<52} {r['count']:>6} {r['percent_of_total']:>8.1f}%")

chart_path = os.path.join(OUT_DIR, "sponsor_breakdown.png")
top15_s = sponsor_counts.head(15)
hbar(top15_s, chart_path, "Top 15 Sponsors Among CNS-Targeted AAV Trials (n=92)")
summary_rows.append({"section": "7. Sponsor", "csv": "sponsor_breakdown.csv",
                     "chart": "sponsor_breakdown.png", "rows": len(sponsor_records)})


# ═══════════════════════════════════════════════════════════════════════════════
# FINAL SUMMARY TABLE
# ═══════════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 65)
print("SUMMARY")
print("=" * 65)
print(f"{'Section':<18} {'CSV':<38} {'Chart':<38} {'Rows':>5}")
print("-" * 102)
for r in summary_rows:
    print(f"{r['section']:<18} {r['csv']:<38} {r['chart']:<38} {r['rows']:>5}")
print(f"\nAll outputs saved to: {OUT_DIR}")
