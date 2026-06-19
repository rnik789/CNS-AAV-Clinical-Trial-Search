#!/usr/bin/env python3
"""
AAV CNS Pipeline v9 - ClinicalTrials.gov REST API v2
Date: 2026-06-18

v9 changes vs v8:
  Fix 8 — Add vesemnogene lantuparvovec to AAV_PROPRIETARY_PRODUCT_NAMES
           (Sharma, Joshi, Kumar (2025), Neurosci.) and add query.intr=vesemnogene
           to QUERY_INTRS. Vesemnogene lantuparvovec (an AAV9-hSMN1 vector) uses
           only its INN in API-searchable fields; no generic AAV term appears in
           intervention_names, intervention_descriptions, or brief_summary.
  Fix 9 — Add bare "AAV" to AAV_GENERIC_TERMS with \\bAAV\\b word-boundary matching.
           Uses the same _compile_aav_generic() pre-compiled pattern approach as
           TIER2_CNS_DISEASE_PATTERNS / _compile_t2(). Catches standalone "AAV"
           tokens and hyphenated notations like "AAV-hSMN1" while preventing
           false-positive substring matches within "rAAV", "scAAV", "AAV9", etc.
           AAV_GENERIC_COMPILED replaces the plain contains_any() call in
           get_aav_confirmation() Step 1.
           Summary includes bare-\\bAAV\\b marginal contribution audit vs v8 baseline.

v8 changes vs v7:
  Fix 7 — Restructured AAV confirmation into sequential two-step logic with
           per-trial method/source/product audit fields.
           AAV_TERMS → AAV_GENERIC_TERMS (Step 1, algorithmic).
           AAV_MANUAL_OVERRIDE_TERMS → AAV_PROPRIETARY_PRODUCT_NAMES as
           list of (code, citation) tuples (Step 2, literature-verified).
           Step 2 only runs if Step 1 finds no match.
           New CSV columns: aav_confirmation_method, aav_confirmation_source,
           matched_product_name. Summary reports breakdown by product code.

v7 changes vs v6:
  Fix 6 — Proprietary product code queries and manual AAV confirmation overrides.
           Trials using VY-AADC01/02, PR001/006, LY3884961/963, NGN-101 as their
           only identifier contain no generic AAV terminology in API-searchable fields.
           These products are independently verified as AAV-based via Nair et al. 2024
           (Neurotherapeutics) and added as an explicit, auditable override list.

v6 changes vs v5:
  Fix 5 — Word-boundary regex (\bTERM\b) for short/acronym Tier 2 disease terms
           (SMA, ALS, SCA1-7, CLN, CNS, MPS I, MPS II) to eliminate false-positive
           substring collisions. Long distinctive terms remain as plain substring matches.

v5 changes vs v4:
  Fix 1 — Expanded Intraparenchymal Tier 1 terms (burr hole, into the brain, into putamen)
  Fix 2 — Explicit AAVrh-numbered queries (API does not prefix-match AAVrh vs AAVrh8/10/74)
  Fix 3 — Post-pipeline Gaucher Type 2/3 manual-review flag in summary
  Fix 4 — Supplementary query.term=adeno-associated+virus with status filter;
           net-new contribution reported in summary
  Doc   — Empty intervention_descriptions count reported as methodology limitation

CNS targeting uses a two-tier route-based logic:
  Tier 1 — Unambiguous CNS routes (intrathecal, intraparenchymal, ICV, cisterna magna,
            explicit CNS delivery statements in intervention text)
  Tier 2 — Intravenous delivery + CNS disease keyword confirmation in conditions/title
Route search corpus: intervention_names, intervention_descriptions, brief_summary,
  eligibilityCriteria, designInfo.interventionModelDescription.
"""

import re
import requests
import csv
import os
import sys
import shutil
from collections import defaultdict, Counter

# ─── Configuration ────────────────────────────────────────────────────────────

TODAY       = "2026-06-18"
OUTPUT_DIR  = r"D:\ClaudeCode\FDA_AAV_Search_061826"
BASE_URL    = "https://clinicaltrials.gov/api/v2/studies"
STATUSES    = "COMPLETED,TERMINATED,ACTIVE_NOT_RECRUITING,RECRUITING,ENROLLING_BY_INVITATION,SUSPENDED"
PAGE_SIZE   = 1000

ALL_CSV       = os.path.join(OUTPUT_DIR, f"aav_all_confirmed_v9_{TODAY}.csv")
CNS_CSV       = os.path.join(OUTPUT_DIR, f"aav_cns_confirmed_v9_{TODAY}.csv")
EXCLUDED_CSV  = os.path.join(OUTPUT_DIR, f"aav_excluded_v9_{TODAY}.csv")
SEROTYPE_CSV  = os.path.join(OUTPUT_DIR, f"aav_serotype_frequency_v9_{TODAY}.csv")
SUMMARY_TXT   = os.path.join(OUTPUT_DIR, f"aav_summary_v9_{TODAY}.txt")
SCRIPT_DEST   = os.path.join(OUTPUT_DIR, f"aav_cns_pipeline_v9_{TODAY}.py")

# ─── Query List ───────────────────────────────────────────────────────────────
# Exhaustive set run in order; per-query new-unique-trial contributions are
# tracked and reported so signal vs. duplicate queries are visible.

QUERY_INTRS = [
    # Vector name queries
    "adeno-associated virus",
    "adeno associated virus",
    "rAAV",
    "AAV1", "AAV2", "AAV3", "AAV4", "AAV5",
    "AAV6", "AAV7", "AAV8", "AAV9",
    "AAVrh", "AAVhu",
    "scAAV", "ssAAV",
    # Disease-specific construct / program name queries
    "AAV2-hAADC",
    "AAV2-GDNF",
    "AAV2-GAD",
    "AAV2-NGF",
    "AAV2-BDNF",
    "AAV9-SMN",
    "AAV-GAD",
    "CERE-110",
    "CERE-120",
    "RGX-121",
    "RGX-314",
    "AMT-130",
    "AMT-162",
    "TSHA-102",
    "FBX-101",
    "AB-1001",
    "SGT-212",
    "GC101",
    "SKG0201",
    "IPS101A",
    "JAG201",
    "EXG001",
    "BBP-812",
    "UX111",
    "AAVhu68",
    "PBGM01",
    "MVX-220",
    "NSR-REP1",
    "ABO-101",
    "ABO-102",
    "ABO-50",
    "ABO-202",
    # Explicit AAVrh-numbered variants — API does not prefix-match query.intr=AAVrh
    # against tokens like AAVrh8, AAVrh10, etc. (confirmed by Diagnostic 3).
    # CSV-scan also confirms AAVrh.10 and AAVrh.74 are present in the dataset.
    "AAVrh8", "AAVrh10", "AAVrh74", "AAVrh20", "AAVrh39",
    # Proprietary product codes verified as AAV-based via Nair et al. 2024
    # (Neurotherapeutics). These products use no generic AAV terminology in
    # API-searchable fields and would be missed by algorithmic queries alone.
    "VY-AADC01", "VY-AADC02", "VY-AADC",
    "PR001", "PR006",
    "LY3884961", "LY3884963",
    "NGN-101",
    # Fix 8 — vesemnogene lantuparvovec (AAV9-hSMN1 for SMA): INN contains no
    # generic AAV terminology; verified AAV-based via Sharma, Joshi, Kumar (2025).
    "vesemnogene",
]

# ─── AAV Confirmation Terms ───────────────────────────────────────────────────
# Two-step sequential confirmation (Fix 7).
# Step 1 — generic terminology (algorithmic, registry text only):
#   aav_confirmation_method = "generic_terminology"
#   aav_confirmation_source = "registry_text"
#   matched_product_name    = ""
# Step 2 — proprietary product name (only if Step 1 fails):
#   aav_confirmation_method = "proprietary_product_name"
#   aav_confirmation_source = <per-entry citation string>
#   matched_product_name    = <matched code>
# Searched in: intervention_names, intervention_descriptions, brief_summary

AAV_GENERIC_TERMS = [
    # Generic / algorithmically discovered AAV terminology
    "adeno-associated virus",
    "adeno-associated viral",
    "rAAV", "scAAV", "ssAAV",
    "AAV1", "AAV2", "AAV3", "AAV4", "AAV5",
    "AAV6", "AAV7", "AAV8", "AAV9",
    "AAVrh", "AAVhu",
    "onasemnogene",
    "voretigene",
    "eladocagene",
    "valoctocogene",
    "etranacogene",
    "delandistrogene",
    # Fix 9 — bare "AAV" with word-boundary protection (see AAV_GENERIC_WB_TERMS).
    # Catches "AAV-hSMN1", "using AAV vector", etc. without matching within
    # "rAAV", "scAAV", "AAV9" (those have adjacent word chars, so \b doesn't fire).
    "AAV",
]

# ─── Fix 9: Word-boundary compilation for AAV_GENERIC_TERMS ──────────────────
# Short bare "AAV" would collide as a substring inside "rAAV", "scAAV", "AAV9",
# etc. if matched plainly. Word-boundary anchoring prevents this. All other
# generic terms are long enough that plain substring matching is safe.
AAV_GENERIC_WB_TERMS = frozenset(["AAV"])


def _compile_aav_generic(term):
    """Mirror of _compile_t2(): apply \\bTERM\\b for terms in AAV_GENERIC_WB_TERMS."""
    pat = (r'\b' + re.escape(term) + r'\b') if term in AAV_GENERIC_WB_TERMS else re.escape(term)
    return re.compile(pat, re.IGNORECASE)


AAV_GENERIC_COMPILED = [(t, _compile_aav_generic(t)) for t in AAV_GENERIC_TERMS]

# Proprietary product codes externally verified as AAV-based. Each entry is a
# (product_code, citation) tuple. Step 2 only runs when Step 1 finds no match,
# so a trial confirmed here is unambiguously proprietary-only in registry text.
AAV_PROPRIETARY_PRODUCT_NAMES = [
    ("VY-AADC01", "Nair et al. 2024, Neurotherapeutics"),
    ("VY-AADC02", "Nair et al. 2024, Neurotherapeutics"),
    ("VY-AADC",   "Nair et al. 2024, Neurotherapeutics"),
    ("PR001",     "Nair et al. 2024, Neurotherapeutics"),
    ("PR006",     "Nair et al. 2024, Neurotherapeutics"),
    ("LY3884961", "Nair et al. 2024, Neurotherapeutics"),
    ("LY3884963", "Nair et al. 2024, Neurotherapeutics"),
    ("NGN-101",   "Nair et al. 2024, Neurotherapeutics"),
    # Fix 8 — vesemnogene lantuparvovec: AAV9-hSMN1 vector for SMA.
    # INN contains no generic AAV terminology in any API-searchable field.
    ("vesemnogene lantuparvovec", "Sharma, Joshi, Kumar (2025), Neurosci."),
]

# ─── Non-CNS Exclusion Terms ─────────────────────────────────────────────────
# Applied after AAV confirmation, before CNS flagging.
# Trials whose conditions/title match any term here are removed from the active
# dataset UNLESS a Tier 1 CNS delivery term is present in the delivery fields
# (the Tier 1 override preserves genuinely novel CNS-targeted trials).

EXCLUSION_TERMS = {
    "Blood disorders": [
        "hemophilia", "haemophilia", "factor VIII", "factor IX", "von Willebrand",
    ],
    "Ocular only": [
        "retinal dystrophy", "retinitis pigmentosa", "macular degeneration",
        "achromatopsia", "leber congenital amaurosis", "RPE65",
        "choroideremia", "retinoschisis", "Stargardt", "RPGR",
    ],
    "Liver/metabolic": [
        "ornithine transcarbamylase", "OTC deficiency",
        "alpha-1 antitrypsin", "AAT deficiency",
        "Crigler-Najjar", "phenylketonuria", "PKU",
        "glycogen storage disease", "Wilson disease",
        "familial hypercholesterolemia",
    ],
    "Muscle": [
        "myotubular myopathy", "Duchenne muscular dystrophy", "DMD",
        "limb girdle muscular dystrophy", "LGMD", "Becker muscular dystrophy",
    ],
    "Cardiac only": [
        "cardiomyopathy of Friedreich", "hypertrophic cardiomyopathy", "heart failure",
    ],
    "Other peripheral": [
        "parotid", "salivary gland", "alpha-1 antitrypsin lung",
        "pulmonary", "arthritis", "joint",
    ],
}

# ─── Tier 1 — Unambiguous CNS Delivery Routes ─────────────────────────────────
# Searched in: build_delivery_corpus() — intervention_names, intervention_descriptions,
#   brief_summary, eligibility_criteria, design_info_text.
# A match here flags CNS regardless of disease indication.

TIER1_ROUTES = {
    "Intrathecal": [
        "intrathecal",
        "intrathecally",
        "lumbar puncture",
        "lumbar intrathecal",
        "injection into the cerebrospinal fluid",
        "CSF administration",
        "IT administration",
    ],
    "Intraparenchymal": [
        "intraparenchymal",
        "stereotaxic",
        "stereotactic",
        "intracerebral",
        "intrathalamic",
        "intrastriatal",
        "intra-striatal",
        "subthalamic nucleus",
        "convection enhanced delivery",
        "direct injection into the brain",
        "direct brain injection",
        "intraputaminal",
        "into the putamen",
        "into putamen",
        "putaminal injection",
        "burr hole",
        "burr holes",
        "into the brain",
        "into brain",
    ],
    "Intracerebroventricular": [
        "intracerebroventricular",
        "intraventricular",
        "ICV",
    ],
    "Cisterna magna": [
        "cisterna magna",
        "intracisternal",
        "intracisternally",
    ],
    "Explicit CNS delivery": [
        "to the central nervous system",
        "to the cns",
        "to the brain and spinal cord",
        "to the brain and peripheral tissues",
        "delivered to the brain",
        "targeting the central nervous system",
        "crosses the blood-brain barrier",
        "crosses the blood brain barrier",
    ],
}

# ─── Tier 2 — Intravenous Delivery Terms ─────────────────────────────────────
# Searched in: intervention_names, intervention_descriptions, brief_summary

TIER2_IV_TERMS = [
    "intravenous",
    "intravenously",
    "IV infusion",
    "systemic administration",
    "systemic delivery",
    "IV administration",
]

# ─── Tier 2 — CNS Disease Confirmation Terms ─────────────────────────────────
# Searched in: conditions, brief_title only.
# An IV trial is CNS-targeted only when at least one of these also matches.

TIER2_CNS_DISEASE_TERMS = [
    "spinal muscular atrophy", "SMA",
    "Huntington", "Parkinson", "Alzheimer",
    "ALS", "amyotrophic lateral sclerosis",
    "Batten", "leukodystrophy", "mucopolysaccharidosis",
    "gangliosidosis", "Rett", "AADC",
    "aromatic L-amino acid", "giant axonal", "Canavan",
    "seizure",
    "spinocerebellar ataxia", "SCA1", "SCA2", "SCA3", "SCA6", "SCA7", "cerebellar ataxia",
    "mesial temporal lobe epilepsy", "CLN",
    "frontotemporal dementia", "FTD-GRN",
    "glioblastoma", "glioma", "cerebellar",
    "motor neuron", "dopaminergic",
    "neuronal ceroid", "metachromatic", "Krabbe",
    "spinal cord", "brain", "central nervous system",
    "CNS", "neurodegenerative", "neurological",
    "Angelman syndrome", "UBE3A",
    "Hunter syndrome", "MPS II", "mucopolysaccharidosis II",
    "Sanfilippo", "Hurler syndrome", "MPS I",
    "Tay-Sachs", "Sandhoff",
]

# Terms that require a CNS-confirming co-term elsewhere in conditions+title+delivery corpus.
# Format: (disease_term, [accepted_co_terms])
TIER2_CNS_DISEASE_TERMS_CONDITIONAL = [
    ("Pompe disease", ["CNS", "neurological"]),
]

# ─── Tier 2 — word-boundary set and pre-compiled patterns (Fix 5) ─────────────
# Short acronyms (<=4 chars) get \b word-boundary anchors to prevent false-positive
# substring collisions (e.g. "SMA" in "plasma", "CLN" in "clinical", "CNS" in
# "AADC/CNS" compound labels, "ALS" in "false"). Long, distinctive terms remain as
# plain substring matches — their length makes collisions negligible.
TIER2_WB_TERMS = frozenset([
    "SMA",
    "ALS",
    "SCA1", "SCA2", "SCA3", "SCA6", "SCA7",
    "CLN",
    "CNS",
    "MPS I", "MPS II",
])


def _compile_t2(term):
    pat = (r'\b' + re.escape(term) + r'\b') if term in TIER2_WB_TERMS else re.escape(term)
    return re.compile(pat, re.IGNORECASE)


TIER2_CNS_DISEASE_PATTERNS = [(t, _compile_t2(t)) for t in TIER2_CNS_DISEASE_TERMS]

# Pre-compiled co-term patterns for conditional checks; "CNS" gets \b anchoring.
_COND_WB = frozenset(["CNS"])
TIER2_CONDITIONAL_COMPILED = [
    (term, [re.compile((r'\b' + re.escape(c) + r'\b') if c in _COND_WB else re.escape(c),
                       re.IGNORECASE)
            for c in co_terms])
    for term, co_terms in TIER2_CNS_DISEASE_TERMS_CONDITIONAL
]

# ─── Helpers ──────────────────────────────────────────────────────────────────

def contains_any(text, terms):
    if not text:
        return False
    tl = text.lower()
    return any(t.lower() in tl for t in terms)


def find_matches(text, terms):
    if not text:
        return []
    tl = text.lower()
    return [t for t in terms if t.lower() in tl]


def find_matches_wb(text, term_patterns):
    """Match Tier 2 disease terms using pre-compiled patterns (some with word boundaries)."""
    if not text:
        return []
    return [term for term, pat in term_patterns if pat.search(text)]


def contains_any_compiled(text, compiled_patterns):
    """Check if text matches any pre-compiled regex pattern."""
    if not text:
        return False
    return any(p.search(text) for p in compiled_patterns)


def dedup_ordered(items):
    seen, out = set(), []
    for x in items:
        k = x.lower()
        if k not in seen:
            seen.add(k)
            out.append(x)
    return out


def build_delivery_corpus(f):
    """Full text corpus searched for delivery route and CNS terms (Fix 3 fields included)."""
    return " ".join([
        f["intervention_names"],
        f["intervention_descriptions"],
        f["brief_summary"],
        f.get("eligibility_criteria", ""),
        f.get("design_info_text", ""),
    ])


def build_query_url(query_intr):
    q = query_intr.replace(" ", "+")
    return (
        f"{BASE_URL}?query.intr={q}"
        f"&filter.overallStatus={STATUSES}"
        f"&pageSize={PAGE_SIZE}&format=json"
    )

def get_nct(study):
    return (study.get("protocolSection", {})
                 .get("identificationModule", {})
                 .get("nctId", ""))

# ─── API Fetch ────────────────────────────────────────────────────────────────

def fetch_all_studies(params):
    """Paginate to exhaustion for any params dict. Caller builds query.intr or query.term."""
    studies, page_num, next_token = [], 0, None
    p = dict(params)
    while True:
        if next_token:
            p["pageToken"] = next_token
        elif "pageToken" in p:
            del p["pageToken"]

        page_num += 1
        print(f"  Page {page_num} ...", flush=True)
        try:
            r = requests.get(BASE_URL, params=p, timeout=90)
            r.raise_for_status()
        except requests.RequestException as exc:
            print(f"  ERROR page {page_num}: {exc}", file=sys.stderr)
            break

        data  = r.json()
        batch = data.get("studies", [])
        studies.extend(batch)
        print(f"    {len(batch)} received  |  running total: {len(studies)}", flush=True)

        next_token = data.get("nextPageToken")
        if not next_token:
            break
    return studies

# ─── Field Extraction ─────────────────────────────────────────────────────────

def extract_fields(study):
    proto  = study.get("protocolSection", {})
    id_m   = proto.get("identificationModule", {})
    st_m   = proto.get("statusModule", {})
    de_m   = proto.get("descriptionModule", {})
    co_m   = proto.get("conditionsModule", {})
    ar_m   = proto.get("armsInterventionsModule", {})
    ds_m   = proto.get("designModule", {})
    el_m   = proto.get("eligibilityModule", {})
    sp_m   = proto.get("sponsorCollaboratorsModule", {})

    phases = ds_m.get("phases", [])
    design_info = ds_m.get("designInfo", {})
    design_info_text = " ".join(v for v in design_info.values() if isinstance(v, str))
    eligibility_criteria = el_m.get("eligibilityCriteria", "")

    iv_names, iv_descs = [], []
    for iv in ar_m.get("interventions", []):
        n = iv.get("name", "").strip()
        d = iv.get("description", "").strip()
        if n: iv_names.append(n)
        if d: iv_descs.append(d)

    return {
        "nct_id":                    id_m.get("nctId", ""),
        "brief_title":               id_m.get("briefTitle", ""),
        "overall_status":            st_m.get("overallStatus", ""),
        "sponsor":                   sp_m.get("leadSponsor", {}).get("name", ""),
        "phase":                     "; ".join(phases) if phases else "N/A",
        "start_date":                st_m.get("startDateStruct", {}).get("date", ""),
        "primary_completion_date":   st_m.get("primaryCompletionDateStruct", {}).get("date", ""),
        "conditions":                "; ".join(co_m.get("conditions", [])),
        "intervention_names":        "; ".join(iv_names),
        "intervention_descriptions": "; ".join(iv_descs),
        "brief_summary":             de_m.get("briefSummary", ""),
        "eligibility_criteria":      eligibility_criteria,
        "design_info_text":          design_info_text,
    }

# ─── AAV Confirmation ─────────────────────────────────────────────────────────

def get_aav_confirmation(f):
    """
    Sequential two-step AAV confirmation.
    Returns (method, source, matched_product_name) or None.

    Step 1: generic terminology — uses AAV_GENERIC_COMPILED (pre-compiled patterns);
            bare "AAV" gets \\bAAV\\b word-boundary anchoring (Fix 9); all other
            terms use plain case-insensitive substring match via compiled re.escape.
    Step 2: proprietary product name — only runs if Step 1 finds no match; plain
            case-insensitive substring match; requires per-entry literature citation.
    """
    corpus = " ".join([f["intervention_names"], f["intervention_descriptions"], f["brief_summary"]])
    # Step 1 — generic terminology (word-boundary-aware via compiled patterns)
    if any(pat.search(corpus) for _, pat in AAV_GENERIC_COMPILED):
        return ("generic_terminology", "registry_text", "")
    # Step 2 — proprietary product name (sequential: only reached if Step 1 failed)
    cl = corpus.lower()
    for code, citation in AAV_PROPRIETARY_PRODUCT_NAMES:
        if code.lower() in cl:
            return ("proprietary_product_name", citation, code)
    return None

# ─── Non-CNS Exclusion ────────────────────────────────────────────────────────

def has_tier1_delivery(f):
    corpus = build_delivery_corpus(f)
    return any(contains_any(corpus, terms) for terms in TIER1_ROUTES.values())


def check_exclusion(f):
    """
    Returns (is_excluded, category, triggering_term).
    Excluded when: (conditions OR brief_title) matches an EXCLUSION_TERMS entry
                   AND no Tier 1 delivery term is present (Tier 1 override).
    """
    title_cond = f["conditions"] + " " + f["brief_title"]
    for category, terms in EXCLUSION_TERMS.items():
        hits = find_matches(title_cond, terms)
        if hits:
            if has_tier1_delivery(f):
                return False, "", ""
            return True, category, hits[0]
    return False, "", ""

# ─── CNS Targeting — Two-Tier Logic ──────────────────────────────────────────

def check_tier1(f):
    corpus = build_delivery_corpus(f)
    matched_routes, all_hits = [], []
    for route, terms in TIER1_ROUTES.items():
        hits = find_matches(corpus, terms)
        if hits:
            matched_routes.append(route)
            all_hits.extend(hits)
    return matched_routes, dedup_ordered(all_hits)


def check_tier2(f):
    deliv_corpus = build_delivery_corpus(f)
    iv_hits = find_matches(deliv_corpus, TIER2_IV_TERMS)
    if not iv_hits:
        return False, [], ""

    disease_corpus = f["conditions"] + " " + f["brief_title"]
    # Fix 5: use word-boundary patterns for short/acronym terms
    disease_hits = find_matches_wb(disease_corpus, TIER2_CNS_DISEASE_PATTERNS)

    for term, co_patterns in TIER2_CONDITIONAL_COMPILED:
        if contains_any(disease_corpus, [term]):
            if contains_any_compiled(disease_corpus + " " + deliv_corpus, co_patterns):
                disease_hits.append(term)

    if not disease_hits:
        return False, iv_hits, ""

    return True, iv_hits, disease_hits[0]


def determine_cns(f):
    t1_routes, t1_terms              = check_tier1(f)
    t2_match, t2_iv_terms, t2_disease = check_tier2(f)

    is_t1 = bool(t1_routes)
    is_t2 = t2_match

    if not is_t1 and not is_t2:
        return None

    if is_t1 and is_t2:
        tier = "Tier 1 + Tier 2"
    elif is_t1:
        tier = "Tier 1"
    else:
        tier = "Tier 2"

    route_cats    = list(t1_routes) + (["Intravenous"] if is_t2 else [])
    trigger_terms = dedup_ordered(t1_terms + (t2_iv_terms if is_t2 else []))

    return {
        "cns_tier":                 tier,
        "cns_route_categories":     "|".join(route_cats),
        "cns_triggering_terms":     "|".join(trigger_terms),
        "iv_cns_confirmation_term": t2_disease if is_t2 else "",
    }


def detect_iv_for_non_cns(f):
    corpus = build_delivery_corpus(f)
    hits   = dedup_ordered(find_matches(corpus, TIER2_IV_TERMS))
    return "|".join(hits) if hits else "not detected"

# ─── Serotype Extraction ──────────────────────────────────────────────────────
# Design decisions:
# - sc/ss prefixes (scAAV9) are kept as distinct serotype entries — they denote
#   a different capsid configuration with different packaging capacity.
# - Leading 'r' recombinant prefix (rAAVhu68) is stripped — it denotes production
#   method, not capsid identity, so rAAVhu68 normalizes to AAVhu68.
# - rh/hu arms now require trailing digits (e.g. [- .]?\d+) to prevent matching
#   'AAVhuman', which the prior '\w*' arm was incorrectly consuming as AAVhu+man.
# - Cargo/transgene suffixes on rh/hu variants (AAVrh.10CUARSA, AAVhu68hFXN) are
#   stripped to a normalized parent for counting; the full string is preserved in
#   construct_variant for audit.
# - Chimeric notations (AAV2/9, AAV2/5) are kept as-is — distinct capsids.
# - All normalized serotypes are uppercased for consistent Counter grouping.

SEROTYPE_RE = re.compile(
    r'\b(?:r|sc|ss)?AAV'
    r'(?:'
    r'[- .]?(?:rh|hu)[- .]?\d+'        # AAVrh10, AAV-rh10, AAVrh.10, AAVhu68,
                                        # AAV-hu68, AAV.hu68, rAAVhu68
    r'|[- .]?php[- .]?\w*'             # AAV-PHP.B, AAV-PHP.eB, AAV-PHP.S
    r'|[- .]?(?:DJ|2i8|B1)(?:/\d+)?'  # AAV-DJ, AAV-DJ/8, AAV2i8, AAVB1
    r'|\d+(?:[./]\d+)?'                # AAV2, AAV9, AAV2/5, AAV2.5, AAV2/9
    r')\b',
    re.IGNORECASE,
)

# Detects cargo/transgene suffixes on rh and hu family variants.
# E.g. AAVRH.10CUARSA → parent AAVRH.10, suffix CUARSA.
CARGO_RE = re.compile(
    r'^((?:SC|SS)?AAV(?:RH\.?\d+|HU\.?\d+))([A-Z][A-Z0-9]*)$'
)

# Diagnostic: catch 'AAVhu<letter>' patterns that would be false positives.
AAVHU_NONDIGIT_RE = re.compile(r'AAVhu([a-zA-Z][^\s]*)', re.IGNORECASE)


def normalize_raw_match(raw):
    """Uppercase, strip internal whitespace, remove leading 'r' recombinant prefix."""
    s = re.sub(r'\s+', '', raw).upper()
    # rAAVhu68 → RAAVHU68: strip the R (production method, not capsid identity).
    # sc/ss prefixes are NOT stripped — they change the vector biology.
    if s.startswith('R') and len(s) > 3 and s[1:3] == 'AA':
        s = s[1:]
    # Canonicalize AAVrh family: insert dot before digit (AAVrh10 → AAVrh.10).
    # The dot notation is the IUPAC/Gao-lab canonical form; both forms appear in text.
    s = re.sub(r'^((?:SC|SS)?AAVRH)(\d)', r'\1.\2', s)
    return s


def split_serotype_cargo(normalized):
    """
    Returns (parent_serotype, construct_variant).
    Strips transgene cargo suffix from rh/hu family strings.
    Non-rh/hu serotypes and chimeric notations (AAV2/9) are returned unchanged.
      'AAVRH.10CUARSA' → ('AAVRH.10', 'AAVRH.10CUARSA')
      'AAV9'           → ('AAV9', '')
      'AAV2/9'         → ('AAV2/9', '')
    """
    m = CARGO_RE.match(normalized)
    if m:
        return m.group(1), normalized
    return normalized, ""


def extract_serotypes_with_variants(f):
    """
    Return set of (parent_serotype, construct_variant) pairs from the trial's
    three delivery fields. Set deduplication ensures each pair counted once per trial.
    """
    corpus = " ".join([
        f["intervention_names"],
        f["intervention_descriptions"],
        f["brief_summary"],
    ])
    pairs = set()
    for raw in SEROTYPE_RE.findall(corpus):
        norm   = normalize_raw_match(raw)
        parent, variant = split_serotype_cargo(norm)
        pairs.add((parent, variant))
    return pairs


def find_aavhu_nondigit_context(f):
    """
    Return list of {match, context} dicts for any 'AAVhu<letter>...' hit in the
    trial corpus. Used to surface AAVHUMAN-type regex false positives for audit.
    """
    corpus = " ".join([
        f["intervention_names"],
        f["intervention_descriptions"],
        f["brief_summary"],
    ])
    hits = []
    for m in AAVHU_NONDIGIT_RE.finditer(corpus):
        s   = max(0, m.start() - 50)
        e   = min(len(corpus), m.end() + 50)
        hits.append({
            "match":   m.group(),
            "context": f"...{corpus[s:m.start()]}[{m.group()}]{corpus[m.end():e]}...",
        })
    return hits


# ─── Main ─────────────────────────────────────────────────────────────────────

def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    SEP = "=" * 70
    print(SEP)
    print("AAV CNS Pipeline v9  |  ClinicalTrials.gov REST API v2")
    print(f"Date of query: {TODAY}")
    print(f"Total queries to run: {len(QUERY_INTRS)}")
    print(SEP)

    # ── 1. Fetch — per-query with contribution tracking ───────────────────────
    nct_map      = {}   # NCT ID -> study dict; built incrementally for dedup
    total_raw    = 0
    contributions = []  # [(query_intr, raw_count, new_unique_count)]

    for q in QUERY_INTRS:
        print(f"\n[{len(contributions)+1}/{len(QUERY_INTRS)}] query.intr={q!r}", flush=True)
        batch     = fetch_all_studies({
            "query.intr": q,
            "filter.overallStatus": STATUSES,
            "pageSize": PAGE_SIZE,
            "format": "json",
        })
        raw_count = len(batch)
        new_count = 0
        for s in batch:
            nct = get_nct(s)
            if nct and nct not in nct_map:
                nct_map[nct] = s
                new_count += 1
        total_raw += raw_count
        contributions.append((q, raw_count, new_count))
        print(f"  => {raw_count} raw  |  {new_count} new unique  |  cumulative unique: {len(nct_map)}")

    # ── 1b. Supplementary query.term fetch (Fix 4) ────────────────────────────
    # Run AFTER all query.intr queries so net-new contribution is measurable.
    print(f"\n[SUPPLEMENTARY] query.term='adeno-associated virus' (status filter applied)",
          flush=True)
    nct_count_pre_supp = len(nct_map)
    supp_batch = fetch_all_studies({
        "query.term": "adeno-associated virus",
        "filter.overallStatus": STATUSES,
        "pageSize": PAGE_SIZE,
        "format": "json",
    })
    supp_raw = len(supp_batch)
    for s in supp_batch:
        nct = get_nct(s)
        if nct and nct not in nct_map:
            nct_map[nct] = s
    supp_new = len(nct_map) - nct_count_pre_supp
    print(f"  => {supp_raw} raw  |  {supp_new} net new unique beyond query.intr set")

    total_deduped = len(nct_map)
    deduped = list(nct_map.values())

    print(f"\n{'─'*50}")
    print(f"Total raw across all queries : {total_raw}")
    print(f"Supplementary query.term raw : {supp_raw}  ({supp_new} net new unique)")
    print(f"Total unique after dedup     : {total_deduped}")

    # ── 2. Extract fields + AAV confirmation ──────────────────────────────────
    confirmed, dropped = [], 0
    # (nct_id, matched_code, citation) for proprietary-only confirmations
    proprietary_confirmed = []
    for s in deduped:
        f = extract_fields(s)
        if not f["nct_id"]:
            dropped += 1
            continue
        conf = get_aav_confirmation(f)
        if conf:
            method, source, matched_code = conf
            f["aav_confirmation_method"] = method
            f["aav_confirmation_source"] = source
            f["matched_product_name"]    = matched_code
            confirmed.append(f)
            if method == "proprietary_product_name":
                proprietary_confirmed.append((f["nct_id"], matched_code, source))
        else:
            dropped += 1

    total_confirmed = len(confirmed)
    print(f"Dropped (AAV filter)         : {dropped}")
    print(f"Confirmed AAV trials         : {total_confirmed}")
    print(f"  generic_terminology        : {total_confirmed - len(proprietary_confirmed)}")
    print(f"  proprietary_product_name   : {len(proprietary_confirmed)}")

    # ── 2b. Fix 8/9 attribution vs v8 baseline ────────────────────────────────
    # Load v8 confirmed NCT IDs for marginal-contribution analysis.
    # aav_all_confirmed_v8 holds only post-exclusion active trials (162);
    # aav_excluded_v8 holds the 156 excluded trials — both needed for all 318.
    v8_nct_ids = set()
    for _v8fn in ["aav_all_confirmed_v8_2026-06-18.csv",
                  "aav_excluded_v8_2026-06-18.csv"]:
        _v8path = os.path.join(OUTPUT_DIR, _v8fn)
        if os.path.exists(_v8path):
            with open(_v8path, newline="", encoding="utf-8") as _fh:
                v8_nct_ids.update(r["nct_id"] for r in csv.DictReader(_fh))

    # Counterfactual: would a trial have been confirmed by v8's logic?
    # v8 used contains_any(corpus, AAV_GENERIC_TERMS_WITHOUT_AAV) and the old
    # proprietary list (without vesemnogene lantuparvovec).
    _V8_GENERIC_TERMS   = [t for t in AAV_GENERIC_TERMS if t != "AAV"]
    _V8_PROPRIETARY     = [(c, cit) for c, cit in AAV_PROPRIETARY_PRODUCT_NAMES
                           if c != "vesemnogene lantuparvovec"]
    _BARE_AAV_PAT       = re.compile(r'\bAAV\b', re.IGNORECASE)

    def _would_pass_v8_aav(f):
        corpus = " ".join([f["intervention_names"], f["intervention_descriptions"],
                           f["brief_summary"]])
        if contains_any(corpus, _V8_GENERIC_TERMS):
            return True
        cl = corpus.lower()
        return any(code.lower() in cl for code, _ in _V8_PROPRIETARY)

    new_trials = [f for f in confirmed if f["nct_id"] not in v8_nct_ids]
    bare_aav_new    = []   # new AND confirmed only because of bare \bAAV\b
    vesemnogene_new = []   # new AND vesemnogene text present in corpus

    for f in new_trials:
        corpus = " ".join([f["intervention_names"], f["intervention_descriptions"],
                           f["brief_summary"]])
        if not _would_pass_v8_aav(f) and _BARE_AAV_PAT.search(corpus):
            bare_aav_new.append(f)
        if "vesemnogene" in corpus.lower():
            vesemnogene_new.append(f)

    print(f"\n  [Fix 8/9 vs v8 baseline]")
    print(f"  New trials total            : {len(new_trials)}")
    print(f"  New via vesemnogene text    : {len(vesemnogene_new)}")
    print(f"  New via bare \\bAAV\\b only   : {len(bare_aav_new)}")

    # ── 3. Non-CNS exclusion pass ─────────────────────────────────────────────
    active, excluded = [], []
    for f in confirmed:
        is_excl, excl_cat, excl_term = check_exclusion(f)
        if is_excl:
            f["exclusion_category"] = excl_cat
            f["exclusion_term"]     = excl_term
            excluded.append(f)
        else:
            active.append(f)

    total_excluded = len(excluded)
    print(f"Excluded (non-CNS filter)    : {total_excluded}  (Tier 1 override applied)")
    print(f"Active after exclusion       : {len(active)}")

    # ── 4. Two-tier CNS determination ─────────────────────────────────────────
    cns_trials, non_cns_trials = [], []

    for f in active:
        cns_info = determine_cns(f)
        if cns_info:
            f.update(cns_info)
            f["cns_flagged"]    = "yes"
            f["delivery_route"] = cns_info["cns_route_categories"]
            cns_trials.append(f)
        else:
            f["cns_flagged"]    = "no"
            f["delivery_route"] = detect_iv_for_non_cns(f)
            f["cns_tier"]                 = ""
            f["cns_route_categories"]     = ""
            f["cns_triggering_terms"]     = ""
            f["iv_cns_confirmation_term"] = ""
            non_cns_trials.append(f)

    total_cns = len(cns_trials)
    print(f"CNS-flagged trials           : {total_cns}")

    # ── 4c. Fix 3 — Gaucher Type 2/3 manual-review flag ──────────────────────
    GAUCHER_CNS_TERMS = [
        "gaucher disease type 2", "gaucher disease type 3",
        "type 2 gaucher", "type 3 gaucher",
        "gaucher type 2", "gaucher type 3",
        "neuronopathic gaucher", "neuronopathic gaucher disease",
    ]
    cns_nct_set  = {t["nct_id"] for t in cns_trials}
    excl_nct_set = {t["nct_id"] for t in excluded}
    gaucher_review = []
    for f in confirmed:
        haystack = (f["conditions"] + " " + f["brief_title"] + " " +
                    f["brief_summary"]).lower()
        if any(term in haystack for term in GAUCHER_CNS_TERMS):
            if f["nct_id"] in cns_nct_set:
                flag = "CNS-flagged"
            elif f["nct_id"] in excl_nct_set:
                flag = "excluded"
            else:
                flag = "active-not-CNS-flagged"
            gaucher_review.append((f["nct_id"], f["brief_title"], flag))

    # ── 4d. Empty intervention_descriptions audit ─────────────────────────────
    empty_iv_desc = [
        (f["nct_id"], f["brief_title"])
        for f in confirmed
        if not f.get("intervention_descriptions", "").strip()
    ]

    # ── 4b. Serotype frequency analysis + AAVHUMAN diagnostic ─────────────────
    # Counts per-trial: each serotype counted once per trial even if mentioned many times.
    # (a) all confirmed AAV trials — full serotype landscape (pre-exclusion)
    # (b) CNS-flagged trials only  — CNS-specific serotype usage
    cns_nct_ids = {t["nct_id"] for t in cns_trials}

    serotype_all          = Counter()
    serotype_cns          = Counter()
    variants_all          = defaultdict(set)  # parent → set of full cargo-suffixed strings
    multi_serotype_all_cnt = 0
    multi_serotype_cns_cnt = 0
    aavhuman_hits         = []               # audit log for AAVhu<letter> false positives

    for f in confirmed:
        pairs   = extract_serotypes_with_variants(f)
        parents = {p for p, _ in pairs}
        if len(parents) > 1:
            multi_serotype_all_cnt += 1
        for parent, variant in pairs:
            serotype_all[parent] += 1
            if variant:
                variants_all[parent].add(variant)
        if f["nct_id"] in cns_nct_ids:
            if len(parents) > 1:
                multi_serotype_cns_cnt += 1
            for parent, _ in pairs:
                serotype_cns[parent] += 1

        for hit in find_aavhu_nondigit_context(f):
            aavhuman_hits.append({
                "nct_id":  f["nct_id"],
                "title":   f["brief_title"],
                "match":   hit["match"],
                "context": hit["context"],
            })

    distinct_serotypes_all = len(serotype_all)
    distinct_serotypes_cns = len(serotype_cns)

    print()
    if aavhuman_hits:
        print(f"[AAVHUMAN DIAGNOSTIC] {len(aavhuman_hits)} 'AAVhu<letter>' hit(s) found:")
        for h in aavhuman_hits:
            print(f"  {h['nct_id']}  {h['title']}")
            print(f"    match  : {h['match']}")
            print(f"    context: {h['context']}")
    else:
        print("[AAVHUMAN DIAGNOSTIC] No 'AAVhu<letter>' false-positive hits — regex fix confirmed.")

    # ── 5. Statistics ─────────────────────────────────────────────────────────
    tier_counts     = Counter(t["cns_tier"] for t in cns_trials)
    route_counts    = defaultdict(int)
    multi_route_cnt = 0
    for t in cns_trials:
        cats = t["cns_route_categories"].split("|") if t["cns_route_categories"] else []
        for c in cats:
            route_counts[c] += 1
        if len(cats) > 1:
            multi_route_cnt += 1

    status_counts    = Counter(f["overall_status"] for f in active)
    phase_counts     = Counter(f["phase"] for f in active)
    excl_by_category = Counter(f["exclusion_category"] for f in excluded)

    cond_counter = Counter()
    for t in cns_trials:
        for cond in t["conditions"].split("; "):
            cond = cond.strip()
            if cond:
                cond_counter[cond] += 1
    top10_conditions = cond_counter.most_common(10)

    # ── 6. Write all-confirmed CSV ────────────────────────────────────────────
    all_cols = [
        "nct_id", "brief_title", "overall_status", "phase",
        "start_date", "primary_completion_date", "sponsor", "conditions",
        "intervention_names", "intervention_descriptions", "brief_summary",
        "aav_confirmation_method", "aav_confirmation_source", "matched_product_name",
        "cns_flagged", "delivery_route",
    ]
    with open(ALL_CSV, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=all_cols, extrasaction="ignore")
        w.writeheader()
        w.writerows(cns_trials + non_cns_trials)
    print(f"\nWrote: {ALL_CSV}")

    # ── 7. Write CNS-confirmed CSV ────────────────────────────────────────────
    cns_cols = [
        "nct_id", "brief_title", "overall_status", "phase",
        "start_date", "primary_completion_date", "sponsor", "conditions",
        "intervention_names", "intervention_descriptions", "brief_summary",
        "aav_confirmation_method", "aav_confirmation_source", "matched_product_name",
        "cns_tier", "cns_route_categories", "cns_triggering_terms",
        "iv_cns_confirmation_term",
    ]
    with open(CNS_CSV, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=cns_cols, extrasaction="ignore")
        w.writeheader()
        w.writerows(cns_trials)
    print(f"Wrote: {CNS_CSV}")

    # ── 7b. Write excluded CSV ────────────────────────────────────────────────
    excl_cols = [
        "nct_id", "brief_title", "overall_status", "phase",
        "start_date", "primary_completion_date", "sponsor", "conditions",
        "intervention_names", "intervention_descriptions", "brief_summary",
        "exclusion_category", "exclusion_term",
    ]
    with open(EXCLUDED_CSV, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=excl_cols, extrasaction="ignore")
        w.writeheader()
        w.writerows(excluded)
    print(f"Wrote: {EXCLUDED_CSV}")

    # ── 7c. Write serotype frequency CSV ─────────────────────────────────────
    all_serotypes_found = set(serotype_all.keys()) | set(serotype_cns.keys())
    serotype_rows = sorted(
        all_serotypes_found,
        key=lambda s: (-serotype_cns.get(s, 0), -serotype_all.get(s, 0)),
    )
    with open(SEROTYPE_CSV, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=[
            "serotype", "count_all_trials", "count_cns_trials",
            "percent_of_cns_trials", "construct_variant",
        ])
        w.writeheader()
        for s in serotype_rows:
            cnt_cns  = serotype_cns.get(s, 0)
            variants = "|".join(sorted(variants_all.get(s, set())))
            w.writerow({
                "serotype":               s,
                "count_all_trials":       serotype_all.get(s, 0),
                "count_cns_trials":       cnt_cns,
                "percent_of_cns_trials":  round(cnt_cns / total_cns * 100, 1) if total_cns else 0,
                "construct_variant":      variants,
            })
    print(f"Wrote: {SEROTYPE_CSV}")

    # ── 8. Build summary ──────────────────────────────────────────────────────
    lines = []
    lines.append(SEP)
    lines.append("AAV CNS PIPELINE v9 — QUERY SUMMARY")
    lines.append(SEP)
    lines.append(f"Date of query: {TODAY}")
    lines.append(f"Total queries run: {len(QUERY_INTRS)}")
    lines.append("")
    lines.append(f"Total raw results (query.intr set, pre-dedup)  : {total_raw}")
    lines.append(f"Supplementary query.term raw                   : {supp_raw}  ({supp_new} net new unique)")
    lines.append(f"Total unique after deduplication               : {total_deduped}")
    lines.append(f"Dropped by AAV confirmation filter         : {dropped}")
    lines.append(f"Total confirmed AAV trials                 : {total_confirmed}")
    lines.append(f"Excluded (non-CNS indication filter)       : {total_excluded}")
    lines.append(f"Active trials (post-exclusion)             : {len(active)}")
    lines.append(f"Total CNS-flagged                          : {total_cns}")
    lines.append(f"  of which proprietary_product_name        : {sum(1 for t in cns_trials if t.get('aav_confirmation_method') == 'proprietary_product_name')}")
    lines.append("")

    lines.append("-" * 50)
    lines.append("AAV Confirmation Method Breakdown (Fix 7):")
    n_generic     = sum(1 for f in confirmed if f.get("aav_confirmation_method") == "generic_terminology")
    n_proprietary = len(proprietary_confirmed)
    lines.append(f"  generic_terminology      : {n_generic}")
    lines.append(f"  proprietary_product_name : {n_proprietary}")
    lines.append("")
    if proprietary_confirmed:
        cns_nct_set_local = {t["nct_id"] for t in cns_trials}
        code_rows = {}
        for nct_id, code, citation in proprietary_confirmed:
            code_rows.setdefault(code, []).append(nct_id)
        lines.append("  Breakdown by product code:")
        for code, ncts in sorted(code_rows.items()):
            cns_cnt = sum(1 for n in ncts if n in cns_nct_set_local)
            lines.append(f"    {code:<15} {len(ncts)} confirmed  ({cns_cnt} CNS-flagged)")
        lines.append("")
        lines.append("  Per-trial detail:")
        for nct_id, code, citation in proprietary_confirmed:
            row = next((f for f in confirmed if f["nct_id"] == nct_id), None)
            cns_flag = "CNS-flagged" if nct_id in cns_nct_set_local else "not-CNS-flagged"
            title = row["brief_title"] if row else ""
            lines.append(f"  [{cns_flag}]  {nct_id}  matched={code}  {title}")
    else:
        lines.append("  None — all confirmed trials matched generic AAV terms algorithmically.")
    lines.append("")

    lines.append("-" * 50)
    lines.append("Fix 8/9 — New Trial Attribution vs v8 Baseline:")
    if v8_nct_ids:
        lines.append(f"  v8 baseline confirmed NCT IDs loaded : {len(v8_nct_ids)}")
    else:
        lines.append(f"  (v8 baseline CSV not found — attribution skipped)")
    lines.append(f"  New trials in v9 (not in v8 confirmed) : {len(new_trials)}")
    lines.append(f"  New via vesemnogene lantuparvovec text  : {len(vesemnogene_new)}")
    lines.append(f"  New via bare \\bAAV\\b only (not v8)     : {len(bare_aav_new)}")
    overlap = sum(1 for f in bare_aav_new if "vesemnogene" in " ".join([
        f["intervention_names"], f["intervention_descriptions"], f["brief_summary"]]).lower())
    lines.append(f"  (overlap — both categories)            : {overlap}")
    lines.append("")

    if vesemnogene_new:
        lines.append("  New vesemnogene trials:")
        for f in vesemnogene_new:
            cns_flag = "CNS" if f["nct_id"] in {t["nct_id"] for t in cns_trials} else "not-CNS"
            lines.append(f"    [{cns_flag}]  {f['nct_id']}  {f['brief_title']}")
        lines.append("")

    if bare_aav_new:
        lines.append(f"  New bare \\bAAV\\b-only trials (total {len(bare_aav_new)}) — first 10 for spot-check:")
        cns_nct_set_tmp = {t["nct_id"] for t in cns_trials}
        for f in bare_aav_new[:10]:
            cns_flag = "CNS" if f["nct_id"] in cns_nct_set_tmp else "not-CNS"
            corpus_snip = " ".join([f["intervention_names"],
                                    f["intervention_descriptions"],
                                    f["brief_summary"]])[:200]
            lines.append(f"    [{cns_flag}]  {f['nct_id']}  {f['brief_title']}")
            lines.append(f"           corpus: {corpus_snip!r}")
        if len(bare_aav_new) > 10:
            lines.append(f"    ... and {len(bare_aav_new) - 10} more (see CSV for full list)")
    else:
        lines.append("  No trials newly confirmed exclusively via bare \\bAAV\\b.")
    lines.append("")

    lines.append("-" * 50)
    lines.append("Per-Query Contribution  (raw | new unique):")
    for q, raw_c, new_c in contributions:
        marker = "  *" if new_c > 0 else "   "
        lines.append(f"{marker} {q:<35} {raw_c:>5} raw  |  {new_c:>4} new unique")
    lines.append("")

    lines.append("-" * 50)
    lines.append("Exclusion Breakdown by Category:")
    for cat in EXCLUSION_TERMS:
        lines.append(f"  {cat:<30} {excl_by_category.get(cat, 0)}")
    lines.append("")

    lines.append("-" * 50)
    lines.append("CNS Tier Breakdown:")
    for label in ["Tier 1", "Tier 2", "Tier 1 + Tier 2"]:
        lines.append(f"  {label:<30} {tier_counts.get(label, 0)}")
    lines.append("")

    lines.append("-" * 50)
    lines.append("CNS Delivery Route Category Breakdown:")
    for r in list(TIER1_ROUTES.keys()) + ["Intravenous"]:
        lines.append(f"  {r:<30} {route_counts.get(r, 0)}")
    lines.append(f"  {'Trials with multiple CNS routes':<30} {multi_route_cnt}")
    lines.append("")

    lines.append("-" * 50)
    lines.append("Status Breakdown  (active AAV trials, post-exclusion):")
    for s, n in sorted(status_counts.items()):
        lines.append(f"  {s:<40} {n}")
    lines.append("")

    lines.append("-" * 50)
    lines.append("Phase Breakdown  (active AAV trials, post-exclusion):")
    for p, n in sorted(phase_counts.items()):
        lines.append(f"  {p:<40} {n}")
    lines.append("")

    lines.append("-" * 50)
    lines.append("Top 10 Conditions Among CNS-Flagged Trials:")
    for rank, (cond, cnt) in enumerate(top10_conditions, 1):
        lines.append(f"  {rank:>2}. {cond} ({cnt})")
    lines.append("")

    lines.append("-" * 50)
    lines.append("Serotype Frequency Analysis:")
    lines.append(f"  Distinct serotypes found — all confirmed trials : {distinct_serotypes_all}")
    lines.append(f"  Distinct serotypes found — CNS trials only      : {distinct_serotypes_cns}")
    lines.append(f"  Trials using >1 serotype — all confirmed        : {multi_serotype_all_cnt}")
    lines.append(f"  Trials using >1 serotype — CNS only             : {multi_serotype_cns_cnt}")
    lines.append(f"  Notes: sc/ss prefixes kept distinct; rh/hu cargo suffixes normalized to")
    lines.append(f"         parent serotype (full variant in construct_variant CSV column);")
    lines.append(f"         leading 'r' recombinant prefix stripped; chimeric AAV2/9 etc. kept;")
    lines.append(f"         AAVrhXX and AAVrh.XX collapsed to canonical dot form (AAVrh.XX).")
    lines.append("")
    lines.append("  AAVHUMAN diagnostic:")
    if aavhuman_hits:
        for h in aavhuman_hits:
            lines.append(f"    MATCH  {h['nct_id']}  {h['title']}")
            lines.append(f"           match: {h['match']}")
            lines.append(f"           ctx:   {h['context']}")
    else:
        lines.append("    No 'AAVhu<letter>' hits found — regex tightening confirmed effective.")
    lines.append("")
    lines.append("  Top 10 serotypes by CNS trial count:")
    lines.append(f"  {'Serotype':<20} {'All trials':>10}  {'CNS trials':>10}  {'% of CNS':>9}")
    lines.append(f"  {'-'*20}  {'-'*10}  {'-'*10}  {'-'*9}")
    for rank, (s, cnt_cns) in enumerate(serotype_cns.most_common(10), 1):
        cnt_all = serotype_all.get(s, 0)
        pct     = round(cnt_cns / total_cns * 100, 1) if total_cns else 0
        lines.append(f"  {rank:>2}. {s:<18} {cnt_all:>10}  {cnt_cns:>10}  {pct:>8.1f}%")

    lines.append("-" * 50)
    lines.append("Fix 3 — Gaucher Type 2/3 Manual Review (neuronopathic Gaucher):")
    if gaucher_review:
        for nct_id, title, flag in gaucher_review:
            lines.append(f"  [{flag}]  {nct_id}  {title}")
    else:
        lines.append("  No neuronopathic Gaucher trials found in confirmed set.")
    lines.append("")

    lines.append("-" * 50)
    lines.append("Methodology Limitation — Trials With Empty Intervention Descriptions:")
    lines.append(f"  Count: {len(empty_iv_desc)}")
    lines.append("  These trials cannot be route-detected via intervention text.")
    lines.append("  They rely solely on brief_summary and eligibility corpus for CNS detection.")
    for nct_id, title in empty_iv_desc:
        lines.append(f"    {nct_id}  {title}")
    lines.append(SEP)

    summary_text = "\n".join(lines)
    print("\n" + summary_text)

    with open(SUMMARY_TXT, "w", encoding="utf-8") as fh:
        fh.write(summary_text + "\n")
    print(f"\nWrote: {SUMMARY_TXT}")

    # ── 9. Copy script to output folder ───────────────────────────────────────
    try:
        src = os.path.abspath(__file__)
        if os.path.isfile(src) and os.path.abspath(src) != os.path.abspath(SCRIPT_DEST):
            shutil.copy2(src, SCRIPT_DEST)
            print(f"Wrote: {SCRIPT_DEST}")
    except Exception as exc:
        print(f"Note: could not copy script ({exc})", file=sys.stderr)

    print("\nAll done.")


if __name__ == "__main__":
    main()
