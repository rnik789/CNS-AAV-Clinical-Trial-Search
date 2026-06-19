# AAV CNS Clinical Trial Identification Pipeline

## 1. Overview

This pipeline identifies adeno-associated virus (AAV) gene therapy clinical trials that target the central nervous system (brain or spinal cord) using ClinicalTrials.gov as the primary data source. Trials are first confirmed as AAV-based through a two-step keyword confirmation procedure, then classified for CNS targeting through a two-tier route- and disease-based logic. The pipeline was executed against the live ClinicalTrials.gov REST API v2 on **2026-06-18** and identified **340 confirmed AAV trials**, of which **92 are classified as CNS-targeted**. The complete pipeline is implemented in the single script `aav_cns_pipeline_v9_2026-06-18.py`.

---

## 2. Data Source and Query Strategy

### API Endpoint

All data were retrieved from the ClinicalTrials.gov REST API v2:

```
https://clinicaltrials.gov/api/v2/studies
```

Responses were requested in JSON format with a page size of 1,000 records per request.

### Status Filter

A `filter.overallStatus` parameter was applied to all queries, restricting results to trials with one of the following statuses:

```
COMPLETED, TERMINATED, ACTIVE_NOT_RECRUITING, RECRUITING, ENROLLING_BY_INVITATION, SUSPENDED
```

This filter excludes trials with status `NOT_YET_RECRUITING`, `WITHDRAWN`, and `UNKNOWN`. Trials with these statuses may include legitimate historical AAV CNS programs that fall outside the scope of this search.

### Query Strategy

Sixty-two (62) primary queries were issued using the `query.intr` parameter, which searches the intervention name and related fields in the ClinicalTrials.gov index. Queries are organized into four logical groups:

**Generic vector name queries (16):**
`adeno-associated virus`, `adeno associated virus`, `rAAV`, `AAV1`, `AAV2`, `AAV3`, `AAV4`, `AAV5`, `AAV6`, `AAV7`, `AAV8`, `AAV9`, `AAVrh`, `AAVhu`, `scAAV`, `ssAAV`

**Disease-specific construct and program name queries (33):**
`AAV2-hAADC`, `AAV2-GDNF`, `AAV2-GAD`, `AAV2-NGF`, `AAV2-BDNF`, `AAV9-SMN`, `AAV-GAD`, `CERE-110`, `CERE-120`, `RGX-121`, `RGX-314`, `AMT-130`, `AMT-162`, `TSHA-102`, `FBX-101`, `AB-1001`, `SGT-212`, `GC101`, `SKG0201`, `IPS101A`, `JAG201`, `EXG001`, `BBP-812`, `UX111`, `AAVhu68`, `PBGM01`, `MVX-220`, `NSR-REP1`, `ABO-101`, `ABO-102`, `ABO-50`, `ABO-202`

**Explicit AAVrh-numbered variant queries (5):**
`AAVrh8`, `AAVrh10`, `AAVrh74`, `AAVrh20`, `AAVrh39`

These are issued separately because the ClinicalTrials.gov API does not perform prefix-matching: a `query.intr=AAVrh` search does not return records indexed under the token `AAVrh10` or `AAVrh74`. Direct diagnostic testing confirmed this behavior.

**Proprietary product code queries (9):**
`VY-AADC01`, `VY-AADC02`, `VY-AADC`, `PR001`, `PR006`, `LY3884961`, `LY3884963`, `NGN-101`, `vesemnogene`

These product codes were verified as AAV-based via Nair et al. 2024 (*Neurotherapeutics*) or Sharma, Joshi, Kumar (2025) (*Neurosci.*). Trials using these codes contain no generic AAV terminology in API-searchable fields and would be missed by all other query strategies.

### Supplementary Query Pass

After all 62 `query.intr` queries, a supplementary pass was run using `query.term=adeno-associated+virus`, which searches across all indexed trial fields rather than intervention fields only. This pass returned **265 raw results**, contributing **29 net-new unique trials** beyond those already captured by the `query.intr` set.

### Pagination

All queries use token-based pagination via the `nextPageToken` field in the API response. Each query iterates pages until `nextPageToken` is absent, ensuring exhaustive retrieval regardless of result set size.

### Raw Results and Deduplication

The 62 `query.intr` queries returned **904 raw results** in aggregate (counting duplicates). After merging with the supplementary `query.term` results and deduplicating by NCT ID, **438 unique trials** were retained for downstream processing. Per-query contribution tracking was applied: each query reports the count of raw results and the count of NCT IDs that were new to the dataset at the time of that query.

---

## 3. AAV Confirmation Methodology

After deduplication, each trial undergoes a sequential two-step confirmation to verify that it involves an AAV vector. Trials failing both steps are dropped from the dataset. The confirmation logic searches a **corpus constructed from three fields**: `intervention_names`, `intervention_descriptions`, and `brief_summary`.

### Step 1 — Generic Terminology (Algorithmic)

The corpus is searched using pre-compiled patterns (`AAV_GENERIC_COMPILED`). All terms use case-insensitive substring matching except `AAV`, which is matched with a `\bAAV\b` word-boundary anchor to prevent false-positive matches within `rAAV`, `scAAV`, `AAV9`, etc., while still catching standalone tokens such as `AAV-hSMN1` or `AAV-based`:

```
adeno-associated virus, adeno-associated viral,
rAAV, scAAV, ssAAV,
AAV1, AAV2, AAV3, AAV4, AAV5, AAV6, AAV7, AAV8, AAV9,
AAVrh, AAVhu,
onasemnogene, voretigene, eladocagene,
valoctocogene, etranacogene, delandistrogene,
AAV  (word-boundary anchored: \bAAV\b)
```

A match at Step 1 sets:
- `aav_confirmation_method` = `generic_terminology`
- `aav_confirmation_source` = `registry_text`
- `matched_product_name` = *(empty)*

Step 2 is not executed for trials confirmed at Step 1.

### Step 2 — Proprietary Product Name (Literature-Verified)

Step 2 runs only when Step 1 finds no match. The corpus is searched case-insensitively for each of the following proprietary product codes, which have been independently verified as AAV-based via the cited literature:

| Product Code | Source Citation |
|---|---|
| VY-AADC01 | Nair et al. 2024, Neurotherapeutics |
| VY-AADC02 | Nair et al. 2024, Neurotherapeutics |
| VY-AADC   | Nair et al. 2024, Neurotherapeutics |
| PR001     | Nair et al. 2024, Neurotherapeutics |
| PR006     | Nair et al. 2024, Neurotherapeutics |
| LY3884961 | Nair et al. 2024, Neurotherapeutics |
| LY3884963 | Nair et al. 2024, Neurotherapeutics |
| NGN-101   | Nair et al. 2024, Neurotherapeutics |
| vesemnogene lantuparvovec | Sharma, Joshi, Kumar (2025), Neurosci. |

The first matching code sets:
- `aav_confirmation_method` = `proprietary_product_name`
- `aav_confirmation_source` = the citation string for that code
- `matched_product_name` = the matched product code

Trials confirmed exclusively via Step 2 are unambiguously proprietary-only in registry text: by construction, no generic AAV terminology appeared in any API-searchable field.

### Confirmation Results

Of the 438 unique trials entering this step, **340 were confirmed** as AAV trials and **98 were dropped**:

- Confirmed via `generic_terminology`: **332**
- Confirmed via `proprietary_product_name`: **8**
  - VY-AADC01: 3 trials (NCT01973543, NCT03065192, NCT03733496)
  - PR001: 1 trial (NCT04411654)
  - LY3884961: 1 trial (NCT04127578)
  - LY3884963: 1 trial (NCT04408625)
  - NGN-101: 1 trial (NCT05228145)
  - vesemnogene lantuparvovec: 1 trial (NCT07265232)

Both output CSVs (`aav_all_confirmed_v9_2026-06-18.csv` and `aav_cns_confirmed_v9_2026-06-18.csv`) carry the columns `aav_confirmation_method`, `aav_confirmation_source`, and `matched_product_name` for every trial.

---

## 4. Exclusion Criteria

After AAV confirmation, trials whose disease indication matches a non-CNS category are removed from the active dataset. Exclusion is evaluated by searching a corpus of `conditions` concatenated with `brief_title` for any term in the following category lists:

| Category | Exclusion Terms |
|---|---|
| Blood disorders | hemophilia, haemophilia, factor VIII, factor IX, von Willebrand |
| Ocular only | retinal dystrophy, retinitis pigmentosa, macular degeneration, achromatopsia, leber congenital amaurosis, RPE65, choroideremia, retinoschisis, Stargardt, RPGR |
| Liver/metabolic | ornithine transcarbamylase, OTC deficiency, alpha-1 antitrypsin, AAT deficiency, Crigler-Najjar, phenylketonuria, PKU, glycogen storage disease, Wilson disease, familial hypercholesterolemia |
| Muscle | myotubular myopathy, Duchenne muscular dystrophy, DMD, limb girdle muscular dystrophy, LGMD, Becker muscular dystrophy |
| Cardiac only | cardiomyopathy of Friedreich, hypertrophic cardiomyopathy, heart failure |
| Other peripheral | parotid, salivary gland, alpha-1 antitrypsin lung, pulmonary, arthritis, joint |

### Tier 1 Override

Exclusion is suppressed when the trial's delivery corpus (see Section 5) contains a Tier 1 CNS delivery term. This override preserves trials that target both a non-CNS indication and the CNS—for example, a systemic disease program that also includes direct brain delivery. A trial is excluded only when a non-CNS disease term matches **and** no Tier 1 CNS delivery evidence is present.

Of the 340 confirmed AAV trials, **165 were excluded** under this step, leaving **175 active trials** for CNS classification. Exclusion by category: Blood disorders (42), Ocular only (74), Liver/metabolic (12), Muscle (22), Cardiac only (9), Other peripheral (6).

---

## 5. CNS Classification Methodology

CNS targeting is determined by a two-tier route-based logic applied to the 175 active trials. Tiers are not mutually exclusive: 9 trials satisfied both.

### Delivery Corpus

The delivery corpus searched for route-related evidence consists of five concatenated fields: `intervention_names`, `intervention_descriptions`, `brief_summary`, `eligibilityCriteria`, and `designInfo.interventionModelDescription` (the last two extracted from the API's `eligibilityModule` and `designModule` respectively). This expanded corpus ensures that route language appearing in eligibility or design narrative is captured even when absent from intervention-specific fields.

### Tier 1 — Unambiguous CNS Delivery Routes

A trial is Tier 1 CNS-flagged when the delivery corpus contains any term from the following route categories:

**Intrathecal:** intrathecal, intrathecally, lumbar puncture, lumbar intrathecal, injection into the cerebrospinal fluid, CSF administration, IT administration

**Intraparenchymal:** intraparenchymal, stereotaxic, stereotactic, intracerebral, intrathalamic, intrastriatal, intra-striatal, subthalamic nucleus, convection enhanced delivery, direct injection into the brain, direct brain injection, intraputaminal, into the putamen, into putamen, putaminal injection, burr hole, burr holes, into the brain, into brain

**Intracerebroventricular:** intracerebroventricular, intraventricular, ICV

**Cisterna magna:** cisterna magna, intracisternal, intracisternally

A Tier 1 match alone is sufficient for CNS classification. **77 trials** were classified as Tier 1 only, with the following route category distribution: Intrathecal (49), Intraparenchymal (31), Intracerebroventricular (13), Cisterna magna (10). Twenty-one trials matched more than one route category.

### Tier 2 — Intravenous Delivery with CNS Disease Confirmation

A trial is Tier 2 CNS-flagged when two conditions are both met:

**Condition A — IV delivery:** The delivery corpus (same five-field corpus as Tier 1) contains any of the following terms: `intravenous`, `intravenously`, `IV infusion`, `systemic administration`, `systemic delivery`, `IV administration`.

**Condition B — CNS disease:** A separate disease corpus, consisting only of `conditions` concatenated with `brief_title`, contains at least one CNS disease indicator from the following list:

```
spinal muscular atrophy, SMA, Huntington, Parkinson, Alzheimer,
ALS, amyotrophic lateral sclerosis, Batten, leukodystrophy,
mucopolysaccharidosis, gangliosidosis, Rett, AADC,
aromatic L-amino acid, giant axonal, Canavan, seizure,
spinocerebellar ataxia, SCA1, SCA2, SCA3, SCA6, SCA7,
cerebellar ataxia, mesial temporal lobe epilepsy, CLN,
frontotemporal dementia, FTD-GRN, glioblastoma, glioma,
cerebellar, motor neuron, dopaminergic, neuronal ceroid,
metachromatic, Krabbe, spinal cord, brain, central nervous system,
CNS, neurodegenerative, neurological, Angelman syndrome, UBE3A,
Hunter syndrome, MPS II, mucopolysaccharidosis II, Sanfilippo,
Hurler syndrome, MPS I, Tay-Sachs, Sandhoff
```

**6 trials** were classified as Tier 2 only; **8 trials** satisfied both Tier 1 and Tier 2 (classified as `Tier 1 + Tier 2`).

### Word-Boundary Protection for Short Terms

Terms of four characters or fewer, or short acronyms that are vulnerable to substring collision (e.g., `SMA` matching within `plasma`, `CLN` matching within `clinical`, `ALS` matching within `false`), are matched using `\bTERM\b` word-boundary anchors rather than plain substring search. The following terms receive word-boundary treatment:

```
SMA, ALS, SCA1, SCA2, SCA3, SCA6, SCA7, CLN, CNS, MPS I, MPS II
```

All other Tier 2 disease terms are matched as plain case-insensitive substrings. This distinction was established in pipeline version 6 after the trial NCT06533579 (VNX-101, a hematology trial) was incorrectly CNS-flagged due to `SMA` appearing as a substring in hematology-context text.

Word-boundary patterns are pre-compiled at import time (`TIER2_CNS_DISEASE_PATTERNS`) to avoid per-trial recompilation overhead. The `CNS` co-term in the conditional Pompe disease check is similarly anchored with `\bCNS\b`.

### Conditional Disease Terms

Pompe disease is included in the Tier 2 logic conditionally: it triggers CNS classification only when the combined disease corpus and delivery corpus also contain the term `CNS` (word-boundary-anchored) or `neurological`. This prevents Pompe disease trials targeting only cardiac or skeletal muscle from being incorrectly flagged as CNS.

---

## 6. Serotype Extraction

Serotype annotations are extracted from a three-field corpus: `intervention_names`, `intervention_descriptions`, and `brief_summary`. A single regular expression captures the range of AAV serotype notations present in registry text:

```
\b(?:r|sc|ss)?AAV(?:
    [- .]?(?:rh|hu)[- .]?\d+   # AAVrh10, AAVrh.10, AAVhu68, AAV-rh10
  | [- .]?php[- .]?\w*          # AAV-PHP.B, AAV-PHP.eB, AAV-PHP.S
  | [- .]?(?:DJ|2i8|B1)(?:/\d+)?  # AAV-DJ, AAV-DJ/8, AAV2i8
  | \d+(?:[./]\d+)?             # AAV2, AAV9, AAV2/5, AAV2/9
)\b
```

### Normalization Rules

Raw matches are normalized before counting:

- **Recombinant prefix (`r`):** The leading `r` in `rAAVhu68` is stripped. It denotes production method, not capsid identity, so `rAAVhu68` normalizes to `AAVhu68`. The `sc` and `ss` prefixes (self-complementary, single-stranded) are retained as they denote distinct vector configurations with different packaging capacity.
- **AAVrh canonicalization:** The `rh` family appears in registry text as both `AAVrh10` and `AAVrh.10`. Both forms are normalized to the canonical dot notation (`AAVrh.10`) so that, for example, `AAVrh10` and `AAVrh.10` are counted as the same serotype.
- **Cargo/transgene suffixes:** Transgene cargo suffixes on `rh` and `hu` family variants (e.g., `AAVrh.10CUARSA`, `AAVhu68hFXN`) are stripped to produce a normalized parent serotype (`AAVrh.10`, `AAVhu68`) for frequency counting. The full construct string is preserved in the `construct_variant` column of the serotype frequency CSV for audit purposes.
- **Chimeric notations:** Notations such as `AAV2/9` and `AAV2/5` are retained as-is, as they represent distinct dual-capsid constructs rather than typographic variants of a single serotype.
- **Case normalization:** All normalized serotypes are uppercased for consistent grouping.

### False-Positive Guard

A secondary diagnostic pattern (`AAVhu<letter>`) detects cases where the `rh/hu` regex arm might match non-digit suffixes (e.g., `AAVhuman`). All such matches are surfaced to a diagnostic log. In the 2026-06-18 run, no such false-positive hits were found.

### Frequency Counting

Serotypes are counted on a per-trial basis: each distinct (parent serotype, construct variant) pair is counted at most once per trial regardless of how many times it appears in that trial's text. The serotype frequency CSV reports counts across all 340 confirmed trials and separately across the 92 CNS-flagged trials.

The 2026-06-18 run identified 20 distinct serotypes across all confirmed trials and 9 distinct serotypes among CNS-flagged trials. The top serotype in CNS trials was AAV9 (31 trials, 33.7% of CNS-flagged), followed by AAV2 (12 trials, 13.0%) and scAAV9 (8 trials, 8.7%).

---

## 7. Known Limitations

### Empty Intervention Descriptions

The pipeline's delivery route classification depends on text present in registry fields. Seventeen confirmed AAV trials have empty `intervention_descriptions` fields in the API response. These trials cannot be route-classified via intervention text and rely entirely on `brief_summary` and `eligibilityCriteria` for delivery evidence. The following trials were identified in the 2026-06-18 run as having empty `intervention_descriptions`:

NCT04909346, NCT01801709, NCT00004533, NCT00985517, NCT06311708, NCT00430768, NCT00252850, NCT00515710, NCT02302690, NCT00076557, NCT00195143, NCT01416467, NCT01301573, NCT03733496, NCT03602820, NCT07264166, NCT04272554

Three proprietary-confirmed VY-AADC01 trials (NCT01973543, NCT03065192, NCT03733496) are confirmed AAV trials that are not classified as CNS-targeted despite being intraputaminal AADC delivery programs. These trials contain no delivery route language in any API-accessible field sufficient to trigger Tier 1 classification. NCT03733496 also appears in the empty `intervention_descriptions` list above. This represents a structural ceiling of text-based route detection: the pipeline cannot flag CNS delivery when the registry record does not describe it.

### Status Filter Exclusions

Trials with status `NOT_YET_RECRUITING`, `WITHDRAWN`, or `UNKNOWN` are excluded by the API query filter. (`SUSPENDED` is included in the filter and appears in the confirmed set.) This may omit legitimate AAV CNS programs that have not yet initiated enrollment. The total number of such excluded trials was not quantified in this pipeline.

### Keyword Detection Ceiling

The AAV confirmation step depends on the presence of recognizable AAV terminology or verified proprietary product codes in the trial's intervention fields. Trials that use unrecognized proprietary codes with no accompanying generic AAV language, and whose codes do not appear in the verified product lists or any of the 62 `query.intr` search terms, will be missed. The scope of this gap is unknown by definition.

### Cross-Validation Against Nair et al. 2024

The Nair et al. 2024 systematic review (*Neurotherapeutics*) and Sharma, Joshi, Kumar (2025) (*Neurosci.*) were used as independent cross-references for proprietary product verification. The 9 product codes in `AAV_PROPRIETARY_PRODUCT_NAMES` were drawn from these reviews. The cross-validation established that 8 trials in this pipeline's confirmed set were identifiable only through these proprietary codes and would have been dropped by the algorithmic confirmation step alone. The trial NCT04903288 (Kebilidi/eladocagene pivotal trial) was identified as an expected entry based on Nair et al. but was not returned by any of the 62 `query.intr` queries or the supplementary `query.term` pass in the 2026-06-18 run; its absence may reflect a status outside the filter, an intervention name that does not match any query token, or a ClinicalTrials.gov indexing gap.

---

## 8. Reproducibility

The complete pipeline is implemented as a single self-contained Python script:

```
aav_cns_pipeline_v9_2026-06-18.py
```

version-controlled in this git repository. The script requires the `requests` library and Python 3.

Running the script against the live ClinicalTrials.gov API on a different date may yield different results due to new trial registrations, status changes that move trials into or out of the filtered statuses, or updates to registry records that add or remove text relevant to AAV confirmation or CNS classification. The 2026-06-18 results are fixed in the versioned output CSVs and summary text file produced by that run:

- `aav_all_confirmed_v9_2026-06-18.csv` — all 340 confirmed AAV trials
- `aav_cns_confirmed_v9_2026-06-18.csv` — the 92 CNS-flagged trials
- `aav_cns_confirmed_v9_2026-06-18_no_explicit_route.csv` — the 92 CNS-flagged trials with "Explicit CNS delivery" stripped from `cns_route_categories` (produced by `v9_route_cleanup.py`)
- `aav_excluded_v9_2026-06-18.csv` — the 165 non-CNS-excluded trials
- `aav_serotype_frequency_v9_2026-06-18.csv` — serotype frequency table
- `aav_summary_v9_2026-06-18.txt` — full run summary with per-query contribution breakdown

---

## 9. Landscape Analysis

Two additional scripts produce summary breakdowns and charts from the 92 CNS-flagged trials.

### `v9_route_cleanup.py`

Strips `"Explicit CNS delivery"` from the `cns_route_categories` field of `aav_cns_confirmed_v9_2026-06-18.csv` and saves the result as `aav_cns_confirmed_v9_2026-06-18_no_explicit_route.csv`. Seven trials were affected. The script also regenerates `route_breakdown.csv` and `route_breakdown.png` in `cns_landscape_analysis/` from the cleaned file.

### `cns_landscape_analysis.py`

Generates seven breakdowns from the 92 CNS-flagged trials (v9), outputting one CSV and one chart per section into `cns_landscape_analysis/`:

| Section | Output CSV | Chart |
|---|---|---|
| 1. Trial status | `status_breakdown.csv` | `status_breakdown.png` |
| 2. Trial phase | `phase_breakdown.csv` | `phase_breakdown.png` |
| 3. Start year trend | `start_year_breakdown.csv` | `start_year_trend.png` |
| 4. Disease condition | `condition_breakdown.csv` | `condition_breakdown.png` |
| 5. Delivery route | `route_breakdown.csv` | `route_breakdown.png` |
| 6. Serotype | `serotype_breakdown.csv` | `serotype_breakdown.png` |
| 7. Sponsor | `sponsor_breakdown.csv` | `sponsor_breakdown.png` |

A raw condition value-counts file (`raw_condition_counts.csv`) is also written before synonym mapping is applied.

**Condition mapping:** Disease name synonyms are resolved via a manual mapping dict (e.g., "Parkinson Disease", "Idiopathic Parkinson's Disease", and "PD" all map to "Parkinson's Disease"; all SMA registry variants map to "Spinal Muscular Atrophy"; all CLN subtypes map to "Neuronal Ceroid Lipofuscinosis (Batten Disease)"). Unmapped strings pass through unchanged and are printed to console for review.

**Route counting:** Routes are counted independently per trial from the cleaned (`_no_explicit_route`) file; a trial with two routes contributes to both counts. Total route occurrences (118) exceed the trial count (92) by design.

**Sponsor normalization:** Whitespace is collapsed and values are grouped case-insensitively; the most common casing variant is used as the display form. Explicit alias merges beyond whitespace/case (e.g., "Novartis Pharmaceuticals" → "Novartis Gene Therapies") are applied via a separate `SPONSOR_MAP` dict and logged to console.
