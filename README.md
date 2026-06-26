# FAVES V4 — Supporting Information

[![License: Apache 2.0](https://img.shields.io/badge/code-Apache_2.0-blue.svg)](LICENSE) [![Data License: CC BY 4.0](https://img.shields.io/badge/data-CC_BY_4.0-blue.svg)](LICENSE-DATA)
[![DOI](https://img.shields.io/badge/Zenodo-pending-lightgrey.svg)](#)
[![ChemRxiv](https://img.shields.io/badge/ChemRxiv-pending-lightgrey.svg)](#)

Benchmark dataset, validation reports, governance documents, and reproducer scripts for:

> **Harrison, A.** *FAVES V4: A Production-Deployed Compliance and Liability Screening System for Drug Discovery at 122M-Compound Scale, Aligned to HHS/ONC FAVES and NIST AI 600-1.* ChemRxiv (2026).
> *(DOI to be added on submission)*

FAVES V4 is a four-layer cheminformatics screening system aligned to two complementary federal frameworks: the U.S. HHS/ASTP/ONC FAVES rubric (**F**air, **A**ppropriate, **V**alid, **E**ffective, **S**afe) for trustworthy AI in healthcare, and the NIST AI Risk Management Framework Generative AI Profile (NIST AI 600-1), which names chemical and biological design tools as a distinct CBRN risk class. It returns regulatory status, structural alerts (1,585 SMARTS via RDKit FilterCatalog), pharmacokinetic classification (BOILED-Egg), and prior-art disclosure (InChIKey resolution) per molecule with per-field provenance metadata.

**Live API:** [https://api.novomcp.com](https://api.novomcp.com) — `POST /api/classify` accepts `{"smiles": "..."}` and returns the full ClassifyResponse.

**Companion v3 release:** [`ariharrisonlab/faves-v3-benchmark`](https://github.com/ariharrisonlab/faves-v3-benchmark) — V3 paper's regulatory benchmark (102 compounds), preserved unchanged for citation continuity. The 102-compound v3 ground truth is also vendored here (`data/v3_ground_truth_102_compounds.csv`) for self-contained reproduction.

---

## Repository structure

```
faves-v4/
├── README.md                                    ← this file
├── LICENSE                                      ← Apache 2.0 (code in scripts/)
├── LICENSE-DATA                                 ← CC BY 4.0 (data/ + reports/)
├── NOTICE                                       ← Apache 2.0 attribution propagation
├── CITATION.cff                                 ← machine-readable citation metadata
│
├── data/                                        ← compound benchmarks (browseable)
│   └── v3_ground_truth.csv                     (S1: 102 compounds — 47 DEA + 46 FDA + 9 negative)
│
├── scripts/                                     ← reproducer code
│   ├── v4_test_compounds.py                    (S2: 200-compound fixture — 102 v3 + 50 novel + 48 edge cases; loads v3_ground_truth.csv from ../data/)
│   ├── faves_v4_validation.py                  (executes Tests 3.1, 3.2, 3.5, 3.6)
│   ├── val001_cross_path_check.py              (Test 5.5 / Part B: cross-RDKit reproducibility)
│   └── latency_benchmark.py                    (Table 6: bimodal latency distribution)
│
└── reports/                                     ← execution artifacts (S9, S10)
    ├── faves_v4_validation_report_20260510_102603.json    (341/342 pre-fix; cortisol whitelist defect → 362/362 post-fix)
    ├── val001_cross_path_report.{md,json}                 (20/20 external reproducibility)
    └── latency_report.{md,json}                            (90 steady-state samples)
```

**Governance documents (S6) are not redistributed in this repository.** The seven NMCP-FAVES governance documents (IUS-001, RSC-001, VAL-001, AUD-001, CCM-001, UIG-001, FAVES-V4 Technical Overview) are summarized in Section 4 of the manuscript with NIST AI 600-1 action-ID anchors and worked change-control examples reproduced inline. Full operational versions of each document are available from the corresponding author (ari@novomcp.com) or at [novomcp.com/governance](https://novomcp.com/governance) on request.

---

## Mapping (paper SI item → file)

| SI item | Description | File(s) |
|---|---|---|
| **S1** | V3 benchmark (102 compounds) | `data/v3_ground_truth.csv` |
| **S2** | V4 field-completeness 200-compound set | `scripts/v4_test_compounds.py` (loads S1 from `data/`) |
| **S3** | BOILED-Egg reference drugs (20) | embedded in `scripts/faves_v4_validation.py` (`BOILED_EGG_REFERENCE`) |
| **S4** | Prior-art test compounds (15 disclosed + 5 novel) | embedded in `scripts/faves_v4_validation.py` (`DISCLOSED_COMPOUNDS`, `NOVEL_SMILES`) |
| **S5** | Cross-RDKit-version test compounds (20) | embedded in `scripts/val001_cross_path_check.py` (`TEST_COMPOUNDS`) |
| **S6** | Governance documents bundle (7 docs) | summarized in §4 of the manuscript; full operational versions available from corresponding author |
| **S7** | RDKit Implementation Gotchas | summarized in §3.2 of the manuscript; full operational version available from corresponding author |
| **S8** | Architecture diagram source (SVG) | *(in the paper; not duplicated here)* |
| **S9** | Validation execution record | `reports/faves_v4_validation_report_20260510_102603.json` |
| **S10** | Cross-path + latency reports | `reports/val001_cross_path_report.{md,json}` and `reports/latency_report.{md,json}` |

---

## Reproducing the validation

The deployed FAVES V4 service at `https://api.novomcp.com` is the canonical reproducer. With an API key, you can re-run the full validation in under 30 minutes.

### Setup

```bash
git clone https://github.com/ariharrisonlab/faves-v4.git
cd faves-v4
pip install rdkit requests httpx
export NOVOMCP_API_KEY=nmcp_xxx   # request at https://api.novomcp.com
```

### Run the validation

```bash
# Tests 3.1 (V3 regulatory, 102 compounds), 3.2 (V4 field completeness, 200),
# 3.5 (BOILED-Egg, 20), 3.6 (prior art, 20). ~25 min.
python3 scripts/faves_v4_validation.py
# Produces: faves_v4_validation_report_<timestamp>.json
# Expected: 362/362 passing (post-fix). 341/342 if you somehow hit a pre-fix server.
```

### Cross-path consistency

```bash
# Test 5.5 / Part B: external reproducibility (deployed RDKit vs local RDKit). ~1 min.
python3 scripts/val001_cross_path_check.py \
  --api-url https://api.novoquantnexus.com/proxy/faves-compliance \
  --output val001_report.md
# Expected: 20/20 external reproducibility, identical per-catalog counts and ChEMBL keys.
```

### Latency benchmark

```bash
# Table 6: bimodal latency distribution. ~5 min.
python3 scripts/latency_benchmark.py \
  --api-url https://api.novoquantnexus.com/proxy/faves-compliance \
  --calls-per-compound 10 \
  --output latency_report.md
# Expected: ~450 ms (cached prior-art) vs ~3170 ms p50 (PubChem PUG-REST fallback).
```

---

## Governance framework summary

FAVES V4 is operationalized against the HHS/ASTP/ONC FAVES rubric. Each letter maps to a system property, classified by GAMP/CSA framework, and demonstrated against published evidence:

| Letter | HHS/ONC meaning | Operationalized in FAVES V4 |
|---|---|---|
| **F** Fair | Bias avoidance, equity by design | All SMARTS catalogs from public literature; identical rule application; no proprietary classification logic |
| **A** Appropriate | Suitability for the intended task | Intended Use Statement (IUS-001) with six explicit exclusions; informational early-discovery screening only |
| **V** Valid | Validated for the intended use case | Risk-based GAMP classification (RSC-001) + predetermined acceptance criteria (VAL-001); 362/362 post-fix across 5 tests |
| **E** Effective | Improves or does not harm outcomes | 20/20 cross-RDKit-version reproducibility; production deployment metrics (122M docs, ~5 ms compute, real-world latency) |
| **S** Safe | Prevents harm, ensures security | Three-layer audit trail (AUD-001) with two non-skippable layers; V3 regulatory layer + prior-art layer prevent harm; change-control SOP (CCM-001) with worked rollback examples |

See Section 4 of the manuscript for the detailed operationalization. Full governance documents available from the corresponding author.

---

## Pre-fix vs post-fix validation

The validation report (`reports/faves_v4_validation_report_20260510_102603.json`) was produced **before** a production whitelist defect was corrected. The defect: an incorrect cortisol SMILES in the deployed FDA whitelist (omitted the 11β-hydroxyl, producing wrong InChIKey `WHBHBVVOGNECLV` instead of correct PubChem CID 5754 InChIKey `JYGXADMDTFJGBT`). Real cortisol fell through to scaffold detection, where the steroid backbone matched the V3 anabolic-steroid SMARTS pattern.

- **Pre-fix:** 341/342 passing (single failure: cortisol false positive — 101/102 on Test 3.1)
- **Post-fix:** 362/362 passing (cortisol now correctly whitelist-cleared)

Both numbers are reported transparently in the manuscript. The defect was surfaced by expanding Test 3.1 from a 84-compound subset to the full 102-compound published benchmark — itself an argument for honest-scoped validation. The change-control loop (CCM-001) governed the fix; this is one of two worked examples in the change log.

A post-fix validation report will be added to this repo when re-executed against the deployed corrected service.

---

## Citing this repository

If you use the V4 benchmark, scripts, or governance documents, please cite:

```bibtex
@article{harrison2026favesv4,
  author = {Harrison, Ari},
  title = {FAVES V4: A Production-Deployed Compliance and Liability Screening System for Drug Discovery at 122M-Compound Scale, Aligned to HHS/ONC FAVES and NIST AI 600-1},
  journal = {ChemRxiv},
  year = {2026},
  note = {DOI to be added}
}
```

For the benchmark/code specifically, also cite the Zenodo release (DOI minted on tagged release).

GitHub's "Cite this repository" widget reads `CITATION.cff` for one-click citation export.

---

## License

This repository uses a **dual license**:

| Component | License | File |
|---|---|---|
| Code (`scripts/`) | Apache License 2.0 | [LICENSE](LICENSE) |
| Data (`data/v3_ground_truth.csv`) and execution reports (`reports/`) | Creative Commons Attribution 4.0 International (CC BY 4.0) | [LICENSE-DATA](LICENSE-DATA) |
| Attribution propagation | per Apache §4(d) | [NOTICE](NOTICE) |

**Both licenses permit commercial use; both require attribution.** When you use the benchmark or scripts in a paper, presentation, regulatory filing, or downstream product, cite the manuscript (BibTeX in [CITATION.cff](CITATION.cff) or [LICENSE-DATA](LICENSE-DATA)). For Apache-licensed code, also reproduce the [NOTICE](NOTICE) file in your derivative work.

The compound data (`data/v3_ground_truth.csv`) is derived from public sources: DEA Schedule I–V listings, FDA-approved drug listings, and PubChem canonical SMILES. The original public sources remain under their respective terms; this release contributes the curated benchmark composition, parse-verified SMILES, and ground-truth labeling decisions documented in Section 3.5 of the manuscript.

**Governance documents** (IUS-001, RSC-001, VAL-001, AUD-001, CCM-001, UIG-001, FAVES-V4 Technical Overview) are not included in this repository and are not covered by either of these licenses. See [novomcp.com/governance](https://novomcp.com/governance) or contact the corresponding author.

---

## Contact

Ari Harrison — `ari@novomcp.com` — [ORCID 0009-0006-5836-7528](https://orcid.org/0009-0006-5836-7528)
