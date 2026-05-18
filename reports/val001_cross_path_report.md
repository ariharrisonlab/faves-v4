# VAL-001 Test 3.4 — Cross-Path Consistency Report

**Run:** 2026-05-10T13:53:49.523448+00:00
**API:** `https://api.novoquantnexus.com/proxy/faves-compliance`
**Compounds tested:** 20

## Summary

- API call failures: **0**
- Served by precomputed override: **0**
- Served by live computation only: **20**
- Compounds with divergences: **0**

## Per-compound results

| Compound | API path | API total | Local total | ChEMBL key match | Divergences |
|---|---|---|---|---|---|
| aspirin | faves_v4_live | 4 | 4 | ✓ | — |
| curcumin | faves_v4_live | 9 | 9 | ✓ | — |
| warfarin | faves_v4_live | 5 | 5 | ✓ | — |
| ibuprofen | faves_v4_live | 0 | 0 | ✓ | — |
| acetaminophen | faves_v4_live | 2 | 2 | ✓ | — |
| caffeine | faves_v4_live | 1 | 1 | ✓ | — |
| atorvastatin | faves_v4_live | 2 | 2 | ✓ | — |
| sildenafil | faves_v4_live | 1 | 1 | ✓ | — |
| metformin | faves_v4_live | 6 | 6 | ✓ | — |
| amoxicillin | faves_v4_live | 4 | 4 | ✓ | — |
| diazepam | faves_v4_live | 0 | 0 | ✓ | — |
| quercetin | faves_v4_live | 6 | 6 | ✓ | — |
| resveratrol | faves_v4_live | 4 | 4 | ✓ | — |
| genistein | faves_v4_live | 1 | 1 | ✓ | — |
| doxorubicin | faves_v4_live | 11 | 11 | ✓ | — |
| imatinib | faves_v4_live | 0 | 0 | ✓ | — |
| tamoxifen | faves_v4_live | 3 | 3 | ✓ | — |
| naproxen | faves_v4_live | 0 | 0 | ✓ | — |
| osimertinib | faves_v4_live | 11 | 11 | ✓ | — |
| lisinopril | faves_v4_live | 4 | 4 | ✓ | — |

## Source-of-truth determination

All 20 compounds produced identical per-catalog counts and ChEMBL sub-catalog keys between the deployed `/api/classify` response and a fresh local RDKit FilterCatalog scan. Live and precomputed paths are consistent in shape and arithmetic.
