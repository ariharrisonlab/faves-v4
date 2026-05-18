"""VAL-001 Test 3.4 — Cross-path consistency spot-check (FAVES V4)

Verifies that the precomputed the precomputed structural-alerts cache override produces the same
structural_alert_summary as the live RDKit FilterCatalog scan.

Why this exists
---------------
/api/classify always runs the live FilterCatalog scan, then optionally
swaps in precomputed Cosmos data when the SMILES → CID Redis lookup hits
(see the precomputed-override block). The two paths must produce identical
output. The risk surface is the chembl_alerts dict keying — live builds
keys by splitting "ChEMBL23_Dundee" → "dundee", precomputed reads
alerts_doc["chembl_alerts"] verbatim from the Cosmos loader. If the loader
wrote different keys (e.g., "ChEMBL23_Dundee", or omitted empty buckets),
the two response shapes diverge.

What this script does
---------------------
For each test SMILES, this script makes TWO calls and a local computation:

1. **Default call:** POST /api/classify (no flag). Returns either
   `version: faves_v4_live` (cache miss) or `faves_v4_precomputed` (cache
   hit, override applied).
2. **Forced live call:** POST /api/classify?force_live=true. Bypasses the
   precomputed override entirely. Always returns `version: faves_v4_live`.
3. **Local RDKit call:** independent FilterCatalog(ALL) scan via the same
   categorization logic the live-path module uses.

Three diffs are reported per compound:
  - **In-service cross-path** (default vs forced-live): only meaningful when
    the default call returned `faves_v4_precomputed`. This is Test 3.4 proper.
  - **External reproducibility** (forced-live API vs local RDKit): confirms
    the deployed service matches an independent installation. Always meaningful.
  - **End-to-end consistency** (default vs local): the bottom-line check —
    does the response a real consumer receives match an independent reproducer?

Usage
-----
    pip install rdkit requests
    python val001_cross_path_check.py \\
        --api-url https://api.novoquantnexus.com/proxy/faves-compliance \\
        --output val001_cross_path_check_report.md

Requires RDKit installed locally. Production currently runs RDKit 2023.3.1;
local 2025.09.2 is the typical reproducer environment. The 1,585-pattern
catalog is stable across this version range.
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from dataclasses import dataclass, field
from datetime import datetime, timezone
from typing import Any

import requests

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem  # noqa: F401  (ensures FilterCatalog deps load)
    from rdkit.Chem.rdfiltercatalog import FilterCatalog, FilterCatalogParams
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False


# Twenty well-known drug/natural-product SMILES selected to (a) likely resolve
# in the 122M PubChem-derived index (so we exercise the precomputed override)
# and (b) collectively span PAINS, Brenk, NIH, ZINC, and ChEMBL sub-catalogs.
TEST_COMPOUNDS: list[tuple[str, str]] = [
    ("aspirin",         "CC(=O)Oc1ccccc1C(=O)O"),
    ("curcumin",        "COc1cc(/C=C/C(=O)CC(=O)/C=C/c2ccc(O)c(OC)c2)ccc1O"),
    ("warfarin",        "CC(=O)CC(c1ccccc1)c1c(O)c2ccccc2oc1=O"),
    ("ibuprofen",       "CC(C)Cc1ccc(C(C)C(=O)O)cc1"),
    ("acetaminophen",   "CC(=O)Nc1ccc(O)cc1"),
    ("caffeine",        "Cn1cnc2c1c(=O)n(C)c(=O)n2C"),
    ("atorvastatin",    "CC(C)c1c(C(=O)Nc2ccccc2)c(-c2ccccc2)c(-c2ccc(F)cc2)n1CC[C@@H](O)C[C@@H](O)CC(=O)O"),
    ("sildenafil",      "CCCc1nn(C)c2c1nc([nH]c2=O)c1cc(S(=O)(=O)N2CCN(C)CC2)ccc1OCC"),
    ("metformin",       "CN(C)C(=N)NC(=N)N"),
    ("amoxicillin",     "CC1(C)S[C@@H]2[C@H](NC(=O)[C@H](N)c3ccc(O)cc3)C(=O)N2[C@H]1C(=O)O"),
    ("diazepam",        "CN1C(=O)CN=C(c2ccccc2)c2cc(Cl)ccc21"),
    ("quercetin",       "Oc1cc(O)c2c(c1)oc(-c1ccc(O)c(O)c1)c(O)c2=O"),
    ("resveratrol",     "Oc1ccc(/C=C/c2cc(O)cc(O)c2)cc1"),
    ("genistein",       "Oc1ccc(-c2coc3cc(O)cc(O)c3c2=O)cc1"),
    ("doxorubicin",     "COc1cccc2C(=O)c3c(O)c4C[C@@](O)(C(=O)CO)C[C@H](O[C@H]5C[C@H](N)[C@H](O)[C@H](C)O5)c4c(O)c3C(=O)c12"),
    ("imatinib",        "Cc1ccc(NC(=O)c2ccc(CN3CCN(C)CC3)cc2)cc1Nc1nccc(-c2cccnc2)n1"),
    ("tamoxifen",       "CC/C(=C(/c1ccccc1)c1ccc(OCCN(C)C)cc1)c1ccccc1"),
    ("naproxen",        "COc1ccc2cc(C(C)C(=O)O)ccc2c1"),
    ("osimertinib",     "C=CC(=O)Nc1cc(Nc2nccc(-c3cn(C)c4ccccc34)n2)c(OC)cc1N(C)CCN(C)C"),
    ("lisinopril",      "NCCCC[C@@H](N[C@@H](CCc1ccccc1)C(=O)O)C(=O)N1CCC[C@H]1C(=O)O"),
]


CHEMBL_SETS = {"glaxo", "dundee", "bms", "surechembl", "mlsmr", "inpharmatica", "lint"}


@dataclass
class LocalCounts:
    pains: list[str] = field(default_factory=list)
    brenk: list[str] = field(default_factory=list)
    nih: list[str] = field(default_factory=list)
    zinc: list[str] = field(default_factory=list)
    chembl: dict[str, list[str]] = field(default_factory=dict)
    chembl_total: int = 0

    @property
    def total(self) -> int:
        return len(self.pains) + len(self.brenk) + len(self.nih) + len(self.zinc) + self.chembl_total


def categorize_local(smiles: str) -> LocalCounts | None:
    """Mirror of the live-path categorization step live-path logic."""
    if not RDKIT_AVAILABLE:
        return None
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    params = FilterCatalogParams()
    params.AddCatalog(FilterCatalogParams.FilterCatalogs.ALL)
    catalog = FilterCatalog(params)
    counts = LocalCounts()
    for entry in catalog.GetMatches(mol):
        desc = entry.GetDescription()
        filter_set = ""
        try:
            filter_set = entry.GetProp("FilterSet")
        except KeyError:
            pass
        except Exception:
            pass
        fs_lower = filter_set.lower() if filter_set else ""
        if fs_lower.startswith("pains"):
            counts.pains.append(desc)
        elif fs_lower == "brenk":
            counts.brenk.append(desc)
        elif fs_lower == "nih":
            counts.nih.append(desc)
        elif fs_lower == "zinc":
            counts.zinc.append(desc)
        elif fs_lower.startswith("chembl"):
            sub = fs_lower.split("_", 1)[-1] if "_" in fs_lower else "other"
            counts.chembl.setdefault(sub, []).append(desc)
            counts.chembl_total += 1
        elif not filter_set:
            counts.chembl.setdefault("other", []).append(desc)
            counts.chembl_total += 1
    return counts


def call_classify(api_url: str, smiles: str, force_live: bool = False, timeout: float = 30.0) -> dict[str, Any] | None:
    try:
        params = {"force_live": "true"} if force_live else {}
        r = requests.post(
            f"{api_url.rstrip('/')}/api/classify",
            json={"smiles": smiles},
            params=params,
            timeout=timeout,
        )
        r.raise_for_status()
        body = r.json()
        # The quanta-mcp proxy wraps faves-compliance responses in
        # {service, status, data: {...}}. Direct calls return the flat
        # ClassifyResponse. Unwrap if we see the proxy envelope.
        if isinstance(body, dict) and "data" in body and isinstance(body["data"], dict) and "faves_v4" in body["data"]:
            return body["data"]
        return body
    except Exception as e:
        return {"_error": str(e)}


def extract_counts(api_resp: dict[str, Any]) -> dict[str, Any] | None:
    """Pull the comparison fields from a /api/classify response."""
    if not isinstance(api_resp, dict) or "_error" in api_resp:
        return None
    summary = api_resp.get("structural_alert_summary") or {}
    chembl = summary.get("chembl") or {}
    return {
        "version": (api_resp.get("faves_v4") or {}).get("version"),
        "per_catalog": {
            "pains":  (summary.get("pains")  or {}).get("count", 0),
            "brenk":  (summary.get("brenk")  or {}).get("count", 0),
            "nih":    (summary.get("nih")    or {}).get("count", 0),
            "zinc":   (summary.get("zinc")   or {}).get("count", 0),
            "chembl": chembl.get("count", 0),
        },
        "chembl_keys": sorted((chembl.get("catalogs") or {}).keys()),
        "total":       summary.get("total_alert_count"),
    }


def _diff_pair(label_a: str, a: dict[str, Any] | None, label_b: str, b: dict[str, Any] | None) -> list[str]:
    """Return human-readable divergence strings between two extract_counts dicts."""
    if a is None or b is None:
        return [f"{label_a if a is None else label_b} unavailable — cannot compare"]
    out: list[str] = []
    for cat in ("pains", "brenk", "nih", "zinc", "chembl"):
        if a["per_catalog"][cat] != b["per_catalog"][cat]:
            out.append(f"{cat}: {label_a}={a['per_catalog'][cat]} {label_b}={b['per_catalog'][cat]}")
    if a["total"] != b["total"]:
        out.append(f"total_alert_count: {label_a}={a['total']} {label_b}={b['total']}")
    if set(a["chembl_keys"]) != set(b["chembl_keys"]):
        only_a = sorted(set(a["chembl_keys"]) - set(b["chembl_keys"]))
        only_b = sorted(set(b["chembl_keys"]) - set(a["chembl_keys"]))
        out.append(f"chembl keys differ — {label_a}_only={only_a} {label_b}_only={only_b}")
    return out


def diff_compound(
    name: str,
    smiles: str,
    api_default: dict[str, Any],
    api_forced: dict[str, Any],
    local: LocalCounts | None,
) -> dict[str, Any]:
    """Return a per-compound record with three diffs."""
    default_x = extract_counts(api_default)
    forced_x  = extract_counts(api_forced)
    local_x = None
    if local is not None:
        local_x = {
            "version": "local_rdkit",
            "per_catalog": {
                "pains":  len(local.pains),
                "brenk":  len(local.brenk),
                "nih":    len(local.nih),
                "zinc":   len(local.zinc),
                "chembl": local.chembl_total,
            },
            "chembl_keys": sorted(local.chembl.keys()),
            "total":       local.total,
        }

    record: dict[str, Any] = {
        "name": name,
        "smiles": smiles,
        "default_path":   default_x["version"] if default_x else None,
        "forced_path":    forced_x["version"]  if forced_x  else None,
        "default":        default_x,
        "forced_live":    forced_x,
        "local":          local_x,
        "default_error":  api_default.get("_error") if isinstance(api_default, dict) else None,
        "forced_error":   api_forced.get("_error")  if isinstance(api_forced,  dict) else None,
    }

    # In-service cross-path: only meaningful if default returned precomputed
    in_service_meaningful = (
        default_x is not None and forced_x is not None
        and default_x["version"] == "faves_v4_precomputed"
    )
    record["in_service_meaningful"] = in_service_meaningful
    record["in_service_divergences"] = (
        _diff_pair("precomputed", default_x, "live", forced_x) if in_service_meaningful else []
    )

    record["external_divergences"]    = _diff_pair("api_live", forced_x,  "local", local_x)
    record["end_to_end_divergences"]  = _diff_pair("api_default", default_x, "local", local_x)

    return record


def render_report(records: list[dict[str, Any]], api_url: str) -> str:
    now = datetime.now(timezone.utc).isoformat()
    total = len(records)

    failed_default = sum(1 for r in records if r["default_path"] is None)
    failed_forced  = sum(1 for r in records if r["forced_path"]  is None)
    served_precomputed = sum(1 for r in records if r["default_path"] == "faves_v4_precomputed")
    served_live        = sum(1 for r in records if r["default_path"] == "faves_v4_live")

    in_service_runs       = sum(1 for r in records if r["in_service_meaningful"])
    in_service_divergent  = sum(1 for r in records if r["in_service_divergences"])
    external_divergent    = sum(1 for r in records if r["external_divergences"])
    end_to_end_divergent  = sum(1 for r in records if r["end_to_end_divergences"])

    lines: list[str] = []
    lines.append("# VAL-001 Test 3.4 — Cross-Path Consistency Report")
    lines.append("")
    lines.append(f"**Run:** {now}")
    lines.append(f"**API:** `{api_url}`")
    lines.append(f"**Compounds tested:** {total}")
    lines.append("")
    lines.append("## Summary")
    lines.append("")
    lines.append(f"- Default-call failures: **{failed_default}** / forced-live failures: **{failed_forced}**")
    lines.append(f"- Default served by precomputed override: **{served_precomputed}**")
    lines.append(f"- Default served by live computation only: **{served_live}**")
    lines.append("")
    lines.append("**Three diffs:**")
    lines.append("")
    lines.append(f"- **In-service cross-path** (precomputed vs forced-live, the strict Test 3.4): "
                 f"{in_service_runs - in_service_divergent}/{in_service_runs} aligned. "
                 f"({total - in_service_runs} compounds were not eligible — default call returned live.)")
    lines.append(f"- **External reproducibility** (forced-live API vs local RDKit): "
                 f"{total - external_divergent}/{total} aligned.")
    lines.append(f"- **End-to-end consistency** (default API response vs local RDKit): "
                 f"{total - end_to_end_divergent}/{total} aligned.")
    lines.append("")

    lines.append("## Per-compound results")
    lines.append("")
    lines.append("| Compound | Default path | In-service Δ | External Δ | End-to-end Δ |")
    lines.append("|---|---|---|---|---|")
    for r in records:
        in_s = "n/a" if not r["in_service_meaningful"] else (
            "—" if not r["in_service_divergences"] else str(len(r["in_service_divergences"]))
        )
        ext = "—" if not r["external_divergences"] else str(len(r["external_divergences"]))
        e2e = "—" if not r["end_to_end_divergences"] else str(len(r["end_to_end_divergences"]))
        lines.append(f"| {r['name']} | {r['default_path'] or 'error'} | {in_s} | {ext} | {e2e} |")
    lines.append("")

    divergent_records = [
        r for r in records
        if r["in_service_divergences"] or r["external_divergences"] or r["end_to_end_divergences"]
    ]
    if divergent_records:
        lines.append("## Divergence detail")
        lines.append("")
        for r in divergent_records:
            lines.append(f"### {r['name']}")
            lines.append("")
            lines.append(f"- **SMILES:** `{r['smiles']}`")
            lines.append(f"- **Default path:** `{r['default_path']}` / forced-live path: `{r['forced_path']}`")
            if r["default"]:
                lines.append(f"- **Default API counts:** `{r['default']['per_catalog']}`, "
                             f"chembl_keys={r['default']['chembl_keys']}, total={r['default']['total']}")
            if r["forced_live"]:
                lines.append(f"- **Forced-live API counts:** `{r['forced_live']['per_catalog']}`, "
                             f"chembl_keys={r['forced_live']['chembl_keys']}, total={r['forced_live']['total']}")
            if r["local"]:
                lines.append(f"- **Local RDKit counts:** `{r['local']['per_catalog']}`, "
                             f"chembl_keys={r['local']['chembl_keys']}, total={r['local']['total']}")
            for label, divs in (
                ("In-service (precomputed vs forced-live)", r["in_service_divergences"]),
                ("External (forced-live vs local)",        r["external_divergences"]),
                ("End-to-end (default vs local)",          r["end_to_end_divergences"]),
            ):
                if divs:
                    lines.append(f"- **{label}:**")
                    for d in divs:
                        lines.append(f"  - {d}")
            lines.append("")

    lines.append("## Interpretation")
    lines.append("")
    if in_service_runs == 0:
        lines.append("**Caveat:** none of the test compounds resolved to the precomputed cache in the "
                     "default-call path. The strict in-service cross-path test (Test 3.4 proper) was "
                     "therefore not exercised. The external-reproducibility and end-to-end results are "
                     "still meaningful and reported above. To exercise in-service cross-path, choose "
                     "compounds known to have entries in the SMILES → CID Redis index, or warm the "
                     "cache by querying CIDs directly first.")
    elif in_service_divergent == 0:
        lines.append("**In-service cross-path:** all eligible compounds produced identical results "
                     "across precomputed and forced-live paths. The the precomputed structural-alerts cache and "
                     "the live RDKit FilterCatalog are arithmetically and structurally consistent.")
    else:
        lines.append("**In-service cross-path divergences detected.** Likely causes:")
        lines.append("- ChEMBL keys differ → precomputed loader writes keys in a different format than "
                     "the live-path categorization step produces. Fix in the loader.")
        lines.append("- Per-catalog counts diverge → precomputed data is stale relative to the current "
                     "SMARTS set, or RDKit version drift between the batch pipeline and the live endpoint.")
        lines.append("- `total_alert_count` arithmetic mismatch → bug in the override block "
                     "(the precomputed-override block) or in the loader's stored `chembl_total_count`.")

    if external_divergent == 0:
        lines.append("")
        lines.append("**External reproducibility:** the deployed service produces output identical to a "
                     "fresh local RDKit FilterCatalog scan. Pattern set and categorization logic are "
                     "stable across the deployed and local RDKit versions.")
    lines.append("")
    return "\n".join(lines)


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--api-url", default="https://api.novoquantnexus.com/proxy/faves-compliance",
                   help="Base URL of the deployed faves-compliance service")
    p.add_argument("--output", default="val001_cross_path_check_report.md",
                   help="Output markdown report path")
    p.add_argument("--delay", type=float, default=0.5,
                   help="Seconds between API calls (rate-limit cushion)")
    args = p.parse_args()

    if not RDKIT_AVAILABLE:
        print("ERROR: RDKit is not installed locally — `pip install rdkit` and re-run.", file=sys.stderr)
        return 2

    print(f"Hitting {args.api_url} with {len(TEST_COMPOUNDS)} compounds (2 calls each)...", file=sys.stderr)
    records: list[dict[str, Any]] = []
    for name, smiles in TEST_COMPOUNDS:
        print(f"  {name}...", file=sys.stderr)
        api_default = call_classify(args.api_url, smiles, force_live=False)
        time.sleep(args.delay)
        api_forced  = call_classify(args.api_url, smiles, force_live=True)
        local = categorize_local(smiles)
        records.append(diff_compound(name, smiles, api_default or {}, api_forced or {}, local))
        time.sleep(args.delay)

    report = render_report(records, args.api_url)
    with open(args.output, "w") as f:
        f.write(report)

    raw_path = args.output.replace(".md", ".json")
    with open(raw_path, "w") as f:
        json.dump(records, f, indent=2, default=str)

    in_service_div = sum(1 for r in records if r["in_service_divergences"])
    external_div   = sum(1 for r in records if r["external_divergences"])
    end_to_end_div = sum(1 for r in records if r["end_to_end_divergences"])
    print(f"\nReport written: {args.output}", file=sys.stderr)
    print(f"Raw records:    {raw_path}", file=sys.stderr)
    print(f"In-service divergent:  {in_service_div}/{len(records)}", file=sys.stderr)
    print(f"External divergent:    {external_div}/{len(records)}", file=sys.stderr)
    print(f"End-to-end divergent:  {end_to_end_div}/{len(records)}", file=sys.stderr)
    return 0 if (in_service_div + external_div + end_to_end_div) == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
