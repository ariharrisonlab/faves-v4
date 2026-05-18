"""End-to-end latency benchmark for /api/classify.

Hits the deployed faves-compliance service N times across a diverse compound
set. Reports per-compound and overall percentile latencies. Output is a
markdown table suitable for direct inclusion in the V4 paper §5.6.

The first call to a fresh container can incur warm-up overhead. The script
discards the first call per compound from percentile statistics; that
result is reported separately as "warm-up latency."

Usage
-----
    python latency_benchmark.py \\
        --api-url https://api.novoquantnexus.com/proxy/faves-compliance \\
        --calls-per-compound 10 \\
        --output latency_report.md

Default is 10 SMILES × 10 calls = 100 calls.
"""

from __future__ import annotations

import argparse
import json
import statistics
import sys
import time
from datetime import datetime, timezone
from typing import Any

import requests


# Diverse SMILES covering small → large, with and without alerts. Picked from
# the V4 cross-path test set so percentile numbers are comparable across
# benchmarks.
BENCHMARK_COMPOUNDS: list[tuple[str, str]] = [
    ("aspirin",       "CC(=O)Oc1ccccc1C(=O)O"),
    ("ibuprofen",     "CC(C)Cc1ccc(C(C)C(=O)O)cc1"),
    ("caffeine",      "Cn1cnc2c1c(=O)n(C)c(=O)n2C"),
    ("acetaminophen", "CC(=O)Nc1ccc(O)cc1"),
    ("metformin",     "CN(C)C(=N)NC(=N)N"),
    ("warfarin",      "CC(=O)CC(c1ccccc1)c1c(O)c2ccccc2oc1=O"),
    ("imatinib",      "Cc1ccc(NC(=O)c2ccc(CN3CCN(C)CC3)cc2)cc1Nc1nccc(-c2cccnc2)n1"),
    ("atorvastatin",  "CC(C)c1c(C(=O)Nc2ccccc2)c(-c2ccccc2)c(-c2ccc(F)cc2)n1CC[C@@H](O)C[C@@H](O)CC(=O)O"),
    ("doxorubicin",   "COc1cccc2C(=O)c3c(O)c4C[C@@](O)(C(=O)CO)C[C@H](O[C@H]5C[C@H](N)[C@H](O)[C@H](C)O5)c4c(O)c3C(=O)c12"),
    ("osimertinib",   "C=CC(=O)Nc1cc(Nc2nccc(-c3cn(C)c4ccccc34)n2)c(OC)cc1N(C)CCN(C)C"),
]


def call_classify(api_url: str, smiles: str, timeout: float = 60.0) -> tuple[float | None, str | None]:
    """Return (elapsed_seconds, version_tag) or (None, error_message)."""
    t0 = time.perf_counter()
    try:
        r = requests.post(
            f"{api_url.rstrip('/')}/api/classify",
            json={"smiles": smiles},
            timeout=timeout,
        )
        elapsed = time.perf_counter() - t0
        r.raise_for_status()
        body = r.json()
        if isinstance(body, dict) and "data" in body and isinstance(body["data"], dict) and "faves_v4" in body["data"]:
            body = body["data"]
        version = (body.get("faves_v4") or {}).get("version") if isinstance(body, dict) else None
        return elapsed, version
    except Exception as e:
        return None, str(e)


def percentile(values: list[float], p: float) -> float:
    if not values:
        return float("nan")
    s = sorted(values)
    k = (len(s) - 1) * (p / 100.0)
    lo = int(k)
    hi = min(lo + 1, len(s) - 1)
    if lo == hi:
        return s[lo]
    return s[lo] + (s[hi] - s[lo]) * (k - lo)


def render_report(per_compound: list[dict[str, Any]], all_steady_ms: list[float],
                  api_url: str, calls_per_compound: int) -> str:
    now = datetime.now(timezone.utc).isoformat()

    lines: list[str] = []
    lines.append("# FAVES V4 — End-to-End Latency Benchmark")
    lines.append("")
    lines.append(f"**Run:** {now}")
    lines.append(f"**API:** `{api_url}`")
    lines.append(f"**Calls per compound:** {calls_per_compound}  ")
    lines.append(f"**Compounds:** {len(per_compound)}  ")
    lines.append(f"**Total calls (steady-state):** {len(all_steady_ms)} (first call per compound excluded as warm-up)")
    lines.append("")

    if all_steady_ms:
        lines.append("## Overall steady-state latency")
        lines.append("")
        lines.append(f"- min:  **{min(all_steady_ms):.0f} ms**")
        lines.append(f"- p50:  **{percentile(all_steady_ms, 50):.0f} ms**")
        lines.append(f"- p90:  **{percentile(all_steady_ms, 90):.0f} ms**")
        lines.append(f"- p95:  **{percentile(all_steady_ms, 95):.0f} ms**")
        lines.append(f"- p99:  **{percentile(all_steady_ms, 99):.0f} ms**")
        lines.append(f"- max:  **{max(all_steady_ms):.0f} ms**")
        lines.append(f"- mean: **{statistics.fmean(all_steady_ms):.0f} ms**")
        lines.append("")

    lines.append("## Per-compound latencies (ms)")
    lines.append("")
    lines.append("| Compound | Path | Warm-up | Steady min | p50 | p95 | max | n |")
    lines.append("|---|---|---|---|---|---|---|---|")
    for c in per_compound:
        lines.append(
            f"| {c['name']} | {c['path']} | "
            f"{c['warmup_ms']:.0f} | "
            f"{min(c['steady_ms']):.0f} | "
            f"{percentile(c['steady_ms'], 50):.0f} | "
            f"{percentile(c['steady_ms'], 95):.0f} | "
            f"{max(c['steady_ms']):.0f} | "
            f"{len(c['steady_ms'])} |"
        )
    lines.append("")
    return "\n".join(lines)


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--api-url", default="https://api.novoquantnexus.com/proxy/faves-compliance",
                   help="Base URL of the deployed faves-compliance service")
    p.add_argument("--calls-per-compound", type=int, default=10,
                   help="Calls per compound (first call per compound is warm-up; remaining are steady-state)")
    p.add_argument("--delay", type=float, default=0.2, help="Seconds between calls")
    p.add_argument("--output", default="latency_report.md", help="Output markdown report path")
    args = p.parse_args()

    print(f"Latency benchmark: {len(BENCHMARK_COMPOUNDS)} compounds × {args.calls_per_compound} calls",
          file=sys.stderr)

    per_compound: list[dict[str, Any]] = []
    all_steady_ms: list[float] = []

    for name, smiles in BENCHMARK_COMPOUNDS:
        print(f"  {name}", end=" ", flush=True, file=sys.stderr)
        all_calls: list[float | None] = []
        path = None
        for i in range(args.calls_per_compound):
            elapsed, version_or_err = call_classify(args.api_url, smiles)
            all_calls.append(elapsed)
            if elapsed is not None and path is None:
                path = version_or_err  # first non-error response sets path
            if elapsed is None:
                print(f"\n    call {i+1} FAILED: {version_or_err}", file=sys.stderr)
            else:
                print(f"{elapsed*1000:.0f}ms", end=" ", flush=True, file=sys.stderr)
            time.sleep(args.delay)

        successful = [t for t in all_calls if t is not None]
        if not successful:
            per_compound.append({"name": name, "path": "ERROR", "warmup_ms": 0,
                                 "steady_ms": [0]})
            print(file=sys.stderr)
            continue

        warmup_ms = (successful[0] or 0) * 1000
        steady_ms = [(t * 1000) for t in successful[1:]]
        if not steady_ms:
            steady_ms = [warmup_ms]  # only one call, can't separate

        per_compound.append({
            "name":      name,
            "path":      path or "unknown",
            "warmup_ms": warmup_ms,
            "steady_ms": steady_ms,
        })
        all_steady_ms.extend(steady_ms)
        print(file=sys.stderr)

    report = render_report(per_compound, all_steady_ms, args.api_url, args.calls_per_compound)
    with open(args.output, "w") as f:
        f.write(report)

    raw = {
        "run":                   datetime.now(timezone.utc).isoformat(),
        "api_url":               args.api_url,
        "calls_per_compound":    args.calls_per_compound,
        "per_compound":          per_compound,
        "overall_steady_ms":     all_steady_ms,
    }
    with open(args.output.replace(".md", ".json"), "w") as f:
        json.dump(raw, f, indent=2)

    print(f"\nReport: {args.output}", file=sys.stderr)
    if all_steady_ms:
        print(f"Steady-state p50={percentile(all_steady_ms, 50):.0f}ms "
              f"p95={percentile(all_steady_ms, 95):.0f}ms", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
