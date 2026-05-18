# FAVES V4 — End-to-End Latency Benchmark

**Run:** 2026-05-10T17:06:58.980980+00:00
**API:** `https://api.novoquantnexus.com/proxy/faves-compliance`
**Calls per compound:** 10  
**Compounds:** 10  
**Total calls (steady-state):** 90 (first call per compound excluded as warm-up)

## Overall steady-state latency

- min:  **404 ms**
- p50:  **3126 ms**
- p90:  **3280 ms**
- p95:  **3376 ms**
- p99:  **3460 ms**
- max:  **3479 ms**
- mean: **2363 ms**

## Per-compound latencies (ms)

| Compound | Path | Warm-up | Steady min | p50 | p95 | max | n |
|---|---|---|---|---|---|---|---|
| aspirin | faves_v4_live | 3285 | 3104 | 3175 | 3321 | 3347 | 9 |
| ibuprofen | faves_v4_live | 3070 | 2999 | 3136 | 3301 | 3388 | 9 |
| caffeine | faves_v4_live | 3100 | 3071 | 3183 | 3276 | 3280 | 9 |
| acetaminophen | faves_v4_live | 3195 | 3086 | 3161 | 3212 | 3218 | 9 |
| metformin | faves_v4_live | 445 | 404 | 449 | 503 | 511 | 9 |
| warfarin | faves_v4_live | 3200 | 2999 | 3292 | 3470 | 3479 | 9 |
| imatinib | faves_v4_live | 505 | 415 | 445 | 519 | 536 | 9 |
| atorvastatin | faves_v4_live | 436 | 426 | 472 | 525 | 536 | 9 |
| doxorubicin | faves_v4_live | 3094 | 3063 | 3171 | 3315 | 3379 | 9 |
| osimertinib | faves_v4_live | 3273 | 3048 | 3175 | 3250 | 3274 | 9 |
