"""Reference compound sets for FAVES V4 validation (NMCP-FAVES-VAL-001).

Provides:
- load_v3_ground_truth() — 102 compounds (47 DEA + 46 FDA + 9 negative controls)
  loaded from v3_ground_truth.csv, the same dataset published with the v3 paper
  (https://github.com/realariharrison/faves-v3-benchmark).
- TEST_3_2_NOVEL — 50 novel SMILES (not expected to resolve in PubChem) for
  field-completeness stress-testing of the live computation path.
- TEST_3_2_EDGE_CASES — 48 edge-case SMILES spanning charged species, small
  fragments, large/complex molecules, peptides, heterocycle-heavy structures,
  and unusual functional groups.

Test 3.2 set = 102 (v3) + 50 (novel) + 48 (edge) = 200 compounds, matching
NMCP-FAVES-VAL-001 §3.2 specification.
"""

import csv
import os
from typing import List, Tuple


def load_v3_ground_truth() -> List[dict]:
    """Load the v3 benchmark CSV. Returns list of dicts with name/smiles/cid/expected_controlled/category/schedule."""
    csv_path = os.path.join(os.path.dirname(__file__), "..", "data", "v3_ground_truth.csv")
    with open(csv_path) as f:
        rows = list(csv.DictReader(f))
    for r in rows:
        r["expected_controlled"] = r["expected_controlled"].strip().lower() == "true"
    return rows


# --------------------------------------------------------------------------
# Test 3.2: 50 NOVEL SMILES for stress-testing the live computation path
# --------------------------------------------------------------------------
# Constructed by combining drug-like fragments in arrangements unlikely to
# appear in PubChem. Each tested locally via RDKit MolFromSmiles to confirm
# they parse. Spread across alert catalogs (some PAINS-prone, some clean).

TEST_3_2_NOVEL: List[Tuple[str, str]] = [
    # Heteroaromatic amides (5)
    ("novel_pyrazine_amide_1",   "c1cnc(NC(=O)c2ccccn2)cc1"),
    ("novel_thiazole_amide_1",   "c1csc(NC(=O)c2cnccc2F)n1"),
    ("novel_oxazole_amide_1",    "c1ocnc1NC(=O)c1ccc(N)cc1"),
    ("novel_pyrimidine_amide_1", "c1cnc(NC(=O)c2cccc(O)c2)nc1"),
    ("novel_imidazole_amide_1",  "c1ncn(c1)C(=O)Nc1cccc(F)c1"),
    # Sulfonamides (5)
    ("novel_sulfonamide_1", "CC(=O)Nc1ccc(S(=O)(=O)NC2CC2)cc1F"),
    ("novel_sulfonamide_2", "Cc1ccc(S(=O)(=O)Nc2cccc(C(F)(F)F)c2)cc1"),
    ("novel_sulfonamide_3", "O=S(=O)(N1CCOCC1)c1ccc(C(=O)O)cc1"),
    ("novel_sulfonamide_4", "CN(C)S(=O)(=O)c1ccc(NC(=O)c2cccnc2)cc1"),
    ("novel_sulfonamide_5", "Nc1ccc(S(=O)(=O)Nc2nc3ccccc3s2)cc1"),
    # Substituted indoles (5)
    ("novel_indole_1", "c1ccc2c(c1)cc(CC(=O)N(C)C)n2C"),
    ("novel_indole_2", "Fc1ccc2[nH]cc(C(=O)Nc3ccc(F)cc3)c2c1"),
    ("novel_indole_3", "OC(=O)Cc1cn(CC(=O)N2CCOCC2)c2ccccc12"),
    ("novel_indole_4", "Clc1ccc2c(c1)cc(C(=O)NCc1ccncc1)n2C"),
    ("novel_indole_5", "Nc1ccc2[nH]cc(CCN3CCOCC3)c2c1"),
    # Quinolines / isoquinolines (5)
    ("novel_quinoline_1",   "c1ccc2nc(NC(=O)c3cccs3)ccc2c1"),
    ("novel_quinoline_2",   "Cc1ccc2nc(NC(=O)c3ccc(F)cc3)cc(O)c2c1"),
    ("novel_isoquinoline_1","c1ccc2cnc(NC(=O)C3CC3)cc2c1"),
    ("novel_quinoline_3",   "Fc1ccc2nc(N3CCN(C)CC3)ccc2c1"),
    ("novel_isoquinoline_2","COc1cc2cnc(N3CCOCC3)cc2cc1OC"),
    # Triazoles + tetrazoles (5)
    ("novel_triazole_1", "c1nnc(Cc2ccc(F)cc2)n1Cc1ccc(C(=O)O)cc1"),
    ("novel_triazole_2", "Cn1nnc(C(=O)Nc2ccc(C)cc2)n1"),
    ("novel_tetrazole_1","c1nnnn1Cc1ccc(C(=O)NC2CC2)cc1"),
    ("novel_triazole_3", "Cc1nn(c(=O)n1Cc1ccccc1F)CC(=O)O"),
    ("novel_tetrazole_2","OC(=O)c1ccc(-c2nnn[nH]2)cc1F"),
    # Piperidines / piperazines / morpholines (5)
    ("novel_piperidine_1","O=C(NCC1CCN(C(=O)c2ccncc2)CC1)c1ccoc1"),
    ("novel_piperazine_1","Cc1ccc(N2CCN(C(=O)Nc3ccc(F)cc3)CC2)cc1"),
    ("novel_morpholine_1","O=C(N1CCOCC1)c1ccc(NC(=O)c2cccs2)cc1"),
    ("novel_piperazine_2","CN1CCN(c2ccc(NC(=O)c3cccnc3F)cc2)CC1"),
    ("novel_piperidine_2","O=C(c1ccncc1)N1CCN(C(=O)Cc2ccccc2Cl)CC1"),
    # Fluorinated drug-like (5)
    ("novel_fluoro_1","FC(F)(F)c1ccc(NC(=O)c2cccnc2N)cc1"),
    ("novel_fluoro_2","Fc1ccc(NC(=O)C2CCN(c3ncccn3)CC2)cc1F"),
    ("novel_fluoro_3","FC(F)(F)Oc1ccc(NC(=O)Cn2ccnc2)cc1"),
    ("novel_fluoro_4","Fc1ccc2c(c1)C(=O)N(CCO)CC2"),
    ("novel_fluoro_5","Fc1cccc(NC(=O)c2cnc3[nH]ccc3c2)c1"),
    # Phenol / catechol containers (likely PAINS) (5)
    ("novel_phenol_1",     "Oc1ccc(C(=O)Nc2ccccc2)cc1O"),
    ("novel_phenol_2",     "Oc1cccc(/C=C/c2ccc(F)cc2)c1"),
    ("novel_catechol_1",   "Oc1ccc(CCN2CCN(c3ccccc3)CC2)cc1O"),
    ("novel_phenol_3",     "Oc1cc(C(=O)Cc2ccccc2)ccc1O"),
    ("novel_phenol_4",     "Oc1ccc(C(=O)NN=Cc2ccc(O)cc2)cc1"),
    # Heterocycle-heavy / fused ring (10)
    ("novel_pyrrolopyrimidine_1","c1cc2[nH]cnc2nc1NC(=O)c1ccccn1"),
    ("novel_furopyridine_1",     "O=C(NC1CC1)c1cc2cccnc2o1"),
    ("novel_thienopyridine_1",   "c1csc2c1cnc(N3CCOCC3)c2"),
    ("novel_pyrazolopyridine_1", "Cn1nccc1NC(=O)c1cccnc1"),
    ("novel_imidazopyridine_1",  "c1ccc2nc(N3CCN(C)CC3)cnc2c1"),
    ("novel_benzothiazole_1",    "c1ccc2sc(NC(=O)c3ccoc3)nc2c1"),
    ("novel_benzimidazole_1",    "c1ccc2[nH]c(NC(=O)c3cccnc3)nc2c1"),
    ("novel_benzofuran_1",       "Fc1ccc2oc(C(=O)NCC3CC3)cc2c1"),
    ("novel_benzoxazole_1",      "c1ccc2oc(NC(=O)C3CCCO3)nc2c1"),
    ("novel_indazole_1",         "c1ccc2[nH]ncc2c1NC(=O)C1CC1"),
]


# --------------------------------------------------------------------------
# Test 3.2: 48 EDGE-CASE compounds
# --------------------------------------------------------------------------
# Stress the live path with structures the catalog SMARTS may not encounter
# in typical drug-like screening. All parse via RDKit. Goal: confirm V4
# fields are still populated (no None), regardless of whether alerts fire.

TEST_3_2_EDGE_CASES: List[Tuple[str, str]] = [
    # Very small fragments (MW < 50) (8)
    ("edge_methane",      "C"),
    ("edge_ethane",       "CC"),
    ("edge_methanol",     "CO"),
    ("edge_ethylene",     "C=C"),
    ("edge_acetylene",    "C#C"),
    ("edge_ammonia",      "N"),
    ("edge_water",        "O"),
    ("edge_carbon_diox",  "O=C=O"),
    # Small drug-like (50-150 Da) (6)
    ("edge_glycine",      "NCC(=O)O"),
    ("edge_alanine",      "C[C@@H](N)C(=O)O"),
    ("edge_urea",         "NC(=O)N"),
    ("edge_acetamide",    "CC(=O)N"),
    ("edge_acetic_acid",  "CC(=O)O"),
    ("edge_dimethyl_eth", "COCOC"),
    # Charged / zwitterionic species (8)
    ("edge_choline",          "C[N+](C)(C)CCO"),
    ("edge_glutamate",        "[O-]C(=O)CCC(N)C(=O)[O-]"),
    ("edge_betaine",          "C[N+](C)(C)CC(=O)[O-]"),
    ("edge_taurine",          "NCCS(=O)(=O)[O-]"),
    ("edge_carnitine",        "C[N+](C)(C)CC(O)CC(=O)[O-]"),
    ("edge_arginine",         "N=C(N)NCCCC(N)C(=O)O"),
    ("edge_protonated_amine", "[NH3+]CCCN"),
    ("edge_phosphate",        "OP(=O)(O)O"),
    # Large / peptide-like (MW > 700) (6)
    ("edge_cyclosporine_frag",  "CC(C)CC1NC(=O)C(C(C)C)N(C)C(=O)C(CC(C)C)N(C)C(=O)C(C)NC(=O)C(C)NC(=O)C(C(C)C)N(C)C1=O"),
    ("edge_glycopeptide_frag",  "CC1OC(OC2C(O)C(N)C(C)OC2O)C(O)C(O)C1O"),
    ("edge_paclitaxel",         "CC1=C2C(C(=O)C3(C(CC4C(C3C(C(C2(C)C)(CC1OC(=O)C(C(c1ccccc1)NC(=O)c1ccccc1)O)O)OC(=O)C)(CO4)OC(=O)C)O)C)O"),
    ("edge_amphotericin_frag",  "OC1CC(O)CC(O)CC(O)CC(O)CC(O)CC(O)CCC(O)C(C)OC(=O)CC1"),
    ("edge_macrolide_frag",     "CCC1OC(=O)C(C)C(OC2OC(C)CC(N(C)C)C2O)C(C)C(=O)CC(C)CC(C)C(=O)C(C)CC1C"),
    ("edge_long_peptide",       "NC(C)C(=O)NC(CC(N)=O)C(=O)NC(CCC(=O)O)C(=O)NC(C(C)C)C(=O)NC(CC(=O)O)C(=O)NC(C)C(=O)O"),
    # Unusual functional groups (8)
    ("edge_nitro_compound",   "[O-][N+](=O)c1ccc(N)cc1"),
    ("edge_diazonium",        "[N+](#N)c1ccccc1"),
    ("edge_isocyanate",       "O=C=Nc1ccccc1"),
    ("edge_thiocyanate",      "N#CSc1ccccc1"),
    ("edge_azide",            "[N-]=[N+]=NCc1ccccc1"),
    ("edge_hydrazone",        "c1ccc(/C=N\\Nc2ccccc2)cc1"),
    ("edge_aldehyde",         "O=Cc1ccc(O)cc1"),
    ("edge_thiol_aromatic",   "Sc1ccc(C(=O)O)cc1"),
    # Heavy halogens / unusual atoms (6)
    ("edge_iodobenzene",      "Ic1ccccc1"),
    ("edge_bromohexane",      "BrCCCCCC"),
    ("edge_trifluoromethyl",  "FC(F)(F)c1ccc(C(F)(F)F)cc1"),
    ("edge_perfluoro_acid",   "OC(=O)C(F)(F)C(F)(F)C(F)(F)F"),
    ("edge_silicon_compound", "C[Si](C)(C)c1ccccc1"),
    ("edge_boron_compound",   "OB(O)c1ccc(N)cc1"),
    # Fused poly-aromatic (6)
    ("edge_anthracene",       "c1ccc2cc3ccccc3cc2c1"),
    ("edge_pyrene",           "c1cc2ccc3cccc4ccc(c1)c2c34"),
    ("edge_phenanthrene",     "c1ccc2c(c1)ccc1ccccc12"),
    ("edge_naphthacene",      "c1ccc2cc3cc4ccccc4cc3cc2c1"),
    ("edge_chrysene",         "c1ccc2c(c1)ccc1c2ccc2ccccc12"),
    ("edge_acenaphthylene",   "c1ccc2cc3CC=Cc3cc2c1"),
]


def build_test_3_2_set() -> List[Tuple[str, str]]:
    """Compose the 200-compound Test 3.2 set per VAL-001 §3.2.

    Composition: 102 known (from v3 ground truth) + 50 novel + 48 edge cases = 200.
    """
    compounds: List[Tuple[str, str]] = []

    # 102 known compounds (47 DEA + 46 FDA + 9 negative controls)
    for row in load_v3_ground_truth():
        compounds.append((row["name"], row["smiles"]))

    compounds.extend(TEST_3_2_NOVEL)
    compounds.extend(TEST_3_2_EDGE_CASES)

    assert len(compounds) == 200, f"expected 200 compounds, got {len(compounds)}"
    return compounds


if __name__ == "__main__":
    # Quick sanity check
    rows = load_v3_ground_truth()
    print(f"v3 ground truth: {len(rows)} rows")
    print(f"  controlled: {sum(1 for r in rows if r['expected_controlled'])}")
    print(f"  non-controlled: {sum(1 for r in rows if not r['expected_controlled'])}")
    print(f"novel SMILES: {len(TEST_3_2_NOVEL)}")
    print(f"edge cases: {len(TEST_3_2_EDGE_CASES)}")
    print(f"total Test 3.2 set: {len(build_test_3_2_set())}")

    # Verify all SMILES parse
    try:
        from rdkit import Chem
        bad = []
        for name, smi in TEST_3_2_NOVEL + TEST_3_2_EDGE_CASES:
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                bad.append((name, smi))
        if bad:
            print(f"\nWARNING: {len(bad)} SMILES failed to parse:")
            for name, smi in bad:
                print(f"  {name}: {smi}")
        else:
            print(f"\nAll {len(TEST_3_2_NOVEL) + len(TEST_3_2_EDGE_CASES)} new SMILES parse correctly.")
    except ImportError:
        print("\nRDKit not available — skipping SMILES parse check.")
