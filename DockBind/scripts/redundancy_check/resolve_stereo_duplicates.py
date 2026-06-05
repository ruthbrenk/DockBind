#!/usr/bin/env python3
"""
Stereo and salt resolution of duplicated CHEMBL / BINDINGDB entries

Needs the duplication log produced by find_sdf_duplicates.py

- Retrieves ORIGINAL smiles from their parent datasets
    * CHEMBL  -> from a ChEMBL chemreps TXT file (chembl_id, canonical_smiles)
    * BINDINGDB -> from a BindingDB TSV (BindingDB MonomerID, Ligand SMILES)
- Polishes SMILES by removing trailing non-SMILES annotations
- Counts explicitly defined stereocenters
- Removes compounds with fewer stereocenters (unambiguous losers).
- Detects ambiguous winners (same stereo count), and:
    * resolves them when a salt vs desalted pair is present
    * records these in a separate salt CSV
    * logs unresolved ambiguities.

Outputs:
- A detailed log with statistics and unresolved cases.
- A CSV of unambiguous losers (safe removals) with columns:
    template,ID,smiles
- A CSV of salt–desalt resolutions with columns:
    template,desalted,salt

"""

import argparse
import csv
import re
from rdkit import Chem


# =========================
# Helpers
# =========================

def polish_smiles(smiles: str) -> str:
    """
    Remove everything after the first whitespace.
    Preserves salts and disconnected components.
    """
    return smiles.split()[0]


def is_salt(smiles: str) -> bool:
    """
    Detect salts / non-covalent entities by presence of '.'.
    """
    return "." in smiles


def count_defined_stereocenters(smiles: str):
    """
    Return the number of explicitly defined stereocenters in a SMILES.
    Returns None if SMILES cannot be parsed.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    centers = Chem.FindMolChiralCenters(
        mol, includeUnassigned=True, useLegacyImplementation=False
    )
    return sum(1 for _, tag in centers if tag != "?")


# =========================
# Log parsing
# =========================

TEMPLATE_RE = re.compile(r"^Template:\s+(.+)")
GROUP_RE = re.compile(r"^\s*same compound group")
ID_RE = re.compile(r"^\s+(CHEMBL\d+|BINDINGDB\d+)")


def parse_duplication_log(path):
    templates = {}
    current_template = None
    current_group = None

    with open(path, "r") as f:
        for line in f:
            line = line.rstrip()

            m = TEMPLATE_RE.match(line)
            if m:
                current_template = m.group(1)
                templates[current_template] = []
                current_group = None
                continue

            if GROUP_RE.match(line):
                current_group = []
                templates[current_template].append(current_group)
                continue

            m = ID_RE.match(line)
            if m and current_group is not None:
                current_group.append(m.group(1))

    return templates


# =========================
# Data loading
# =========================

def load_bindingdb_tsv(tsv):
    mapping = {}
    with open(tsv) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            mid = row.get("BindingDB MonomerID", "").strip()
            smi = row.get("Ligand SMILES", "").strip() or None
            if mid:
                mapping[mid] = smi
    return mapping


def load_chembl_chemreps(txt):
    mapping = {}
    with open(txt) as f:
        header = f.readline().split()
        idx = {k: i for i, k in enumerate(header)}
        for line in f:
            parts = line.split()
            if len(parts) <= max(idx["chembl_id"], idx["canonical_smiles"]):
                continue
            cid = parts[idx["chembl_id"]]
            smi = parts[idx["canonical_smiles"]] or None
            mapping[cid] = smi
    return mapping


# =========================
# Main logic
# =========================

def main(log_path, bindingdb_tsv, chemreps_txt, out_log, losers_csv, salts_csv):
    templates = parse_duplication_log(log_path)
    bdb = load_bindingdb_tsv(bindingdb_tsv)
    chembl = load_chembl_chemreps(chemreps_txt)

    # losers: list of (template, ID, smiles)
    losers = []
    # salt_pairs: list of (template, desalted_ID, salt_ID)
    salt_pairs = []

    ambiguous_total = 0
    ambiguous_with_salt = 0
    ambiguous_resolved = 0
    ambiguous_unresolved = []

    with open(out_log, "w") as log:

        for template, groups in templates.items():
            for group in groups:
                for db, ids in [
                    ("CHEMBL", [x for x in group if x.startswith("CHEMBL")]),
                    ("BINDINGDB", [x for x in group if x.startswith("BINDINGDB")]),
                ]:
                    if len(ids) < 2:
                        continue

                    log.write(f"Template: {template}\n")
                    log.write(f"Database: {db}\n")

                    scored = []
                    missing = []

                    for cid in ids:
                        raw = (
                            chembl.get(cid)
                            if db == "CHEMBL"
                            else bdb.get(cid.replace("BINDINGDB", ""))
                        )

                        if raw is None:
                            log.write(f"  {cid}: MISSING\n")
                            missing.append(cid)
                            continue

                        polished = polish_smiles(raw)
                        stereo = count_defined_stereocenters(polished)

                        if stereo is None:
                            log.write(f"  {cid}: INVALID_SMILES\n")
                            missing.append(cid)
                            continue

                        log.write(
                            f"  {cid}: stereo_centers={stereo} "
                            f"SMILES={polished}\n"
                        )
                        scored.append((cid, polished, stereo))

                    if not scored:
                        log.write("  No valid SMILES; skipping group\n\n")
                        continue

                    max_stereo = max(x[2] for x in scored)
                    winners = [x for x in scored if x[2] == max_stereo]
                    losers_local = [x for x in scored if x[2] < max_stereo]

                    # unambiguous losers (template-aware)
                    for cid, smi, _ in losers_local:
                        losers.append((template, cid, smi))

                    if len(winners) == 1:
                        log.write("\n")
                        continue

                    # ambiguous winners
                    ambiguous_total += 1
                    log.write("  WARNING: multiple winners\n")

                    salted = [w for w in winners if is_salt(w[1])]
                    desalted = [w for w in winners if not is_salt(w[1])]

                    if len(salted) == 1 and len(desalted) == 1:
                        ambiguous_with_salt += 1
                        ambiguous_resolved += 1
                        # record template + desalted + salt
                        salt_pairs.append((template, desalted[0][0], salted[0][0]))
                        log.write(
                            f"  Salt resolved: {desalted[0][0]} (desalted) "
                            f"{salted[0][0]} (salt)\n\n"
                        )
                    else:
                        ambiguous_unresolved.append(
                            (template, db, [w[0] for w in winners])
                        )
                        log.write("  Ambiguity not resolvable by salt logic\n\n")

        # summary
        log.write("SUMMARY\n=======\n")
        log.write(f"Ambiguous winner groups: {ambiguous_total}\n")
        log.write(f"With salts detected: {ambiguous_with_salt}\n")
        log.write(f"Resolved by salt logic: {ambiguous_resolved}\n")
        log.write(
            f"Unresolved ambiguities: {len(ambiguous_unresolved)}\n"
        )

        if ambiguous_unresolved:
            log.write("\nUNRESOLVED CASES\n")
            for t, db, ids in ambiguous_unresolved:
                log.write(f"{t} {db}: {', '.join(ids)}\n")

    # write losers csv: template,ID,smiles
    with open(losers_csv, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["template", "ID", "smiles"])
        for template, cid, smi in losers:
            writer.writerow([template, cid, smi])

    # write salt resolutions csv: template,desalted,salt
    with open(salts_csv, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["template", "desalted", "salt"])
        for template, desalted_id, salt_id in salt_pairs:
            writer.writerow([template, desalted_id, salt_id])


# =========================
# Entry point
# =========================

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Resolve stereo-duplicate CHEMBL/BINDINGDB entries and list losers and salt resolutions."
    )
    parser.add_argument("log", help="Path to duplication log file.")
    parser.add_argument(
        "bindingdb_tsv",
        help="Path to BindingDB TSV with 'BindingDB MonomerID' and 'Ligand SMILES'.",
    )
    parser.add_argument(
        "chemreps_txt",
        help="Path to ChEMBL chemreps TXT with 'chembl_id' and 'canonical_smiles'.",
    )
    parser.add_argument(
        "--out-log",
        default="stereo_resolution.log",
        help="Path to detailed resolution log (default: stereo_resolution.log).",
    )
    parser.add_argument(
        "--losers-csv",
        default="losers_to_remove.csv",
        help="Path to CSV with unambiguous losers (default: losers_to_remove.csv).",
    )
    parser.add_argument(
        "--salts-csv",
        default="salt_resolutions.csv",
        help="Path to CSV with salt/desalted pairs (default: salt_resolutions.csv).",
    )
    args = parser.parse_args()

    main(
        args.log,
        args.bindingdb_tsv,
        args.chemreps_txt,
        args.out_log,
        args.losers_csv,
        args.salts_csv,
    )



