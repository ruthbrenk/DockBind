#!/usr/bin/env python3
"""
DockBind final cleaning and deduplication script.

1. Optionally remove specific (template, ID) pairs from:
   - affinities CSV (BM_Template + Docked)
   - pose folders (complexes/<template>/<ID>)
   using:
     - manual CSV (template, ID)
     - stereo CSV  (template, ID)
2. Optionally resolve salt/desalted pairs using:
   - salts CSV (template, desalted, salt)
   Removes salt entries but rescues any unique affinity types
   from salts into the desalted compound.
3. Recompute intra-template duplicates using Standard InChIKey
   on the cleaned pose tree.
4. For each duplicate group, use mcs_rmsd from the affinities CSV
   to keep only the lowest-RMSD compound, but keeps always CHEMBL affinity values
"""

import argparse
import csv
import os
import shutil
from collections import defaultdict

from rdkit import Chem
from rdkit.Chem import inchi


# =========================
# Utility helpers
# =========================

def polish_smiles(smiles: str) -> str:
    """Remove everything after the first whitespace; keep salts/components."""
    return smiles.split()[0]


def is_salt_smiles(smiles: str) -> bool:
    """Heuristic: salt / multiple non-covalent components if '.' present."""
    return "." in smiles


def count_defined_stereocenters(smiles: str):
    """Count explicitly defined stereocenters (None if SMILES invalid)."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    centers = Chem.FindMolChiralCenters(
        mol, includeUnassigned=True, useLegacyImplementation=False
    )
    return sum(1 for _, tag in centers if tag != "?")


def ensure_dir_deleted(path: str) -> bool:
    """Delete directory if it exists; return True if deleted."""
    if os.path.isdir(path):
        shutil.rmtree(path)
        return True
    return False


def remove_template_if_empty(complexes_root: str, template: str) -> bool:
    """
    If template folder has no subdirectories (no compound folders),
    delete the template folder. Return True if deleted.
    """
    tpath = os.path.join(complexes_root, template)
    if not os.path.isdir(tpath):
        return False
    entries = [
        e for e in os.listdir(tpath)
        if os.path.isdir(os.path.join(tpath, e))
    ]
    if entries:
        return False
    shutil.rmtree(tpath)
    return True


# =========================
# Affinities CSV helpers
# =========================

def load_affinities_csv(path):
    with open(path, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        rows = list(reader)
        header = reader.fieldnames
    return rows, header


def write_affinities_csv(rows, header, path):
    with open(path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=header)
        writer.writeheader()
        for r in rows:
            writer.writerow(r)


def filter_rows(rows, template, cid):
    """Return (new_rows, removed_count) filtering out rows with template+cid."""
    kept = []
    removed = 0
    for r in rows:
        if r.get("BM_Template") == template and r.get("Docked") == cid:
            removed += 1
        else:
            kept.append(r)
    return kept, removed


def find_rows(rows, template, cid):
    """Return list of rows matching template+cid."""
    return [
        r for r in rows
        if r.get("BM_Template") == template and r.get("Docked") == cid
    ]


# =========================
# Phases
# =========================

def phase_manual_or_stereo(
    rows,
    complexes_root,
    pairs_csv,
    phase_name,
    log,
    stats,
):
    """
    Generic handler for manual and stereo removal CSVs.
    pairs_csv: path with columns template, ID
    phase_name: "manual" or "stereo"
    """
    if not pairs_csv:
        return rows

    log.write(f"\n=== Phase: {phase_name.upper()} deletions ===\n")
    total_pairs = 0
    total_rows_removed = 0
    total_compounds_removed = 0
    total_templates_removed = 0

    with open(pairs_csv, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for rec in reader:
            template = rec.get("template", "").strip()
            cid = rec.get("ID", "").strip()
            if not template or not cid:
                continue
            total_pairs += 1

            # Affinities removal
            matches = find_rows(rows, template, cid)
            log.write(
                f"[{phase_name}] template={template}, ID={cid}, "
                f"rows_found={len(matches)}\n"
            )
            rows, removed = filter_rows(rows, template, cid)
            total_rows_removed += removed

            # Pose removal
            template_path = os.path.join(complexes_root, template)
            compound_path = os.path.join(template_path, cid)
            if ensure_dir_deleted(compound_path):
                total_compounds_removed += 1
                log.write(
                    f"  Pose folder removed: {compound_path}\n"
                )
                if remove_template_if_empty(complexes_root, template):
                    total_templates_removed += 1
                    log.write(
                        f"  Template folder removed (now empty): {template_path}\n"
                    )
            else:
                log.write(
                    f"  Pose folder not found: {compound_path}\n"
                )

    log.write(
        f"\n[{phase_name}] summary: pairs={total_pairs}, "
        f"rows_removed={total_rows_removed}, "
        f"compounds_removed={total_compounds_removed}, "
        f"templates_removed={total_templates_removed}\n"
    )

    stats[f"{phase_name}_rows_removed"] += total_rows_removed
    stats[f"{phase_name}_compounds_removed"] += total_compounds_removed
    stats[f"{phase_name}_templates_removed"] += total_templates_removed

    return rows


def phase_salts(rows, complexes_root, salts_csv, log, stats):
    """
    Handle salts resolution:
    - remove salts
    - rescue unique affinity types into desalted compound
    - delete salt pose folders
    """
    if not salts_csv:
        return rows

    log.write("\n=== Phase: SALTS resolution ===\n")
    total_pairs = 0
    total_salt_rows_removed = 0
    total_rescued_rows = 0
    total_compounds_removed = 0
    total_templates_removed = 0
    total_pairs_without_desalted_rows = 0

    with open(salts_csv, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for rec in reader:
            template = rec.get("template", "").strip()
            desalted = rec.get("desalted", "").strip()
            salt = rec.get("salt", "").strip()
            if not template or not desalted or not salt:
                continue
            total_pairs += 1
            log.write(
                f"[salts] template={template}, desalted={desalted}, salt={salt}\n"
            )

            D_set = find_rows(rows, template, desalted)
            S_set = find_rows(rows, template, salt)

            log.write(
                f"  desalted_rows={len(D_set)}, salt_rows={len(S_set)}\n"
            )

            # If no desalted rows, we cannot rescue anything
            if not D_set:
                total_pairs_without_desalted_rows += 1
                # just delete salt rows
                rows, removed = filter_rows(rows, template, salt)
                total_salt_rows_removed += removed
                log.write(
                    "  No desalted rows found; deleting salt rows without rescue.\n"
                )
            else:
                # Clone rows for rescue before we delete salts
                # Build set of affinity types present in desalted
                desalted_types = {r.get("Affinity_Type") for r in D_set}
                # Remove salt rows from main table
                rows, removed = filter_rows(rows, template, salt)
                total_salt_rows_removed += removed

                # For each salt row, if its Affinity_Type is not
                # present in D_set, create a new row for desalted
                for srow in S_set:
                    T_salt = srow.get("Affinity_Type")
                    if not T_salt:
                        continue
                    if T_salt in desalted_types:
                        continue
                    base = D_set[0].copy()
                    base["Affinity_Type"] = T_salt
                    base["Binding_Affinity"] = srow.get("Binding_Affinity")
                    base["log_affinity"] = srow.get("log_affinity")
                    rows.append(base)
                    total_rescued_rows += 1
                    desalted_types.add(T_salt)
                    log.write(
                        f"  Rescued affinity_type={T_salt} from salt to desalted.\n"
                    )

            # Pose deletion for salt
            template_path = os.path.join(complexes_root, template)
            salt_path = os.path.join(template_path, salt)
            if ensure_dir_deleted(salt_path):
                total_compounds_removed += 1
                log.write(
                    f"  Pose folder removed (salt): {salt_path}\n"
                )
                if remove_template_if_empty(complexes_root, template):
                    total_templates_removed += 1
                    log.write(
                        f"  Template folder removed (now empty): {template_path}\n"
                    )
            else:
                log.write(
                    f"  Pose folder not found (salt): {salt_path}\n"
                )

    log.write(
        f"\n[salts] summary: pairs={total_pairs}, "
        f"pairs_without_desalted_rows={total_pairs_without_desalted_rows}, "
        f"salt_rows_removed={total_salt_rows_removed}, "
        f"rescued_rows_added={total_rescued_rows}, "
        f"compounds_removed={total_compounds_removed}, "
        f"templates_removed={total_templates_removed}\n"
    )

    stats["salts_rows_removed"] += total_salt_rows_removed
    stats["salts_rows_rescued"] += total_rescued_rows
    stats["salts_compounds_removed"] += total_compounds_removed
    stats["salts_templates_removed"] += total_templates_removed

    return rows


def compute_inchikey_for_sdf(path, log=None):
    """Load SDF and compute Standard InChIKey; return None if fails."""
    try:
        mol = Chem.MolFromMolFile(path, sanitize=True)
        if mol is None:
            if log:
                log.write(f"  Failed to load SDF: {path}\n")
            return None
        return inchi.MolToInchiKey(mol)
    except Exception as e:
        if log:
            log.write(f"  Exception reading SDF {path}: {e}\n")
        return None


def phase_inchikey_dedup(rows, complexes_root, log, stats):
    """
    Recompute intra-template duplicates using Standard InChIKey and 
    remove all but the lowest-RMSD compound in each group. 
    Keep always CHEMBL affinity values when possible
    """
    log.write("\n=== Phase: InChIKey-based RMSD deduplication ===\n")

    total_templates = 0
    templates_with_duplicates = 0
    total_duplicate_groups = 0
    total_rows_removed = 0
    total_compounds_removed = 0
    total_templates_removed = 0
    total_chembl_affinity_transfers = 0

    for template in sorted(os.listdir(complexes_root)):
        tpath = os.path.join(complexes_root, template)
        if not os.path.isdir(tpath):
            continue
        total_templates += 1

        compounds = [
            d for d in os.listdir(tpath)
            if os.path.isdir(os.path.join(tpath, d))
        ]
        if len(compounds) < 2:
            continue

        ik_groups = defaultdict(list)
        for cid in compounds:
            sdf_path = os.path.join(tpath, cid, f"{cid}.sdf")
            if not os.path.isfile(sdf_path):
                log.write(
                    f"Template {template}: SDF not found for {cid} at {sdf_path}\n"
                )
                continue
            ik = compute_inchikey_for_sdf(sdf_path, log=log)
            if ik is None:
                continue
            ik_groups[ik].append(cid)

        num_compounds = len(compounds)
        unique_compounds = len(ik_groups)
        dup_to_remove = num_compounds - unique_compounds
        if dup_to_remove <= 0:
            continue

        templates_with_duplicates += 1
        log.write(
            f"\nTemplate: {template}\n"
            f"Number of compounds: {num_compounds}\n"
            f"Number of unique compounds: {unique_compounds}\n"
            f"Compounds to remove for uniqueness: {dup_to_remove}\n\n"
        )

        for idx, (ik, cids) in enumerate(sorted(ik_groups.items()), start=1):
            if len(cids) < 2:
                continue

            total_duplicate_groups += 1
            log.write(f"  same compound group {idx} (InChIKey = {ik}):\n")
            for cid in cids:
                log.write(f"    {cid}\n")
            log.write("\n")

            candidate_data = []
            rows_by_cid = {}

            for cid in cids:
                matches = find_rows(rows, template, cid)
                rows_by_cid[cid] = matches

                if not matches:
                    candidate_data.append((cid, None))
                    continue

                rmsd_str = matches[0].get("mcs_rmsd")
                try:
                    rmsd_val = float(rmsd_str) if rmsd_str is not None else None
                except ValueError:
                    rmsd_val = None

                candidate_data.append((cid, rmsd_val))

            valid_candidates = [
                (cid, rmsd) for cid, rmsd in candidate_data
                if rmsd is not None
            ]

            if valid_candidates:
                winner_cid, winner_rmsd = min(valid_candidates, key=lambda x: x[1])
            else:
                winner_cid, winner_rmsd = candidate_data[0]

            log.write(
                f"  RMSD selection: winner={winner_cid}, mcs_rmsd={winner_rmsd}\n"
            )

            loser_cids = {cid for cid, _ in candidate_data if cid != winner_cid}

            # CHEMBL affinity priority rule:
            # If winner is BINDINGDB, use deleted CHEMBL affinities for matching Affinity_Type.
            if winner_cid.startswith("BINDINGDB"):
                winner_rows = rows_by_cid.get(winner_cid, [])
                chembl_loser_rows = []

                for lc in loser_cids:
                    if lc.startswith("CHEMBL"):
                        chembl_loser_rows.extend(rows_by_cid.get(lc, []))

                chembl_by_type = {}
                for r in chembl_loser_rows:
                    affinity_type = r.get("Affinity_Type")
                    if affinity_type and affinity_type not in chembl_by_type:
                        chembl_by_type[affinity_type] = r

                for wrow in winner_rows:
                    affinity_type = wrow.get("Affinity_Type")
                    if affinity_type in chembl_by_type:
                        source = chembl_by_type[affinity_type]
                        old_aff = wrow.get("Binding_Affinity")
                        old_log = wrow.get("log_affinity")

                        wrow["Binding_Affinity"] = source.get("Binding_Affinity")
                        wrow["log_affinity"] = source.get("log_affinity")

                        total_chembl_affinity_transfers += 1
                        log.write(
                            f"    CHEMBL affinity transfer for {winner_cid}, "
                            f"Affinity_Type={affinity_type}: "
                            f"Binding_Affinity {old_aff} -> {wrow.get('Binding_Affinity')}, "
                            f"log_affinity {old_log} -> {wrow.get('log_affinity')}\n"
                        )

            for lc in loser_cids:
                rows, removed = filter_rows(rows, template, lc)
                total_rows_removed += removed

                if removed > 0:
                    log.write(
                        f"    Removed {removed} rows for loser {lc} in affinities.\n"
                    )

                lpath = os.path.join(tpath, lc)
                if ensure_dir_deleted(lpath):
                    total_compounds_removed += 1
                    log.write(f"    Pose folder removed for loser {lc}: {lpath}\n")

            if remove_template_if_empty(complexes_root, template):
                total_templates_removed += 1
                log.write(
                    f"  Template folder removed after dedup (now empty): {tpath}\n"
                )
                break

    log.write(
        f"\n[InChIKey dedup] summary: templates_scanned={total_templates}, "
        f"templates_with_duplicates={templates_with_duplicates}, "
        f"duplicate_groups={total_duplicate_groups}, "
        f"rows_removed={total_rows_removed}, "
        f"compounds_removed={total_compounds_removed}, "
        f"templates_removed={total_templates_removed}, "
        f"chembl_affinity_transfers={total_chembl_affinity_transfers}\n"
    )

    stats["inchikey_rows_removed"] += total_rows_removed
    stats["inchikey_compounds_removed"] += total_compounds_removed
    stats["inchikey_templates_removed"] += total_templates_removed
    stats["inchikey_chembl_affinity_transfers"] += total_chembl_affinity_transfers

    return rows


# =========================
# Main
# =========================

def main():
    parser = argparse.ArgumentParser(
        description="DockBind final cleaning & deduplication."
    )
    parser.add_argument("affinities_csv", help="Original DockBind affinities CSV.")
    parser.add_argument("complexes_root", help="Original complexes root folder.")
    parser.add_argument(
        "--manual-csv",
        help="CSV with columns 'template','ID' for manual deletions.",
    )
    parser.add_argument(
        "--salts-csv",
        help="CSV with columns 'template','desalted','salt' for salt resolution.",
    )
    parser.add_argument(
        "--stereo-csv",
        help="CSV with columns 'template','ID' for stereo-based deletions.",
    )
    parser.add_argument(
        "--out-affinities",
        help="Output affinities CSV path (default: basename + _deduped.csv in CWD).",
    )
    parser.add_argument(
        "--out-complexes",
        help="Output complexes root folder (default: basename + _deduped in CWD).",
    )
    parser.add_argument(
        "--log",
        default="dockbind_dedup.log",
        help="Log file path (default: dockbind_dedup.log).",
    )
    args = parser.parse_args()

    # Determine outputs
    aff_base = os.path.splitext(os.path.basename(args.affinities_csv))[0]
    out_aff = args.out_affinities or f"{aff_base}_deduped.csv"

    comp_base = os.path.basename(os.path.normpath(args.complexes_root))
    out_comp = args.out_complexes or f"{comp_base}_deduped"

    # Prepare output complexes tree
    if os.path.exists(out_comp):
        raise RuntimeError(f"Output complexes dir already exists: {out_comp}")
    shutil.copytree(args.complexes_root, out_comp)

    # Load affinities
    rows, header = load_affinities_csv(args.affinities_csv)

    # Stats accumulator
    stats = defaultdict(int)

    with open(args.log, "w", encoding="utf-8") as log:
        log.write("DockBind final cleaning & deduplication\n")
        log.write("======================================\n\n")
        log.write(f"Input affinities: {args.affinities_csv}\n")
        log.write(f"Input complexes:  {args.complexes_root}\n")
        log.write(f"Output affinities: {out_aff}\n")
        log.write(f"Output complexes:  {out_comp}\n\n")

        # Phase 1: manual
        rows = phase_manual_or_stereo(
            rows, out_comp, args.manual_csv, "manual", log, stats
        )

        # Phase 2: salts
        rows = phase_salts(rows, out_comp, args.salts_csv, log, stats)

        # Phase 3: stereo
        rows = phase_manual_or_stereo(
            rows, out_comp, args.stereo_csv, "stereo", log, stats
        )

        # Phase 4: InChIKey dedup + RMSD choice
        rows = phase_inchikey_dedup(rows, out_comp, log, stats)

        # Global summary
        total_rows_removed = (
            stats["manual_rows_removed"]
            + stats["salts_rows_removed"]
            + stats["stereo_rows_removed"]
            + stats["inchikey_rows_removed"]
        )
        total_compounds_removed = (
            stats["manual_compounds_removed"]
            + stats["salts_compounds_removed"]
            + stats["stereo_compounds_removed"]
            + stats["inchikey_compounds_removed"]
        )
        total_templates_removed = (
            stats["manual_templates_removed"]
            + stats["salts_templates_removed"]
            + stats["stereo_templates_removed"]
            + stats["inchikey_templates_removed"]
        )

        log.write("\n=== GLOBAL SUMMARY ===\n")
        log.write(f"Total rows removed: {total_rows_removed}\n")
        log.write(f"  - manual:  {stats['manual_rows_removed']}\n")
        log.write(f"  - salts:   {stats['salts_rows_removed']}\n")
        log.write(f"  - stereo:  {stats['stereo_rows_removed']}\n")
        log.write(f"  - InChIKey:{stats['inchikey_rows_removed']}\n")

        log.write(f"\nTotal compounds removed: {total_compounds_removed}\n")
        log.write(f"  - manual:  {stats['manual_compounds_removed']}\n")
        log.write(f"  - salts:   {stats['salts_compounds_removed']}\n")
        log.write(f"  - stereo:  {stats['stereo_compounds_removed']}\n")
        log.write(f"  - InChIKey:{stats['inchikey_compounds_removed']}\n")

        log.write(f"\nTotal templates removed: {total_templates_removed}\n")
        log.write(f"  - manual:  {stats['manual_templates_removed']}\n")
        log.write(f"  - salts:   {stats['salts_templates_removed']}\n")
        log.write(f"  - stereo:  {stats['stereo_templates_removed']}\n")
        log.write(f"  - InChIKey:{stats['inchikey_templates_removed']}\n")

    # Write final affinities CSV
    write_affinities_csv(rows, header, out_aff)


if __name__ == "__main__":
    main()

