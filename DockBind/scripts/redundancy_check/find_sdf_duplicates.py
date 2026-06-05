#!/usr/bin/env python3
"""
Find intra-template duplicate compounds using Standard InChIKey.
"""


import argparse
import os
import sys
from collections import defaultdict

from rdkit import Chem
from rdkit.Chem import inchi


def find_templates(root_dir: str):
    """Yield absolute paths to template folders (top-level directories in root_dir)."""
    with os.scandir(root_dir) as it:
        for entry in it:
            if entry.is_dir():
                yield entry.path


def find_compound_folders(template_dir: str):
    """Yield (compound_name, absolute_path) for each compound folder in a template."""
    with os.scandir(template_dir) as it:
        for entry in it:
            if entry.is_dir():
                yield entry.name, entry.path


def get_sdf_path_for_compound(compound_name: str, compound_dir: str):
    sdf_path = os.path.join(compound_dir, f"{compound_name}.sdf")
    return sdf_path if os.path.isfile(sdf_path) else None


def mol_to_inchikey(sdf_path: str):
    try:
        mol = Chem.MolFromMolFile(sdf_path, removeHs=False)
    except Exception as exc:
        sys.stderr.write(f"[WARN] RDKit error reading SDF '{sdf_path}': {exc}\n")
        return None

    if mol is None:
        sys.stderr.write(f"[WARN] Failed to parse SDF '{sdf_path}'\n")
        return None

    try:
        return inchi.MolToInchiKey(mol)
    except Exception as exc:
        sys.stderr.write(f"[WARN] InChIKey generation failed for '{sdf_path}': {exc}\n")
        return None


def process_template(template_dir: str):
    template_name = os.path.basename(template_dir)
    inchikey_to_compounds = defaultdict(list)
    valid_compounds = 0

    for compound_name, compound_dir in find_compound_folders(template_dir):
        sdf_path = get_sdf_path_for_compound(compound_name, compound_dir)
        if sdf_path is None:
            sys.stderr.write(
                f"[WARN] Missing SDF '{compound_name}.sdf' in '{compound_dir}'\n"
            )
            continue

        key = mol_to_inchikey(sdf_path)
        if key is None:
            continue

        valid_compounds += 1
        inchikey_to_compounds[key].append(compound_name)

    if valid_compounds == 0:
        return {
            "template_name": template_name,
            "n_compounds": 0,
            "n_unique": 0,
            "to_remove": 0,
            "has_duplication": False,
            "duplicate_groups": [],
        }

    n_unique = len(inchikey_to_compounds)
    to_remove = valid_compounds - n_unique

    duplicate_groups = [
        {"inchi_key": k, "compounds": sorted(v)}
        for k, v in inchikey_to_compounds.items()
        if len(v) > 1
    ]

    return {
        "template_name": template_name,
        "n_compounds": valid_compounds,
        "n_unique": n_unique,
        "to_remove": to_remove,
        "has_duplication": bool(duplicate_groups),
        "duplicate_groups": duplicate_groups,
    }


def write_log(output_path, global_stats, templates_with_duplicates):
    lines = []

    lines.append("GLOBAL SUMMARY")
    lines.append("==============")
    for k, v in global_stats.items():
        lines.append(f"{k.replace('_', ' ').capitalize()}: {v}")
    lines.append("")
    lines.append("TEMPLATE DUPLICATION DETAILS")
    lines.append("============================\n")

    if not templates_with_duplicates:
        lines.append("No templates with intra-template duplications were found.")
    else:
        for info in templates_with_duplicates:
            lines.append(f"Template: {info['template_name']}")
            lines.append(f"Number of compounds: {info['n_compounds']}")
            lines.append(f"Number of unique compounds: {info['n_unique']}")
            lines.append(f"Compounds to remove for uniqueness: {info['to_remove']}\n")
            for i, group in enumerate(info["duplicate_groups"], 1):
                lines.append(
                    f"  same compound group {i} (InChIKey = {group['inchi_key']}):"
                )
                for name in group["compounds"]:
                    lines.append(f"    {name}")
                lines.append("")
            lines.append("-" * 40 + "\n")

    with open(output_path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines))

    print(f"\nLog written to: {output_path}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("root_dir")
    parser.add_argument("-o", "--output", default="dedup_report.log")
    args = parser.parse_args()

    root_dir = os.path.abspath(args.root_dir)
    templates = list(find_templates(root_dir))   # NEW
    total_templates = len(templates)             # NEW

    stats = {
        "total_templates": total_templates,
        "templates_with_single_compound": 0,
        "templates_with_more_than_one_compound": 0,
        "templates_with_at_least_one_duplication": 0,
        "total_sdf_poses": 0,
        "total_unique_poses": 0,
        "total_to_remove": 0,
    }

    templates_with_duplicates = []

    for idx, template_dir in enumerate(templates, start=1):   # NEW
        info = process_template(template_dir)

        # Progress reporting (NEW)
        percent = (idx / total_templates) * 100
        print(f"Progress: {percent:6.2f}% ({idx}/{total_templates})", end="\r")

        n = info["n_compounds"]
        if n == 0:
            continue

        stats["total_sdf_poses"] += n
        stats["total_unique_poses"] += info["n_unique"]
        stats["total_to_remove"] += info["to_remove"]

        if n == 1:
            stats["templates_with_single_compound"] += 1
        else:
            stats["templates_with_more_than_one_compound"] += 1

        if info["has_duplication"]:
            stats["templates_with_at_least_one_duplication"] += 1
            templates_with_duplicates.append(info)

    write_log(args.output, stats, templates_with_duplicates)


if __name__ == "__main__":
    main()

