#!/usr/bin/env python3

import csv
import importlib.util
import os
import sys


def reverse_complement_rna(seq):
    comp = str.maketrans({"A": "U", "U": "A", "G": "C", "C": "G"})
    return seq.translate(comp)[::-1]


def overlapping_positions(seq, motif):
    positions = []
    start = 0
    while True:
        pos = seq.find(motif, start)
        if pos == -1:
            break
        positions.append(pos + 1)
        start = pos + 1
    return positions


def read_tsv(path):
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path, rows, fieldnames):
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def parse_int(value):
    try:
        return int(value)
    except Exception:
        return 0


def main():
    out_dir = os.path.abspath(sys.argv[1])
    utr_table = read_tsv(os.path.join(out_dir, "transcript_utr_reconstruction.tsv"))
    script_path = os.path.join(os.path.dirname(__file__), "sirna_specificity.py")
    spec = importlib.util.spec_from_file_location("sirna_specificity_module", script_path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    if len(sys.argv) > 2:
        module.configure_seed_length(int(sys.argv[2]))

    strand_rows = []
    for sirna, seq in module.SIRNAS.items():
      derived = module.derive_strands(seq)
      derived["sirna"] = sirna
      strand_rows.append(derived)

    derivation_fieldnames = list(strand_rows[0].keys())
    with open(os.path.join(out_dir, "sirna_strand_derivation.csv"), "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=derivation_fieldnames)
        writer.writeheader()
        writer.writerows(strand_rows)

    seed_label = module.SEED_LABEL

    seeds = {}
    for row in strand_rows:
        seeds[f"{row['sirna']}_guide"] = row[f"guide_seed_{seed_label}"]
        seeds[f"{row['sirna']}_passenger"] = row[f"passenger_seed_{seed_label}"]

    tx_rows = []
    gene_rows = {}

    for row in utr_table:
        if parse_int(row.get("utr_reconstructable", 0)) != 1:
            continue
        utr_seq = row.get("utr_seq", "")
        entry = {
            "transcript_id": row["transcript_id"],
            "transcript_id_stable": row["transcript_id_stable"],
            "gene_id": row["gene_id"],
            "gene_id_stable": row["gene_id_stable"],
            "gene_name": row["gene_name"],
            "transcript_biotype": row["transcript_biotype"],
            "gene_biotype": row["gene_biotype"],
            "utr_length": parse_int(row.get("computed_utr_len", 0)),
        }
        for seed_name, seed in seeds.items():
            target = reverse_complement_rna(seed)
            positions = overlapping_positions(utr_seq, target)
            entry[f"{seed_name}_seed"] = seed
            entry[f"{seed_name}_target_motif"] = target
            entry[f"{seed_name}_match_count"] = len(positions)
            entry[f"{seed_name}_has_match"] = 1 if positions else 0
        tx_rows.append(entry)

        gkey = row["gene_id_stable"]
        if gkey not in gene_rows:
            gene_rows[gkey] = {
                "gene_id": row["gene_id"],
                "gene_id_stable": row["gene_id_stable"],
                "gene_name": row["gene_name"],
                "gene_biotype": row["gene_biotype"],
                "max_utr_length": parse_int(row.get("computed_utr_len", 0)),
                "n_utr_transcripts": 1,
            }
            for seed_name in seeds:
                gene_rows[gkey][f"{seed_name}_has_match"] = entry[f"{seed_name}_has_match"]
                gene_rows[gkey][f"{seed_name}_match_count"] = entry[f"{seed_name}_match_count"]
        else:
            gene_rows[gkey]["max_utr_length"] = max(gene_rows[gkey]["max_utr_length"], parse_int(row.get("computed_utr_len", 0)))
            gene_rows[gkey]["n_utr_transcripts"] += 1
            for seed_name in seeds:
                gene_rows[gkey][f"{seed_name}_has_match"] = max(
                    gene_rows[gkey][f"{seed_name}_has_match"],
                    entry[f"{seed_name}_has_match"]
                )
                gene_rows[gkey][f"{seed_name}_match_count"] += entry[f"{seed_name}_match_count"]

    tx_fieldnames = list(tx_rows[0].keys())
    gene_fieldnames = list(next(iter(gene_rows.values())).keys())
    write_tsv(os.path.join(out_dir, "sirna_seed_annotations.tsv"), tx_rows, tx_fieldnames)
    write_tsv(os.path.join(out_dir, "sirna_seed_annotations_genes.tsv"), list(gene_rows.values()), gene_fieldnames)

    summary_rows = []
    near_counts = {}
    near_path = os.path.join(out_dir, "candidate_offtargets_nearperfect.tsv")
    if os.path.exists(near_path):
        near_rows = read_tsv(near_path)
        for row in near_rows:
            key = (row["sirna"], row["strand_role"], row["mismatches"])
            near_counts[key] = near_counts.get(key, 0) + 1
    for row in strand_rows:
        for strand_role in ("guide", "passenger"):
            key = f"{row['sirna']}_{strand_role}_has_match"
            summary_rows.append({
                "sirna": row["sirna"],
                "strand_role": strand_role,
                "candidate_sequence": row[f"{strand_role}_candidate"],
                f"seed_{seed_label}": row[f"{strand_role}_seed_{seed_label}"],
                "n_nearperfect_hits_0mm": near_counts.get((row["sirna"], strand_role, "0"), 0),
                "n_nearperfect_hits_1mm": near_counts.get((row["sirna"], strand_role, "1"), 0),
                "n_genes_with_3utr_seed_match": sum(parse_int(g[key]) for g in gene_rows.values()),
            })

    csv_fieldnames = list(summary_rows[0].keys())
    with open(os.path.join(out_dir, "sirna_specificity_summary.csv"), "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=csv_fieldnames)
        writer.writeheader()
        writer.writerows(summary_rows)
    write_tsv(os.path.join(out_dir, "supplemental_sirna_specificity_summary.tsv"), summary_rows, csv_fieldnames)


if __name__ == "__main__":
    main()
