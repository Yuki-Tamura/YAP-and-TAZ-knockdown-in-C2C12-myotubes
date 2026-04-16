#!/usr/bin/env python3

import csv
import gzip
import math
import os
import re
import statistics
import subprocess
import sys
from collections import Counter, defaultdict


SIRNAS = {
    "YAP": "GAUACUUCUUAAAUCACAAtt",
    "TAZ": "AUACUUCCUUAAUCACAUAtt",
}

SEED_START = 1


def configure_seed_length(seed_len=None):
    global SEED_END, SEED_LABEL, SEED_LEN
    if seed_len is None:
        seed_len = int(os.environ.get("SIRNA_SEED_LEN", "6"))
    if seed_len < 1:
        raise ValueError("Seed length must be at least 1.")
    SEED_LEN = seed_len
    SEED_END = SEED_START + seed_len
    SEED_LABEL = f"{SEED_START + 1}_{SEED_END}"


configure_seed_length()


def normalize_rna(seq):
    seq = seq.strip()
    seq = re.sub(r"[^AUGCaugctT]", "", seq)
    return seq.upper().replace("T", "U")


def reverse_complement_rna(seq):
    comp = str.maketrans({"A": "U", "U": "A", "G": "C", "C": "G"})
    return seq.translate(comp)[::-1]


def strip_version(identifier):
    return identifier.split(".")[0]


def derive_strands(sense_seq):
    sense = normalize_rna(sense_seq)
    core = sense[:-2] if len(sense) >= 2 else sense
    passenger = sense
    guide_core = reverse_complement_rna(core)
    guide = guide_core + "UU"
    return {
        "input_sense": sense,
        "core_19nt": core,
        "passenger_candidate": passenger,
        "guide_candidate": guide,
        f"guide_seed_{SEED_LABEL}": guide[SEED_START:SEED_END] if len(guide) >= SEED_END else "",
        f"passenger_seed_{SEED_LABEL}": passenger[SEED_START:SEED_END] if len(passenger) >= SEED_END else "",
    }


def find_local_resources(project_root):
    resources = {
        "transcriptome_fasta": [],
        "utr_fasta": [],
        "annotation_gtf_gff": [],
    }
    for root, _, files in os.walk(project_root):
        for name in files:
            lower = name.lower()
            path = os.path.join(root, name)
            if lower.endswith((".fa", ".fasta", ".fna", ".fa.gz", ".fasta.gz", ".fna.gz")):
                if "utr" in lower or "3utr" in lower or "3'utr" in lower:
                    resources["utr_fasta"].append(path)
                else:
                    resources["transcriptome_fasta"].append(path)
            elif lower.endswith((".gtf", ".gtf.gz", ".gff", ".gff3")):
                resources["annotation_gtf_gff"].append(path)
    return resources


def open_text(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "rt")


def parse_attributes(attr_text):
    attrs = {}
    for key, value in re.findall(r'([A-Za-z0-9_]+) "([^"]*)"', attr_text):
        attrs[key] = value
    return attrs


def load_fasta_with_metadata(path):
    seqs = {}
    meta = {}
    current_id = None
    current_seq = []
    with open_text(path) as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_id is not None:
                    seqs[current_id] = "".join(current_seq).upper().replace("T", "U")
                header = line[1:]
                current_id = header.split()[0]
                current_seq = []
                header_meta = {"header": header, "transcript_id": current_id, "transcript_id_stable": strip_version(current_id)}
                gene_match = re.search(r"gene:([A-Za-z0-9_.-]+)", header)
                gene_symbol_match = re.search(r"gene_symbol:([^\s]+)", header)
                transcript_biotype_match = re.search(r"transcript_biotype:([^\s]+)", header)
                gene_biotype_match = re.search(r"gene_biotype:([^\s]+)", header)
                if gene_match:
                    header_meta["gene_id"] = gene_match.group(1)
                    header_meta["gene_id_stable"] = strip_version(gene_match.group(1))
                if gene_symbol_match:
                    header_meta["gene_name"] = gene_symbol_match.group(1)
                if transcript_biotype_match:
                    header_meta["transcript_biotype"] = transcript_biotype_match.group(1)
                if gene_biotype_match:
                    header_meta["gene_biotype"] = gene_biotype_match.group(1)
                meta[current_id] = header_meta
            else:
                current_seq.append(line)
    if current_id is not None:
        seqs[current_id] = "".join(current_seq).upper().replace("T", "U")
    return seqs, meta


def init_transcript_record(attrs, strand):
    transcript_id = attrs.get("transcript_id")
    transcript_version = attrs.get("transcript_version")
    transcript_full = transcript_id if not transcript_version else f"{transcript_id}.{transcript_version}"
    gene_id = attrs.get("gene_id")
    gene_version = attrs.get("gene_version")
    gene_full = gene_id if not gene_version else f"{gene_id}.{gene_version}"
    return {
        "transcript_id": transcript_full,
        "transcript_id_stable": transcript_id,
        "gene_id": gene_full,
        "gene_id_stable": gene_id,
        "gene_name": attrs.get("gene_name", ""),
        "gene_biotype": attrs.get("gene_biotype", ""),
        "transcript_name": attrs.get("transcript_name", ""),
        "transcript_biotype": attrs.get("transcript_biotype", ""),
        "strand": strand,
        "exons": [],
        "cds": [],
        "stop_codons": [],
        "three_prime_utr_features": [],
    }


def parse_gtf(path):
    transcripts = {}
    feature_counts = Counter()
    with open_text(path) as handle:
        for raw in handle:
            if raw.startswith("#"):
                continue
            fields = raw.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            feature = fields[2]
            start = int(fields[3])
            end = int(fields[4])
            strand = fields[6]
            attrs = parse_attributes(fields[8])
            transcript_id = attrs.get("transcript_id")
            if not transcript_id:
                continue
            feature_counts[feature] += 1
            if transcript_id not in transcripts:
                transcripts[transcript_id] = init_transcript_record(attrs, strand)
            rec = transcripts[transcript_id]
            if feature == "exon":
                rec["exons"].append((start, end))
            elif feature == "CDS":
                rec["cds"].append((start, end))
            elif feature == "stop_codon":
                rec["stop_codons"].append((start, end))
            elif feature == "three_prime_utr":
                rec["three_prime_utr_features"].append((start, end))
    return transcripts, feature_counts


def transcript_model_metrics(rec, tx_sequence):
    exons = list(rec["exons"])
    cds = list(rec["cds"])
    stop_codons = list(rec["stop_codons"])
    if not exons:
      return None
    strand = rec["strand"]
    ordered_exons = sorted(exons, key=lambda x: x[0], reverse=(strand == "-"))
    transcript_len = 0
    cds_t_start = None
    cds_t_end = None
    stop_t_start = None
    stop_t_end = None
    explicit_utr_t_coords = []
    for exon_start, exon_end in ordered_exons:
        exon_len = exon_end - exon_start + 1
        exon_t_start = transcript_len + 1
        exon_t_end = transcript_len + exon_len
        transcript_len += exon_len
        for feature_intervals, kind in ((cds, "cds"), (stop_codons, "stop")):
            for cds_start, cds_end in feature_intervals:
                ov_start = max(exon_start, cds_start)
                ov_end = min(exon_end, cds_end)
                if ov_start > ov_end:
                    continue
                if strand == "+":
                    local_start = ov_start - exon_start
                    local_end = ov_end - exon_start
                else:
                    local_start = exon_end - ov_end
                    local_end = exon_end - ov_start
                this_t_start = exon_t_start + local_start
                this_t_end = exon_t_start + local_end
                if kind == "cds":
                    cds_t_start = this_t_start if cds_t_start is None else min(cds_t_start, this_t_start)
                    cds_t_end = this_t_end if cds_t_end is None else max(cds_t_end, this_t_end)
                else:
                    stop_t_start = this_t_start if stop_t_start is None else min(stop_t_start, this_t_start)
                    stop_t_end = this_t_end if stop_t_end is None else max(stop_t_end, this_t_end)
        for utr_start_g, utr_end_g in rec["three_prime_utr_features"]:
            ov_start = max(exon_start, utr_start_g)
            ov_end = min(exon_end, utr_end_g)
            if ov_start > ov_end:
                continue
            if strand == "+":
                local_start = ov_start - exon_start
                local_end = ov_end - exon_start
            else:
                local_start = exon_end - ov_end
                local_end = exon_end - ov_start
            explicit_utr_t_coords.append((exon_t_start + local_start, exon_t_start + local_end))

    fasta_len = len(tx_sequence) if tx_sequence else None
    has_length_match = fasta_len == transcript_len if fasta_len is not None else False

    utr_start = None
    utr_end = None
    if cds_t_start is not None and cds_t_end is not None:
        coding_t_start = cds_t_start if stop_t_start is None else min(cds_t_start, stop_t_start)
        coding_t_end = cds_t_end if stop_t_end is None else max(cds_t_end, stop_t_end)
        if strand == "+":
            if coding_t_end < transcript_len:
                utr_start = coding_t_end + 1
                utr_end = transcript_len
        else:
            if coding_t_start > 1:
                utr_start = 1
                utr_end = coding_t_start - 1

    inferred_utr_seq = ""
    if utr_start is not None and utr_end is not None and has_length_match:
        inferred_utr_seq = tx_sequence[utr_start - 1:utr_end]

    explicit_three_prime_utr_len = sum(end - start + 1 for start, end in rec["three_prime_utr_features"])
    explicit_utr_seq = ""
    if explicit_utr_t_coords and has_length_match:
        pieces = [tx_sequence[start - 1:end] for start, end in sorted(explicit_utr_t_coords)]
        explicit_utr_seq = "".join(pieces)

    utr_seq = explicit_utr_seq if explicit_utr_seq else inferred_utr_seq
    computed_utr_len = len(utr_seq)
    return {
        "transcript_len_from_exons": transcript_len,
        "fasta_len": fasta_len,
        "has_length_match": has_length_match,
        "cds_t_start": cds_t_start,
        "cds_t_end": cds_t_end,
        "stop_t_start": stop_t_start,
        "stop_t_end": stop_t_end,
        "computed_utr_start": utr_start,
        "computed_utr_end": utr_end,
        "inferred_utr_len": len(inferred_utr_seq),
        "computed_utr_len": computed_utr_len,
        "explicit_three_prime_utr_len": explicit_three_prime_utr_len,
        "explicit_utr_seq_len": len(explicit_utr_seq),
        "utr_seq": utr_seq,
        "utr_source": "explicit_three_prime_utr" if explicit_utr_seq else ("inferred_from_cds_stop" if inferred_utr_seq else ""),
        "coding_with_cds": cds_t_start is not None and cds_t_end is not None,
        "utr_reconstructable": bool(utr_seq),
        "inferred_matches_explicit_len": int(bool(explicit_three_prime_utr_len) and len(inferred_utr_seq) == explicit_three_prime_utr_len),
    }


def generate_variants(query):
    bases = "AUCG"
    variants = {query: 0}
    for i, base in enumerate(query):
        for alt in bases:
            if alt != base:
                variants[query[:i] + alt + query[i + 1:]] = 1
    return variants


def overlapping_positions(seq, motif):
    start = 0
    positions = []
    while True:
        pos = seq.find(motif, start)
        if pos == -1:
            break
        positions.append(pos + 1)
        start = pos + 1
    return positions


def search_near_perfect(transcript_seqs, transcript_meta, queries):
    rows = []
    for query_name, query_seq in queries.items():
        target_motif = reverse_complement_rna(query_seq)
        variants = generate_variants(target_motif)
        for tx_id, seq in transcript_seqs.items():
            for motif, mismatches in variants.items():
                positions = overlapping_positions(seq, motif)
                if not positions:
                    continue
                txm = transcript_meta.get(tx_id, {})
                for pos in positions:
                    rows.append({
                        "query_name": query_name,
                        "sirna": query_name.split("_")[0],
                        "strand_role": query_name.split("_")[1],
                        "query_core_rna": query_seq,
                        "target_motif_rna": target_motif,
                        "transcript_id": tx_id,
                        "transcript_id_stable": strip_version(tx_id),
                        "gene_id": txm.get("gene_id", ""),
                        "gene_id_stable": txm.get("gene_id_stable", ""),
                        "gene_name": txm.get("gene_name", ""),
                        "transcript_biotype": txm.get("transcript_biotype", ""),
                        "position_1based": pos,
                        "mismatches": mismatches,
                        "matched_sequence_rna": motif,
                        "is_same_gene_symbol_as_intended": int(txm.get("gene_name", "").upper() == query_name.split("_")[0]),
                    })
    return rows


def build_seed_annotations(utr_table, derived_rows):
    seeds = {}
    for row in derived_rows:
        seeds[f"{row['sirna']}_guide"] = row[f"guide_seed_{SEED_LABEL}"]
        seeds[f"{row['sirna']}_passenger"] = row[f"passenger_seed_{SEED_LABEL}"]

    out_rows = []
    gene_collapse = {}
    for row in utr_table:
        if not row["utr_reconstructable"]:
            continue
        utr_seq = row["utr_seq"]
        entry = {
            "transcript_id": row["transcript_id"],
            "transcript_id_stable": row["transcript_id_stable"],
            "gene_id": row["gene_id"],
            "gene_id_stable": row["gene_id_stable"],
            "gene_name": row["gene_name"],
            "transcript_biotype": row["transcript_biotype"],
            "gene_biotype": row["gene_biotype"],
            "utr_length": int(row["computed_utr_len"]),
        }
        for seed_name, seed in seeds.items():
            target = reverse_complement_rna(seed)
            positions = overlapping_positions(utr_seq, target)
            entry[f"{seed_name}_seed"] = seed
            entry[f"{seed_name}_target_motif"] = target
            entry[f"{seed_name}_match_count"] = len(positions)
            entry[f"{seed_name}_has_match"] = 1 if positions else 0
        out_rows.append(entry)

        gkey = row["gene_id_stable"]
        if gkey not in gene_collapse:
            gene_collapse[gkey] = {
                "gene_id": row["gene_id"],
                "gene_id_stable": row["gene_id_stable"],
                "gene_name": row["gene_name"],
                "gene_biotype": row["gene_biotype"],
                "max_utr_length": int(row["computed_utr_len"]),
                "n_utr_transcripts": 1,
            }
            for seed_name in seeds:
                gene_collapse[gkey][f"{seed_name}_has_match"] = entry[f"{seed_name}_has_match"]
                gene_collapse[gkey][f"{seed_name}_match_count"] = entry[f"{seed_name}_match_count"]
        else:
            gene_collapse[gkey]["max_utr_length"] = max(gene_collapse[gkey]["max_utr_length"], int(row["computed_utr_len"]))
            gene_collapse[gkey]["n_utr_transcripts"] += 1
            for seed_name in seeds:
                gene_collapse[gkey][f"{seed_name}_has_match"] = max(gene_collapse[gkey][f"{seed_name}_has_match"], entry[f"{seed_name}_has_match"])
                gene_collapse[gkey][f"{seed_name}_match_count"] += entry[f"{seed_name}_match_count"]
    return out_rows, list(gene_collapse.values())


def write_tsv(path, rows, fieldnames):
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def write_csv(path, rows, fieldnames):
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def write_resource_status(path, resources):
    with open(path, "w") as handle:
        handle.write("Local annotation resource scan\n\n")
        for key, paths in resources.items():
            handle.write(f"{key}: {'FOUND' if paths else 'MISSING'}\n")
            for p in paths:
                handle.write(f"  - {p}\n")


def call_r_stats(project_root, out_dir):
    script = os.path.join(project_root, "analysis", "R", "sirna_seed_stats.R")
    cmd = ["Rscript", script, project_root, out_dir]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(f"R seed statistics failed.\nSTDOUT:\n{proc.stdout}\nSTDERR:\n{proc.stderr}")
    return proc.stdout, proc.stderr


def main():
    out_dir = sys.argv[1] if len(sys.argv) > 1 else "."
    out_dir = os.path.abspath(out_dir)
    analysis_dir = os.path.dirname(out_dir)
    project_root = os.path.abspath(sys.argv[2]) if len(sys.argv) > 2 else os.path.dirname(analysis_dir)
    os.makedirs(out_dir, exist_ok=True)

    derived_rows = []
    for sirna, seq in SIRNAS.items():
        derived = derive_strands(seq)
        derived["sirna"] = sirna
        derived_rows.append(derived)

    write_csv(
        os.path.join(out_dir, "sirna_strand_derivation.csv"),
        derived_rows,
        list(derived_rows[0].keys())
    )

    derivation_notes = [
        "siRNA strand derivation notes",
        "",
        "Input sequences were treated as sense/passenger strands with canonical 3' dinucleotide overhangs.",
        "The final two nucleotides were treated as overhangs and removed to define a 19-nt duplex core.",
        "The candidate guide strand was derived as the RNA reverse complement of the 19-nt core, then assigned a 3' UU overhang for a duplex-style representation.",
        "This is a supportive approximation and does not prove the strand-loading preference used in cells.",
    ]
    with open(os.path.join(out_dir, "sirna_derivation_notes.txt"), "w") as handle:
        handle.write("\n".join(derivation_notes) + "\n")

    resources = find_local_resources(project_root)
    write_resource_status(os.path.join(out_dir, "required_inputs_status.txt"), resources)

    transcriptome_candidates = [p for p in resources["transcriptome_fasta"] if "cdna" in os.path.basename(p).lower()]
    gtf_candidates = [p for p in resources["annotation_gtf_gff"] if os.path.basename(p).lower().endswith((".gtf", ".gtf.gz"))]
    if not transcriptome_candidates or not gtf_candidates:
        raise RuntimeError("Required Ensembl cDNA FASTA and GTF were not found locally.")

    transcriptome_path = transcriptome_candidates[0]
    gtf_path = gtf_candidates[0]
    transcript_seqs, transcript_meta = load_fasta_with_metadata(transcriptome_path)
    gtf_records, feature_counts = parse_gtf(gtf_path)

    utr_table = []
    length_match_count = 0
    coding_with_cds = 0
    reconstructable = 0
    explicit_utr_with_match = 0
    explicit_utr_checked = 0
    explicit_utr_used = 0
    inferred_utr_used = 0

    for stable_tx_id, rec in gtf_records.items():
        versioned_id = rec["transcript_id"]
        tx_sequence = transcript_seqs.get(versioned_id)
        if tx_sequence is None:
            continue
        if versioned_id not in transcript_meta:
            transcript_meta[versioned_id] = {}
        transcript_meta[versioned_id].update({
            "gene_id": rec["gene_id"],
            "gene_id_stable": rec["gene_id_stable"],
            "gene_name": rec["gene_name"],
            "transcript_biotype": rec["transcript_biotype"],
            "gene_biotype": rec["gene_biotype"],
        })
        metrics = transcript_model_metrics(rec, tx_sequence)
        if metrics is None:
            continue
        if metrics["has_length_match"]:
            length_match_count += 1
        if metrics["coding_with_cds"]:
            coding_with_cds += 1
        if metrics["utr_reconstructable"]:
            reconstructable += 1
        if rec["three_prime_utr_features"]:
            explicit_utr_checked += 1
            if metrics["inferred_matches_explicit_len"]:
                explicit_utr_with_match += 1
        if metrics["utr_source"] == "explicit_three_prime_utr":
            explicit_utr_used += 1
        elif metrics["utr_source"] == "inferred_from_cds_stop":
            inferred_utr_used += 1
        utr_table.append({
            "transcript_id": rec["transcript_id"],
            "transcript_id_stable": rec["transcript_id_stable"],
            "gene_id": rec["gene_id"],
            "gene_id_stable": rec["gene_id_stable"],
            "gene_name": rec["gene_name"],
            "gene_biotype": rec["gene_biotype"],
            "transcript_biotype": rec["transcript_biotype"],
            "strand": rec["strand"],
            "transcript_len_from_exons": metrics["transcript_len_from_exons"],
            "fasta_len": metrics["fasta_len"],
            "has_length_match": int(metrics["has_length_match"]),
            "coding_with_cds": int(metrics["coding_with_cds"]),
            "computed_utr_start": metrics["computed_utr_start"] or "",
            "computed_utr_end": metrics["computed_utr_end"] or "",
            "utr_source": metrics["utr_source"],
            "inferred_utr_len": metrics["inferred_utr_len"],
            "computed_utr_len": metrics["computed_utr_len"],
            "explicit_three_prime_utr_len": metrics["explicit_three_prime_utr_len"],
            "explicit_utr_seq_len": metrics["explicit_utr_seq_len"],
            "inferred_matches_explicit_len": metrics["inferred_matches_explicit_len"],
            "utr_reconstructable": int(metrics["utr_reconstructable"]),
            "utr_seq": metrics["utr_seq"],
        })

    write_tsv(
        os.path.join(out_dir, "transcript_utr_reconstruction.tsv"),
        utr_table,
        list(utr_table[0].keys())
    )

    with open(os.path.join(out_dir, "utr_reconstruction_status.txt"), "w") as handle:
        handle.write("3'UTR reconstruction status\n\n")
        handle.write(f"Transcriptome FASTA used: {transcriptome_path}\n")
        handle.write(f"GTF used: {gtf_path}\n")
        handle.write(f"Transcript sequences loaded: {len(transcript_seqs)}\n")
        handle.write(f"GTF transcript records loaded: {len(gtf_records)}\n")
        handle.write(f"Transcripts with exon-length = cDNA-length agreement: {length_match_count}\n")
        handle.write(f"Transcripts with CDS annotation among cDNA transcripts: {coding_with_cds}\n")
        handle.write(f"Transcripts with reconstructable 3'UTRs: {reconstructable}\n")
        handle.write(f"Reconstructed using explicit three_prime_utr features: {explicit_utr_used}\n")
        handle.write(f"Reconstructed using CDS/stop-codon inference fallback: {inferred_utr_used}\n")
        handle.write(f"Transcripts with explicit three_prime_utr feature checked: {explicit_utr_checked}\n")
        if explicit_utr_checked > 0:
            handle.write(f"Exact match between CDS/stop-derived and explicit 3'UTR length: {explicit_utr_with_match} ({100.0 * explicit_utr_with_match / explicit_utr_checked:.1f}%)\n")
        if reconstructable > 0:
            handle.write("\nConclusion: transcript-level 3'UTR reconstruction is feasible. Explicit three_prime_utr annotations were used whenever present, and CDS/stop-codon inference was used only as a fallback.\n")
        else:
            handle.write("\nConclusion: transcript-level 3'UTR reconstruction was not feasible from the available files alone.\n")

    queries = {}
    for row in derived_rows:
        queries[f"{row['sirna']}_guide"] = row["guide_candidate"][:-2]
        queries[f"{row['sirna']}_passenger"] = row["passenger_candidate"][:-2]

    near_perfect_rows = search_near_perfect(transcript_seqs, transcript_meta, queries)
    near_perfect_rows.sort(key=lambda x: (x["query_name"], x["mismatches"], x["gene_name"], x["transcript_id"], x["position_1based"]))
    write_tsv(
        os.path.join(out_dir, "candidate_offtargets_nearperfect.tsv"),
        near_perfect_rows,
        list(near_perfect_rows[0].keys()) if near_perfect_rows else [
            "query_name", "sirna", "strand_role", "query_core_rna", "target_motif_rna",
            "transcript_id", "transcript_id_stable", "gene_id", "gene_id_stable", "gene_name",
            "transcript_biotype", "position_1based", "mismatches", "matched_sequence_rna",
            "is_same_gene_symbol_as_intended"
        ]
    )

    seed_tx_rows, seed_gene_rows = build_seed_annotations(utr_table, derived_rows)
    write_tsv(
        os.path.join(out_dir, "sirna_seed_annotations.tsv"),
        seed_tx_rows,
        list(seed_tx_rows[0].keys()) if seed_tx_rows else []
    )
    write_tsv(
        os.path.join(out_dir, "sirna_seed_annotations_genes.tsv"),
        seed_gene_rows,
        list(seed_gene_rows[0].keys()) if seed_gene_rows else []
    )

    summary_rows = []
    near_counts = Counter((row["sirna"], row["strand_role"], row["mismatches"]) for row in near_perfect_rows)
    gene_level_match_counts = {}
    for row in seed_gene_rows:
        for sirna in ("YAP", "TAZ"):
            for strand_role in ("guide", "passenger"):
                key = f"{sirna}_{strand_role}_has_match"
                gene_level_match_counts.setdefault((sirna, strand_role), 0)
                gene_level_match_counts[(sirna, strand_role)] += int(row.get(key, 0))

    for row in derived_rows:
        for strand_role in ("guide", "passenger"):
            summary_rows.append({
                "sirna": row["sirna"],
                "strand_role": strand_role,
                "candidate_sequence": row[f"{strand_role}_candidate"],
                f"seed_{SEED_LABEL}": row[f"{strand_role}_seed_{SEED_LABEL}"],
                "n_nearperfect_hits_0mm": near_counts[(row["sirna"], strand_role, 0)],
                "n_nearperfect_hits_1mm": near_counts[(row["sirna"], strand_role, 1)],
                "n_genes_with_3utr_seed_match": gene_level_match_counts.get((row["sirna"], strand_role), 0),
            })
    write_tsv(
        os.path.join(out_dir, "supplemental_sirna_specificity_summary.tsv"),
        summary_rows,
        list(summary_rows[0].keys())
    )
    write_csv(
        os.path.join(out_dir, "sirna_specificity_summary.csv"),
        summary_rows,
        list(summary_rows[0].keys())
    )

    r_stdout, r_stderr = call_r_stats(project_root, out_dir)
    with open(os.path.join(out_dir, "seed_stats_r_log.txt"), "w") as handle:
        handle.write(r_stdout)
        handle.write("\n--- STDERR ---\n")
        handle.write(r_stderr)

    with open(os.path.join(out_dir, "required_input_specification.txt"), "w") as handle:
        handle.write("Required input specification for full siRNA off-target support analysis\n\n")
        handle.write("Used in this run:\n")
        handle.write(f"- Transcriptome FASTA: {transcriptome_path}\n")
        handle.write(f"- GTF annotation: {gtf_path}\n\n")
        handle.write("Optional future input:\n")
        handle.write("- Dedicated mouse 3'UTR FASTA retaining transcript identifiers.\n")
        handle.write("  This would allow direct confirmation of seed-match results without reconstructing 3'UTRs from exon/CDS structure.\n")


if __name__ == "__main__":
    main()
