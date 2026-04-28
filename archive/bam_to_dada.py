#!/usr/bin/env python3
"""
BAM → DADA pipeline for nano-COP data.

Steps:
  1. Parse GENCODE GTF → gene boundaries + annotated intron donors/acceptors
  2. Extract junctions from BAM CIGAR strings (N operations = spliced introns)
  3. Cluster junctions (±1bp tolerance for GT-AG alignment wobble)
  4. Assign reads to genes by TSS (5'-most aligned position on coding strand)
  5. For each gene: build observed junction set, construct DADA, compute IPD via EM
  6. Report landscape: intron counts, state counts, read counts per gene

Output: pickle with per-gene DADA data ready for rate recovery.
"""

import sys
import os
import re
import pickle
import time
import bisect
import numpy as np
import pandas as pd
from collections import defaultdict, Counter
from dataclasses import dataclass, field

import pysam

sys.path.insert(0, os.path.dirname(__file__))
# KineticDADA not needed — states are derived from read data, not enumerated


# ══════════════════════════════════════════════════════════════
#  Step 1: Parse GTF for gene annotations
# ══════════════════════════════════════════════════════════════

def parse_gtf_genes(gtf_path, gene_type="protein_coding"):
    """
    Parse GENCODE GTF → DataFrame of gene boundaries + intron coordinates.

    Returns genes_df with columns:
      gene_id, gene_name, chrom, strand, start, end, tss,
      intron_donors (list), intron_acceptors (list)
    """
    print(f"  Parsing GTF: {gtf_path}")

    genes = {}       # gene_id → {chrom, strand, start, end, gene_name}
    transcripts = {} # transcript_id → {gene_id, exons: [(start,end), ...]}

    with open(gtf_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            if len(fields) < 9:
                continue

            chrom, source, feature, start, end, _, strand, _, attrs = fields
            start, end = int(start) - 1, int(end)  # 0-based half-open

            # Parse attributes
            attr_dict = {}
            for attr in attrs.strip().split(";"):
                attr = attr.strip()
                if not attr:
                    continue
                parts = attr.split(" ", 1)
                if len(parts) == 2:
                    key = parts[0]
                    val = parts[1].strip('"')
                    attr_dict[key] = val

            if feature == "gene":
                gt = attr_dict.get("gene_type", "")
                if gene_type and gt != gene_type:
                    continue
                gid = attr_dict.get("gene_id", "")
                gname = attr_dict.get("gene_name", "")
                genes[gid] = {
                    "gene_id": gid, "gene_name": gname,
                    "chrom": chrom, "strand": strand,
                    "start": start, "end": end,
                }

            elif feature == "exon":
                tid = attr_dict.get("transcript_id", "")
                gid = attr_dict.get("gene_id", "")
                if gid not in genes:
                    continue
                if tid not in transcripts:
                    transcripts[tid] = {"gene_id": gid, "exons": []}
                transcripts[tid]["exons"].append((start, end))

    # Compute intron donors/acceptors per gene (union across transcripts)
    gene_introns = defaultdict(set)  # gene_id → set of (donor, acceptor)

    for tid, tdata in transcripts.items():
        gid = tdata["gene_id"]
        exons = sorted(tdata["exons"], key=lambda x: x[0])
        for i in range(len(exons) - 1):
            donor = exons[i][1]      # end of upstream exon = donor position
            acceptor = exons[i+1][0] # start of downstream exon = acceptor position
            if acceptor > donor:     # sanity
                gene_introns[gid].add((donor, acceptor))

    # Build DataFrame
    records = []
    for gid, ginfo in genes.items():
        introns = sorted(gene_introns.get(gid, []))
        donors = [d for d, a in introns]
        acceptors = [a for d, a in introns]

        # TSS
        if ginfo["strand"] == "+":
            tss = ginfo["start"]
        else:
            tss = ginfo["end"]

        records.append({
            **ginfo,
            "tss": tss,
            "intron_donors": donors,
            "intron_acceptors": acceptors,
            "n_annotated_introns": len(introns),
        })

    df = pd.DataFrame(records)
    print(f"  Parsed {len(df)} protein-coding genes, "
          f"{df['n_annotated_introns'].sum()} annotated introns")
    return df


# ══════════════════════════════════════════════════════════════
#  Step 2: Extract junctions from BAM
# ══════════════════════════════════════════════════════════════

@dataclass(slots=True)
class ReadJunctions:
    """One read's junction information (memory-optimized)."""
    chrom: str
    strand: str          # from read orientation
    read_start: int      # 5'-most aligned position (reference coords)
    read_end: int        # 3'-most aligned position
    junctions: tuple     # ((donor, acceptor), ...) from CIGAR N operations


def extract_junctions_from_bam(bam_path, min_intron=50, max_intron=500000):
    """
    Extract splice junctions from BAM using CIGAR N operations.

    Each N in the CIGAR string represents a spliced intron.
    Returns list of ReadJunctions.
    """
    print(f"  Extracting junctions from: {bam_path}")

    reads = []
    n_total = 0
    n_spliced = 0
    n_unspliced = 0

    bam = pysam.AlignmentFile(bam_path, "rb")

    for read in bam.fetch():
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue
        n_total += 1

        # Get reference positions of N (skip/intron) operations
        cigar = read.cigartuples
        if cigar is None:
            continue

        ref_pos = read.reference_start
        junctions = []
        for op, length in cigar:
            if op == 3:  # N = reference skip (intron)
                donor = ref_pos
                acceptor = ref_pos + length
                if min_intron <= length <= max_intron:
                    junctions.append((donor, acceptor))
                ref_pos += length
            elif op in (0, 2, 7, 8):  # M, D, =, X consume reference
                ref_pos += length
            # S, I, H, P don't consume reference

        # Strand from read: nano-COP is direct RNA, so read strand = transcript strand
        strand = "-" if read.is_reverse else "+"

        # Read boundaries (avoid get_blocks() — just use reference positions)
        read_start = read.reference_start
        read_end = read.reference_end
        if read_end is None:
            read_end = read_start  # fallback for edge cases

        rj = ReadJunctions(
            chrom=read.reference_name,
            strand=strand,
            read_start=read_start,
            read_end=read_end,
            junctions=tuple(junctions),
        )
        reads.append(rj)

        if junctions:
            n_spliced += 1
        else:
            n_unspliced += 1

    bam.close()

    print(f"  Total reads: {n_total}, spliced: {n_spliced}, "
          f"unspliced: {n_unspliced}")
    return reads


# ══════════════════════════════════════════════════════════════
#  Step 3: Cluster junctions (±1bp tolerance)
# ══════════════════════════════════════════════════════════════

def cluster_junctions(all_junctions, tolerance=1):
    """
    Cluster nearby junctions (within ±tolerance bp) to handle
    GT-AG alignment wobble.

    Input: list of (donor, acceptor) tuples
    Output: dict mapping each (donor, acceptor) → canonical (donor, acceptor)
    """
    if not all_junctions:
        return {}

    # Count occurrences
    counts = defaultdict(int)
    for j in all_junctions:
        counts[j] += 1

    # Sort by count (most frequent first = canonical representative)
    sorted_juncs = sorted(counts.keys(), key=lambda j: -counts[j])

    cluster_map = {}
    # Use a dict keyed by (donor // (tolerance+1)) for O(1) lookup
    # instead of O(n) linear scan through canonical_list
    donor_bins = defaultdict(list)  # binned_donor → list of canonical junctions

    for junc in sorted_juncs:
        d, a = junc
        # Check nearby bins for merge candidates
        d_bin = d // (tolerance + 1)
        merged = False
        for db in range(d_bin - 1, d_bin + 2):
            if merged:
                break
            for canon in donor_bins.get(db, []):
                cd, ca = canon
                if abs(d - cd) <= tolerance and abs(a - ca) <= tolerance:
                    cluster_map[junc] = canon
                    merged = True
                    break
        if not merged:
            cluster_map[junc] = junc
            donor_bins[d // (tolerance + 1)].append(junc)

    n_raw = len(counts)
    n_clustered = len(donor_bins)  # approximate
    n_canonical = sum(1 for j, c in cluster_map.items() if j == c)
    if n_raw != n_canonical:
        print(f"  Clustered {n_raw} → {n_canonical} junctions "
              f"(±{tolerance}bp)")
    sys.stdout.flush()

    return cluster_map


# ══════════════════════════════════════════════════════════════
#  Step 4: Assign reads to genes by TSS
# ══════════════════════════════════════════════════════════════

def assign_reads_to_genes(reads, genes_df):
    """
    Assign each read to a gene by TSS using bisect for O(log n) lookup.

    Returns dict: gene_name → list of ReadJunctions
    """
    print(f"  Assigning {len(reads)} reads to genes...")
    t0 = time.time()

    # Build sorted interval arrays per (chrom, strand) for bisect lookup
    # For each (chrom, strand): parallel arrays of (start-500, end+500, gene_name, tss)
    intervals = defaultdict(lambda: {"starts": [], "ends": [], "names": [], "tss": []})

    for _, row in genes_df.iterrows():
        key = (row["chrom"], row["strand"])
        intervals[key]["starts"].append(row["start"] - 500)
        intervals[key]["ends"].append(row["end"] + 500)
        intervals[key]["names"].append(row["gene_name"])
        tss = row["start"] if row["strand"] == "+" else row["end"]
        intervals[key]["tss"].append(tss)

    # Sort each group by start position
    for key in intervals:
        iv = intervals[key]
        order = sorted(range(len(iv["starts"])), key=lambda i: iv["starts"][i])
        iv["starts"] = [iv["starts"][i] for i in order]
        iv["ends"] = [iv["ends"][i] for i in order]
        iv["names"] = [iv["names"][i] for i in order]
        iv["tss"] = [iv["tss"][i] for i in order]

    gene_reads = defaultdict(list)
    n_assigned = 0
    n_ambiguous = 0
    n_unassigned = 0

    for ri, read in enumerate(reads):
        if ri % 500000 == 0 and ri > 0:
            print(f"    ... processed {ri}/{len(reads)} reads "
                  f"({time.time()-t0:.0f}s)")

        key = (read.chrom, read.strand)
        if key not in intervals:
            n_unassigned += 1
            continue

        iv = intervals[key]
        tss_pos = read.read_start if read.strand == "+" else read.read_end

        # bisect: find genes whose start <= tss_pos
        idx = bisect.bisect_right(iv["starts"], tss_pos) - 1

        # Check a window around the bisect point (genes can overlap)
        candidates = []
        for i in range(max(0, idx - 5), min(len(iv["starts"]), idx + 10)):
            if iv["starts"][i] <= tss_pos <= iv["ends"][i]:
                candidates.append(i)

        if len(candidates) == 1:
            gene_reads[iv["names"][candidates[0]]].append(read)
            n_assigned += 1
        elif len(candidates) > 1:
            best = min(candidates, key=lambda i: abs(tss_pos - iv["tss"][i]))
            gene_reads[iv["names"][best]].append(read)
            n_assigned += 1
            n_ambiguous += 1
        else:
            n_unassigned += 1

    elapsed = time.time() - t0
    print(f"  Assigned: {n_assigned} ({n_ambiguous} ambiguous), "
          f"unassigned: {n_unassigned} ({elapsed:.1f}s)")
    print(f"  Genes with reads: {len(gene_reads)}")
    sys.stdout.flush()

    return gene_reads


# ══════════════════════════════════════════════════════════════
#  Step 5: Build per-gene DADA + IPD
# ══════════════════════════════════════════════════════════════

@dataclass
class GeneDadaResult:
    gene_name: str
    chrom: str
    strand: str
    n_reads: int
    n_spliced_reads: int
    observed_junctions: list   # [(donor, acceptor), ...] genomic coords
    site_types: list           # ["D", "A", "D", "A", ...]
    junction_map: dict         # genomic (d,a) → DADA index (di, ai)
    dada_junctions: list       # [(di, ai), ...] in DADA index space
    n_states: int
    n_edges: int
    state_counts: np.ndarray
    ipd: np.ndarray
    isoform_labels: list
    dada_string: str           # e.g. "DADADA"


def build_gene_dada(gene_name, gene_reads, genes_df, cluster_map,
                    min_reads=50, min_junctions=2):
    """
    Build DADA graph for one gene from its assigned reads.

    DATA-DRIVEN state construction: instead of enumerating all 2^J possible
    junction subsets (exponential blowup), we derive states from the junction
    patterns actually observed in reads. This makes even 50-junction genes
    tractable — the state count equals the number of unique read patterns,
    not 2^50.

    1. Collect observed junctions from reads (clustered)
    2. For each read, determine its junction pattern (spliced/retained/unknown)
    3. Collect unique observed spliced-junction-sets → these ARE the states
    4. Build DAG edges: states differing by one junction
    5. EM to estimate IPD from partial reads
    """
    if isinstance(genes_df, dict):
        gene_row = genes_df[gene_name]
    else:
        gene_row = genes_df[genes_df["gene_name"] == gene_name].iloc[0]
    chrom = gene_row["chrom"]
    strand = gene_row["strand"]

    reads = gene_reads[gene_name]
    if len(reads) < min_reads:
        return None

    # ── Collect all observed junctions from reads, apply clustering ──
    junc_counts = Counter()
    for read in reads:
        for j in read.junctions:
            canonical = cluster_map.get(j, j)
            junc_counts[canonical] += 1

    observed_juncs = sorted(junc_counts.keys())

    if len(observed_juncs) < min_junctions:
        return None

    # ── Build site_types and junction index map ──
    donor_set = set(d for d, a in observed_juncs)
    acceptor_set = set(a for d, a in observed_juncs)
    all_donors = sorted(donor_set)
    all_acceptors = sorted(acceptor_set)
    all_sites = sorted(donor_set | acceptor_set)

    site_entries = []
    for pos in all_sites:
        is_d = pos in donor_set
        is_a = pos in acceptor_set
        if is_d and is_a:
            site_entries.append((pos, "A"))
            site_entries.append((pos, "D"))
        elif is_d:
            site_entries.append((pos, "D"))
        else:
            site_entries.append((pos, "A"))

    if strand == "-":
        site_entries = list(reversed(site_entries))

    site_types = [s[1] for s in site_entries]

    # Map genomic junctions to index pairs
    pos_to_idx = {}
    for idx, (pos, stype) in enumerate(site_entries):
        pos_to_idx[(pos, stype)] = idx

    junction_map = {}    # genomic (d,a) → (di, ai)
    valid_juncs = []     # list of (di, ai) index pairs
    for d, a in observed_juncs:
        di = pos_to_idx.get((d, "D"))
        ai = pos_to_idx.get((a, "A"))
        if di is not None and ai is not None:
            # Use (min, max) for consistent DADA ordering.
            # On minus strand, di > ai after reversal — that's fine,
            # just normalize so the smaller index comes first.
            lo, hi = min(di, ai), max(di, ai)
            junction_map[(d, a)] = (lo, hi)
            valid_juncs.append((lo, hi))

    if len(valid_juncs) < min_junctions:
        return None

    valid_juncs_set = set(valid_juncs)
    genomic_for_idx = {v: k for k, v in junction_map.items()}

    # ── Classify each read's junction observations ──
    read_patterns = []   # per-read: (spliced_frozenset, retained_frozenset)
    observed_states = set()
    observed_states.add(frozenset())  # null state (nothing spliced) always exists
    n_spliced_reads = 0

    for read in reads:
        read_juncs_clustered = set(
            cluster_map.get(rj, rj) for rj in read.junctions)

        spliced = set()
        retained = set()
        for dj in valid_juncs:
            gj = genomic_for_idx.get(dj)
            if gj is None:
                continue  # unknown
            gd, ga = gj
            if read.read_start <= gd and read.read_end >= ga:
                # Read spans this junction — we know its status
                if gj in read_juncs_clustered:
                    spliced.add(dj)
                else:
                    retained.add(dj)
            # else: junction outside read window → unknown

        spliced_fs = frozenset(spliced)
        retained_fs = frozenset(retained)
        read_patterns.append((spliced_fs, retained_fs))

        # The spliced set is a directly observed state
        if spliced_fs:
            observed_states.add(spliced_fs)
            n_spliced_reads += 1

    # ── Build data-driven state space ──
    # States = unique spliced-junction-sets seen in reads.
    # Also add intermediate states: if we see {A,B,C} spliced, the
    # intermediates {A}, {A,B}, {B}, {A,C}, etc. that are subsets
    # ONLY if they're also seen. We already collected all observed patterns.
    #
    # Additionally, for each pair of observed states that differ by one
    # junction, the connection is a DADA edge. States not connected to
    # anything are still valid — they just have no edges (terminal or orphan).

    # Cap state count to prevent OOM on genes with many singleton patterns.
    # If too many states, keep only those observed by ≥2 reads, plus null state.
    MAX_STATES = 500
    if len(observed_states) > MAX_STATES:
        state_read_counts = Counter()
        for spliced_fs, _ in read_patterns:
            state_read_counts[spliced_fs] += 1
        # Keep null state + states with ≥2 observations
        observed_states = {s for s in observed_states
                          if s == frozenset() or state_read_counts.get(s, 0) >= 2}
        if len(observed_states) > MAX_STATES:
            # Still too many — keep top MAX_STATES by read count
            ranked = sorted(observed_states,
                           key=lambda s: -state_read_counts.get(s, 0))
            observed_states = set(ranked[:MAX_STATES])
            observed_states.add(frozenset())  # always keep null

    states = sorted(observed_states, key=lambda s: (len(s), sorted(s)))
    state_idx = {s: i for i, s in enumerate(states)}
    n_states = len(states)

    # Build DAG edges: states differing by exactly one junction
    edges = []
    for si, src in enumerate(states):
        for junc in valid_juncs:
            if junc not in src:
                candidate = frozenset(src | {junc})
                if candidate in state_idx:
                    edges.append((si, state_idx[candidate], junc))

    # ── EM for IPD ──
    state_counts, ipd = _compute_ipd_em_datadriven(
        read_patterns, states, state_idx, valid_juncs)

    # ── Build isoform labels ──
    isoform_labels = []
    for s in states:
        if not s:
            isoform_labels.append("()")
        else:
            parts = [f"({di},{ai})" for di, ai in sorted(s)]
            isoform_labels.append(",".join(parts))

    dada_string = "".join(site_types)

    return GeneDadaResult(
        gene_name=gene_name,
        chrom=chrom,
        strand=strand,
        n_reads=len(reads),
        n_spliced_reads=n_spliced_reads,
        observed_junctions=valid_juncs,
        site_types=site_types,
        junction_map=junction_map,
        dada_junctions=valid_juncs,
        n_states=n_states,
        n_edges=len(edges),
        state_counts=state_counts,
        ipd=ipd,
        isoform_labels=isoform_labels,
        dada_string=dada_string,
    )


def _compute_ipd_em_datadriven(read_patterns, states, state_idx, all_junctions,
                                max_iter=100, tol=1e-8):
    """
    EM for IPD from read observations against data-driven states.

    Each read has (spliced_frozenset, retained_frozenset). A state is
    compatible with a read if:
      - all spliced junctions are IN the state
      - no retained junctions are IN the state
    """
    n_states = len(states)

    # Precompute compatible states per read
    read_compat = []
    for spliced, retained in read_patterns:
        compatible = []
        for idx, state in enumerate(states):
            if not spliced.issubset(state):
                continue
            if retained & state:
                continue
            compatible.append(idx)
        read_compat.append(compatible)

    ipd = np.ones(n_states) / n_states
    for iteration in range(max_iter):
        state_counts = np.zeros(n_states)
        for compatible in read_compat:
            if not compatible:
                continue
            if len(compatible) == 1:
                state_counts[compatible[0]] += 1.0
            else:
                weights = ipd[compatible]
                total = weights.sum()
                if total > 0:
                    weights = weights / total
                else:
                    weights = np.ones(len(compatible)) / len(compatible)
                for i, idx in enumerate(compatible):
                    state_counts[idx] += weights[i]

        total = state_counts.sum()
        new_ipd = state_counts / total if total > 0 else ipd
        delta = np.abs(new_ipd - ipd).max()
        ipd = new_ipd
        if delta < tol:
            break

    return state_counts, ipd


# ══════════════════════════════════════════════════════════════
#  Step 6: Run full pipeline
# ══════════════════════════════════════════════════════════════

def run_pipeline(bam_paths, gtf_path, output_path,
                 min_reads=50, min_junctions=2):
    """Full pipeline: BAMs + GTF → per-gene DADA results."""
    # Force line-buffered stdout so SLURM .out files update in real time
    sys.stdout.reconfigure(line_buffering=True)
    sys.stderr.reconfigure(line_buffering=True)

    # 1. Parse GTF
    genes_df = parse_gtf_genes(gtf_path)

    # 2. Extract junctions from all BAMs
    all_reads = []
    for bam_path in bam_paths:
        reads = extract_junctions_from_bam(bam_path)
        all_reads.extend(reads)
    print(f"\n  Total reads across BAMs: {len(all_reads)}")
    sys.stdout.flush()

    # 3. Cluster junctions
    print(f"\n  Collecting raw junctions...")
    sys.stdout.flush()
    all_raw_juncs = []
    for read in all_reads:
        all_raw_juncs.extend(read.junctions)
    print(f"  {len(all_raw_juncs)} raw junction instances")
    sys.stdout.flush()
    cluster_map = cluster_junctions(all_raw_juncs, tolerance=1)
    print(f"  Clustering done")
    sys.stdout.flush()

    # 4. Assign reads to genes
    gene_reads = assign_reads_to_genes(all_reads, genes_df)

    # Pre-build gene lookup dict to avoid repeated DataFrame filtering
    # Use first occurrence per gene_name to avoid cross-chromosome contamination
    gene_lookup = {}
    for _, row in genes_df.iterrows():
        if row["gene_name"] not in gene_lookup:
            gene_lookup[row["gene_name"]] = row

    # 5. Build per-gene DADA + IPD
    eligible_genes = [g for g in sorted(gene_reads.keys())
                      if len(gene_reads[g]) >= min_reads]
    print(f"\n  Building DADA graphs for {len(eligible_genes)} eligible genes "
          f"(of {len(gene_reads)} with reads)...")
    t0 = time.time()
    results = []
    n_success = 0
    n_fail = 0

    for gi, gene_name in enumerate(eligible_genes):
        if gi % 200 == 0 and gi > 0:
            print(f"    ... {gi}/{len(eligible_genes)} genes, "
                  f"{n_success} built, {n_fail} filtered ({time.time()-t0:.0f}s)")
            sys.stdout.flush()

        result = build_gene_dada(gene_name, gene_reads, gene_lookup,
                                 cluster_map, min_reads, min_junctions)
        if result is not None:
            results.append(result)
            n_success += 1
        else:
            n_fail += 1

    elapsed = time.time() - t0
    print(f"\n  Built {n_success} DADA graphs ({n_fail} failed filters) "
          f"in {elapsed:.1f}s")

    # 6. Report landscape
    print(f"\n{'=' * 70}")
    print(f"  Data Landscape")
    print(f"{'=' * 70}")

    if results:
        dada_strings = defaultdict(list)
        for r in results:
            dada_strings[r.dada_string].append(r)

        print(f"\n  {'DADA':>12s} {'Genes':>6s} {'States':>7s} "
              f"{'Edges':>6s} {'MedReads':>9s} {'MedSpliced':>11s}")
        print(f"  {'-'*60}")

        for ds in sorted(dada_strings.keys(), key=lambda s: (len(s), s)):
            genes = dada_strings[ds]
            states = [g.n_states for g in genes]
            edges = [g.n_edges for g in genes]
            nreads = [g.n_reads for g in genes]
            nspliced = [g.n_spliced_reads for g in genes]
            print(f"  {ds:>12s} {len(genes):>6d} {np.median(states):>7.0f} "
                  f"{np.median(edges):>6.0f} {np.median(nreads):>9.0f} "
                  f"{np.median(nspliced):>11.0f}")

        # Overall stats
        all_nreads = [r.n_reads for r in results]
        all_states = [r.n_states for r in results]
        all_juncs = [len(r.observed_junctions) for r in results]

        print(f"\n  Overall: {len(results)} genes")
        print(f"  Reads/gene:    median={np.median(all_nreads):.0f}, "
              f"max={np.max(all_nreads):.0f}")
        print(f"  Junctions/gene: median={np.median(all_juncs):.0f}, "
              f"max={np.max(all_juncs):.0f}")
        print(f"  States/gene:   median={np.median(all_states):.0f}, "
              f"max={np.max(all_states):.0f}")

    # Save
    with open(output_path, "wb") as f:
        pickle.dump(results, f)
    print(f"\n  Saved to {output_path}")

    return results


# ══════════════════════════════════════════════════════════════
#  Main
# ══════════════════════════════════════════════════════════════

if __name__ == "__main__":
    base = "/grid/wsbs/home_norepl/smittal/KOO/projects/Alphagenome_splice_isoform_prediction"

    bam_paths = [
        f"{base}/kinetic_splicing/data/bam/SRR20215246.sorted.bam",
        f"{base}/kinetic_splicing/data/bam/SRR20215247.sorted.bam",
    ]
    gtf_path = f"{base}/refs/gencode.v49.annotation.gtf"
    output_path = f"{base}/kinetic_splicing/data/processed/gene_dada_bam.pkl"

    os.makedirs(os.path.dirname(output_path), exist_ok=True)

    results = run_pipeline(bam_paths, gtf_path, output_path,
                           min_reads=50, min_junctions=2)
