#!/usr/bin/env python3
"""
Run bam_to_dada pipeline on ALL 10 BAM samples.

Combines reads from:
  - Choquet 2023 (SRR20215246, SRR20215247) — already processed individually
  - Drexler 2020 untreated K562 chromatin (SRR8268942-44, SRR10097603-05)
  - Drexler 2020 DMSO controls (SRR8932660-61)

All are K562 chromatin-associated RNA, Oxford Nanopore long-read.
More reads per gene = better IPD estimation.

Output: gene_dada_bam_all.pkl (separate from gene_dada_bam.pkl)
"""
import os
import sys

sys.path.insert(0, os.path.dirname(__file__))
from bam_to_dada import run_pipeline

base = "/grid/wsbs/home_norepl/smittal/KOO/projects/Alphagenome_splice_isoform_prediction"
bam_dir = f"{base}/kinetic_splicing/data/bam"

bam_paths = [
    # Choquet 2023 (K562 chromatin, Nanopore)
    f"{bam_dir}/SRR20215246.sorted.bam",
    f"{bam_dir}/SRR20215247.sorted.bam",
    # Drexler 2020 untreated K562 chromatin (MinION reps 1-3)
    f"{bam_dir}/SRR8268942.sorted.bam",
    f"{bam_dir}/SRR8268943.sorted.bam",
    f"{bam_dir}/SRR8268944.sorted.bam",
    # Drexler 2020 untreated K562 chromatin (PromethION reps 4-5)
    f"{bam_dir}/SRR10097603.sorted.bam",
    f"{bam_dir}/SRR10097604.sorted.bam",
    f"{bam_dir}/SRR10097605.sorted.bam",
    # Drexler 2020 DMSO vehicle controls
    f"{bam_dir}/SRR8932660.sorted.bam",
    f"{bam_dir}/SRR8932661.sorted.bam",
]

gtf_path = f"{base}/refs/gencode.v49.annotation.gtf"
output_path = f"{base}/kinetic_splicing/data/processed/gene_dada_bam_all.pkl"

os.makedirs(os.path.dirname(output_path), exist_ok=True)

print(f"Processing {len(bam_paths)} BAM files:")
for p in bam_paths:
    name = os.path.basename(p)
    size_mb = os.path.getsize(p) / 1e6
    print(f"  {name} ({size_mb:.0f} MB)")
print()

results = run_pipeline(bam_paths, gtf_path, output_path,
                       min_reads=50, min_junctions=2)
