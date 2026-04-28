#!/usr/bin/env python3
"""
TRUE DADA Graphs from Churchman Data
=====================================
Validates predetermined splicing patterns using nano-COP data
analyzed with DADA graphs.

Key insight: predetermined splicing = context-dependent (CD) rates.
If intron 1 must splice before intron 2, the rate of junction (2,3)
depends on whether (0,1) is already spliced. The DADA graph encodes
all possible orders; the observed read distribution reveals which
splice programs are actually played out.

Output: kinetic_splicing/results/true_dada_churchman/
"""
import sys
import os
import re
import csv
import json
import time
import pickle
import numpy as np
from itertools import combinations
from dataclasses import dataclass
from collections import defaultdict
from types import ModuleType
from pathlib import Path

from radial_dada_html import write_radial_html

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

sys.stdout.reconfigure(line_buffering=True)
sys.stderr.reconfigure(line_buffering=True)

# ══════════════════════════════════════════════════════════════
#  Config
# ══════════════════════════════════════════════════════════════
BASE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(BASE, "results", "true_dada_churchman")
os.makedirs(OUT, exist_ok=True)

DATA_PKL = os.path.join(BASE, "data", "processed", "gene_dada_bam_all.pkl")
CI_CSV = os.path.join(BASE, "results", "ci_export_bam", "ci_test_bam.csv")

plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.size': 10,
    'axes.spines.top': False,
    'axes.spines.right': False,
    'figure.dpi': 150,
    'savefig.bbox': 'tight',
    'savefig.pad_inches': 0.15,
})

# ══════════════════════════════════════════════════════════════
#  Section A: Data Loading
# ══════════════════════════════════════════════════════════════

@dataclass
class GeneDadaResult:
    gene_name: str; chrom: str; strand: str; n_reads: int; n_spliced_reads: int
    observed_junctions: list; site_types: list; junction_map: dict
    dada_junctions: list; n_states: int; n_edges: int
    state_counts: object; ipd: object; isoform_labels: list; dada_string: str

_fake = ModuleType('bam_to_dada')
_fake.GeneDadaResult = GeneDadaResult
sys.modules['bam_to_dada'] = _fake


def load_genes(pkl_path):
    """Load GeneDadaResult list from pickle."""
    print(f"  Loading {pkl_path}...")
    with open(pkl_path, 'rb') as f:
        genes = pickle.load(f)
    print(f"  Loaded {len(genes)} genes")
    return genes


def load_ci_metrics(csv_path):
    """Load CI fit metrics from CSV. Returns dict gene_name → row dict."""
    if not os.path.exists(csv_path):
        print(f"  WARNING: CI metrics not found at {csv_path}")
        return {}
    ci = {}
    with open(csv_path, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            ci[row['gene']] = row
    print(f"  Loaded CI metrics for {len(ci)} genes")
    return ci


# ══════════════════════════════════════════════════════════════
#  Section B: DADA Reconstruction + Flow Estimation
# ══════════════════════════════════════════════════════════════

def parse_label(label):
    """Parse isoform label like '(0,1),(2,3)' → frozenset of junction tuples."""
    if label == "()":
        return frozenset()
    pairs = re.findall(r'\((\d+),(\d+)\)', label)
    return frozenset((int(a), int(b)) for a, b in pairs)


def reconstruct_dada(gene):
    """Reconstruct states, edges, junction list from GeneDadaResult."""
    states = [parse_label(lab) for lab in gene.isoform_labels]
    state_idx = {s: i for i, s in enumerate(states)}

    all_juncs = sorted(set().union(*states)) if any(states) else []
    junc_idx = {j: i for i, j in enumerate(all_juncs)}

    edges = []
    for si, src in enumerate(states):
        for j in all_juncs:
            if j not in src:
                candidate = frozenset(src | {j})
                if candidate in state_idx:
                    di = state_idx[candidate]
                    edges.append((si, di, j, junc_idx[j]))

    layers = [len(s) for s in states]
    return states, edges, all_juncs, junc_idx, layers


def estimate_edge_flows(edges, layers, observed_ipd):
    """
    Estimate flow on each edge from observed node occupancies.
    Top-down layer-by-layer allocation using SUBTREE DEMAND: at each node,
    distribute outgoing flow proportionally to the total observed IPD in
    each child's subtree (not just the direct child). This prevents
    high-flow grandchildren from being starved when intermediate nodes
    have zero reads.
    """
    n_states = len(layers)
    children_of = defaultdict(list)
    for k, (src, dst, j, ji) in enumerate(edges):
        children_of[src].append((k, dst))

    # Compute subtree demand: sum of observed_ipd for node + all descendants
    subtree_demand = np.array(observed_ipd, dtype=float)
    max_layer = max(layers) if layers else 0
    # Bottom-up pass: accumulate child demands into parents
    for layer in range(max_layer, -1, -1):
        for si in [i for i, l in enumerate(layers) if l == layer]:
            for _, ci in children_of[si]:
                subtree_demand[si] += subtree_demand[ci]

    edge_flows = np.zeros(len(edges))
    node_inflow = np.zeros(n_states)
    node_inflow[0] = 1.0  # production into null state

    for layer in range(max_layer + 1):
        for si in [i for i, l in enumerate(layers) if l == layer]:
            total_in = node_inflow[si]
            child_list = children_of[si]
            if not child_list:
                continue

            # Use subtree demand for allocation weights
            child_weights = np.array([max(subtree_demand[ci], 1e-30)
                                      for _, ci in child_list])
            total_w = child_weights.sum()
            if total_w < 1e-30:
                continue

            for (ek, ci), w in zip(child_list, child_weights):
                edge_flows[ek] = total_in * (w / total_w)
                node_inflow[ci] += edge_flows[ek]

    return edge_flows


# ══════════════════════════════════════════════════════════════
#  Section C: Predeterminedness Metrics
# ══════════════════════════════════════════════════════════════

def compute_pairwise_ordering(states, layers, observed_ipd):
    """
    Pairwise ordering consistency: for each junction pair, which is spliced first?
    Returns (mean_consistency, ordering_matrix, junction_list).
    """
    all_juncs = sorted(set().union(*states)) if any(s for s in states) else []
    if len(all_juncs) < 2:
        return 1.0, np.ones((1, 1)), all_juncs

    n_juncs = len(all_juncs)
    ordering = np.full((n_juncs, n_juncs), 0.5)  # P(j_i before j_j)

    for i, j1 in enumerate(all_juncs):
        for j, j2 in enumerate(all_juncs):
            if i >= j:
                continue
            # Find states containing j1 but not j2, and vice versa
            ipd_j1_only = 0.0
            ipd_j2_only = 0.0
            for si, s in enumerate(states):
                if j1 in s and j2 not in s:
                    ipd_j1_only += observed_ipd[si]
                elif j2 in s and j1 not in s:
                    ipd_j2_only += observed_ipd[si]

            total = ipd_j1_only + ipd_j2_only
            if total > 1e-12:
                p_j1_first = ipd_j1_only / total
                ordering[i, j] = p_j1_first
                ordering[j, i] = 1.0 - p_j1_first

    # Mean consistency = average max(P, 1-P) over all pairs
    consistencies = []
    for i in range(n_juncs):
        for j in range(i + 1, n_juncs):
            consistencies.append(max(ordering[i, j], ordering[j, i]))

    mean_consistency = np.mean(consistencies) if consistencies else 1.0
    return mean_consistency, ordering, all_juncs


def enumerate_paths(edges, layers, n_states, max_paths=10000):
    """
    Enumerate root→terminal paths in the observed DAG.
    Returns list of paths, each a list of state indices.
    """
    children_of = defaultdict(list)
    for k, (src, dst, j, ji) in enumerate(edges):
        children_of[src].append(dst)

    # Terminal states = no children
    terminal = set(range(n_states)) - set(children_of.keys())
    if not terminal:
        terminal = {i for i, l in enumerate(layers) if l == max(layers)}

    paths = []
    stack = [(0, [0])]  # (current_node, path_so_far)
    while stack and len(paths) < max_paths:
        node, path = stack.pop()
        if node in terminal:
            paths.append(path)
        else:
            for child in children_of[node]:
                stack.append((child, path + [child]))
    return paths


def compute_path_metrics(edges, layers, observed_ipd, edge_flows):
    """
    Compute path-based metrics: dominant path fraction and path entropy.
    """
    n_states = len(layers)
    paths = enumerate_paths(edges, layers, n_states)
    if not paths:
        return 1.0, 0.0, []

    # Build edge lookup for flows
    edge_flow_map = {}
    for k, (src, dst, j, ji) in enumerate(edges):
        edge_flow_map[(src, dst)] = edge_flows[k]

    # Compute branching probabilities at each node
    children_of = defaultdict(list)
    for k, (src, dst, j, ji) in enumerate(edges):
        children_of[src].append((dst, edge_flows[k]))

    node_total_outflow = {}
    for src, children in children_of.items():
        node_total_outflow[src] = sum(f for _, f in children)

    # Compute path probabilities
    path_probs = []
    for path in paths:
        prob = 1.0
        for i in range(len(path) - 1):
            src, dst = path[i], path[i + 1]
            total_out = node_total_outflow.get(src, 0)
            edge_f = edge_flow_map.get((src, dst), 0)
            if total_out > 1e-30:
                prob *= edge_f / total_out
            else:
                prob = 0
                break
        path_probs.append(prob)

    # Normalize
    total_prob = sum(path_probs)
    if total_prob > 1e-30:
        path_probs = [p / total_prob for p in path_probs]

    # Metrics
    dominant_fraction = max(path_probs) if path_probs else 1.0

    # Path entropy (normalized by log of number of paths)
    entropy = 0.0
    for p in path_probs:
        if p > 1e-30:
            entropy -= p * np.log(p)
    max_entropy = np.log(len(paths)) if len(paths) > 1 else 1.0
    norm_entropy = entropy / max_entropy if max_entropy > 0 else 0.0

    # Find the dominant path
    best_idx = np.argmax(path_probs)
    dominant_path = paths[best_idx]

    return dominant_fraction, norm_entropy, dominant_path


def analyze_gene(gene, ci_row=None):
    """Full analysis of a single gene. Returns dict of metrics or None."""
    try:
        states, edges, all_juncs, junc_idx, layers = reconstruct_dada(gene)
    except Exception:
        return None

    n_states = len(states)
    if n_states < 3 or not edges:
        return None

    observed_ipd = np.array(gene.state_counts, dtype=float)
    total = observed_ipd.sum()
    if total < 1:
        return None
    observed_ipd = observed_ipd / total

    n_juncs = len(all_juncs)
    edge_flows = estimate_edge_flows(edges, layers, observed_ipd)

    # Pairwise ordering
    mean_consistency, ordering_matrix, ord_juncs = compute_pairwise_ordering(
        states, layers, observed_ipd)

    # Path-based (for genes with manageable complexity)
    dominant_fraction = np.nan
    norm_entropy = np.nan
    dominant_path = []
    if n_juncs <= 10:
        dominant_fraction, norm_entropy, dominant_path = compute_path_metrics(
            edges, layers, observed_ipd, edge_flows)

    # CI fit from CSV
    ci_r2 = np.nan
    ci_kl = np.nan
    if ci_row:
        try:
            ci_r2 = float(ci_row.get('ci_r2_log', 'nan'))
            ci_kl = float(ci_row.get('ci_kl', 'nan'))
        except (ValueError, TypeError):
            pass

    return {
        'gene_name': gene.gene_name,
        'chrom': gene.chrom,
        'strand': gene.strand,
        'n_reads': gene.n_reads,
        'n_spliced_reads': gene.n_spliced_reads,
        'n_junctions': n_juncs,
        'n_states': n_states,
        'n_edges': len(edges),
        'max_layer': max(layers) if layers else 0,
        'dada_string': gene.dada_string,
        'ordering_consistency': mean_consistency,
        'dominant_path_fraction': dominant_fraction,
        'path_entropy': norm_entropy,
        'ci_r2_log': ci_r2,
        'ci_kl': ci_kl,
        # For showcase genes — not written to CSV
        '_states': states,
        '_edges': edges,
        '_all_juncs': all_juncs,
        '_junc_idx': junc_idx,
        '_layers': layers,
        '_observed_ipd': observed_ipd,
        '_edge_flows': edge_flows,
        '_ordering_matrix': ordering_matrix,
        '_dominant_path': dominant_path,
        '_isoform_labels': gene.isoform_labels,
        '_state_counts': np.array(gene.state_counts),
        '_site_types': gene.site_types if hasattr(gene, 'site_types') else list(gene.dada_string),
    }


# ══════════════════════════════════════════════════════════════
#  Section D: Gene Selection
# ══════════════════════════════════════════════════════════════

def select_genes(results, n_per_category=3):
    """
    Select showcase genes:
    - Top N by ordering consistency (most predetermined / "tree")
    - Bottom N by ordering consistency (least predetermined / "bushy")
    - A few from the middle
    """
    # Filter: 3-8 junctions, >50 reads, >5 states
    filtered = [r for r in results
                if r is not None
                and 3 <= r['n_junctions'] <= 8
                and r['n_spliced_reads'] >= 50
                and r['n_states'] >= 5
                and not np.isnan(r['ordering_consistency'])]

    if not filtered:
        print("  WARNING: No genes pass showcase filters")
        return []

    sorted_by_oc = sorted(filtered, key=lambda x: x['ordering_consistency'],
                          reverse=True)

    showcase = []
    seen = set()

    # Most predetermined (tree-like)
    for r in sorted_by_oc:
        if r['gene_name'] not in seen:
            showcase.append(r)
            seen.add(r['gene_name'])
        if len(showcase) >= n_per_category:
            break

    # Least predetermined (bushy)
    for r in reversed(sorted_by_oc):
        if r['gene_name'] not in seen:
            showcase.append(r)
            seen.add(r['gene_name'])
        if len(showcase) >= 2 * n_per_category:
            break

    # Middle range
    mid = len(sorted_by_oc) // 2
    for r in sorted_by_oc[mid:mid+20]:
        if r['gene_name'] not in seen:
            showcase.append(r)
            seen.add(r['gene_name'])
        if len(showcase) >= 3 * n_per_category:
            break

    print(f"  Selected {len(showcase)} showcase genes")
    for r in showcase:
        tag = "TREE" if r['ordering_consistency'] > 0.85 else (
              "BUSHY" if r['ordering_consistency'] < 0.65 else "MID")
        print(f"    [{tag}] {r['gene_name']:15s}  OC={r['ordering_consistency']:.3f}  "
              f"J={r['n_junctions']}  S={r['n_states']}  "
              f"reads={r['n_spliced_reads']}")

    return showcase


# ══════════════════════════════════════════════════════════════
#  Section E: Interactive HTML Gallery
# ══════════════════════════════════════════════════════════════

def build_html_data(result):
    """Convert an analysis result dict to the JSON-serializable format for HTML."""
    layers_dict = defaultdict(list)
    for i, l in enumerate(result['_layers']):
        layers_dict[str(l)].append(i)

    return {
        'geneName': result['gene_name'],
        'chrom': result['chrom'],
        'strand': result['strand'],
        'nReads': result['n_reads'],
        'nSplicedReads': result['n_spliced_reads'],
        'dadaString': result['dada_string'],
        'nStates': result['n_states'],
        'nEdges': result['n_edges'],
        'nJunctions': result['n_junctions'],
        'maxLayer': result['max_layer'],
        'labels': result['_isoform_labels'],
        'observedIPD': result['_observed_ipd'].tolist(),
        'stateCounts': result['_state_counts'].tolist(),
        'edges': [[s, d, [j[0], j[1]]] for s, d, j, ji in result['_edges']],
        'edgeFlows': result['_edge_flows'].tolist(),
        'dominantPath': result['_dominant_path'],
        'layers': dict(layers_dict),
        'orderingConsistency': round(result['ordering_consistency'], 4),
        'dominantPathFraction': round(result['dominant_path_fraction'], 4)
            if not np.isnan(result['dominant_path_fraction']) else None,
        'pathEntropy': round(result['path_entropy'], 4)
            if not np.isnan(result['path_entropy']) else None,
        'ciR2': round(result['ci_r2_log'], 4)
            if not np.isnan(result['ci_r2_log']) else None,
        'orderingMatrix': result['_ordering_matrix'].tolist(),
        'junctionLabels': [f"({j[0]},{j[1]})" for j in result['_all_juncs']],
    }


def write_gallery_html(showcase_results, output_path, title="TRUE DADA: Churchman Data"):
    """Write multi-gene interactive HTML gallery."""
    genes_data = [build_html_data(r) for r in showcase_results]
    model_json = json.dumps(genes_data)

    html = _build_gallery_html(model_json, title)
    Path(output_path).write_text(html)
    print(f"  Saved {output_path}")


def _build_gallery_html(genes_json, title):
    return f'''<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>{title}</title>
<style>
*, *::before, *::after {{ box-sizing: border-box; margin: 0; padding: 0; }}
:root {{
  --pri: #335E95; --sec: #ADC9E9; --hi: #F1C045;
  --bg: #F4F7FB; --white: #FFFFFF; --text: #1e2d40;
  --tree: #2d8a4e; --bushy: #c0392b;
}}
body {{ font-family: 'Segoe UI', system-ui, sans-serif; background: var(--bg); color: var(--text); font-size: 14px; }}
.header {{
  background: var(--pri); color: var(--white);
  padding: 14px 28px; display: flex; align-items: center; gap: 14px;
  box-shadow: 0 2px 6px rgba(0,0,0,0.18);
}}
.header h1 {{ font-size: 1.25em; font-weight: 700; }}
.header .sub {{ font-size: 0.8em; opacity: 0.75; margin-top: 3px; font-style: italic; }}
.header select {{
  padding: 6px 12px; border-radius: 5px; border: none;
  font-size: 0.9em; background: rgba(255,255,255,0.9); color: var(--text);
  cursor: pointer; margin-left: auto;
}}
.main-grid {{
  display: grid; grid-template-columns: 240px 1fr 280px;
  gap: 12px; padding: 12px 14px; min-height: calc(100vh - 60px);
}}
@media (max-width: 900px) {{ .main-grid {{ grid-template-columns: 1fr; }} }}
.card {{
  background: var(--white); border-radius: 8px;
  box-shadow: 0 1px 5px rgba(51,94,149,0.10); padding: 14px;
}}
.card-title {{
  font-size: 0.78em; font-weight: 700; text-transform: uppercase;
  letter-spacing: 0.08em; color: var(--pri);
  border-bottom: 2px solid var(--sec); padding-bottom: 6px; margin-bottom: 12px;
}}
.info-row {{
  display: flex; justify-content: space-between; font-size: 0.82em;
  margin-bottom: 4px; padding: 2px 0;
}}
.info-row .label {{ color: #666; }}
.info-row .value {{ font-weight: 600; font-family: monospace; color: var(--pri); }}
.metric-row {{
  display: flex; justify-content: space-between; font-size: 0.82em;
  margin-bottom: 6px; padding: 4px 6px; background: #f0f5fc; border-radius: 4px;
}}
.metric-row .val {{ font-weight: 700; font-family: monospace; }}
.metric-row .val.high {{ color: var(--tree); }}
.metric-row .val.low {{ color: var(--bushy); }}
.bar-block {{ margin-bottom: 5px; }}
.bar-name {{
  font-size: 0.72em; font-weight: 600; color: var(--text);
  display: flex; justify-content: space-between; margin-bottom: 1px;
}}
.bar-name .pct {{ color: var(--pri); font-family: monospace; }}
.bar-track {{ height: 12px; background: #edf2f8; border-radius: 4px; overflow: hidden; }}
.bar-fill {{ height: 100%; border-radius: 4px; transition: width 0.2s ease; }}
.legend {{
  display: flex; gap: 14px; font-size: 0.73em; color: #666;
  margin-top: 8px; padding-top: 6px; border-top: 1px solid #eee;
}}
.leg-dot {{
  width: 10px; height: 10px; border-radius: 2px;
  display: inline-block; margin-right: 4px; vertical-align: middle;
}}
svg {{ width: 100%; height: 100%; }}
#tooltip {{
  position: fixed; pointer-events: none; padding: 8px 12px;
  background: rgba(30,45,64,0.92); color: #fff; border-radius: 6px;
  font-size: 12px; font-family: monospace; z-index: 999;
  display: none; max-width: 260px; line-height: 1.5;
}}
.ordering-grid {{
  display: grid; gap: 1px; margin-top: 8px; font-size: 0.7em;
  font-family: monospace; text-align: center;
}}
.ordering-cell {{
  padding: 3px 2px; border-radius: 2px; color: #fff; font-weight: 600;
}}
.ordering-label {{ color: #666; padding: 3px 2px; font-weight: 700; }}
</style>
</head>
<body>

<div class="header">
  <div>
    <h1>{title}</h1>
    <div class="sub">Validating predetermined splicing patterns from nano-COP chromatin RNA</div>
  </div>
  <select id="gene-select"></select>
</div>

<div class="main-grid">
  <!-- Left: Gene info + Metrics -->
  <div>
    <div class="card" id="info-card">
      <div class="card-title">Gene Info</div>
      <div id="gene-info"></div>
    </div>
    <div class="card" style="margin-top:12px" id="metrics-card">
      <div class="card-title">Predeterminedness</div>
      <div id="metrics-panel"></div>
    </div>
    <div class="card" style="margin-top:12px" id="ordering-card">
      <div class="card-title">Junction Ordering</div>
      <div id="ordering-panel"></div>
    </div>
  </div>

  <!-- Center: DADA graph -->
  <div class="card" style="padding:8px; display:flex; flex-direction:column;">
    <div class="card-title">DADA Graph — Observed Splice Programs</div>
    <svg id="state-graph" style="flex:1; min-height:400px;">
      <defs>
        <marker id="arr" markerWidth="8" markerHeight="6" refX="7" refY="3"
                orient="auto" markerUnits="strokeWidth">
          <polygon points="0 0, 8 3, 0 6" fill="#335E95" opacity="0.7"/>
        </marker>
        <marker id="arr-dom" markerWidth="8" markerHeight="6" refX="7" refY="3"
                orient="auto" markerUnits="strokeWidth">
          <polygon points="0 0, 8 3, 0 6" fill="#F1C045" opacity="0.9"/>
        </marker>
      </defs>
      <g id="graph-arrows"></g>
      <g id="graph-nodes"></g>
    </svg>
  </div>

  <!-- Right: IPD bars -->
  <div class="card" id="ipd-card" style="overflow-y:auto; max-height:calc(100vh - 90px);">
    <div class="card-title">Isoform Probability (IPD)</div>
    <div id="ipd-bars"></div>
    <div class="legend" style="margin-top:12px;">
      <div><span class="leg-dot" style="background:#335E95;"></span>Observed IPD</div>
    </div>
  </div>
</div>

<div id="tooltip"></div>

<script>
const GENES = ''' + genes_json + ''';

const C = {
  pri: '#335E95', sec: '#ADC9E9', hi: '#F1C045',
  tree: '#2d8a4e', bushy: '#c0392b',
};

// Viridis-like layer colors
const LAYER_COLORS = ['#8b949e','#440154','#31688e','#35b779','#fde725','#d4aa00','#b07800','#8b0000'];
function layerColor(l) {
  return LAYER_COLORS[Math.min(l, LAYER_COLORS.length-1)];
}

function viridis(t) {
  // Simplified viridis colormap
  t = Math.max(0, Math.min(1, t));
  const r = Math.round(68 + t * (253 - 68));
  const g = Math.round(1 + t * (231 - 1));
  const b = Math.round(84 + t * (37 - 84));
  return `rgb(${r},${g},${b})`;
}

let G; // current gene data

// ── Populate gene selector ──
const sel = document.getElementById('gene-select');
GENES.forEach((g, i) => {
  const opt = document.createElement('option');
  const oc = g.orderingConsistency;
  const tag = oc > 0.85 ? 'TREE' : oc < 0.65 ? 'BUSHY' : 'MID';
  opt.value = i;
  opt.textContent = `${g.geneName} [${tag}] — OC=${oc.toFixed(3)}, J=${g.nJunctions}`;
  sel.appendChild(opt);
});
sel.addEventListener('change', () => loadGene(parseInt(sel.value)));

function loadGene(idx) {
  G = GENES[idx];
  drawInfo();
  drawMetrics();
  drawOrdering();
  drawGraph();
  drawIPD();
}

// ── Info panel ──
function drawInfo() {
  const el = document.getElementById('gene-info');
  el.innerHTML = [
    ['Gene', G.geneName],
    ['Location', `${G.chrom} (${G.strand})`],
    ['DADA sites', G.dadaString],
    ['Junctions', G.nJunctions],
    ['States', G.nStates],
    ['Edges', G.nEdges],
    ['Reads', G.nReads.toLocaleString()],
    ['Spliced reads', G.nSplicedReads.toLocaleString()],
    ['Max layer', G.maxLayer],
  ].map(([l,v]) => `<div class="info-row"><span class="label">${l}</span><span class="value">${v}</span></div>`).join('');
}

// ── Metrics panel ──
function drawMetrics() {
  const el = document.getElementById('metrics-panel');
  const oc = G.orderingConsistency;
  const dp = G.dominantPathFraction;
  const pe = G.pathEntropy;
  const ci = G.ciR2;

  function cls(v, highIsGood, thresh) {
    if (v === null) return '';
    return highIsGood ? (v > thresh ? 'high' : v < (1-thresh+0.5) ? 'low' : '') :
                        (v < (1-thresh) ? 'high' : v > thresh ? 'low' : '');
  }

  el.innerHTML = `
    <div class="metric-row"><span>Ordering consistency</span>
      <span class="val ${oc > 0.8 ? 'high' : oc < 0.6 ? 'low' : ''}">${oc.toFixed(3)}</span></div>
    <div class="metric-row"><span>Dominant path fraction</span>
      <span class="val ${dp !== null && dp > 0.7 ? 'high' : dp !== null && dp < 0.3 ? 'low' : ''}">${dp !== null ? dp.toFixed(3) : 'N/A'}</span></div>
    <div class="metric-row"><span>Path entropy (norm)</span>
      <span class="val ${pe !== null && pe < 0.3 ? 'high' : pe !== null && pe > 0.7 ? 'low' : ''}">${pe !== null ? pe.toFixed(3) : 'N/A'}</span></div>
    <div class="metric-row"><span>CI fit (R² log)</span>
      <span class="val">${ci !== null ? ci.toFixed(3) : 'N/A'}</span></div>
    <div style="margin-top:10px; font-size:0.75em; color:#666; line-height:1.5;">
      <b style="color:var(--tree)">Green</b> = strongly predetermined<br>
      <b style="color:var(--bushy)">Red</b> = diverse / alternative splicing
    </div>
  `;
}

// ── Ordering matrix ──
function drawOrdering() {
  const el = document.getElementById('ordering-panel');
  const mat = G.orderingMatrix;
  const labs = G.junctionLabels;
  const n = labs.length;
  if (n < 2 || n > 12) {
    el.innerHTML = '<div style="font-size:0.8em;color:#888;">Too few/many junctions for matrix</div>';
    return;
  }
  const sz = Math.max(28, Math.min(40, 240 / (n + 1)));
  let html = `<div class="ordering-grid" style="grid-template-columns:${sz}px repeat(${n}, ${sz}px);">`;
  html += `<div></div>`;
  for (const l of labs) html += `<div class="ordering-label">${l}</div>`;
  for (let i = 0; i < n; i++) {
    html += `<div class="ordering-label">${labs[i]}</div>`;
    for (let j = 0; j < n; j++) {
      if (i === j) {
        html += `<div style="background:#eee;border-radius:2px;">—</div>`;
      } else {
        const v = mat[i][j];
        const intensity = Math.abs(v - 0.5) * 2;
        const bg = v > 0.5
          ? `rgba(45,138,78,${0.2 + intensity * 0.7})`
          : `rgba(192,57,43,${0.2 + intensity * 0.7})`;
        html += `<div class="ordering-cell" style="background:${bg}">${v.toFixed(2)}</div>`;
      }
    }
  }
  html += '</div>';
  html += '<div style="font-size:0.68em;color:#888;margin-top:6px;">P(row junction spliced before column junction)</div>';
  el.innerHTML = html;
}

// ── SVG graph ──
const NS = 'http://www.w3.org/2000/svg';
function mkEl(tag, attrs) {
  const el = document.createElementNS(NS, tag);
  for (const [k,v] of Object.entries(attrs)) el.setAttribute(k, v);
  return el;
}

function drawGraph() {
  const svg = document.getElementById('state-graph');
  const rect = svg.getBoundingClientRect();
  const W = Math.max(rect.width, 400);
  const H = Math.max(rect.height, 350);
  svg.setAttribute('viewBox', `0 0 ${W} ${H}`);

  // Layout
  const pad = 50;
  const nLayers = G.maxLayer + 1;
  const layerY = l => pad + l * ((H - 2*pad) / Math.max(nLayers-1, 1));
  const positions = new Array(G.nStates);
  for (let l = 0; l <= G.maxLayer; l++) {
    const nodes = G.layers[String(l)] || [];
    const n = nodes.length;
    for (let xi = 0; xi < n; xi++) {
      const x = n === 1 ? W/2 : pad + xi * ((W - 2*pad) / (n-1));
      positions[nodes[xi]] = { x, y: layerY(l) };
    }
  }

  // Dominant path set for highlighting
  const domSet = new Set();
  if (G.dominantPath && G.dominantPath.length > 1) {
    for (let i = 0; i < G.dominantPath.length - 1; i++) {
      domSet.add(G.dominantPath[i] + ',' + G.dominantPath[i+1]);
    }
  }

  // Radii from state counts
  const maxCount = Math.max(...G.stateCounts, 1);
  const minR = 12, maxR = 34;
  const radii = G.stateCounts.map(c => minR + (c/maxCount)*(maxR-minR));

  // Max edge flow for scaling
  const maxFlow = Math.max(...G.edgeFlows, 1e-9);

  const arrowsG = document.getElementById('graph-arrows');
  const nodesG = document.getElementById('graph-nodes');
  arrowsG.innerHTML = '';
  nodesG.innerHTML = '';

  // Draw edges
  G.edges.forEach(([si, di, junc], ei) => {
    const flow = G.edgeFlows[ei];
    if (flow < 1e-10) return;
    const src = positions[si], dst = positions[di];
    if (!src || !dst) return;
    const dx = dst.x-src.x, dy = dst.y-src.y;
    const dist = Math.sqrt(dx*dx+dy*dy) || 1;
    const ux = dx/dist, uy = dy/dist;
    const isDom = domSet.has(si+','+di);
    const sw = isDom ? 2.5 + (flow/maxFlow)*5 : 1 + (flow/maxFlow)*3.5;

    const srcL = G.labels[si] === '()' ? 0 : (G.labels[si].match(/[(]/g)||[]).length;
    const dstL = G.labels[di] === '()' ? 0 : (G.labels[di].match(/[(]/g)||[]).length;
    const layerDiff = dstL - srcL;
    const color = isDom ? C.hi : C.pri;
    const opacity = isDom ? 0.9 : 0.5;
    const marker = isDom ? 'url(#arr-dom)' : 'url(#arr)';

    if (layerDiff === 1 && Math.abs(dx) < W*0.6) {
      arrowsG.appendChild(mkEl('line', {
        x1: src.x+ux*radii[si], y1: src.y+uy*radii[si],
        x2: dst.x-ux*(radii[di]+8), y2: dst.y-uy*(radii[di]+8),
        stroke: color, 'stroke-width': sw, 'marker-end': marker, opacity,
      }));
    } else {
      const mx = (src.x+dst.x)/2 + (src.y-dst.y)*0.25;
      const my = (src.y+dst.y)/2 + (dst.x-src.x)*0.15;
      arrowsG.appendChild(mkEl('path', {
        d: `M${src.x+ux*radii[si]},${src.y+uy*radii[si]} Q${mx},${my} ${dst.x-ux*(radii[di]+8)},${dst.y-uy*(radii[di]+8)}`,
        fill: 'none', stroke: color, 'stroke-width': sw,
        'stroke-dasharray': isDom ? 'none' : '4,2',
        'marker-end': marker, opacity,
      }));
    }
  });

  // Draw nodes
  const tooltip = document.getElementById('tooltip');
  const domNodeSet = new Set(G.dominantPath || []);
  for (let si = 0; si < G.nStates; si++) {
    const p = positions[si];
    if (!p) continue;
    const r = radii[si];
    const layer = G.labels[si] === '()' ? 0 : (G.labels[si].match(/[(]/g)||[]).length;
    const isDom = domNodeSet.has(si);
    const col = isDom ? C.hi : layerColor(layer);
    const g = mkEl('g', { cursor: 'pointer' });

    if (isDom) {
      g.appendChild(mkEl('circle', {
        cx: p.x, cy: p.y, r: r+3.5,
        fill: 'none', stroke: C.hi, 'stroke-width': '2.5', opacity: '0.6',
      }));
    }

    g.appendChild(mkEl('circle', {
      cx: p.x, cy: p.y, r, fill: col, stroke: '#fff', 'stroke-width': '1.5',
    }));

    // Label
    let lbl = G.labels[si];
    if (lbl.length > 14) lbl = lbl.substring(0, 12) + '..';
    const txt = mkEl('text', {
      x: p.x, y: p.y + 1, 'text-anchor': 'middle', 'dominant-baseline': 'central',
      fill: layer === 0 ? '#333' : '#fff', 'font-size': Math.max(7, 10 - G.nStates*0.1) + 'px',
      'font-family': 'monospace', 'font-weight': isDom ? '700' : '500',
    });
    txt.textContent = lbl;
    g.appendChild(txt);

    // Tooltip
    g.addEventListener('mouseenter', e => {
      tooltip.style.display = 'block';
      tooltip.innerHTML = `<b>${G.labels[si]}</b><br>` +
        `Reads: ${G.stateCounts[si].toFixed(1)}<br>` +
        `IPD: ${(G.observedIPD[si]*100).toFixed(2)}%<br>` +
        `Layer: ${layer}` +
        (isDom ? '<br><b style="color:#F1C045">★ Dominant path</b>' : '');
    });
    g.addEventListener('mousemove', e => {
      tooltip.style.left = (e.clientX+14)+'px';
      tooltip.style.top = (e.clientY-10)+'px';
    });
    g.addEventListener('mouseleave', () => { tooltip.style.display = 'none'; });

    nodesG.appendChild(g);
  }
}

// ── IPD bar chart ──
function drawIPD() {
  const el = document.getElementById('ipd-bars');
  const maxIPD = Math.max(...G.observedIPD, 0.01);
  const domSet = new Set(G.dominantPath || []);

  el.innerHTML = G.observedIPD.map((ipd, i) => {
    const pct = (ipd * 100).toFixed(2);
    const w = (ipd / maxIPD * 100).toFixed(1);
    const layer = G.labels[i] === '()' ? 0 : (G.labels[i].match(/[(]/g)||[]).length;
    const col = domSet.has(i) ? C.hi : layerColor(layer);
    const star = domSet.has(i) ? ' ★' : '';
    let lbl = G.labels[i];
    if (lbl.length > 20) lbl = lbl.substring(0, 18) + '..';
    return `<div class="bar-block">
      <div class="bar-name"><span title="${G.labels[i]}">${lbl}${star}</span><span class="pct">${pct}%</span></div>
      <div class="bar-track"><div class="bar-fill" style="width:${w}%; background:${col};"></div></div>
    </div>`;
  }).join('');
}

// Load first gene
loadGene(0);
</script>
</body>
</html>'''


# ══════════════════════════════════════════════════════════════
#  Section F: Summary Figures
# ══════════════════════════════════════════════════════════════

def save_metrics_csv(results, output_path):
    """Save genome-wide metrics to CSV."""
    fields = ['gene_name', 'chrom', 'strand', 'n_reads', 'n_spliced_reads',
              'n_junctions', 'n_states', 'n_edges', 'max_layer', 'dada_string',
              'ordering_consistency', 'dominant_path_fraction', 'path_entropy',
              'ci_r2_log', 'ci_kl']

    with open(output_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        for r in results:
            if r is not None:
                row = {k: r[k] for k in fields}
                writer.writerow(row)

    print(f"  Saved {output_path} ({sum(1 for r in results if r is not None)} genes)")


def fig1_predeterminedness_dist(results, output_path):
    """Histograms of predeterminedness metrics across all genes."""
    valid = [r for r in results if r is not None]
    oc = [r['ordering_consistency'] for r in valid if not np.isnan(r['ordering_consistency'])]
    dp = [r['dominant_path_fraction'] for r in valid if not np.isnan(r['dominant_path_fraction'])]
    pe = [r['path_entropy'] for r in valid if not np.isnan(r['path_entropy'])]

    fig, axes = plt.subplots(1, 3, figsize=(14, 4))
    fig.suptitle('Predeterminedness of Splicing Order Across the Genome',
                 fontsize=13, fontweight='bold')

    # Ordering consistency
    ax = axes[0]
    ax.hist(oc, bins=30, color='#335E95', alpha=0.8, edgecolor='white')
    high_pct = 100 * sum(1 for x in oc if x > 0.8) / len(oc) if oc else 0
    ax.axvline(0.8, color='#2d8a4e', ls='--', lw=1.5, label=f'>0.8: {high_pct:.0f}%')
    ax.set_xlabel('Ordering Consistency')
    ax.set_ylabel('Number of genes')
    ax.set_title('Pairwise Junction Ordering')
    ax.legend(fontsize=9)

    # Dominant path fraction
    ax = axes[1]
    ax.hist(dp, bins=30, color='#2d8a4e', alpha=0.8, edgecolor='white')
    high_pct = 100 * sum(1 for x in dp if x > 0.5) / len(dp) if dp else 0
    ax.axvline(0.5, color='#c0392b', ls='--', lw=1.5, label=f'>0.5: {high_pct:.0f}%')
    ax.set_xlabel('Dominant Path Fraction')
    ax.set_title('Path Concentration')
    ax.legend(fontsize=9)

    # Path entropy
    ax = axes[2]
    ax.hist(pe, bins=30, color='#F1C045', alpha=0.8, edgecolor='white')
    low_pct = 100 * sum(1 for x in pe if x < 0.3) / len(pe) if pe else 0
    ax.axvline(0.3, color='#2d8a4e', ls='--', lw=1.5, label=f'<0.3: {low_pct:.0f}%')
    ax.set_xlabel('Normalized Path Entropy')
    ax.set_title('Splice Program Diversity')
    ax.legend(fontsize=9)

    plt.tight_layout(rect=[0, 0, 1, 0.92])
    fig.savefig(output_path)
    plt.close(fig)
    print(f"  Saved {output_path}")


def fig2_ordering_heatmaps(showcase, output_path):
    """Grid of ordering heatmaps for showcase genes."""
    # Select top 3 (tree) and bottom 3 (bushy)
    sorted_sc = sorted(showcase, key=lambda x: x['ordering_consistency'], reverse=True)
    top3 = sorted_sc[:3]
    bot3 = sorted_sc[-3:]
    genes = top3 + bot3

    fig, axes = plt.subplots(2, 3, figsize=(14, 8))
    fig.suptitle('Junction Splicing Order: Predetermined (top) vs Alternative (bottom)',
                 fontsize=13, fontweight='bold')

    for idx, (r, ax) in enumerate(zip(genes, axes.flat)):
        mat = r['_ordering_matrix']
        n = mat.shape[0]
        jlabs = [f"({j[0]},{j[1]})" for j in r['_all_juncs']]

        # Custom colormap: red (0) → white (0.5) → green (1)
        cmap = mcolors.LinearSegmentedColormap.from_list(
            'order', ['#c0392b', '#ffffff', '#2d8a4e'])

        im = ax.imshow(mat, cmap=cmap, vmin=0, vmax=1, interpolation='nearest')
        ax.set_xticks(range(n))
        ax.set_yticks(range(n))
        ax.set_xticklabels(jlabs, fontsize=6, rotation=45, ha='right')
        ax.set_yticklabels(jlabs, fontsize=6)

        # Annotate cells
        for i in range(n):
            for j in range(n):
                if i != j:
                    color = 'white' if abs(mat[i, j] - 0.5) > 0.3 else 'black'
                    ax.text(j, i, f'{mat[i,j]:.2f}', ha='center', va='center',
                            fontsize=5.5, color=color, fontweight='bold')

        tag = "TREE" if idx < 3 else "BUSHY"
        ax.set_title(f"{r['gene_name']} [{tag}]\nOC={r['ordering_consistency']:.3f}",
                     fontsize=9, fontweight='bold')

    fig.colorbar(im, ax=axes, fraction=0.02, pad=0.04,
                 label='P(row spliced before column)')
    plt.tight_layout(rect=[0, 0, 0.95, 0.92])
    fig.savefig(output_path)
    plt.close(fig)
    print(f"  Saved {output_path}")


def fig3_tree_vs_bushy(showcase, output_path):
    """Side-by-side DADA graph examples using networkx-style layout."""
    sorted_sc = sorted(showcase, key=lambda x: x['ordering_consistency'], reverse=True)
    tree_gene = sorted_sc[0]
    bushy_gene = sorted_sc[-1]

    fig, axes = plt.subplots(1, 2, figsize=(14, 7))
    fig.suptitle('Splice Program Structure: Predetermined vs Alternative',
                 fontsize=13, fontweight='bold')

    for ax, gene, label in [(axes[0], tree_gene, 'Predetermined (Tree)'),
                             (axes[1], bushy_gene, 'Alternative (Bushy)')]:
        layers = gene['_layers']
        ipd = gene['_observed_ipd']
        edges = gene['_edges']
        edge_flows = gene['_edge_flows']
        labels = gene['_isoform_labels']
        dom_path = set(gene['_dominant_path']) if gene['_dominant_path'] else set()

        max_layer = max(layers)
        n_states = len(layers)

        # Compute positions (layered)
        layer_nodes = defaultdict(list)
        for i, l in enumerate(layers):
            layer_nodes[l].append(i)

        pos = {}
        for l in range(max_layer + 1):
            nodes = layer_nodes[l]
            for xi, ni in enumerate(nodes):
                x = (xi + 0.5) / max(len(nodes), 1)
                y = 1.0 - l / max(max_layer, 1)
                pos[ni] = (x, y)

        # Draw edges
        max_flow = max(edge_flows) if len(edge_flows) > 0 else 1
        dom_edges = set()
        dp = gene['_dominant_path']
        if dp and len(dp) > 1:
            for i in range(len(dp) - 1):
                dom_edges.add((dp[i], dp[i+1]))

        for k, (src, dst, j, ji) in enumerate(edges):
            if src not in pos or dst not in pos:
                continue
            flow = edge_flows[k]
            lw = 0.5 + (flow / max(max_flow, 1e-9)) * 4
            is_dom = (src, dst) in dom_edges
            color = '#F1C045' if is_dom else '#335E95'
            alpha = 0.9 if is_dom else 0.3
            ax.annotate('', xy=pos[dst], xytext=pos[src],
                        arrowprops=dict(arrowstyle='->', color=color,
                                       lw=lw, alpha=alpha))

        # Draw nodes
        max_ipd = max(ipd) if len(ipd) > 0 else 1
        for si in range(n_states):
            if si not in pos:
                continue
            x, y = pos[si]
            size = 100 + (ipd[si] / max(max_ipd, 1e-9)) * 600
            is_dom = si in dom_path
            color = '#F1C045' if is_dom else layerColor(layers[si])
            ax.scatter(x, y, s=size, c=color, edgecolors='white',
                      linewidth=1.5, zorder=5)

            if ipd[si] > 0.03:
                lbl = labels[si] if len(labels[si]) <= 12 else labels[si][:10] + '..'
                ax.annotate(lbl, (x, y), fontsize=5, ha='center', va='center',
                           color='white' if not is_dom else '#333',
                           fontweight='bold', zorder=6)

        ax.set_xlim(-0.05, 1.05)
        ax.set_ylim(-0.1, 1.1)
        ax.set_title(f"{label}\n{gene['gene_name']} — "
                     f"OC={gene['ordering_consistency']:.3f}, "
                     f"J={gene['n_junctions']}, S={gene['n_states']}",
                     fontsize=10, fontweight='bold')
        ax.set_xlabel('← States at same layer →', fontsize=9)
        ax.set_ylabel('Splicing progression (layer) →', fontsize=9)
        ax.set_xticks([])
        ax.set_yticks([])

    plt.tight_layout(rect=[0, 0, 1, 0.92])
    fig.savefig(output_path)
    plt.close(fig)
    print(f"  Saved {output_path}")


# Use existing layerColor function for matplotlib
def layerColor(l):
    colors = ['#8b949e', '#440154', '#31688e', '#35b779', '#fde725', '#d4aa00', '#b07800', '#8b0000']
    return colors[min(l, len(colors) - 1)]


# ══════════════════════════════════════════════════════════════
#  Main
# ══════════════════════════════════════════════════════════════

def main():
    t0 = time.time()
    print("=" * 60)
    print("  TRUE DADA: Validating Predetermined Splicing")
    print("=" * 60)

    # A. Load data
    print("\n[A] Loading data...")
    genes = load_genes(DATA_PKL)
    ci_metrics = load_ci_metrics(CI_CSV)

    # B+C. Analyze each gene
    print("\n[B+C] Analyzing genes (reconstruction + metrics)...")
    results = []
    n_pass = 0
    for i, gene in enumerate(genes):
        if gene.n_spliced_reads < 10:
            results.append(None)
            continue
        try:
            states, edges, all_juncs, junc_idx, layers = reconstruct_dada(gene)
        except Exception:
            results.append(None)
            continue

        if len(all_juncs) < 2:
            results.append(None)
            continue

        ci_row = ci_metrics.get(gene.gene_name)
        result = analyze_gene(gene, ci_row)
        results.append(result)
        if result is not None:
            n_pass += 1

        if (i + 1) % 500 == 0:
            print(f"    Processed {i+1}/{len(genes)} ({n_pass} pass)")

    print(f"  Total: {n_pass} genes analyzed")

    # D. Select showcase genes
    print("\n[D] Selecting showcase genes...")
    showcase = select_genes(results)

    if not showcase:
        print("  ERROR: No showcase genes selected. Exiting.")
        return

    # E. Write interactive HTML
    print("\n[E] Writing interactive HTML gallery...")
    html_path = os.path.join(OUT, "showcase_genes.html")
    write_gallery_html(showcase, html_path)

    # E2. Write radial domtree HTML
    print("\n[E2] Writing radial domtree visualization...")
    radial_path = os.path.join(OUT, "radial_showcase.html")
    write_radial_html(showcase, radial_path,
                      title="TRUE DADA — Radial Splice Programs (Churchman nano-COP)")
    print(f"  Radial HTML: {radial_path}")

    # F. Save metrics + figures
    print("\n[F] Saving metrics and figures...")
    csv_path = os.path.join(OUT, "predeterminedness_metrics.csv")
    save_metrics_csv(results, csv_path)

    fig1_predeterminedness_dist(
        results, os.path.join(OUT, "fig1_predeterminedness_dist.png"))

    if len(showcase) >= 6:
        fig2_ordering_heatmaps(
            showcase, os.path.join(OUT, "fig2_ordering_heatmaps.png"))

    if len(showcase) >= 2:
        fig3_tree_vs_bushy(
            showcase, os.path.join(OUT, "fig3_tree_vs_bushy.png"))

    elapsed = time.time() - t0
    print(f"\n{'=' * 60}")
    print(f"  Done in {elapsed:.1f}s. All outputs in:")
    print(f"  {OUT}/")
    print(f"{'=' * 60}")


if __name__ == '__main__':
    main()
