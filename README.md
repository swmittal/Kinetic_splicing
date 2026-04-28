# Kinetic Splicing — the DADA Rate Recovery Sandbox

An interactive browser tool for exploring when the splicing-rate inverse problem is identifiable.

Given the steady-state isoform-proportion data (IPD) of a gene's alternative-splicing isoforms, when can we recover the underlying kinetic rates by non-negative least squares? The answer depends on the topological structure of the splicing graph and on which prior we place on the rate vector. This sandbox is designed to make the dependency visible, knob by knob.

Open [results/nnls_test/dada_rate_sandbox.html](results/nnls_test/dada_rate_sandbox.html) in a browser. No backend needed.

---

## Problem setup

A gene is represented by a **DADA string** of splice sites, letters D (donor) and A (acceptor) along the transcript (e.g. `DADADA` for three canonical introns). The isoform state graph `G = (V, E)` is built by enumerating all legal junctions `(D_i, A_j)` and all sets of pairwise-compatible junctions:

- **states** `V`: sets of already-spliced junctions; `s_0 = ∅` is the null (unspliced) state.
- **edges** `E`: transitions `s → s ∪ {j}` that fire junction `j` compatible with state `s`.
- `N = |V|`, `M = |E|`.

A continuous-time Markov chain on `G` models splicing + export. New pre-mRNA enters `s_0` at production rate 1; each state `i` exports at rate `λ_i`; edges fire at rates `r_k > 0`. At steady state:

    for every state i:    Σ_{k: dst(k)=i} r_k p_{src(k)}  =  (Σ_{k: src(k)=i} r_k + λ_i) p_i  −  b_i

where `b_0 = 1` and `b_i = 0` otherwise. Rearranging gives a linear system in `r`:

    A r = y ,   A ∈ R^{N × M} ,   y_i = λ_i p_i − b_i

`A` is a **weighted incidence matrix**: column `k` has `+p_{src(k)}` at row `dst(k)` and `−p_{src(k)}` at row `src(k)`. Its column sums are zero, so `rank(A) ≤ N − 1`. The inverse problem is: given `p` and `λ`, recover `r`.

We solve `A r = y` subject to `r ≥ 0` using non-negative least squares (NNLS). The sandbox offers two NNLS variants and four prior regimes on the ground-truth rate vector.

---

## The four prior regimes

The sandbox samples ground-truth rates `r*` under exactly one of four structural priors. They differ in (i) which edges of `G` are active and (ii) how rates are distributed across active edges. One shared log-normal scale `r ~ R_0 · exp(σ · N(0, 1))`, with `R_0 = 75` and `σ` = Spread slider, controls rate magnitude.

**Null.** Full graph active. Every edge is drawn independently. No pooling. Unknowns `= M`.

**Pure kinetic.** Full graph active. One rate `r_j` per junction, shared across every edge that uses it. Unknowns `= K_junc` (the number of junctions).

**Purely predetermined.** A branching tree rooted at `s_0`, grown by breadth-first search; at each state, each outgoing edge is independently included with probability `β` (the Branching slider), with a floor of one child per node. Per-edge rates on the tree; non-tree edges are inactive. Active subgraph has cyclomatic number `c(G*) = 0` (tree). Unknowns `= |tree| − 1`.

**Mixed.** The Predetermined tree + a controlled number of cycle-closing extra edges ("diamonds") whose both endpoints are already in the tree. Each diamond adds one cycle; `c(G*) = n_diamonds`. Per-edge rates on every active edge. Unknowns `= |tree| − 1 + n_diamonds`.

## The two NNLS solvers

**Per-edge NNLS** treats each active edge as an independent unknown. `A` has one column per active edge. Unknowns `= M_active`.

**Pooled NNLS** sums the columns of `A` that share a junction into a single column, producing a reduced design matrix over junctions. Unknowns `= K_junc_active` (junctions with at least one active edge). After solving, the recovered per-edge rate is `r̂_k = r̂_{j(k)}` for every edge `k`.

Both are Lawson–Hanson NNLS. The difference is the column grouping: the pooled solver assumes the truth is junction-wise context-independent (one rate per junction), the per-edge solver assumes no such pooling.

---

## The identifiability matrix

NNLS on `A r = y` uniquely recovers `r*` iff the effective number of unknowns is at most the number of independent balance rows (`N − 1`) **and** the imposed rate-structure prior matches the truth. The 4 × 2 matrix:

| prior \ solver | per-edge NNLS | pooled NNLS |
|---|---|---|
| **Null** | fails (M unknowns ≫ N) | fails (projects pooled structure onto per-edge truth) |
| **Pure kinetic** | fails (M unknowns ≫ N) | **exact** (K_junc ≪ N, structure matches) |
| **Purely predetermined** | **exact** (tree, M* = |tree|−1 ≤ N−1) | exact iff no junction repeats in the tree; fails as soon as branching reuses any junction |
| **Mixed** | fails once `n_diamonds > 0` (cycle-space null direction) | fails (per-edge truth + junction reuse) |

The sandbox displays this live. Pick a prior, pick a solver, hit Randomize; the diagnostics bar reports `rank(A)`, `Nullity`, `log₁₀ κ(A)`, `R² on observable edges`, and the residual `‖A r̂ − y‖₂`.

## Why NNLS fails in the underdetermined case

When `rank(A) < M_active`, the solution set `{r : A r = y, r ≥ 0}` is a polytope with multiple vertices. NNLS (Lawson–Hanson) returns one vertex — the one its greedy active-set walk reaches first. The cyclomatic number `c(G*) = M_active − N_active + 1` is exactly the dimension of the null space of `A`; every added cycle gives NNLS one more direction in which it can shuffle flux without changing the observed IPD. Generically the returned vertex has `c(G*)` edges pinned to zero, and the remaining edges carry all the redistributed flux.

Non-negativity does not rescue this in our setting. The null-space directions of `A` are cycle flows (positive on forward edges, negative on reverse), so non-negativity constrains but does not eliminate the ambiguity. The condition under which non-negativity alone yields exact recovery — Slawski & Hein's "M+ property" — does not hold for CTMC balance matrices.

See [notes/cyclomatic_primer.html](notes/cyclomatic_primer.html) for the full derivation (textbook style, opens in browser, math via KaTeX).

---

## Using the sandbox

Open [results/nnls_test/dada_rate_sandbox.html](results/nnls_test/dada_rate_sandbox.html) in a browser.

**Controls.**
- **DADA string.** Any sequence of `D` and `A`, e.g. `DADADA`, `DADDAADA`. Hit Apply.
- **Prior.** One of Null / Pure kinetic / Purely predetermined / Mixed. Picks the ground-truth generator.
- **Solver.** Per-edge or Pooled. Picks the NNLS variant used for inversion.
- **Spread slider** (`σ`). Log-normal standard deviation on generated rates. 0 means all rates equal.
- **Structure slider** (label varies by prior: Branching for Predetermined, Diamonds for Mixed; disabled for Null / Pure kinetic).
- **Export profile + null fraction.** Shape of `λ(i)` as a function of splice completion. Default "parabolic" with null fraction 0 works as a reasonable assumption.
- **Randomize** draws a fresh instance under the current prior. **Reset** returns to uniform rates.

**Diagnostics bar.** Nine metrics:
1. **Balance rows**  — observable balance equations used in the solver.
2. **Unknowns**  — columns of `A` (edges or pooled junctions).
3. **rank(A)**.
4. **Nullity**  — `Unknowns − rank(A)`. Zero when identifiable; positive means cycle-space ambiguity.
5. **log₁₀ κ(A)**  — condition number of `A`. Above 6 means noise-sensitive even if formally full rank.
6. **R² (observable edges)**  — coefficient of determination on edges whose endpoints both have above-threshold abundance and whose true rate is above trace level. `N/A` if fewer than two such edges or their true rates have zero variance.
7. **Observable edges** — `n_obs / n_total`.
8. `‖A r̂ − y‖₂` — balance residual. Machine-precision (`~10⁻¹⁶`) when NNLS has a feasible exact solution.
9. **p(s_0)** — unspliced-fraction, a diagnostic on how much mass sits at the null state.

**DAG view.** Nodes are isoform states laid out by layer; edge width ∝ rate. Click any edge or True-rate cell in the recovery table to edit by hand. Toggles for overlay recovered-rate overlay, colour-by-identifiability, and hide-unobservable are just below the diagnostics bar.

**Sashimi panel.** Hover a node to see its splice pattern (the arcs of the mature mRNA at that state).

**About / Math modal.** Click the button in the top-right for the self-contained math writeup of everything above. Uses KaTeX to render equations.

---

## Scaling to 8-intron genes

The full DADA graph at `K = 8` introns has 1,597 states and 5,911 edges — too large for the dense JS solver in the browser. Two Python companions run the same math at scale:

- [enforce/forced_order.py](enforce/forced_order.py) — samples a random total order over junctions, walks the DAG greedily, and runs forward-solve + NNLS on the resulting single path. Under the Predetermined prior at `Branching = 0`.
- [enforce/random_tree.py](enforce/random_tree.py) — samples a uniform random spanning tree rooted at the null state, runs forward-solve + NNLS on it. The richer branching limit of the Predetermined prior. At `K = 8`, NNLS on the spanning tree is in principle exact (tree) but numerically marginal — condition number growth is the practical ceiling.

Both scripts use `scipy.optimize.nnls` and `matplotlib` for visualization, and run in seconds on a laptop.

### Parameter sweeps

The `rate_sweep.js` + `rate_sweep_plot.py` + `run_rate_sweep.sh` triple sweeps `(sparsity, spread, coupling)` grids via SLURM and produces heatmaps and marginals of R² across the parameter space. This was useful for characterising NNLS behaviour in the earlier 3-knob model and is retained for reference. It is not directly relevant to the current four-prior formulation.

---

## File layout

```
kinetic_splicing/
  README.md                                  this file
  results/nnls_test/dada_rate_sandbox.html   the interactive sandbox (open in browser)
  notes/
    cyclomatic_primer.html                   textbook-style theory primer (KaTeX)
    slide_identifiability.html               one-pager identifiability slide
    textbook.css                             shared stylesheet for the notes
  enforce/
    forced_order.py                          single-path predetermined (Python + scipy)
    random_tree.py                           spanning-tree predetermined (Python + scipy)
  rate_sweep.js                              parameter-sweep driver (headless Node)
  rate_sweep_plot.py                         sweep-figure generator (matplotlib)
  run_rate_sweep.sh                          SLURM wrapper for the sweep
```

---

## Notation reference

| symbol | meaning |
|---|---|
| `N` | number of reachable states in the DADA graph |
| `M` | number of edges in the DADA graph |
| `K_junc` | number of distinct junctions |
| `c(G) = M − N + 1` | cyclomatic number (active subgraph's cycle count) |
| `r` | edge rate vector, `r_k ≥ 0` |
| `r*` | ground-truth rate vector generated by the prior |
| `r̂` | NNLS-recovered rate vector |
| `p` | steady-state state-occupancy vector (the IPD) |
| `λ_i` | export rate from state `i` |
| `A`, `y` | balance matrix and RHS of the inverse problem |
| `κ(A)` | condition number of `A`: `σ_max / σ_min` |
