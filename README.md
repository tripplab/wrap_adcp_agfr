# wrap_adcp_agfr

`wrap_adcp_agfr` is a Python script suite for peptide-control generation and ADCP/AGFR campaign execution and analysis.

## Quick install (clone + conda/micromamba + pip)

```bash
# 1) Clone the repository
git clone https://github.com/tripplab/wrap_adcp_agfr.git
cd wrap_adcp_agfr

# 2) Create and activate an environment (choose one)
micromamba create -n wrap_adcp_agfr python=3.11 -y && micromamba activate wrap_adcp_agfr
# or:
conda create -n wrap_adcp_agfr python=3.11 -y && conda activate wrap_adcp_agfr

# 3) Upgrade pip in that environment
python -m pip install --upgrade pip

# 4) Verify the main CLI
python pepctrl.py --help
```

> This project uses Python standard-library modules only for normal operation (no required third-party Python packages).

## Suggested workflow order

1. `pepctrl.py`
2. `tools_generate_run_manifests.py`
3. `run_adcp_agfr_replicas_campaign.py`
4. `anal_adcp_agfr_replicas_campaign.py`

## Runtime and dependencies

- Recommended runtime: Python 3.8+
- Primary modules used are from the Python standard library (`argparse`, `csv`, `json`, `math`, `os`, `pathlib`, `re`, `subprocess`, `typing`, etc.).

---

## Script 1: `pepctrl.py`

### Purpose
Generate peptide control sequences for a peptide of interest (POI).

### Modes
- `decoy`: matched decoys with composition/physicochemical constraints.
- `scramble`: residue permutations of the POI.
- `random`: random AA20 peptides with the same length as the POI.

### Core usage
```bash
python3 pepctrl.py --help
python3 pepctrl.py --poi ACDEFGHIKLMNPQRSTVWY --type decoy --n 25
```

### Input rules
- Provide POI with `--poi` or `--poi-file`.
- Allowed residues: `ACDEFGHIKLMNPQRSTVWY`.
- `--type` and `--n` are required.

### Output formats
- `--format`: `csv` (default), `tsv`, `fasta`, `jsonl`
- `--out -` writes to stdout; otherwise provide a file path.
- CSV/TSV include the POI row by default (disable with `--no-poi-row`).

### Key options
- Reproducibility: `--seed`
- Motif exclusion: `--exclude-motif`, `--exclude-motif-file`
- Diversity: `--max-identity-between-decoys`, `--max-identity`
- Generation methods: `--method` (`mutate`, `compose`, `hybrid`)
- Runtime controls: `--batch-size`, `--max-attempts`, `--allow-partial`, `-v/-vv/-vvv`, `--timing`

### Example commands
```bash
# Decoys
python3 pepctrl.py --poi ACDEFGHIKLMNPQRSTVWY --type decoy --n 25

# Scrambles to TSV
python3 pepctrl.py --poi ACDEFGHIKLMNPQRSTVWY --type scramble --n 20 --format tsv --out scrambles.tsv

# Random controls with fixed seed
python3 pepctrl.py --poi ACDEFGHIKLMNPQRSTVWY --type random --n 50 --seed 42 --out random.csv
```

### Exit codes
- `0`: success
- `2`: partial/insufficient result conditions
- `3`: configuration or input error

---

## Script 2: `tools_generate_run_manifests.py`

### Purpose
Generate campaign directory skeletons and per-condition `run_manifest.json` files for all protein × peptide combinations.

### Required inputs (in current working directory)
- `proteins.csv`
- `peptides.csv`
- `experiment_manifest.json`

### What it creates
For each `<protein_id>/<peptide_id>` condition:
- `runs/<protein_id>/<peptide_id>/run_manifest.json`
- Subdirectories: `targets/`, `replicas/`, `results/`, `results_parsed/`, `logs/`

Global output:
- `analysis/campaigns_matrix.csv`

### Behavior summary
- Reads proteins and peptides.
- Computes Cartesian product (all campaigns).
- Writes one run manifest per campaign.
- Builds a campaign matrix for downstream execution and analysis.

### Usage
```bash
# Show CLI help
python3 tools_generate_run_manifests.py -h

# Run with defaults (expects inputs in current directory)
python3 tools_generate_run_manifests.py

# Run using a specific experiment root and custom paths
python3 tools_generate_run_manifests.py \
  --exp-dir /path/to/experiment \
  --proteins-csv proteins.csv \
  --peptides-csv peptides.csv \
  --experiment-manifest experiment_manifest.json \
  --runs-root runs \
  --analysis-dir analysis
```

### CLI options (summary)
- `--exp-dir`: experiment root directory (default: current directory)
- `--proteins-csv`: proteins table path (default: `proteins.csv`)
- `--peptides-csv`: peptides table path (default: `peptides.csv`)
- `--experiment-manifest`: manifest JSON path (default: `experiment_manifest.json`)
- `--runs-root`: runs output directory (default: `runs`)
- `--analysis-dir`: analysis output directory (default: `analysis`)
- `--matrix-csv`: campaign matrix filename inside analysis dir (default: `campaigns_matrix.csv`)

### Notes
- Replica defaults are read from `defaults.replicas` (for example `planned_count_default`, `index_start`, `suffix_width`, `seed_base`) with compatibility fallbacks for older layouts.
- Config defaults support both `defaults.config.<key>` and legacy `defaults.<key>` fields.
- The script validates that every `receptor_pdb_file` referenced by `proteins.csv` exists as a file under `<exp_dir>/<inputs_dir>/receptors` and exits early with a detailed error if any are missing.
- Campaign count equals `N_proteins × N_peptides`.

---

## Script 3: `run_adcp_agfr_replicas_campaign.py`

### Purpose
Run a single campaign condition end-to-end: preflight checks, target initialization (`prepare_receptor`/`agfr`), and ADCP replicas.

### Expected context
Run against a condition manifest, typically:
`runs/<protein_id>/<peptide_id>/run_manifest.json`

The script resolves relative manifest paths against the experiment root inferred from manifest location.

### Typical experiment layout
```text
<exp_root>/
  inputs/
    receptors/
    peptides/
  runs/
    <protein_id>/<peptide_id>/
      run_manifest.json
      targets/
      replicas/
      results/
      logs/
```

### Typical usage
```bash
# Validate only
python3 run_adcp_agfr_replicas_campaign.py --run-manifest runs/P001/POI/run_manifest.json --validate-only --verbosity 2

# Initialize target only
python3 run_adcp_agfr_replicas_campaign.py --run-manifest runs/P001/POI/run_manifest.json --init-only --real-time --verbosity 2

# Run replicas
python3 run_adcp_agfr_replicas_campaign.py --run-manifest runs/P001/POI/run_manifest.json --real-time --verbosity 3
```

### Expected tools in PATH
- `prepare_receptor`
- `agfr`
- `adcp`

### Outputs
- Condition targets: `targets/receptor.pdbqt`, `targets/receptor.trg`
- Condition logs: `logs/preflight.log`, `logs/prepare_receptor.log`, `logs/agfr.log`, `logs/timing.json`, `logs/timing.csv`
- Replica artifacts under `replicas/replica_###/`:
  - `logs/adcp.log`
  - `docking/work/...`
  - `results/*_summary.dlg`, `results/*_out.pdb`, `results/README.txt`
- Campaign summary pointer in `results/README.txt`
- Runtime updates written back into `run_manifest.json`

---

## Script 4: `anal_adcp_agfr_replicas_campaign.py`

### Purpose
Parse campaign replica outputs and generate analysis tables for comparison, convergence, and quality checks.

### Analysis intent
Use replica-level distributions to compare peptides within each protein context (for example, POI vs controls), rather than relying on a single best run.

### Main data sources per replica
- Canonical: `results/*_summary.dlg`
- Optional: `results/*_omm_rescored_out.pdb`
- Optional: `docking/work/*_out.pdb`

### Winner definition per replica
- Winner is the `RankOpenMM = 1` model from the OMM Ranking section in `*_summary.dlg`.
- Primary score for comparisons: `winner_omm_dE_interaction`.

### Typical outputs (under `analysis/`)
- `replicas_parsed.csv`: one row per replica winner
- `topk_parsed.csv`: top-k models per replica
- `clusters_parsed.csv`: ADCP cluster rows per replica
- `group_summary.csv`: per protein-peptide summary statistics
- `discrimination.csv`: within-protein peptide comparisons

### How to use outputs
- Compare peptide groups using medians/IQR and effect sizes.
- Inspect convergence using winner cluster mode frequency/entropy.
- Track missing/incomplete data via QA/status fields.

---

## Minimal validation commands

```bash
# CLI sanity check
python3 pepctrl.py --help

# Smoke generation check
python3 pepctrl.py --poi ACDEFGHIKLMNPQRSTVWY --type scramble --n 3 --seed 1 --format csv | head
```

If these commands run successfully and output looks correct, the local setup is functional.

---

Author: trippm@tripplab.com [Feb 2026]
