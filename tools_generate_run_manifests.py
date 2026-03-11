#!/usr/bin/env python3
import argparse
import csv
import json
from pathlib import Path
from datetime import datetime

DEFAULT_EXP_DIR = Path(".")


def read_csv(path: Path):
    with path.open(newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def safe_mkdir(p: Path):
    p.mkdir(parents=True, exist_ok=True)


def pick_first(*values):
    for v in values:
        if v is not None:
            return v
    return None


def parse_int_or_keep(value):
    if isinstance(value, bool):
        return value
    if isinstance(value, int):
        return value
    if isinstance(value, str) and value.strip().isdigit():
        return int(value.strip())
    return value


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="tools_generate_run_manifests.py",
        description=(
            "Generate runs/<protein_id>/<peptide_id>/run_manifest.json for every "
            "protein×peptide combination and write analysis/campaigns_matrix.csv."
        ),
    )
    p.add_argument(
        "--exp-dir",
        default=str(DEFAULT_EXP_DIR),
        help="Experiment root directory (default: current directory).",
    )
    p.add_argument("--proteins-csv", default="proteins.csv", help="Proteins CSV path (relative to --exp-dir if not absolute).")
    p.add_argument("--peptides-csv", default="peptides.csv", help="Peptides CSV path (relative to --exp-dir if not absolute).")
    p.add_argument(
        "--experiment-manifest",
        default="experiment_manifest.json",
        help="Experiment manifest JSON path (relative to --exp-dir if not absolute).",
    )
    p.add_argument("--runs-root", default="runs", help="Runs output directory (relative to --exp-dir if not absolute).")
    p.add_argument("--analysis-dir", default="analysis", help="Analysis output directory (relative to --exp-dir if not absolute).")
    p.add_argument("--matrix-csv", default="campaigns_matrix.csv", help="Campaign matrix filename under --analysis-dir.")
    return p


def resolve_path(exp_dir: Path, candidate: str) -> Path:
    p = Path(candidate)
    return p if p.is_absolute() else (exp_dir / p)


def main(argv=None):
    args = build_parser().parse_args(argv)

    exp_dir = Path(args.exp_dir).resolve()
    proteins_csv = resolve_path(exp_dir, args.proteins_csv)
    peptides_csv = resolve_path(exp_dir, args.peptides_csv)
    exp_manifest_json = resolve_path(exp_dir, args.experiment_manifest)
    runs_root = resolve_path(exp_dir, args.runs_root)
    analysis_dir = resolve_path(exp_dir, args.analysis_dir)
    matrix_csv = analysis_dir / args.matrix_csv

    if not proteins_csv.exists():
        raise SystemExit(f"Missing {proteins_csv}")
    if not peptides_csv.exists():
        raise SystemExit(f"Missing {peptides_csv}")
    if not exp_manifest_json.exists():
        raise SystemExit(f"Missing {exp_manifest_json}")

    proteins = read_csv(proteins_csv)
    peptides = read_csv(peptides_csv)

    expm = json.loads(exp_manifest_json.read_text(encoding="utf-8"))
    experiment_id = expm.get("experiment_id", "expXXX")
    layout_version = expm.get("layout_version", "1.0")
    defaults = expm.get("defaults", {}) or {}
    defaults_config = defaults.get("config", {}) if isinstance(defaults.get("config", {}), dict) else {}
    defaults_replicas = defaults.get("replicas", {}) if isinstance(defaults.get("replicas", {}), dict) else {}
    inputs_dir = Path(expm.get("inputs_dir", "inputs"))
    receptors_dir = exp_dir / inputs_dir / "receptors"

    default_reps = pick_first(
        defaults_replicas.get("planned_count_default"),
        defaults_replicas.get("planned_count"),
        defaults.get("planned_count_default"),
        defaults.get("planned_count"),
        expm.get("planned_count_default"),
    )
    if default_reps is None:
        default_reps = 30

    default_index_start = pick_first(defaults_replicas.get("index_start"), 1)
    default_suffix_width = pick_first(defaults_replicas.get("suffix_width"), 3)
    default_seed_base = defaults_replicas.get("seed_base")

    safe_mkdir(runs_root)
    safe_mkdir(analysis_dir)

    rows = []
    now = datetime.now().isoformat(timespec="seconds")

    for p in proteins:
        protein_id = p["protein_id"].strip()
        protein_name = p.get("protein_name", "").strip()
        protein_pdb = p.get("protein_pdb", "").strip()
        receptor_pdb_file = p["receptor_pdb_file"].strip()

        receptor_pdb_path = str((receptors_dir / receptor_pdb_file).as_posix())

        for pep in peptides:
            peptide_id = pep["peptide_id"].strip()
            peptide_set = pep.get("peptide_set", "").strip()
            peptide_seq = pep["peptide_seq"].strip()

            run_dir = runs_root / protein_id / peptide_id
            targets_dir = run_dir / "targets"
            replicas_dir = run_dir / "replicas"
            results_dir = run_dir / "results"
            results_parsed_dir = run_dir / "results_parsed"
            logs_dir = run_dir / "logs"

            for d in [targets_dir, replicas_dir, results_dir, results_parsed_dir, logs_dir]:
                safe_mkdir(d)

            target_trg = str((targets_dir / "receptor.trg").as_posix())

            run_manifest = {
                "experiment_id": experiment_id,
                "layout_version": layout_version,
                "created_at": now,
                "engine": expm.get("engine", "ADCP+OpenMM"),
                "directory_layout": expm.get("directory_layout", {}),
                "campaign": {
                    "campaign_id": f"{experiment_id}_{protein_id}_{peptide_id}",
                    "protein_id": protein_id,
                    "protein_name": protein_name,
                    "protein_pdb": protein_pdb,
                    "peptide_id": peptide_id,
                    "peptide_set": peptide_set,
                },
                "inputs": {
                    "receptor_pdb_file": receptor_pdb_file,
                    "receptor_pdb_path": receptor_pdb_path,
                    "peptide_seq": peptide_seq,
                },
                "targets": {
                    "target_trg": target_trg,
                    "target_trg_mode": "per_condition",
                },
                "config": {
                    "adcp_numSteps": pick_first(defaults_config.get("adcp_numSteps"), defaults.get("adcp_numSteps")),
                    "adcp_nbRuns": pick_first(defaults_config.get("adcp_nbRuns"), defaults.get("adcp_nbRuns")),
                    "adcp_nmodes": pick_first(defaults_config.get("adcp_nmodes"), defaults.get("adcp_nmodes")),
                    "adcp_maxCores": pick_first(defaults_config.get("adcp_maxCores"), defaults.get("adcp_maxCores")),
                    "omm_nmin": pick_first(defaults_config.get("omm_nmin"), defaults.get("omm_nmin")),
                    "omm_max_itr": pick_first(defaults_config.get("omm_max_itr"), defaults.get("omm_max_itr")),
                    "omm_environment": pick_first(defaults_config.get("omm_environment"), defaults.get("omm_environment")),
                    "rerank_mode": pick_first(defaults_config.get("rerank_mode"), defaults.get("rerank_mode")),
                    "cluster_cutoff": pick_first(defaults_config.get("cluster_cutoff"), defaults.get("cluster_cutoff")),
                },
                "replicas": {
                    "planned_count": parse_int_or_keep(default_reps),
                    "index_start": parse_int_or_keep(default_index_start),
                    "suffix_width": parse_int_or_keep(default_suffix_width),
                    "seed_base": parse_int_or_keep(default_seed_base),
                    "note": "Adjust planned_count per campaign if needed (replicas may vary across campaigns).",
                },
                "paths": {
                    "run_dir": str(run_dir.as_posix()),
                    "targets_dir": str(targets_dir.as_posix()),
                    "replicas_dir": str(replicas_dir.as_posix()),
                    "results_dir": str(results_dir.as_posix()),
                    "results_parsed_dir": str(results_parsed_dir.as_posix()),
                    "logs_dir": str(logs_dir.as_posix()),
                },
                "status": "planned",
            }

            out_path = run_dir / "run_manifest.json"
            out_path.write_text(json.dumps(run_manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

            rows.append(
                {
                    "campaign_id": run_manifest["campaign"]["campaign_id"],
                    "protein_id": protein_id,
                    "protein_name": protein_name,
                    "protein_pdb": protein_pdb,
                    "peptide_id": peptide_id,
                    "peptide_set": peptide_set,
                    "peptide_seq": peptide_seq,
                    "receptor_pdb_file": receptor_pdb_file,
                    "receptor_pdb_path": receptor_pdb_path,
                    "run_dir": str(run_dir.as_posix()),
                    "target_trg": target_trg,
                    "planned_replicas": run_manifest["replicas"]["planned_count"],
                }
            )

    fieldnames = [
        "campaign_id",
        "protein_id",
        "protein_name",
        "protein_pdb",
        "peptide_id",
        "peptide_set",
        "peptide_seq",
        "receptor_pdb_file",
        "receptor_pdb_path",
        "run_dir",
        "target_trg",
        "planned_replicas",
    ]
    with matrix_csv.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for r in rows:
            w.writerow(r)

    print(f"[ok] Created {len(rows)} campaigns.")
    print(f"[ok] Wrote campaign matrix: {matrix_csv}")
    if rows:
        print(f"[ok] Example run_manifest: {rows[0]['run_dir']}/run_manifest.json")


if __name__ == "__main__":
    main()
