#!/usr/bin/env python3
import csv, json, os
from pathlib import Path
from datetime import datetime

EXP_DIR = Path(".")

proteins_csv = EXP_DIR / "proteins.csv"
peptides_csv = EXP_DIR / "peptides.csv"
exp_manifest_json = EXP_DIR / "experiment_manifest.json"

runs_root = EXP_DIR / "runs"
analysis_dir = EXP_DIR / "analysis"
matrix_csv = analysis_dir / "campaigns_matrix.csv"

def read_csv(path: Path):
    with path.open(newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))

def safe_mkdir(p: Path):
    p.mkdir(parents=True, exist_ok=True)

def main():
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
    inputs_dir = Path(expm.get("inputs_dir", "inputs"))
    receptors_dir = EXP_DIR / inputs_dir / "receptors"

    # réplicas por defecto (editable por campaña)
    #default_reps = expm.get("defaults_replicas_planned_count", None)
    default_reps = expm.get("planned_count_default", None)
    if default_reps is None:
        # fallback conservador
        default_reps = 30

    safe_mkdir(runs_root)
    safe_mkdir(analysis_dir)

    rows = []
    now = datetime.now().isoformat(timespec="seconds")

    for p in proteins:
        protein_id = p["protein_id"].strip()
        protein_name = p.get("protein_name","").strip()
        protein_pdb = p.get("protein_pdb","").strip()
        receptor_pdb_file = p["receptor_pdb_file"].strip()

        # ruta absoluta/relativa coherente dentro del experimento
        receptor_pdb_path = str((receptors_dir / receptor_pdb_file).as_posix())

        for pep in peptides:
            peptide_id = pep["peptide_id"].strip()
            peptide_set = pep.get("peptide_set","").strip()
            peptide_seq = pep["peptide_seq"].strip()

            run_dir = runs_root / protein_id / peptide_id
            targets_dir = run_dir / "targets"
            replicas_dir = run_dir / "replicas"
            results_dir = run_dir / "results"
            results_parsed_dir = run_dir / "results_parsed"
            logs_dir = run_dir / "logs"

            for d in [targets_dir, replicas_dir, results_dir, results_parsed_dir, logs_dir]:
                safe_mkdir(d)

            # En opción 4A cada condición tiene su propio trg
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
                    "adcp_numSteps": defaults.get("adcp_numSteps", None),
                    "adcp_nbRuns": defaults.get("adcp_nbRuns", None),
                    "adcp_nmodes": defaults.get("adcp_nmodes", None),
                    "adcp_maxCores": defaults.get("adcp_maxCores", None),
                    "omm_nmin": defaults.get("omm_nmin", None),
                    "omm_max_itr": defaults.get("omm_max_itr", None),
                    "omm_environment": defaults.get("omm_environment", None),
                    "rerank_mode": defaults.get("rerank_mode", None),
                    "cluster_cutoff": defaults.get("cluster_cutoff", None),
                },
                "replicas": {
                    "planned_count": int(default_reps) if str(default_reps).isdigit() else default_reps,
                    "index_start": 1,
                    "suffix_width": 3,
                    "seed_base": None,
                    "note": "Adjust planned_count per campaign if needed (replicas may vary across campaigns)."
                },
                "paths": {
                    "run_dir": str(run_dir.as_posix()),
                    "targets_dir": str(targets_dir.as_posix()),
                    "replicas_dir": str(replicas_dir.as_posix()),
                    "results_dir": str(results_dir.as_posix()),
                    "results_parsed_dir": str(results_parsed_dir.as_posix()),
                    "logs_dir": str(logs_dir.as_posix()),
                },
                "status": "planned"
            }

            out_path = run_dir / "run_manifest.json"
            out_path.write_text(json.dumps(run_manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

            rows.append({
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
            })

    # escribe matriz
    fieldnames = [
        "campaign_id","protein_id","protein_name","protein_pdb",
        "peptide_id","peptide_set","peptide_seq",
        "receptor_pdb_file","receptor_pdb_path",
        "run_dir","target_trg","planned_replicas"
    ]
    with matrix_csv.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for r in rows:
            w.writerow(r)

    print(f"[ok] Created {len(rows)} campaigns.")
    print(f"[ok] Wrote campaign matrix: {matrix_csv}")
    print(f"[ok] Example run_manifest: {rows[0]['run_dir']}/run_manifest.json")

if __name__ == "__main__":
    main()
