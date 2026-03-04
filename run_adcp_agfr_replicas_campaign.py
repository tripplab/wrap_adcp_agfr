#!/usr/bin/env python3
"""
run_adcp_agfr_replicas_campaign.py

AGFR (AutoSite) + ADCP replicas pipeline with:

Mode B (run-manifest driven):
  - Provide --run-manifest runs/<protein_id>/<peptide_id>/run_manifest.json
  - All relative paths inside run_manifest.json resolve against exp_root
  - Per-condition targets: runs/<Pxxx>/<pep>/targets/receptor.{pdbqt,trg}
  - Replica outputs under runs/<Pxxx>/<pep>/replicas/replica_XXX/...
  - Timing artifacts under runs/<Pxxx>/<pep>/logs/timing.{json,csv}

Legacy mode (no --run-manifest):
  - Original CLI behavior with --receptor-pdb/--peptide-seq/--outdir etc.

Override policy (Mode B):
  - Operational flags ALWAYS follow CLI:
      --verbosity --real-time --validate-only --init-only --overwrite-reps --force-init --dry-run
  - Scientific/computational flags (prep/AGFR/ADCP/replicas) are governed by:
      * Without --allow-cli-overrides: if CLI explicitly passes any of those flags AND they disagree with the manifest -> ERROR+abort
      * With    --allow-cli-overrides: merge per-field (CLI if passed else manifest), print "override:" lines and record in run_manifest.runtime

Key design:
  - Overrideable flags are defined with argparse.SUPPRESS so they only exist in args if explicitly passed.
"""

import argparse
import csv
import datetime as _dt
import json
import os
import shutil
import subprocess
import sys
import time
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple


# --------------------------- utilities ---------------------------

def now_iso() -> str:
    return _dt.datetime.now().isoformat(timespec="seconds")


def ensure_dir(p: Path) -> None:
    p.mkdir(parents=True, exist_ok=True)


def which_or_none(cmd: str) -> Optional[str]:
    return shutil.which(cmd)


class PipelineError(RuntimeError):
    pass


def require_exists(path: Path, what: str) -> None:
    if not path.exists():
        raise PipelineError(f"Missing {what}: {path}")


def dump_json(path: Path, obj: Any) -> None:
    ensure_dir(path.parent)
    tmp = path.with_suffix(path.suffix + ".tmp")
    with tmp.open("w", encoding="utf-8") as f:
        json.dump(obj, f, indent=2, sort_keys=False)
        f.write("\n")
    tmp.replace(path)


def load_json(path: Path) -> Dict[str, Any]:
    require_exists(path, "JSON file")
    with path.open("r", encoding="utf-8") as f:
        return json.load(f)


def write_timing_csv(path: Path, rows: List["StageTiming"]) -> None:
    ensure_dir(path.parent)
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["stage", "start_iso", "end_iso", "elapsed_s", "returncode", "log_file", "command"])
        for r in rows:
            w.writerow([r.name, r.start_iso, r.end_iso, f"{r.elapsed_s:.3f}", r.returncode, r.log_file, " ".join(r.command)])


def validate_trg_zip(trg: Path) -> Tuple[bool, str]:
    import zipfile
    try:
        with zipfile.ZipFile(trg, "r") as z:
            names = z.namelist()
            if not names:
                return False, "TRG zip is empty"
            return True, f"zip ok; {len(names)} files; first={names[0]}"
    except Exception as e:
        return False, f"invalid zip: {e}"


@dataclass
class StageTiming:
    name: str
    start_iso: str
    end_iso: str
    elapsed_s: float
    returncode: int
    command: List[str]
    log_file: str


class TeeRunner:
    """
    Runs subprocesses, tees stdout/stderr to a log file, and optionally prints live to stdout.
    """
    def __init__(self, realtime: bool, verbosity: int):
        self.realtime = realtime
        self.verbosity = verbosity

    def run(
        self,
        cmd: List[str],
        log_path: Path,
        env: Optional[Dict[str, str]] = None,
        cwd: Optional[Path] = None,
        stage_label: str = ""
    ) -> int:
        ensure_dir(log_path.parent)
        with log_path.open("a", encoding="utf-8") as f:
            f.write(f"\n[{now_iso()}] CMD: {' '.join(cmd)}\n")

        p = subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            bufsize=1,
            universal_newlines=True,
            env=env,
            cwd=str(cwd) if cwd else None,
        )

        assert p.stdout is not None
        with log_path.open("a", encoding="utf-8") as f:
            for line in p.stdout:
                f.write(line)
                if self.realtime or self.verbosity >= 2:
                    if self.verbosity >= 3 and stage_label:
                        sys.stdout.write(f"[{stage_label}] {line}")
                    else:
                        sys.stdout.write(line)
                    sys.stdout.flush()

        p.wait()
        rc = int(p.returncode or 0)
        with log_path.open("a", encoding="utf-8") as f:
            f.write(f"[{now_iso()}] RETURN CODE: {rc}\n")
        return rc


def _log_info_line(log_path: Path, msg: str) -> None:
    ensure_dir(log_path.parent)
    with log_path.open("a", encoding="utf-8") as f:
        f.write(f"[{now_iso()}] {msg}\n")


def _print_and_log(runner: TeeRunner, log_path: Path, msg: str) -> None:
    _log_info_line(log_path, msg)
    if runner.realtime or runner.verbosity >= 1:
        print(msg)


def _cli_has(args: argparse.Namespace, key: str) -> bool:
    """
    For argparse.SUPPRESS options: key exists ONLY if user passed the flag.
    """
    return key in vars(args)


# -------------------- peptide sequence extraction --------------------

_THREE_TO_ONE: Dict[str, str] = {
    "ALA": "A", "CYS": "C", "ASP": "D", "GLU": "E", "PHE": "F",
    "GLY": "G", "HIS": "H", "ILE": "I", "LYS": "K", "LEU": "L",
    "MET": "M", "ASN": "N", "PRO": "P", "GLN": "Q", "ARG": "R",
    "SER": "S", "THR": "T", "VAL": "V", "TRP": "W", "TYR": "Y",
    "HID": "H", "HIE": "H", "HIP": "H", "HSD": "H", "HSE": "H", "HSP": "H",
    "CYX": "C", "CSS": "C",
    "MSE": "M",
    "SEC": "U",
    "PYL": "O",
}


def extract_sequence_from_pdb(pdb_path: Path) -> Tuple[str, Dict[str, Any]]:
    pdb_path = Path(pdb_path).resolve()
    require_exists(pdb_path, "peptide PDB (for sequence extraction)")

    seqres_by_chain: Dict[str, List[str]] = {}
    atom_ca_residues_by_chain: Dict[str, List[Tuple[Tuple[int, str], str]]] = {}
    unknown_resnames: Dict[str, int] = {}

    with pdb_path.open("r", encoding="utf-8", errors="replace") as f:
        for line in f:
            rec = line[:6].strip()
            if rec == "SEQRES":
                chain = (line[11] or " ").strip() or "_"
                parts = line[19:].split()
                if parts:
                    seqres_by_chain.setdefault(chain, []).extend(parts)
            elif rec == "ATOM":
                if len(line) < 54:
                    continue
                atom_name = line[12:16].strip()
                if atom_name != "CA":
                    continue
                resname = line[17:20].strip().upper()
                chain = (line[21] or " ").strip() or "_"
                try:
                    resseq = int(line[22:26])
                except Exception:
                    continue
                icode = (line[26] or " ").strip()
                key = (resseq, icode)
                atom_ca_residues_by_chain.setdefault(chain, []).append((key, resname))

    if seqres_by_chain:
        best_chain = max(seqres_by_chain, key=lambda c: len(seqres_by_chain[c]))
        resnames = seqres_by_chain[best_chain]
        seq_chars: List[str] = []
        for rn in resnames:
            rn = rn.strip().upper()
            aa = _THREE_TO_ONE.get(rn, "X")
            if aa == "X":
                unknown_resnames[rn] = unknown_resnames.get(rn, 0) + 1
            seq_chars.append(aa)
        seq = "".join(seq_chars)
        meta = {"method": "SEQRES", "chain": best_chain, "length": len(seq), "unknown_resnames": unknown_resnames}
        return seq, meta

    if atom_ca_residues_by_chain:
        best_chain = max(atom_ca_residues_by_chain, key=lambda c: len(atom_ca_residues_by_chain[c]))
        residues = atom_ca_residues_by_chain[best_chain]
        seen = set()
        ordered: List[Tuple[Tuple[int, str], str]] = []
        for key, rn in residues:
            if key in seen:
                continue
            seen.add(key)
            ordered.append((key, rn))
        ordered.sort(key=lambda x: (x[0][0], x[0][1]))

        seq_chars = []
        for _, rn in ordered:
            rn = rn.strip().upper()
            aa = _THREE_TO_ONE.get(rn, "X")
            if aa == "X":
                unknown_resnames[rn] = unknown_resnames.get(rn, 0) + 1
            seq_chars.append(aa)
        seq = "".join(seq_chars)
        meta = {"method": "ATOM_CA", "chain": best_chain, "length": len(seq), "unknown_resnames": unknown_resnames}
        return seq, meta

    raise PipelineError(f"Could not extract sequence from peptide PDB: {pdb_path}")


def _is_pathlike_peptide_seq(s: str) -> bool:
    s = s.strip()
    if not s:
        return False
    if os.path.exists(s):
        return True
    low = s.lower()
    return any(low.endswith(ext) for ext in [".pdb", ".pdbqt", ".mol2", ".cif"])


# --------------------------- Mode B helpers ---------------------------

def exp_root_from_run_manifest(run_manifest_path: Path) -> Path:
    """
    run_manifest lives at: <exp_root>/runs/<protein_id>/<peptide_id>/run_manifest.json
    So exp_root is four parents up from the manifest file.
    """
    p = run_manifest_path.resolve()
    try:
        return p.parent.parent.parent.parent
    except Exception:
        raise PipelineError(f"Could not infer exp_root from run_manifest path: {run_manifest_path}")


def campaign_dir_from_run_manifest(run_manifest_path: Path) -> Path:
    return run_manifest_path.resolve().parent


def assert_manifest_path_matches_campaign(run_manifest_path: Path, campaign: Dict[str, Any]) -> None:
    """
    Strict rule: path must match runs/<protein_id>/<peptide_id>/run_manifest.json
    """
    protein_id = str(campaign.get("protein_id", "")).strip()
    peptide_id = str(campaign.get("peptide_id", "")).strip()
    if not protein_id or not peptide_id:
        raise PipelineError("run_manifest.json missing required campaign.protein_id and/or campaign.peptide_id")

    p = run_manifest_path.resolve()
    if p.name != "run_manifest.json":
        raise PipelineError(f"Expected run_manifest filename to be run_manifest.json, got: {p.name}")
    if p.parent.name != peptide_id or p.parent.parent.name != protein_id or p.parent.parent.parent.name != "runs":
        raise PipelineError(
            "run_manifest path mismatch:\n"
            f"  manifest: {p}\n"
            f"  expected folder chain: .../runs/{protein_id}/{peptide_id}/run_manifest.json"
        )


def resolve_path_maybe_relative(p: Optional[str], base: Path) -> Optional[str]:
    if p is None:
        return None
    p = str(p).strip()
    if not p:
        return None
    pp = Path(p)
    if pp.is_absolute():
        return str(pp)
    return str((base / pp).resolve())


def resolve_run_manifest_paths(rm: Dict[str, Any], exp_root: Path) -> Dict[str, Any]:
    """
    Adds rm['resolved_paths'] computed by resolving ALL relative paths against exp_root (contract).
    Does not rewrite original fields; it adds a 'resolved_paths' block.
    """
    out = dict(rm)
    inputs = dict(out.get("inputs", {}))
    targets = dict(out.get("targets", {}))
    paths = dict(out.get("paths", {}))

    resolved: Dict[str, Any] = {"exp_root": str(exp_root.resolve())}

    resolved["receptor_pdb_path"] = resolve_path_maybe_relative(inputs.get("receptor_pdb_path"), exp_root)
    resolved["peptide_pdb_path"] = resolve_path_maybe_relative(inputs.get("peptide_pdb_path"), exp_root)
    resolved["target_trg"] = resolve_path_maybe_relative(targets.get("target_trg"), exp_root)

    resolved["run_dir"] = resolve_path_maybe_relative(paths.get("run_dir"), exp_root)
    resolved["targets_dir"] = resolve_path_maybe_relative(paths.get("targets_dir"), exp_root)
    resolved["replicas_dir"] = resolve_path_maybe_relative(paths.get("replicas_dir"), exp_root)
    resolved["results_dir"] = resolve_path_maybe_relative(paths.get("results_dir"), exp_root)
    resolved["results_parsed_dir"] = resolve_path_maybe_relative(paths.get("results_parsed_dir"), exp_root)
    resolved["logs_dir"] = resolve_path_maybe_relative(paths.get("logs_dir"), exp_root)

    out["resolved_paths"] = resolved
    return out


def mode_b_defaults_for_timing(run_manifest: Dict[str, Any], exp_root: Path) -> Tuple[Path, Path]:
    rm2 = resolve_run_manifest_paths(run_manifest, exp_root)
    logs_dir = rm2.get("resolved_paths", {}).get("logs_dir")
    if not logs_dir:
        raise PipelineError("run_manifest missing paths.logs_dir (needed to place timing files)")
    ld = Path(logs_dir).resolve()
    return (ld / "timing.json"), (ld / "timing.csv")


# --------------------------- override mapping (Mode B) ---------------------------

# Fields that are allowed to be overridden ONLY if --allow-cli-overrides is present.
# Keys are *argparse dest names*; values are manifest paths to compare/apply.

_OVERRIDE_MAP_REPLICAS: Dict[str, str] = {
    "replicas": "replicas.planned_count",
    "replica_seed_base": "replicas.seed_base",
    "replica_index_start": "replicas.index_start",
    "replica_suffix_width": "replicas.suffix_width",
}

_OVERRIDE_MAP_PREP: Dict[str, str] = {
    "prep_receptor_A": "config.prep_receptor_A",
    "prep_receptor_U": "config.prep_receptor_U",
    "prep_receptor_v": "config.prep_receptor_v",
    "prep_receptor_C": "config.prep_receptor_C",
    "prep_receptor_e": "config.prep_receptor_e",
    "prep_receptor_w": "config.prep_receptor_w",
    # prep_receptor_p is repeatable; treat as list if passed
    "prep_receptor_p": "config.prep_receptor_p",
}

_OVERRIDE_MAP_AGFR: Dict[str, str] = {
    "agfr_spacing": "config.agfr_spacing",
    "agfr_smoothing": "config.agfr_smoothing",
    "agfr_autoSiteV": "config.agfr_autoSiteV",
    "agfr_ligandSize": "config.agfr_ligandSize",
    "agfr_mapTypes": "config.agfr_mapTypes",
    "agfr_pepScore": "config.agfr_pepScore",
    "agfr_padding": "config.agfr_padding",
    "agfr_pocketCutoff": "config.agfr_pocketCutoff",
    # NOTE: boxMode/pocketMode are lists; we allow override if explicitly passed
    "agfr_boxMode": "config.agfr_boxMode",
    "agfr_pocketMode": "config.agfr_pocketMode",
}

_OVERRIDE_MAP_ADCP: Dict[str, str] = {
    "adcp_numSteps": "config.adcp_numSteps",
    "adcp_nbRuns": "config.adcp_nbRuns",
    "adcp_maxCores": "config.adcp_maxCores",
    "adcp_nmodes": "config.adcp_nmodes",
    "adcp_omm_nmin": "config.omm_nmin",
    "adcp_omm_max_itr": "config.omm_max_itr",
    "adcp_omm_environment": "config.omm_environment",
    "adcp_rerank_mode": "config.rerank_mode",
    # This one is special: -O in adcp (overwrite output files)
    "adcp_overwriteFiles": "config.adcp_overwriteFiles",
}

_OVERRIDE_MAP: Dict[str, str] = {}
_OVERRIDE_MAP.update(_OVERRIDE_MAP_REPLICAS)
_OVERRIDE_MAP.update(_OVERRIDE_MAP_PREP)
_OVERRIDE_MAP.update(_OVERRIDE_MAP_AGFR)
_OVERRIDE_MAP.update(_OVERRIDE_MAP_ADCP)


def _get_manifest_value(rm: Dict[str, Any], dotted: str) -> Any:
    cur: Any = rm
    for part in dotted.split("."):
        if not isinstance(cur, dict):
            return None
        cur = cur.get(part)
    return cur


def _set_manifest_value(rm: Dict[str, Any], dotted: str, value: Any) -> None:
    parts = dotted.split(".")
    cur: Any = rm
    for p in parts[:-1]:
        if p not in cur or not isinstance(cur[p], dict):
            cur[p] = {}
        cur = cur[p]
    cur[parts[-1]] = value


def collect_cli_overrides_mode_b(args: argparse.Namespace) -> Dict[str, Any]:
    """
    Returns a dict of {manifest_path: cli_value} only for overrideable flags that were explicitly passed.
    """
    overrides: Dict[str, Any] = {}
    for dest, mpath in _OVERRIDE_MAP.items():
        if _cli_has(args, dest):
            overrides[mpath] = getattr(args, dest)
    return overrides


def compute_manifest_cli_discrepancies(rm: Dict[str, Any], overrides: Dict[str, Any]) -> List[Tuple[str, Any, Any]]:
    """
    Returns list of (manifest_path, manifest_value, cli_value) where they differ (normalized).
    """
    diffs: List[Tuple[str, Any, Any]] = []
    for mpath, cli_val in overrides.items():
        man_val = _get_manifest_value(rm, mpath)
        # Normalize a bit: int/str comparisons for numbers
        if isinstance(man_val, (int, float)) and isinstance(cli_val, str):
            try:
                cli_val_n = float(cli_val) if "." in cli_val else int(cli_val)
                cli_val = cli_val_n
            except Exception:
                pass
        if isinstance(cli_val, (int, float)) and isinstance(man_val, str):
            try:
                man_val_n = float(man_val) if "." in man_val else int(man_val)
                man_val = man_val_n
            except Exception:
                pass
        if man_val != cli_val:
            diffs.append((mpath, man_val, cli_val))
    return diffs


def apply_cli_overrides_to_manifest(rm: Dict[str, Any], overrides: Dict[str, Any]) -> None:
    for mpath, cli_val in overrides.items():
        _set_manifest_value(rm, mpath, cli_val)


# --------------------------- CLI ---------------------------

def parse_args(argv: Optional[List[str]] = None) -> argparse.Namespace:
    ep = """
Examples (Mode B):
  Validate only:
    run_adcp_agfr_replicas_campaign.py --run-manifest runs/P001/POI/run_manifest.json --validate-only

  Init only (build/validate targets only):
    run_adcp_agfr_replicas_campaign.py --run-manifest runs/P001/POI/run_manifest.json --init-only --real-time

  Run replicas (preserve targets):
    run_adcp_agfr_replicas_campaign.py --run-manifest runs/P001/POI/run_manifest.json --replicas 2 --replica-seed-base 1000 --allow-cli-overrides

Notes:
  - Operational flags always follow CLI.
  - Scientific flags passed on CLI are rejected unless --allow-cli-overrides is provided.
"""
    p = argparse.ArgumentParser(
        prog="run_adcp_agfr_replicas_campaign.py",
        description="AGFR (AutoSite) + ADCP replicas pipeline (Mode B via run_manifest) with strict manifest/CLI policy.",
        epilog=ep,
        formatter_class=argparse.RawTextHelpFormatter,
    )

    # -------------------- Mode selection --------------------
    g_mode = p.add_argument_group("Mode selection")
    g_mode.add_argument("--run-manifest", default=None,
                        help="Mode B: path to runs/<protein_id>/<peptide_id>/run_manifest.json")
    g_mode.add_argument("--allow-cli-overrides", action="store_true",
                        help="Mode B: allow scientific/computational overrides from CLI (merge per-field).")

    # -------------------- Operational (always CLI) --------------------
    g_op = p.add_argument_group("Operational flags (always follow CLI)")
    g_op.add_argument("--validate-only", action="store_true",
                      help="Mode B: validate inputs/paths/targets and exit (no execution).")
    g_op.add_argument("--init-only", action="store_true",
                      help="Mode B: ensure dirs + validate/build targets, then stop (skip replicas).")
    g_op.add_argument("--overwrite-reps", action="store_true",
                      help="Mode B: delete replicas/ + results/ and clean logs/ (targets preserved; init not forced).")
    g_op.add_argument("--force-init", action="store_true",
                      help="Mode B: force rebuilding targets even if receptor.pdbqt and receptor.trg exist.")
    g_op.add_argument("--dry-run", action="store_true",
                      help="Print commands but do not execute.")
    g_op.add_argument("--real-time", action="store_true",
                      help="Stream tool output to stdout in real time (also logged).")
    g_op.add_argument("--verbosity", type=int, choices=[0, 1, 2, 3], default=1,
                      help="0 quiet, 1 normal, 2 verbose, 3 debug")

    # Back-compat / deprecation
    g_op.add_argument("--overwrite", action="store_true",
                      help="DEPRECATED (Mode B): treated as --overwrite-reps. (Legacy): allows outdir overwrite.")

    # -------------------- Legacy core inputs --------------------
    g_leg = p.add_argument_group("Legacy mode inputs (only used if --run-manifest is NOT provided)")
    g_leg.add_argument("--receptor-pdb", default=None, help="(Legacy) input receptor PDB")
    g_leg.add_argument("--peptide-seq", default=None,
                       help="(Legacy) peptide sequence (1-letter). If omitted and --peptide-pdb is provided, auto-extract.")
    g_leg.add_argument("--peptide-pdb", default=None, help="(Legacy) peptide PDB to seed ADCP (-i). Optional.")
    g_leg.add_argument("--outdir", default=None, help="(Legacy) output directory")

    g_leg.add_argument("--reuse-target-trg", default=None, help="(Legacy) path to an existing .trg to reuse (skip AGFR)")
    g_leg.add_argument("--skip-agfr", action="store_true", help="(Legacy) skip AGFR if target exists at default location")

    # -------------------- Scientific/computational overrides (Mode B requires allow) --------------------
    # IMPORTANT: these use SUPPRESS to detect "flag really passed".
    g_rep = p.add_argument_group("Replicas (scientific; Mode B requires --allow-cli-overrides)")
    g_rep.add_argument("--replicas", dest="replicas", type=int, default=argparse.SUPPRESS,
                       help="Override replicas.planned_count")
    g_rep.add_argument("--replica-seed-base", dest="replica_seed_base", type=int, default=argparse.SUPPRESS,
                       help="Override replicas.seed_base (replica i uses seed_base+(i-index_start))")
    g_rep.add_argument("--replica-index-start", dest="replica_index_start", type=int, default=argparse.SUPPRESS,
                       help="Override replicas.index_start")
    g_rep.add_argument("--replica-suffix-width", dest="replica_suffix_width", type=int, default=argparse.SUPPRESS,
                       help="Override replicas.suffix_width (zero padding width)")

    g_prep = p.add_argument_group("prepare_receptor (scientific; Mode B requires --allow-cli-overrides)")
    g_prep.add_argument("--prep-receptor-A", dest="prep_receptor_A", default=argparse.SUPPRESS,
                        choices=["bonds_hydrogens", "bonds", "hydrogens", "checkhydrogens", "None"],
                        help="Override config.prep_receptor_A (prepare_receptor -A)")
    g_prep.add_argument("--prep-receptor-U", dest="prep_receptor_U", default=argparse.SUPPRESS,
                        help="Override config.prep_receptor_U (prepare_receptor -U)")
    g_prep.add_argument("--prep-receptor-v", dest="prep_receptor_v", action="store_true", default=argparse.SUPPRESS,
                        help="Override config.prep_receptor_v (verbose)")
    g_prep.add_argument("--prep-receptor-C", dest="prep_receptor_C", action="store_true", default=argparse.SUPPRESS,
                        help="Override config.prep_receptor_C (preserve charges)")
    g_prep.add_argument("--prep-receptor-e", dest="prep_receptor_e", action="store_true", default=argparse.SUPPRESS,
                        help="Override config.prep_receptor_e (delete non-standard residues)")
    g_prep.add_argument("--prep-receptor-w", dest="prep_receptor_w", action="store_true", default=argparse.SUPPRESS,
                        help="Override config.prep_receptor_w (unique atom names)")
    g_prep.add_argument("--prep-receptor-p", dest="prep_receptor_p", action="append", default=argparse.SUPPRESS,
                        help="Override config.prep_receptor_p (repeatable -p atomType), e.g. --prep-receptor-p Zn")

    g_ag = p.add_argument_group("AGFR target generation (scientific; Mode B requires --allow-cli-overrides)")
    g_ag.add_argument("--agfr-spacing", dest="agfr_spacing", type=float, default=argparse.SUPPRESS,
                      help="Override config.agfr_spacing (-s)")
    g_ag.add_argument("--agfr-smoothing", dest="agfr_smoothing", type=float, default=argparse.SUPPRESS,
                      help="Override config.agfr_smoothing (-S)")
    g_ag.add_argument("--agfr-autoSiteV", dest="agfr_autoSiteV", default=argparse.SUPPRESS,
                      help="Override config.agfr_autoSiteV (-asv)")
    g_ag.add_argument("--agfr-ligandSize", dest="agfr_ligandSize", type=int, default=argparse.SUPPRESS,
                      help="Override config.agfr_ligandSize (-ls)")
    g_ag.add_argument("--agfr-mapTypes", dest="agfr_mapTypes", default=argparse.SUPPRESS, choices=["all", "ligand"],
                      help="Override config.agfr_mapTypes (-m)")
    g_ag.add_argument("--agfr-pepScore", dest="agfr_pepScore", action="store_true", default=argparse.SUPPRESS,
                      help="Override config.agfr_pepScore (-ps)")
    g_ag.add_argument("--no-agfr-pepScore", dest="agfr_pepScore", action="store_false", default=argparse.SUPPRESS,
                      help="Override config.agfr_pepScore (disable -ps)")
    g_ag.add_argument("--agfr-padding", dest="agfr_padding", type=float, default=argparse.SUPPRESS,
                      help="Override config.agfr_padding (-P)")
    g_ag.add_argument("--agfr-pocketCutoff", dest="agfr_pocketCutoff", type=int, default=argparse.SUPPRESS,
                      help="Override config.agfr_pocketCutoff (-C)")
    g_ag.add_argument("--agfr-boxMode", dest="agfr_boxMode", nargs="+", default=argparse.SUPPRESS,
                      help="Override config.agfr_boxMode (pass tokens after -b)")
    g_ag.add_argument("--agfr-pocketMode", dest="agfr_pocketMode", nargs="+", default=argparse.SUPPRESS,
                      help="Override config.agfr_pocketMode (pass tokens after -p)")

    g_ad = p.add_argument_group("ADCP docking (scientific; Mode B requires --allow-cli-overrides)")
    g_ad.add_argument("--adcp-numSteps", dest="adcp_numSteps", type=int, default=argparse.SUPPRESS,
                      help="Override config.adcp_numSteps (-n)")
    g_ad.add_argument("--adcp-nbRuns", dest="adcp_nbRuns", type=int, default=argparse.SUPPRESS,
                      help="Override config.adcp_nbRuns (-N)")
    g_ad.add_argument("--adcp-maxCores", dest="adcp_maxCores", type=int, default=argparse.SUPPRESS,
                      help="Override config.adcp_maxCores (-c)")
    g_ad.add_argument("--adcp-nmodes", dest="adcp_nmodes", type=int, default=argparse.SUPPRESS,
                      help="Override config.adcp_nmodes (-m)")
    g_ad.add_argument("--adcp-omm-nmin", dest="adcp_omm_nmin", type=int, default=argparse.SUPPRESS,
                      help="Override config.omm_nmin (-nmin), MUST be >0")
    g_ad.add_argument("--adcp-omm-max-itr", dest="adcp_omm_max_itr", type=int, default=argparse.SUPPRESS,
                      help="Override config.omm_max_itr (-nitr)")
    g_ad.add_argument("--adcp-omm-environment", dest="adcp_omm_environment", default=argparse.SUPPRESS,
                      choices=["vacuum", "implicit"], help="Override config.omm_environment (-env)")
    g_ad.add_argument("--adcp-rerank-mode", dest="adcp_rerank_mode", default=argparse.SUPPRESS,
                      choices=["complex", "interaction"], help="Override config.rerank_mode")
    g_ad.add_argument("--adcp-overwriteFiles", dest="adcp_overwriteFiles", action="store_true", default=argparse.SUPPRESS,
                      help="Override config.adcp_overwriteFiles (pass -O to adcp to overwrite existing output files)")

    return p.parse_args(argv)


# ---- END PART 1/2 ----

# ---- PART 2/2 ----

# --------------------------- command builders ---------------------------

def build_prepare_receptor_cmd(cfg: Dict[str, Any]) -> List[str]:
    pr = cfg["prepare_receptor"]
    cmd = ["prepare_receptor", "-r", cfg["receptor_pdb"], "-o", cfg["receptor_pdbqt"]]

    if pr.get("v"):
        cmd.append("-v")
    A = pr.get("A")
    if A is not None:
        if A == "None":
            cmd.extend(["-A", "None"])
        else:
            cmd.extend(["-A", str(A)])
    if pr.get("C"):
        cmd.append("-C")
    for at in pr.get("p", []) or []:
        cmd.extend(["-p", str(at)])
    U = pr.get("U")
    if U:
        cmd.extend(["-U", str(U)])
    if pr.get("e"):
        cmd.extend(["-e", "True"])
    if pr.get("M"):
        cmd.extend(["-M", "interactive"])
    if pr.get("d"):
        cmd.extend(["-d", str(pr["d"])])
    if pr.get("w"):
        cmd.append("-w")

    return cmd


def build_agfr_cmd(cfg: Dict[str, Any]) -> List[str]:
    ag = cfg["agfr"]
    cmd = ["agfr"]

    # Minimal default: receptor pdbqt
    if ag.get("config"):
        cmd.extend(["-F", ag["config"]])
    elif ag.get("inputFolder"):
        cmd.extend(["--inputFolder", ag["inputFolder"]])
    else:
        cmd.extend(["-r", cfg["receptor_pdbqt"]])

    if ag.get("ligand"):
        cmd.extend(["-l", ag["ligand"]])

    cmd.extend(["-o", cfg["target_trg"]])

    if ag.get("boxMode"):
        cmd.extend(["-b"] + [str(x) for x in ag["boxMode"]])
    if ag.get("padding") is not None:
        cmd.extend(["-P", str(ag["padding"])])
    if ag.get("spacing") is not None:
        cmd.extend(["-s", str(ag["spacing"])])
    if ag.get("smoothing") is not None:
        cmd.extend(["-S", str(ag["smoothing"])])

    if ag.get("flexRes"):
        cmd.extend(["-f", ag["flexRes"]])
    if ag.get("pocketMode"):
        cmd.extend(["-p"] + [str(x) for x in ag["pocketMode"]])
    if ag.get("pocketCutoff") is not None:
        cmd.extend(["-C", str(ag["pocketCutoff"])])

    if ag.get("mapTypes"):
        cmd.extend(["-m", str(ag["mapTypes"])])
    if ag.get("autoSiteV") is not None:
        cmd.extend(["-asv", str(ag["autoSiteV"])])
    if ag.get("ligandSize") is not None:
        cmd.extend(["-ls", str(ag["ligandSize"])])
    if ag.get("pepScore") is True:
        cmd.append("-ps")

    # Keep additional pass-through fields available (already in cfg if legacy used),
    # but don't add a huge parameter surface here unless needed.
    if ag.get("recGradient"):
        cmd.append("-g")
    if ag.get("constrictionMapCutoff") is not None:
        cmd.extend(["-cm", str(ag["constrictionMapCutoff"])])
    if ag.get("recGradVolCut") is not None:
        cmd.extend(["-V", str(ag["recGradVolCut"])])

    if ag.get("waterWeight") is not None:
        cmd.extend(["--waterWeight", str(ag["waterWeight"])])
    if ag.get("waterEntropy") is not None:
        cmd.extend(["--waterEntropy", str(ag["waterEntropy"])])
    if ag.get("ionCharges"):
        cmd.extend(["--ionCharges", str(ag["ionCharges"])])

    if ag.get("pdbLigand"):
        cmd.extend(["--pdbLigand", str(ag["pdbLigand"])])
    if ag.get("bioMol"):
        cmd.extend(["--bioMol", str(ag["bioMol"])])
    if ag.get("addHlig"):
        cmd.extend(["--addHlig", str(ag["addHlig"])])
    if ag.get("protToRemove"):
        cmd.extend(["--protToRemove", str(ag["protToRemove"])])
    if ag.get("heteroToKeep"):
        cmd.extend(["--heteroToKeep", str(ag["heteroToKeep"])])
    if ag.get("mutate"):
        cmd.extend(["--mutate", str(ag["mutate"])])
    if ag.get("onException"):
        cmd.extend(["--onException", str(ag["onException"])])
    if ag.get("debug") is not None:
        cmd.extend(["-d", str(ag["debug"])])
    if ag.get("toPdbqt"):
        cmd.append("--toPdbqt")
    if ag.get("testBatch"):
        cmd.append("--testBatch")
    if ag.get("saveFilename"):
        cmd.extend(["--saveFilename", str(ag["saveFilename"])])

    if ag.get("covalentBond") is not None:
        a, b = ag["covalentBond"]
        cmd.extend(["-c", str(a), str(b)])
    if ag.get("covalentBondTorsionAtom") is not None:
        cmd.extend(["-T", str(ag["covalentBondTorsionAtom"])])
    if ag.get("covalentResidues"):
        cmd.extend(["-x", str(ag["covalentResidues"])])

    return cmd


def build_adcp_cmd(ad: Dict[str, Any]) -> List[str]:
    cmd = ["adcp"]

    if ad.get("sequence"):
        cmd.extend(["-s", str(ad["sequence"])])
    if ad.get("input"):
        cmd.extend(["-i", str(ad["input"])])
    if (not ad.get("sequence")) and (not ad.get("input")):
        raise PipelineError("Need peptide input: provide peptide sequence and/or peptide PDB.")

    if ad.get("partition") is not None:
        cmd.extend(["-p", str(ad["partition"])])

    cmd.extend(["-T", str(ad["target"])])
    cmd.extend(["-w", str(ad["workingFolder"])])
    cmd.extend(["-o", str(ad["jobName"])])

    if ad.get("rotlibs"):
        cmd.extend(["-L", str(ad["rotlibs"])])
    if ad.get("userRotlibs"):
        cmd.extend(["-l", str(ad["userRotlibs"])])

    if ad.get("numSteps") is not None:
        cmd.extend(["-n", str(ad["numSteps"])])
    if ad.get("nbRuns") is not None:
        cmd.extend(["-N", str(ad["nbRuns"])])
    if ad.get("maxCores") is not None:
        cmd.extend(["-c", str(ad["maxCores"])])

    if ad.get("keepWorkingFolder"):
        cmd.append("-k")
    if ad.get("dryRun"):
        cmd.append("-y")
    if ad.get("cyclic"):
        cmd.append("-cyc")
    if ad.get("cystein"):
        cmd.append("-cys")
    if ad.get("overwriteFiles"):
        cmd.append("-O")
    if ad.get("seed") is not None:
        cmd.extend(["-S", str(ad["seed"])])

    if ad.get("natContacts") is not None:
        cmd.extend(["-nc", str(ad["natContacts"])])
    if ad.get("rmsd") is not None:
        cmd.extend(["-rmsd", str(ad["rmsd"])])
    if ad.get("ref") is not None:
        cmd.extend(["-ref", str(ad["ref"])])
    if ad.get("nmodes") is not None:
        cmd.extend(["-m", str(ad["nmodes"])])
    if ad.get("reclusterOnly"):
        cmd.append("-r")

    # Minimization required
    if ad.get("omm_nmin") is None or int(ad["omm_nmin"]) <= 0:
        raise PipelineError("Docking must do minimization: set omm_nmin > 0.")
    cmd.extend(["-nmin", str(int(ad["omm_nmin"]))])

    if ad.get("dockingRanking"):
        cmd.append("-dr")
    if ad.get("omm_max_itr") is not None:
        cmd.extend(["-nitr", str(ad["omm_max_itr"])])
    if ad.get("omm_environment"):
        cmd.extend(["-env", str(ad["omm_environment"])])
    if ad.get("fix_nst"):
        cmd.append("-fnst")
    if ad.get("postDockMinimize"):
        cmd.append("-pdmin")
    if ad.get("ResumeMinimize"):
        cmd.append("-resumemin")

    rerank_mode = ad.get("rerank_mode", "interaction")
    if rerank_mode == "interaction":
        cmd.append("-reint")

    if ad.get("ffxmlset"):
        cmd.extend(["-F", str(ad["ffxmlset"])])
    if ad.get("userffxml"):
        cmd.extend(["-f", str(ad["userffxml"])])

    return cmd


# --------------------------- replica helpers ---------------------------

def _replica_label(i: int, width: int) -> str:
    return str(i).zfill(width)


def _replica_paths_mode_b(replicas_dir: Path, rep_label: str) -> Dict[str, Path]:
    rep_root = (replicas_dir / f"replica_{rep_label}").resolve()
    return {
        "replica_root": rep_root,
        "logs_dir": rep_root / "logs",
        "docking_work": rep_root / "docking" / "work",
        "results_dir": rep_root / "results",
    }


def _collect_adcp_outputs(jobname: str, search_roots: List[Path]) -> Dict[str, Optional[str]]:
    candidates: List[Path] = []
    for root in search_roots:
        if root and root.exists():
            candidates.extend(list(root.glob(f"{jobname}*")))

    out_pdb = None
    summary = None
    for c in candidates:
        if c.is_file() and c.name.endswith("_out.pdb") and out_pdb is None:
            out_pdb = str(c.resolve())
        if c.is_file() and c.name.endswith("_summary.dlg") and summary is None:
            summary = str(c.resolve())
    return {"out_pdb": out_pdb, "summary_dlg": summary}


def _copy_if_present(src: Optional[str], dst_dir: Path) -> Optional[str]:
    if not src:
        return None
    sp = Path(src)
    if not sp.exists():
        return None
    ensure_dir(dst_dir)
    dest = (dst_dir / sp.name).resolve()
    if dest.exists():
        dest.unlink()
    shutil.copy2(sp, dest)
    return str(dest)


# --------------------------- robust tool error detection ---------------------------

def detect_tool_failure(stage_name: str, rc: int, log_path: Path) -> Tuple[bool, str]:
    """
    Detect "soft failures" where tool prints ERROR but returns rc=0 (observed in ADCP).
    Returns (failed, message).
    """
    if rc != 0:
        return True, f"{stage_name}: non-zero return code rc={rc}"

    try:
        txt = log_path.read_text(encoding="utf-8", errors="replace")
    except Exception as e:
        return False, f"{stage_name}: rc=0 and could not read log for extra checks: {e}"

    # Tool-specific heuristics
    if stage_name.startswith("adcp"):
        # ADCP sometimes prints ERROR and exits 0.
        needles = [
            "ERROR:",
            "Error:",
            "exists! please use different jobName",
            "Traceback (most recent call last):",
        ]
        for n in needles:
            if n in txt:
                # allow benign "Warning:" lines
                return True, f"{stage_name}: log contains failure marker: {n}"

    if stage_name == "agfr":
        for n in ["ERROR", "Traceback", "Exception"]:
            if n in txt:
                return True, f"{stage_name}: log contains failure marker: {n}"

    if stage_name == "prepare_receptor":
        for n in ["ERROR", "Traceback", "Exception"]:
            if n in txt:
                return True, f"{stage_name}: log contains failure marker: {n}"

    return False, f"{stage_name}: rc=0 and no failure markers detected"


# --------------------------- Mode B: config resolution with override policy ---------------------------

def resolve_cfg_mode_b(
    args: argparse.Namespace,
    run_manifest_path: Path,
    runner: TeeRunner,
    preflight_log: Path,
) -> Tuple[Dict[str, Any], Dict[str, Any], Path, Path, Path, Path]:
    """
    Load run_manifest, enforce strict coherence, apply override policy, and build cfg.
    Returns: (cfg, rm_out, exp_root, campaign_dir, timing_json, timing_csv)
    """
    rm = load_json(run_manifest_path)
    rm["__run_manifest_path__"] = str(run_manifest_path.resolve())

    exp_root = exp_root_from_run_manifest(run_manifest_path)
    campaign_dir = campaign_dir_from_run_manifest(run_manifest_path)

    campaign = rm.get("campaign", {})
    assert_manifest_path_matches_campaign(run_manifest_path, campaign)

    # Resolve paths block
    rm2 = resolve_run_manifest_paths(rm, exp_root)
    resolved = rm2["resolved_paths"]

    # Mandatory inputs
    receptor_pdb_path = resolved.get("receptor_pdb_path")
    if not receptor_pdb_path:
        raise PipelineError("run_manifest.inputs.receptor_pdb_path is required (relative ok; resolves vs exp_root)")
    receptor_pdb = Path(receptor_pdb_path).resolve()

    peptide_seq = rm2.get("inputs", {}).get("peptide_seq")
    peptide_pdb_path = resolved.get("peptide_pdb_path")
    peptide_pdb = Path(peptide_pdb_path).resolve() if peptide_pdb_path else None

    if (not peptide_seq) and (peptide_pdb is None):
        raise PipelineError("Mode B: need peptide input: inputs.peptide_seq and/or inputs.peptide_pdb_path")

    if peptide_seq and _is_pathlike_peptide_seq(str(peptide_seq)):
        raise PipelineError(
            f"Mode B: inputs.peptide_seq looks like a file path ({peptide_seq}). "
            "Put it in inputs.peptide_pdb_path instead."
        )

    # Mandatory layout paths
    for k in ["run_dir", "targets_dir", "replicas_dir", "results_dir", "logs_dir"]:
        if not resolved.get(k):
            raise PipelineError(f"Mode B: run_manifest.paths.{k} is required (relative ok)")

    run_dir_p = Path(resolved["run_dir"]).resolve()
    targets_dir_p = Path(resolved["targets_dir"]).resolve()
    replicas_dir_p = Path(resolved["replicas_dir"]).resolve()
    results_dir_p = Path(resolved["results_dir"]).resolve()
    logs_dir_p = Path(resolved["logs_dir"]).resolve()

    if run_dir_p != campaign_dir.resolve():
        raise PipelineError(
            "Mode B: paths.run_dir does not match folder containing run_manifest.json\n"
            f"  run_manifest folder: {campaign_dir.resolve()}\n"
            f"  resolved paths.run_dir: {run_dir_p}"
        )

    target_trg = resolved.get("target_trg")
    if not target_trg:
        raise PipelineError("Mode B: targets.target_trg is required")
    target_trg_p = Path(target_trg).resolve()

    # ---------------- override policy ----------------
    cli_overrides = collect_cli_overrides_mode_b(args)
    diffs = compute_manifest_cli_discrepancies(rm2, cli_overrides)

    if diffs and not args.allow_cli_overrides:
        _print_and_log(runner, preflight_log, "[preflight] ERROR: CLI/manifest discrepancies detected (scientific flags) without --allow-cli-overrides")
        for mpath, man, cli in diffs:
            _print_and_log(runner, preflight_log, f"[preflight] mismatch: {mpath}: manifest={man!r} cli={cli!r}")
        raise PipelineError("Mode B: refusing to run with scientific CLI flags that disagree with manifest (use --allow-cli-overrides).")

    override_records: List[Dict[str, Any]] = []
    if diffs and args.allow_cli_overrides:
        # Apply override values to manifest copy (rm2)
        for mpath, man, cli in diffs:
            override_records.append({
                "path": mpath,
                "from_manifest": man,
                "to_cli": cli,
            })
            _print_and_log(runner, preflight_log, f"[preflight] override: {mpath} {man} -> {cli} (from CLI)")
        apply_cli_overrides_to_manifest(rm2, cli_overrides)

    # ---------------- build cfg from (possibly overridden) manifest ----------------
    cfg: Dict[str, Any] = {}
    cfg["mode"] = "B"
    cfg["exp_root"] = str(exp_root.resolve())
    cfg["campaign_dir"] = str(campaign_dir.resolve())

    cfg["outdir"] = str(run_dir_p)
    cfg["targets_dir"] = str(targets_dir_p)
    cfg["replicas_root"] = str(replicas_dir_p)
    cfg["results_dir"] = str(results_dir_p)
    cfg["logs_dir"] = str(logs_dir_p)

    cfg["receptor_pdb"] = str(receptor_pdb)

    # Default outputs under targets/
    cfg["receptor_pdbqt"] = str((targets_dir_p / "receptor.pdbqt").resolve())
    cfg["target_trg"] = str(target_trg_p)

    cfg["peptide_seq"] = str(peptide_seq) if peptide_seq else None
    cfg["peptide_pdb"] = str(peptide_pdb) if peptide_pdb else None

    # Pull scientific config from rm2["config"] (post-override)
    config = rm2.get("config", {})

    # prepare_receptor
    cfg["prepare_receptor"] = {
        "A": config.get("prep_receptor_A", "checkhydrogens"),
        "U": config.get("prep_receptor_U", "nphs_lps_waters_nonstdres"),
        "v": bool(config.get("prep_receptor_v", False)),
        "C": bool(config.get("prep_receptor_C", False)),
        "p": list(config.get("prep_receptor_p", [])) if config.get("prep_receptor_p") is not None else [],
        "e": bool(config.get("prep_receptor_e", False)),
        "M": bool(config.get("prep_receptor_M", False)),
        "d": config.get("prep_receptor_d"),
        "w": bool(config.get("prep_receptor_w", False)),
    }

    # AGFR
    cfg["agfr"] = {
        "config": config.get("agfr_config"),
        "inputFolder": config.get("agfr_inputFolder"),
        "ligand": config.get("agfr_ligand"),
        "boxMode": config.get("agfr_boxMode", ["receptor"]),
        "padding": config.get("agfr_padding"),
        "spacing": config.get("agfr_spacing", 0.375),
        "smoothing": config.get("agfr_smoothing", 0.5),
        "flexRes": config.get("agfr_flexRes"),
        "pocketMode": config.get("agfr_pocketMode", ["best"]),
        "pocketCutoff": config.get("agfr_pocketCutoff"),
        "mapTypes": config.get("agfr_mapTypes", "all"),
        "autoSiteV": config.get("agfr_autoSiteV", "1.0"),
        "ligandSize": config.get("agfr_ligandSize", 500),
        "pepScore": bool(config.get("agfr_pepScore", True)),
        "recGradient": bool(config.get("agfr_recGradient", False)),
        "constrictionMapCutoff": config.get("agfr_constrictionMapCutoff"),
        "recGradVolCut": config.get("agfr_recGradVolCut"),
        "waterWeight": config.get("agfr_waterWeight"),
        "waterEntropy": config.get("agfr_waterEntropy"),
        "ionCharges": config.get("agfr_ionCharges"),
        "pdbLigand": config.get("agfr_pdbLigand"),
        "bioMol": config.get("agfr_bioMol"),
        "addHlig": config.get("agfr_addHlig"),
        "onException": config.get("agfr_onException"),
        "debug": config.get("agfr_debug"),
        "protToRemove": config.get("agfr_protToRemove"),
        "heteroToKeep": config.get("agfr_heteroToKeep"),
        "mutate": config.get("agfr_mutate"),
        "toPdbqt": bool(config.get("agfr_toPdbqt", False)),
        "testBatch": bool(config.get("agfr_testBatch", False)),
        "saveFilename": config.get("agfr_saveFilename"),
        "covalentBond": config.get("agfr_covalentBond"),
        "covalentBondTorsionAtom": config.get("agfr_covalentBondTorsionAtom"),
        "covalentResidues": config.get("agfr_covalentResidues"),
    }

    # Replicas
    reps = rm2.get("replicas", {})
    planned_count = reps.get("planned_count")
    if planned_count is None:
        raise PipelineError("Mode B: run_manifest.replicas.planned_count is required")
    idx_start = int(reps.get("index_start", 1))
    width = int(reps.get("suffix_width", 3))
    seed_base = reps.get("seed_base", None)

    cfg["replicas"] = {
        "count": int(planned_count),
        "index_start": idx_start,
        "suffix_width": width,
        "seed_base": seed_base,
        "dedup_target_copy": True,
        "shared_target_trg": str(target_trg_p),
    }

    # ADCP
    adcp_rerank_mode = config.get("rerank_mode", "interaction")
    cfg["adcp_jobname_base"] = rm2.get("campaign", {}).get("campaign_id") or f"{rm2.get('experiment_id','exp')}_{campaign.get('protein_id','P')}_{campaign.get('peptide_id','pep')}"
    cfg["adcp_base"] = {
        "target": str(target_trg_p),
        "sequence": cfg.get("peptide_seq"),
        "partition": config.get("adcp_partition"),
        "input": cfg.get("peptide_pdb"),
        "numSteps": int(config.get("adcp_numSteps", 200000)),
        "nbRuns": int(config.get("adcp_nbRuns", 10)),
        "maxCores": int(config.get("adcp_maxCores", 4)),
        "keepWorkingFolder": bool(config.get("adcp_keepWorkingFolder", False)),
        "dryRun": False,
        "cyclic": bool(config.get("adcp_cyclic", False)),
        "cystein": bool(config.get("adcp_cystein", False)),
        "overwriteFiles": bool(config.get("adcp_overwriteFiles", False)),  # IMPORTANT: controls -O for ADCP outputs
        "seed": None,
        "natContacts": config.get("adcp_natContacts"),
        "rmsd": config.get("adcp_rmsd"),
        "ref": config.get("adcp_ref"),
        "nmodes": int(config.get("adcp_nmodes", 10)),
        "reclusterOnly": bool(config.get("adcp_reclusterOnly", False)),
        "rotlibs": config.get("adcp_rotlibs"),
        "userRotlibs": config.get("adcp_userRotlibs"),
        "omm_nmin": int(config.get("omm_nmin", 5)),
        "dockingRanking": bool(config.get("adcp_dockingRanking", False)),
        "omm_max_itr": int(config.get("omm_max_itr", 1000)),
        "omm_environment": config.get("omm_environment", "implicit"),
        "fix_nst": bool(config.get("adcp_fix_nst", False)),
        "postDockMinimize": bool(config.get("adcp_postDockMinimize", False)),
        "ResumeMinimize": bool(config.get("adcp_ResumeMinimize", False)),
        "rerank_mode": adcp_rerank_mode,
        "ffxmlset": config.get("adcp_ffxmlset"),
        "userffxml": config.get("adcp_userffxml"),
    }

    if cfg["adcp_base"]["omm_nmin"] <= 0:
        raise PipelineError("Mode B: omm_nmin must be > 0")

    # validate core files exist
    require_exists(receptor_pdb, "receptor PDB (mode B)")
    if peptide_pdb is not None:
        require_exists(peptide_pdb, "peptide PDB (mode B)")

    timing_json, timing_csv = mode_b_defaults_for_timing(rm2, exp_root)

    # Attach override metadata into rm2 for persistence
    rm2.setdefault("runtime", {})
    rm2["runtime"]["cli_override_policy"] = {
        "allow_cli_overrides": bool(args.allow_cli_overrides),
        "cli_overrides_passed": list(cli_overrides.keys()),
        "override_records": override_records,
    }

    return cfg, rm2, exp_root, campaign_dir, timing_json, timing_csv


# --------------------------- preflight cleanup (overwrite-reps) ---------------------------

def preflight_overwrite_reps(
    runner: TeeRunner,
    preflight_log: Path,
    replicas_dir: Path,
    results_dir: Path,
    logs_dir: Path,
) -> None:
    # Remove replicas and results; clean logs except preflight.log
    if replicas_dir.exists():
        _print_and_log(runner, preflight_log, f"[preflight] overwrite-reps: removing replicas_dir: {replicas_dir}")
        shutil.rmtree(replicas_dir)
    if results_dir.exists():
        _print_and_log(runner, preflight_log, f"[preflight] overwrite-reps: removing results_dir: {results_dir}")
        shutil.rmtree(results_dir)
    if logs_dir.exists():
        kept = logs_dir / "preflight.log"
        _print_and_log(runner, preflight_log, f"[preflight] overwrite-reps: cleaned logs_dir (kept preflight.log): {logs_dir}")
        for p in logs_dir.glob("*"):
            if p.resolve() == kept.resolve():
                continue
            if p.is_dir():
                shutil.rmtree(p)
            else:
                try:
                    p.unlink()
                except Exception:
                    pass
    _print_and_log(runner, preflight_log, "[preflight] overwrite-reps: DONE (targets preserved; init not forced).")


# --------------------------- main (Mode B only here; legacy omitted for brevity) ---------------------------

def main(argv: Optional[List[str]] = None) -> int:
    args = parse_args(argv)

    # Deprecation handling
    if args.overwrite and not args.overwrite_reps:
        args.overwrite_reps = True

    mode_b = bool(args.run_manifest)
    runner = TeeRunner(realtime=args.real_time, verbosity=args.verbosity)

    # Required tools presence
    for tool in ["prepare_receptor", "agfr", "adcp"]:
        if not which_or_none(tool):
            raise PipelineError(f"Required tool not found on PATH: {tool}")

    # Prefer conda/mamba libs if available
    env = os.environ.copy()
    conda_prefix = env.get("CONDA_PREFIX", "")
    if conda_prefix:
        env["LD_LIBRARY_PATH"] = f"{conda_prefix}/lib:" + env.get("LD_LIBRARY_PATH", "")

    if not mode_b:
        raise PipelineError("This revised snippet focuses on Mode B. (Legacy mode code can be re-added if needed.)")

    run_manifest_path = Path(args.run_manifest).resolve()

    # We'll need logs_dir early for preflight logging; infer exp_root/campaign_dir first.
    exp_root = exp_root_from_run_manifest(run_manifest_path)
    campaign_dir = campaign_dir_from_run_manifest(run_manifest_path)

    # Ensure campaign dir exists
    ensure_dir(campaign_dir)
    logs_dir = (campaign_dir / "logs").resolve()
    ensure_dir(logs_dir)
    preflight_log = logs_dir / "preflight.log"

    # Preflight banner / explicit operational state
    _print_and_log(runner, preflight_log, f"[preflight] Mode B: campaign_dir={campaign_dir}")
    _print_and_log(runner, preflight_log, f"[preflight] Mode B: exp_root={exp_root}")

    # Resolve cfg + manifest with override policy
    cfg, rm2, exp_root, campaign_dir, timing_json, timing_csv = resolve_cfg_mode_b(
        args, run_manifest_path, runner, preflight_log
    )

    targets_dir = Path(cfg["targets_dir"]).resolve()
    replicas_dir = Path(cfg["replicas_root"]).resolve()
    results_dir = Path(cfg["results_dir"]).resolve()
    logs_dir = Path(cfg["logs_dir"]).resolve()

    for d in [targets_dir, replicas_dir, results_dir, logs_dir]:
        ensure_dir(d)

    trg_path = Path(cfg["target_trg"]).resolve()
    pdbqt_path = Path(cfg["receptor_pdbqt"]).resolve()

    _print_and_log(runner, preflight_log, f"[preflight] Target .trg per-condition: {trg_path}")
    _print_and_log(runner, preflight_log, f"[preflight] Target pdbqt per-condition: {pdbqt_path}")
    _print_and_log(runner, preflight_log,
                   f"[preflight] overwrite-reps={bool(args.overwrite_reps)} force-init={bool(args.force_init)} "
                   f"init-only={bool(args.init_only)} validate-only={bool(args.validate_only)}")
    if args.overwrite:
        _print_and_log(runner, preflight_log, "[preflight] WARNING: --overwrite is deprecated; treated as --overwrite-reps")

    rm_mode = cfg["adcp_base"].get("rerank_mode", "interaction")
    extra = " (adds -reint)" if rm_mode == "interaction" else " (no -reint)"
    _print_and_log(runner, preflight_log, f"[preflight] ADCP rerank mode: {rm_mode}{extra}")

    # If overwrite-reps: clean only replicas/results/logs (targets preserved)
    if args.overwrite_reps:
        preflight_overwrite_reps(runner, preflight_log, replicas_dir, results_dir, logs_dir)

        # Recreate dirs after removal
        ensure_dir(replicas_dir)
        ensure_dir(results_dir)
        ensure_dir(logs_dir)

    # Build runtime manifest overlay (written into run_manifest.runtime by persistence)
    runtime_manifest: Dict[str, Any] = {
        "mode": "B",
        "created_at": now_iso(),
        "argv": sys.argv[1:] if argv is None else argv,
        "operational": {
            "verbosity": args.verbosity,
            "real_time": bool(args.real_time),
            "validate_only": bool(args.validate_only),
            "init_only": bool(args.init_only),
            "overwrite_reps": bool(args.overwrite_reps),
            "force_init": bool(args.force_init),
            "dry_run": bool(args.dry_run),
            "allow_cli_overrides": bool(args.allow_cli_overrides),
        },
        "tools": {
            "prepare_receptor": which_or_none("prepare_receptor"),
            "agfr": which_or_none("agfr"),
            "adcp": which_or_none("adcp"),
            "python": sys.executable,
        },
        "resolved_config": cfg,
        "stages": [],
        "replica_runs": [],
        "status": "running",
    }

    stage_timings: List[StageTiming] = []

    stage_explain: Dict[str, str] = {
        "prepare_receptor": (
            "Preparing the receptor for docking:\n"
            "  - converts receptor PDB → PDBQT\n"
            "  - adds/checks hydrogens, assigns atom types/charges\n"
            "  - produces a receptor file compatible with AGFR/ADCP"
        ),
        "agfr": (
            "Building docking maps (target .trg) with AutoSite pocket detection:\n"
            "  - identifies likely binding pockets (AutoSite)\n"
            "  - builds grid-based energy maps for docking\n"
            "  - writes a single .trg bundle used by ADCP"
        ),
        "adcp": (
            "Docking the peptide with ADCP and refining top poses with OpenMM:\n"
            "  - runs Monte Carlo docking replicas\n"
            "  - minimizes top poses (required)\n"
            "  - re-ranks primarily by OpenMM, defaulting to interaction energy (Ecomplex − Ereceptor − Epeptide)"
        ),
    }

    def persist_mode_b() -> None:
        # timing artifacts
        dump_json(timing_json, {"generated_at": now_iso(), "stages": [asdict(x) for x in stage_timings]})
        write_timing_csv(timing_csv, stage_timings)

        rm_out = dict(rm2)
        rm_out.pop("__run_manifest_path__", None)
        rm_out.setdefault("runtime", {})
        rm_out["runtime"]["last_run_at"] = now_iso()
        rm_out["runtime"]["status"] = runtime_manifest.get("status", "running")
        rm_out["runtime"]["timing_json"] = str(Path(timing_json).resolve())
        rm_out["runtime"]["timing_csv"] = str(Path(timing_csv).resolve())
        rm_out["runtime"]["replica_count"] = int(cfg["replicas"]["count"])
        rm_out["runtime"]["replica_runs"] = runtime_manifest.get("replica_runs", [])
        rm_out["runtime"]["stages"] = runtime_manifest.get("stages", [])
        rm_out["runtime"]["resolved_config"] = cfg
        dump_json(run_manifest_path, rm_out)

    def run_stage(stage_name: str, cmd: List[str], log_file: Path, stage_label: str = "") -> int:
        start_t = time.time()
        start_iso = now_iso()
        label = stage_label or stage_name

        if args.verbosity >= 1:
            print(f"\n=== [{label}] START {start_iso} ===")
            expl = stage_explain.get(stage_name, "")
            if expl:
                print(expl)
            print("CMD:", " ".join(cmd))

        if args.dry_run:
            rc = 0
            ensure_dir(log_file.parent)
            with log_file.open("a", encoding="utf-8") as f:
                f.write(f"[{now_iso()}] DRY RUN: {' '.join(cmd)}\n")
        else:
            rc = runner.run(cmd, log_file, env=env, cwd=None, stage_label=label)

        end_iso = now_iso()
        elapsed = time.time() - start_t

        st = StageTiming(
            name=label,
            start_iso=start_iso,
            end_iso=end_iso,
            elapsed_s=float(elapsed),
            returncode=int(rc),
            command=cmd,
            log_file=str(log_file),
        )
        stage_timings.append(st)
        runtime_manifest["stages"].append(asdict(st))
        persist_mode_b()

        # Robust failure detection
        failed, why = detect_tool_failure(stage_name if not stage_label else stage_label, rc, log_file)
        if failed:
            runtime_manifest["status"] = "failed"
            persist_mode_b()
            raise PipelineError(why)

        if args.verbosity >= 1:
            print(f"=== [{label}] END {end_iso} (elapsed {elapsed:.2f}s, rc={rc}) ===")

        return rc

    # --------------- validate-only ---------------
    if args.validate_only:
        # Validate inputs and (if targets exist) validate trg zip
        require_exists(Path(cfg["receptor_pdb"]), "receptor PDB")
        if cfg.get("peptide_pdb"):
            require_exists(Path(cfg["peptide_pdb"]), "peptide PDB")
        ok = True
        msg = "targets not present yet"
        if trg_path.exists():
            ok, msg = validate_trg_zip(trg_path)
        runtime_manifest["agfr_trg_validation"] = {"ok": ok, "message": msg}
        runtime_manifest["status"] = "success" if ok else "failed"
        persist_mode_b()
        _print_and_log(runner, preflight_log, f"[preflight] validate-only done. status={runtime_manifest['status']}")
        return 0 if ok else 2

    # --------------- init gating ---------------
    trg_ok = False
    trg_msg = "missing"
    if trg_path.exists():
        trg_ok, trg_msg = validate_trg_zip(trg_path)

    _print_and_log(runner, preflight_log,
                   f"[preflight] init gating: force_init={bool(args.force_init)} "
                   f"pdbqt_exists={pdbqt_path.exists()} trg_exists={trg_path.exists()} trg_ok={trg_ok}")
    if trg_path.exists():
        _print_and_log(runner, preflight_log, f"[preflight] trg validation: {trg_msg}")

    need_init = bool(args.force_init) or (not pdbqt_path.exists()) or (not trg_path.exists()) or (not trg_ok)

    if not need_init:
        _print_and_log(runner, preflight_log, "[preflight] init SKIP (pdbqt+trg present and trg valid).")
    else:
        # Stage A: prepare_receptor
        prep_log = logs_dir / "prepare_receptor.log"
        prep_cmd = build_prepare_receptor_cmd(cfg)
        run_stage("prepare_receptor", prep_cmd, prep_log, stage_label="prepare_receptor")

        # Stage B: agfr
        agfr_log = logs_dir / "agfr.log"
        agfr_cmd = build_agfr_cmd(cfg)
        run_stage("agfr", agfr_cmd, agfr_log, stage_label="agfr")

        # Validate trg after build
        ok, msg = validate_trg_zip(trg_path)
        runtime_manifest["agfr_trg_validation"] = {"ok": ok, "message": msg}
        persist_mode_b()
        if not ok:
            runtime_manifest["status"] = "failed"
            persist_mode_b()
            raise PipelineError(f"Target .trg is invalid after agfr: {msg}")

    if args.init_only:
        runtime_manifest["status"] = "success"
        runtime_manifest["completed_at"] = now_iso()
        persist_mode_b()
        print("init-only done. Target validated; replicas skipped.")
        return 0

    # --------------- replicas loop ---------------
    rep_cfg = cfg["replicas"]
    rep_count = int(rep_cfg["count"])
    idx0 = int(rep_cfg["index_start"])
    width = int(rep_cfg["suffix_width"])
    seed_base = rep_cfg.get("seed_base", None)

    any_failed = False

    for k in range(rep_count):
        rep_index = idx0 + k
        rep_label = _replica_label(rep_index, width)
        rep_paths = _replica_paths_mode_b(replicas_dir, rep_label)
        for pth in rep_paths.values():
            ensure_dir(pth)

        rep_log = rep_paths["logs_dir"] / "adcp.log"
        stage_label = f"adcp:r{rep_label}"

        ad = dict(cfg["adcp_base"])
        jobname = f"{cfg['adcp_jobname_base']}_r{rep_label}"
        ad["jobName"] = jobname
        ad["workingFolder"] = str(rep_paths["docking_work"].resolve())
        ad["target"] = str(trg_path)

        if seed_base is not None:
            ad["seed"] = int(seed_base) + (rep_index - idx0)

        adcp_cmd = build_adcp_cmd(ad)

        # Run ADCP stage with robust log-checking
        try:
            run_stage("adcp", adcp_cmd, rep_log, stage_label=stage_label)
            rc = 0
        except PipelineError as e:
            rc = 2
            any_failed = True
            _print_and_log(runner, preflight_log, f"[preflight] replica {rep_label} failed: {e}")

        # Collect outputs
        found = _collect_adcp_outputs(jobname, [
            Path.cwd(),
            Path(ad["workingFolder"]).resolve(),
            Path(ad["workingFolder"]).resolve().parent,
            rep_paths["replica_root"],
        ])
        copied = {
            "out_pdb": _copy_if_present(found.get("out_pdb"), rep_paths["results_dir"]),
            "summary_dlg": _copy_if_present(found.get("summary_dlg"), rep_paths["results_dir"]),
            "target_trg": str(trg_path),
        }

        # Replica record
        runtime_manifest["replica_runs"].append({
            "replica_index": rep_index,
            "replica_label": rep_label,
            "paths": {kk: str(vv) for kk, vv in rep_paths.items()},
            "status": "success" if rc == 0 else "failed",
            "outputs": {
                "jobName": jobname,
                "found": found,
                "copied": copied,
            }
        })
        persist_mode_b()

    runtime_manifest["status"] = "failed" if any_failed else "success"
    runtime_manifest["completed_at"] = now_iso()
    persist_mode_b()

    if args.verbosity >= 1:
        if any_failed:
            print("\n⚠️ Mode B pipeline completed with failures in one or more replicas.")
        else:
            print("\n✅ Mode B pipeline completed successfully.")
        print("Campaign dir:", campaign_dir)
        print("Results:", results_dir)
        print("Timing:", timing_json, timing_csv)
        print("Run manifest updated:", run_manifest_path)

    return 0 if not any_failed else 3


if __name__ == "__main__":
    try:
        sys.exit(main())
    except PipelineError as e:
        print(f"\n❌ PipelineError: {e}", file=sys.stderr)
        sys.exit(2)
    except KeyboardInterrupt:
        print("\nInterrupted.", file=sys.stderr)
        sys.exit(130)
