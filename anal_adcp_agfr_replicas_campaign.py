#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
anal_adcp_agfr_replicas_campaign.py

Parser + consolidación + export PDB + QA para réplicas ADCP + OpenMM (re-ranking -reint).

Invocación sugerida:
  python anal_adcp_agfr_replicas_campaign.py --exp-root exp001 --outdir analysis

(Después, si deseas, puedes instalarlo como comando `anal_adcp_agfr_replicas_campaign` con un entry-point.)
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import re
import sys
from dataclasses import dataclass, asdict
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Iterable


SCRIPT_NAME = "anal_adcp_agfr_replicas_campaign.py"
SCRIPT_VERSION = "v2"



# ----------------------------
# Utilidades generales
# ----------------------------

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)


def vlog(level: int, verbose: int, msg: str) -> None:
    if verbose >= level:
        eprint(msg)


def safe_mkdir(p: Path) -> None:
    p.mkdir(parents=True, exist_ok=True)


def read_text(path: Path) -> str:
    with path.open("r", encoding="utf-8", errors="replace") as f:
        return f.read()


def write_text(path: Path, text: str) -> None:
    safe_mkdir(path.parent)
    with path.open("w", encoding="utf-8") as f:
        f.write(text)


def now_iso() -> str:
    return datetime.now().isoformat(timespec="seconds")


def median(xs: List[float]) -> float:
    ys = sorted(xs)
    n = len(ys)
    if n == 0:
        return float("nan")
    mid = n // 2
    if n % 2 == 1:
        return ys[mid]
    return 0.5 * (ys[mid - 1] + ys[mid])


def percentile(xs: List[float], p: float) -> float:
    """p in [0,100]. Linear interpolation."""
    ys = sorted(xs)
    n = len(ys)
    if n == 0:
        return float("nan")
    if n == 1:
        return ys[0]
    rank = (p / 100.0) * (n - 1)
    lo = int(math.floor(rank))
    hi = int(math.ceil(rank))
    if lo == hi:
        return ys[lo]
    w = rank - lo
    return (1 - w) * ys[lo] + w * ys[hi]


def iqr(xs: List[float]) -> float:
    if len(xs) == 0:
        return float("nan")
    return percentile(xs, 75) - percentile(xs, 25)


def shannon_entropy_from_counts(counts: Dict[str, int]) -> float:
    n = sum(counts.values())
    if n <= 0:
        return float("nan")
    h = 0.0
    for c in counts.values():
        if c <= 0:
            continue
        p = c / float(n)
        h -= p * math.log(p)
    return h


# ----------------------------
# Descubrimiento de réplicas
# ----------------------------

@dataclass
class ReplicaKey:
    experiment_id: str
    protein_id: str
    peptide_id: str
    replica_id: str

    def group_key(self) -> Tuple[str, str, str]:
        return (self.experiment_id, self.protein_id, self.peptide_id)


@dataclass
class ReplicaPaths:
    key: ReplicaKey
    replica_dir: Path
    summary_path: Optional[Path]
    rescored_pdb_path: Optional[Path]
    out_pdb_path: Optional[Path]


def parse_replica_key_from_path(exp_root: Path, replica_dir: Path) -> Optional[ReplicaKey]:
    """
    Espera rutas tipo:
      <exp_root>/runs/<Protein>/<Peptide>/replicas/replica_XXX
    """
    try:
        rel = replica_dir.resolve().relative_to(exp_root.resolve())
    except Exception:
        return None

    parts = rel.parts
    # runs / P001 / polya / replicas / replica_001
    if len(parts) < 5:
        return None
    if parts[0] != "runs":
        return None
    if parts[3] != "replicas":
        return None

    protein_id = parts[1]
    peptide_id = parts[2]
    replica_leaf = parts[4]
    m = re.match(r"replica_(\d+)$", replica_leaf)
    if not m:
        return None
    replica_id = m.group(1)

    experiment_id = exp_root.name  # e.g. exp001
    return ReplicaKey(experiment_id=experiment_id, protein_id=protein_id, peptide_id=peptide_id, replica_id=replica_id)


def pick_preferred_file(candidates: List[Path]) -> Optional[Path]:
    """
    Preferimos resultados finales en results/; si hay varios, tomamos el más reciente.
    """
    if not candidates:
        return None

    def score(p: Path) -> Tuple[int, float]:
        # results/ preferido
        in_results = 1 if ("results" in p.parts) else 0
        try:
            mtime = p.stat().st_mtime
        except Exception:
            mtime = 0.0
        return (in_results, mtime)

    return sorted(candidates, key=score, reverse=True)[0]


def discover_replicas(exp_root: Path) -> Tuple[List[ReplicaPaths], Dict[str, Dict[str, int]]]:
    """
    Recorre runs/**/replicas/replica_*/ y arma rutas por réplica.
    También produce un resumen por (protein, peptide) de cuántas réplicas se detectaron/parseables.
    """
    runs_dir = exp_root / "runs"
    discovered: List[ReplicaPaths] = []
    group_counts: Dict[str, Dict[str, int]] = {}  # key "P001/polya" -> counters

    if not runs_dir.exists():
        return [], {}

    for replica_dir in runs_dir.glob("*/*/replicas/replica_*"):
        if not replica_dir.is_dir():
            continue

        key = parse_replica_key_from_path(exp_root, replica_dir)
        if key is None:
            continue

        # localizar summary/rescored/out
        summary_candidates = []
        summary_candidates += list(replica_dir.glob("results/*_summary.dlg"))
        summary_candidates += list(replica_dir.glob("docking/work/*_summary.dlg"))
        # algunos wrappers podrían anidar más profundo:
        summary_candidates += list(replica_dir.glob("docking/work/**/*_summary.dlg"))

        rescored_candidates = []
        rescored_candidates += list(replica_dir.glob("results/*_omm_rescored_out.pdb"))
        rescored_candidates += list(replica_dir.glob("docking/work/*_omm_rescored_out.pdb"))
        rescored_candidates += list(replica_dir.glob("docking/work/**/*_omm_rescored_out.pdb"))

        out_candidates = []
        out_candidates += list(replica_dir.glob("docking/work/*_out.pdb"))
        out_candidates += list(replica_dir.glob("docking/work/**/*_out.pdb"))

        summary_path = pick_preferred_file(summary_candidates)
        rescored_path = pick_preferred_file(rescored_candidates)
        out_path = pick_preferred_file(out_candidates)

        discovered.append(
            ReplicaPaths(
                key=key,
                replica_dir=replica_dir,
                summary_path=summary_path,
                rescored_pdb_path=rescored_path,
                out_pdb_path=out_path,
            )
        )

        gk = f"{key.protein_id}/{key.peptide_id}"
        if gk not in group_counts:
            group_counts[gk] = {"replicas_detected": 0, "replicas_with_summary": 0, "replicas_with_rescored": 0, "replicas_with_out": 0}
        group_counts[gk]["replicas_detected"] += 1
        if summary_path is not None:
            group_counts[gk]["replicas_with_summary"] += 1
        if rescored_path is not None:
            group_counts[gk]["replicas_with_rescored"] += 1
        if out_path is not None:
            group_counts[gk]["replicas_with_out"] += 1

    return discovered, group_counts


# ----------------------------
# Parsing de *_summary.dlg
# ----------------------------

@dataclass
class ParsedClusterRow:
    mode: int
    adcp_affinity_kcalmol: float
    cluster_size: int
    best_run_pose_id: str  # keep as string because can include leading zeros


@dataclass
class ParsedTopKRow:
    model_id: int
    rank_openmm: int
    rank_adcp: int
    omm_score_complex_minus_receptor: float
    omm_dE_interaction: float
    adcp_affinity_kcalmol: float
    cluster_size: int
    best_run_pose_id: str
    cluster_mode_id: Optional[int] = None


@dataclass
class ParsedReplica:
    # IDs
    experiment_id: str
    protein_id: str
    peptide_id: str
    replica_id: str

    # Paths
    summary_path: str
    rescored_pdb_path: str
    out_pdb_path: str

    # Params/metadata
    sequence: str
    N_runs: int
    n_evals: int
    nmin: int
    env: str
    nitr: int

    # ADCP raw
    adcp_bestEnergy: float
    adcp_bestEnergy_run: int
    adcp_bestEnergies_list: str

    # Winner
    winner_rank_openmm: int
    winner_model_id: int
    winner_omm_dE_interaction: float
    winner_omm_dE_complex_minus_receptor: float
    winner_adcp_affinity: float
    winner_rank_adcp: int
    winner_cluster_size: int
    winner_best_run_pose_id: str
    winner_cluster_mode_id: Optional[int]

    # QA
    qa_status: str
    qa_message: str


_RE_BESTENERGIES = re.compile(r"bestEnergies\s+\[(.*?)\]\s*", re.IGNORECASE)
_RE_BESTENERGY_IN_RUN = re.compile(r"bestEnergy in run\s+(\d+)\s+([-0-9\.Ee+]+)", re.IGNORECASE)
_RE_PERFORMING = re.compile(r"Performing\s+(\d+)\s+MC searches using\s+(\d+)\s+evals each", re.IGNORECASE)
_RE_MC_CMD_SEQ = re.compile(r'MC search command:.*?-t\s+\d+\s+"([^"]+)"', re.IGNORECASE)
_RE_OMM_SETTINGS = re.compile(r'OpenMM minimization settings:\s+Environment="([^"]+)";\s+Max_itr=(\d+)', re.IGNORECASE)
_RE_FROM_TOTAL = re.compile(r"From total\s+(\d+)\s+models,\s+minimizing top\s+(\d+)", re.IGNORECASE)

# cluster table rows, e.g.: "   1        -14.1      0.0      35      NA      NA    488"
_RE_CLUSTER_ROW = re.compile(
    r"^\s*(\d+)\s+([-0-9\.Ee+]+)\s+([-0-9\.Ee+]+)\s+(\d+)\s+\S+\s+\S+\s+(\S+)\s*$"
)

# OMM Energy block:
# OMM Energy: E_Complex =  -8849.81; E_Receptor =  -8676.93; E_Peptide  =    -97.06
# OMM Energy: dE_Interaction =    -75.82; dE_Complex-Receptor =   -172.88
_RE_OMM_ECOMP = re.compile(r"OMM Energy:\s+E_Complex\s*=\s*([-0-9\.Ee+]+);\s*E_Receptor\s*=\s*([-0-9\.Ee+]+);\s*E_Peptide\s*=\s*([-0-9\.Ee+]+)", re.IGNORECASE)
_RE_OMM_DE = re.compile(r"OMM Energy:\s+dE_Interaction\s*=\s*([-0-9\.Ee+]+);\s*dE_Complex-Receptor\s*=\s*([-0-9\.Ee+]+)", re.IGNORECASE)

# OMM Ranking rows:
# OMM Ranking:      1      1      2     -172.9        -75.8         -11.8      0.0      15      NA      NA    450
_RE_OMM_RANK_ROW = re.compile(
    r"^OMM Ranking:\s+(\d+)\s+(\d+)\s+(\d+)\s+([-0-9\.Ee+]+)\s+([-0-9\.Ee+]+)\s+([-0-9\.Ee+]+)\s+([-0-9\.Ee+]+)\s+(\d+)\s+\S+\s+\S+\s+(\S+)\s*$"
)


def parse_summary_dlg(text: str) -> Tuple[Dict[str, object], List[ParsedClusterRow], List[ParsedTopKRow], List[str]]:
    """
    Extrae:
      - meta dict
      - clusters (mode 1..10)
      - topk (OMM Ranking, típicamente 5 filas)
      - warnings (strings)
    """
    warnings: List[str] = []

    # Defaults
    meta = {
        "sequence": "",
        "N_runs": 0,
        "n_evals": 0,
        "nmin": 0,
        "env": "",
        "nitr": 0,
        "adcp_bestEnergy": float("nan"),
        "adcp_bestEnergy_run": 0,
        "adcp_bestEnergies_list": "",
    }

    m = _RE_PERFORMING.search(text)
    if m:
        meta["N_runs"] = int(m.group(1))
        meta["n_evals"] = int(m.group(2))
    else:
        warnings.append("Missing 'Performing ... MC searches' line (N_runs/n_evals).")

    m = _RE_MC_CMD_SEQ.search(text)
    if m:
        meta["sequence"] = m.group(1).strip()
    else:
        warnings.append("Missing sequence in 'MC search command' line.")

    m = _RE_OMM_SETTINGS.search(text)
    if m:
        meta["env"] = m.group(1).strip()
        meta["nitr"] = int(m.group(2))
    else:
        warnings.append("Missing 'OpenMM minimization settings' (env/nitr).")

    m = _RE_FROM_TOTAL.search(text)
    if m:
        meta["nmin"] = int(m.group(2))
    else:
        # fallback: count OMM Ranking lines later
        pass

    m = _RE_BESTENERGIES.search(text)
    if m:
        meta["adcp_bestEnergies_list"] = m.group(1).strip()
    else:
        warnings.append("Missing 'bestEnergies [...]' list.")

    m = _RE_BESTENERGY_IN_RUN.search(text)
    if m:
        meta["adcp_bestEnergy_run"] = int(m.group(1))
        meta["adcp_bestEnergy"] = float(m.group(2))
    else:
        warnings.append("Missing 'bestEnergy in run ...' line.")

    # Parse cluster rows
    clusters: List[ParsedClusterRow] = []
    in_cluster_table = False
    for line in text.splitlines():
        if line.strip().startswith("mode |") and "affinity" in line and "clust." in line:
            in_cluster_table = True
            continue
        if in_cluster_table:
            if line.strip().startswith("-----+"):
                continue
            if line.strip().startswith("Calculations completed"):
                in_cluster_table = False
                continue
            mm = _RE_CLUSTER_ROW.match(line)
            if mm:
                mode = int(mm.group(1))
                affinity = float(mm.group(2))
                cluster_size = int(mm.group(4))
                best_run_pose_id = str(mm.group(5))
                clusters.append(ParsedClusterRow(mode=mode, adcp_affinity_kcalmol=affinity, cluster_size=cluster_size, best_run_pose_id=best_run_pose_id))

    if not clusters:
        warnings.append("No cluster rows parsed from ADCP cluster table.")

    # Parse OMM Ranking topk rows
    topk: List[ParsedTopKRow] = []
    for line in text.splitlines():
        mm = _RE_OMM_RANK_ROW.match(line)
        if mm:
            model_id = int(mm.group(1))
            rank_openmm = int(mm.group(2))
            rank_adcp = int(mm.group(3))
            e_complex_minus_receptor = float(mm.group(4))
            de_interaction = float(mm.group(5))
            affinity = float(mm.group(6))
            # mm.group(7) = ref.fnc (often 0.0)
            cluster_size = int(mm.group(8))
            best_run_pose_id = str(mm.group(9))
            topk.append(
                ParsedTopKRow(
                    model_id=model_id,
                    rank_openmm=rank_openmm,
                    rank_adcp=rank_adcp,
                    omm_score_complex_minus_receptor=e_complex_minus_receptor,
                    omm_dE_interaction=de_interaction,
                    adcp_affinity_kcalmol=affinity,
                    cluster_size=cluster_size,
                    best_run_pose_id=best_run_pose_id,
                )
            )

    if not topk:
        warnings.append("No OMM Ranking rows parsed (top-k).")
    else:
        # If nmin missing, infer it
        if meta["nmin"] == 0:
            meta["nmin"] = len(topk)

    # Map best_run_pose_id -> mode for Definition B
    best_to_mode: Dict[str, int] = {}
    for c in clusters:
        # if duplicate best_run_pose_id appears, keep first (stable)
        if c.best_run_pose_id not in best_to_mode:
            best_to_mode[c.best_run_pose_id] = c.mode

    for row in topk:
        row.cluster_mode_id = best_to_mode.get(row.best_run_pose_id)

    # winner selection
    # We expect rank_openmm==1 is the winner row. If not, pick smallest rank_openmm.
    return meta, clusters, topk, warnings


# ----------------------------
# Parsing PDB multi-model & export
# ----------------------------

def iter_pdb_models(pdb_text: str) -> Iterable[Tuple[int, str]]:
    """
    Yields (model_number_in_file, model_block_text_including_MODEL_ENDMDL).
    Expects MODEL/ENDMDL records.
    """
    lines = pdb_text.splitlines(True)  # keep line endings
    cur: List[str] = []
    cur_id: Optional[int] = None

    for ln in lines:
        if ln.startswith("MODEL"):
            # flush previous if any
            if cur and cur_id is not None:
                yield (cur_id, "".join(cur))
                cur = []
                cur_id = None
            # new model
            cur = [ln]
            m = re.match(r"MODEL\s+(\d+)", ln.strip())
            cur_id = int(m.group(1)) if m else -1
        else:
            if cur:
                cur.append(ln)
                if ln.startswith("ENDMDL"):
                    if cur_id is None:
                        cur_id = -1
                    yield (cur_id, "".join(cur))
                    cur = []
                    cur_id = None

    # flush if file lacks ENDMDL at end
    if cur and cur_id is not None:
        yield (cur_id, "".join(cur))


def extract_first_k_models(pdb_path: Path, k: int) -> List[str]:
    text = read_text(pdb_path)
    blocks: List[str] = []
    for idx, (mid, block) in enumerate(iter_pdb_models(text), start=1):
        blocks.append(block)
        if idx >= k:
            break
    return blocks


def extract_model1(pdb_path: Path) -> Optional[str]:
    blocks = extract_first_k_models(pdb_path, 1)
    return blocks[0] if blocks else None


def renumber_models(blocks: List[str], start_at: int = 1) -> str:
    """
    Given list of MODEL...ENDMDL blocks, renumber sequentially in output.
    """
    out_lines: List[str] = []
    mnum = start_at
    for block in blocks:
        # Replace first line MODEL and keep rest
        lines = block.splitlines(True)
        if not lines:
            continue
        # Make sure we output MODEL with correct formatting
        lines[0] = f"MODEL     {mnum:4d}\n"
        out_lines.extend(lines)
        mnum += 1
    return "".join(out_lines)


def split_model_target_and_pose(block: str) -> Tuple[List[str], List[str], str]:
    """
    Split MODEL block body lines into target-protein and pose sections.

    Heuristic:
      - Prefer the FIRST TER record as boundary:
          target = lines up to (and including) that TER
          pose   = lines after that TER
      - If no TER is present, try chain-based split:
          identify first atom chain as target chain,
          split before first atom that belongs to a different chain.
      - If chain-based split cannot be inferred, leave target empty and
        mark whole model as pose.
    """
    lines = block.splitlines(True)
    if not lines:
        return [], [], "empty"

    # MODEL ... ENDMDL excluded from split payload.
    body = lines[1:]
    if body and body[-1].startswith("ENDMDL"):
        body = body[:-1]

    ter_idx = -1
    for i, ln in enumerate(body):
        if ln.startswith("TER"):
            ter_idx = i
            break

    if ter_idx >= 0:
        return body[:ter_idx + 1], body[ter_idx + 1:], "first_ter"

    def _chain_id(pdb_line: str) -> str:
        # PDB chain ID is column 22 (1-indexed), index 21 in Python.
        return pdb_line[21].strip() if len(pdb_line) > 21 else ""

    first_chain = ""
    split_idx = -1
    for i, ln in enumerate(body):
        if not (ln.startswith("ATOM") or ln.startswith("HETATM")):
            continue
        cid = _chain_id(ln)
        if not first_chain and cid:
            first_chain = cid
            continue
        if first_chain and cid and cid != first_chain:
            split_idx = i
            break

    if split_idx >= 0:
        return body[:split_idx], body[split_idx:], "chain_fallback"

    return [], body, "unsplit"


def extract_unrestrained_receptor_residues(target_lines: List[str]) -> List[Tuple[str, str, int]]:
    """
    Parse residues listed after:
      USER: RECEPTOR RESIDUES NOT RESTRAINED DURING MINIMIZATION:

    Expected token pattern inside USER lines: _<CHAIN>_<RESNAME>_<RESNUM>
    e.g. _A_ALA_52
    """
    capture = False
    chunks: List[str] = []
    for ln in target_lines:
        if not ln.startswith("USER"):
            continue

        payload = ln[4:].strip()
        if "RECEPTOR RESIDUES NOT RESTRAINED DURING MINIMIZATION:" in payload:
            capture = True
            _, _, suffix = payload.partition("RECEPTOR RESIDUES NOT RESTRAINED DURING MINIMIZATION:")
            if suffix:
                chunks.append(suffix)
            continue

        if capture:
            chunks.append(payload)

    text = " ".join(chunks)
    out: List[Tuple[str, str, int]] = []
    for m in re.finditer(r"_([A-Za-z0-9])_([A-Za-z]{3})_(-?\d+)", text):
        out.append((m.group(1), m.group(2).upper(), int(m.group(3))))
    return out


def format_unrestrained_receptor_residues_dat(residues: List[Tuple[str, str, int]]) -> str:
    return "\n".join(f"Chain  {chain}  {resname:>3}    {resnum}" for chain, resname, resnum in residues) + "\n"


# ----------------------------
# Sanity filter
# ----------------------------

@dataclass
class SanityConfig:
    enable: bool = True
    # If value is outside [min, max], mark QA. Defaults are permissive.
    eint_min: float = -1.0e6
    eint_max: float =  1.0e6


def sanity_check_eint(val: float, cfg: SanityConfig) -> Tuple[bool, str]:
    if not cfg.enable:
        return True, ""
    if val is None or (isinstance(val, float) and math.isnan(val)):
        return False, "winner_omm_dE_interaction is NaN"
    if val < cfg.eint_min or val > cfg.eint_max:
        return False, f"winner_omm_dE_interaction out of range [{cfg.eint_min}, {cfg.eint_max}]: {val}"
    return True, ""


# ----------------------------
# Escritura CSV
# ----------------------------

def write_csv(path: Path, rows: List[Dict[str, object]], fieldnames: List[str]) -> None:
    safe_mkdir(path.parent)
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)



def _import_phase2_deps():
    try:
        import numpy as np  # type: ignore
        import pandas as pd  # type: ignore
        from scipy.stats import mannwhitneyu  # type: ignore
        return np, pd, mannwhitneyu
    except Exception as ex:
        raise RuntimeError(
            "Phase 2 requires numpy, pandas, and scipy. Install them or run with --no-phase2."
        ) from ex


# ----------------------------
# PHASE 2: metrics + significance
# ----------------------------

PHASE2_GROUP_SUMMARY_FIELDS = [
    "protein_id", "peptide_id",
    "n_detected", "n_used", "fail_rate",
    "n_excluded_openmm_missing", "n_excluded_legacy", "n_excluded_other",
    "Eint_median", "Eint_IQR", "Eint_mean", "Eint_std", "Eint_p10", "Eint_p90",
    "cluster_mode_top", "f1_cluster", "H_cluster", "Neff_cluster",
    "flag_low_n", "flag_no_data",
]

PHASE2_DISCRIMINATION_FIELDS = [
    "protein_id", "poi_id", "control_id",
    "n_poi", "n_control",
    "median_poi", "median_control", "delta_median",
    "auc", "cliffs_delta",
    "p_mwu", "q_fdr", "p_perm",
    "f1_cluster_poi", "neff_cluster_poi",
    "f1_cluster_control", "neff_cluster_control",
    "delta_f1", "delta_neff",
    "flag_low_n", "flag_skipped",
]

PHASE2_GLOBAL_RANK_FIELDS = [
    "protein_id", "best_control_type", "auc_best", "q_best", "delta_f1_best",
    "fail_rate_poi", "fail_rate_controls_median", "verdict",
]

SET_SUMMARY_FIELDS = [
    "protein_id", "peptide_set",
    "n_peptides_detected", "n_peptides_used", "n_replicas_used",
    "member_peptides", "flag_low_n_set",
    "Eint_seq_median", "Eint_seq_IQR", "Eint_seq_mean", "Eint_seq_std",
    "Eint_raw_median", "Eint_raw_IQR", "Eint_raw_mean", "Eint_raw_std",
    "f1_seq_median", "f1_seq_IQR", "Neff_seq_median", "Neff_seq_IQR",
    "f1_raw_pool", "Neff_raw_pool",
]

PROTEIN_SET_DISCRIMINATION_FIELDS = [
    "protein_id", "poi_id", "control_set",
    "n_poi_replicas", "n_control_peptides", "n_control_replicas_total",
    "poi_Eint_seq_summary", "control_Eint_seq_median", "control_Eint_seq_IQR", "delta_seq_median",
    "poi_rank_among_control_seq", "poi_percentile_in_control_seq", "auc_seq", "p_empirical_seq",
    "poi_f1", "control_f1_seq_median", "delta_f1_seq",
    "poi_Neff", "control_Neff_seq_median", "delta_Neff_seq",
    "n_poi_replicas_raw", "n_control_replicas_raw",
    "poi_raw_median", "control_raw_median", "delta_raw_median",
    "auc_raw", "cliffs_delta_raw", "p_mwu_raw", "q_fdr_raw", "p_perm_raw",
    "delta_f1_rawpool", "delta_Neff_rawpool",
    "flag_low_n_set", "flag_skipped",
]

PROTEIN_SET_RANK_FIELDS = [
    "protein_id", "best_control_set", "primary_percentile", "primary_empirical_p",
    "secondary_auc_raw", "secondary_q_fdr_raw", "delta_f1_seq", "delta_f1_rawpool",
    "verdict_set_level",
]

SCORE_SYNONYMS = ["winner_omm_dE_interaction", "winner_omm_de_interaction", "winner_omm_deint", "winner_eint"]
CLUSTER_SYNONYMS = ["winner_cluster_mode_id", "winner_cluster_id", "cluster_mode_id"]
QA_SYNONYMS = ["qa_status", "status"]

VERDICT_THRESHOLDS = {
    "strong_auc": 0.8,
    "moderate_auc": 0.7,
    "weak_auc_low": 0.4,
    "weak_auc_high": 0.6,
    "strong_delta_f1": 0.2,
    "moderate_delta_f1": 0.15,
    "poi_fail_rate_max": 0.2,
    "polya_auc_confounded": 0.3,
    "confounded_delta_median": 5.0,
}

SET_VERDICT_THRESHOLDS = {
    "strong_percentile": 0.8,
    "moderate_percentile": 0.65,
    "strong_empirical_p": 0.1,
    "moderate_empirical_p": 0.25,
    "strong_auc_raw": 0.8,
    "moderate_auc_raw": 0.7,
    "strong_delta_f1": 0.2,
    "moderate_delta_f1": 0.1,
}


def _resolve_col(df: pd.DataFrame, candidates: List[str], required: bool = True) -> Optional[str]:
    lower_map = {c.lower(): c for c in df.columns}
    for c in candidates:
        if c in df.columns:
            return c
        if c.lower() in lower_map:
            return lower_map[c.lower()]
    if required:
        raise ValueError(f"Required column not found. Candidates={candidates}; available={list(df.columns)}")
    return None


def _is_finite_number(x: object) -> bool:
    try:
        xv = float(x)
    except Exception:
        return False
    return math.isfinite(xv)


def _qa_indicates_bad(qa_status: str) -> bool:
    s = (qa_status or "").upper()
    bad_tokens = ["OPENMM_MISSING", "INCOMPLETE_OPENMM", "MISSING_SUMMARY", "PARSE_ERROR", "SANITY_FAIL", "INCOMPLETE"]
    return any(t in s for t in bad_tokens)


def _excluded_reason(row: pd.Series, score_col: str, cluster_col: str, qa_col: str,
                     exclude_suffix: str, sanity_cfg: SanityConfig) -> str:
    pep = str(row.get("peptide_id", ""))
    if exclude_suffix and pep.endswith(exclude_suffix):
        return "legacy"

    score = row.get(score_col, float("nan"))
    cluster = row.get(cluster_col, float("nan"))
    qa_status = str(row.get(qa_col, ""))

    valid_score = _is_finite_number(score)
    if valid_score and sanity_cfg.enable:
        sval = float(score)
        valid_score = sanity_cfg.eint_min <= sval <= sanity_cfg.eint_max

    valid_cluster = _is_finite_number(cluster)

    if not valid_score or ("OPENMM" in qa_status.upper()) or _qa_indicates_bad(qa_status):
        return "openmm_missing"
    if not valid_cluster:
        return "other"
    if qa_status and qa_status.upper() not in {"OK", "PARSE_WARN"}:
        return "other"
    return "used"


def _bh_fdr(pvals: List[float]) -> List[float]:
    n = len(pvals)
    if n == 0:
        return []
    indexed = sorted(list(enumerate(pvals)), key=lambda t: t[1])
    q = [float("nan")] * n
    prev = 1.0
    for rank in range(n, 0, -1):
        i, p = indexed[rank - 1]
        val = min(prev, p * n / rank)
        prev = val
        q[i] = max(0.0, min(1.0, val))
    return q


def _compute_auc_and_cliffs(poi, ctrl) -> Tuple[float, float, float]:
    _, _, mannwhitneyu = _import_phase2_deps()
    # "better" means lower (more negative) energy
    res = mannwhitneyu(poi, ctrl, alternative="two-sided", method="auto")
    u_raw = float(res.statistic)
    n1, n2 = len(poi), len(ctrl)
    u_min = n1 * n2 - u_raw  # POI wins for smaller values
    auc = u_min / (n1 * n2)
    cliffs = 2.0 * auc - 1.0
    return u_raw, auc, cliffs


def _perm_pvalue_delta_median(poi, ctrl, n_perm: int, rng) -> float:
    np, _, _ = _import_phase2_deps()
    obs = float(np.median(poi) - np.median(ctrl))
    pooled = np.concatenate([poi, ctrl])
    n1 = len(poi)
    cnt = 0
    for _ in range(n_perm):
        idx = rng.permutation(len(pooled))
        a = pooled[idx[:n1]]
        b = pooled[idx[n1:]]
        stat = float(np.median(a) - np.median(b))
        if abs(stat) >= abs(obs):
            cnt += 1
    return (cnt + 1.0) / (n_perm + 1.0)


def _cluster_metrics_from_values(values) -> Tuple[float, float]:
    np, _, _ = _import_phase2_deps()
    arr = np.asarray(values)
    arr = arr[np.isfinite(arr)]
    if arr.size == 0:
        return float("nan"), float("nan")
    vals, cnts = np.unique(arr.astype(int), return_counts=True)
    if cnts.size == 0:
        return float("nan"), float("nan")
    f1 = float(cnts.max() / cnts.sum())
    probs = cnts / cnts.sum()
    h = float(-np.sum(probs * np.log(probs)))
    neff = float(np.exp(h))
    return f1, neff


def _prefer_exact_poi_set(peptide_set: str, poi_id: str) -> bool:
    normalized = str(peptide_set or "").strip().lower()
    return normalized == str(poi_id).strip().lower()


def _control_set_priority(peptide_set: str) -> Tuple[int, str]:
    name = str(peptide_set or "")
    lname = name.lower()
    if "matched" in lname:
        return (0, lname)
    if "decoy" in lname:
        return (1, lname)
    return (2, lname)


def _load_peptide_metadata(exp_root: Path, poi_id: str):
    _, pd, _ = _import_phase2_deps()
    peptides_csv = exp_root / "peptides.csv"
    if not peptides_csv.exists():
        raise FileNotFoundError(
            f"Set-level analysis requires metadata file: {peptides_csv}. "
            "Expected columns include peptide_id and peptide_set."
        )
    meta = pd.read_csv(peptides_csv)
    required = ["peptide_id", "peptide_set"]
    missing = [c for c in required if c not in meta.columns]
    if missing:
        raise ValueError(f"peptides.csv missing required columns: {missing}; available={list(meta.columns)}")
    meta["peptide_id"] = meta["peptide_id"].astype(str)
    meta["peptide_set"] = meta["peptide_set"].astype(str)

    dup = (
        meta[["peptide_id", "peptide_set"]]
        .drop_duplicates()
        .groupby("peptide_id")["peptide_set"]
        .nunique()
    )
    bad = dup[dup > 1]
    if not bad.empty:
        offenders = sorted(bad.index.tolist())
        raise ValueError(
            "Inconsistent peptide_id -> peptide_set mapping in peptides.csv for: "
            + ", ".join(offenders[:10])
            + (" ..." if len(offenders) > 10 else "")
        )

    mapping = (
        meta[["peptide_id", "peptide_set"]]
        .drop_duplicates(subset=["peptide_id"])
        .set_index("peptide_id")["peptide_set"]
        .to_dict()
    )
    poi_metadata_ids = sorted(meta.loc[meta["peptide_set"].map(lambda s: _prefer_exact_poi_set(s, poi_id)), "peptide_id"].unique().tolist())
    return peptides_csv, meta, mapping, poi_metadata_ids


def run_phase2_from_replicas_csv(
    replicas_csv: Path,
    outdir: Path,
    parse_report_json: Path,
    run_info_json: Path,
    sanity_cfg: SanityConfig,
    run_phase2: bool,
    poi_id: str,
    exclude_peptide_suffix: str,
    min_replicas: int,
    allow_low_n_comparisons: bool,
    alpha: float,
    fdr_method: str,
    permutation_tests: bool,
    n_perm: int,
    seed: Optional[int],
    exp_root: Path,
    min_set_peptides: int,
    min_set_replicas: int,
    verbose: int = 0,
) -> None:
    if not run_phase2:
        vlog(1, verbose, "[INFO] Phase 2 disabled by CLI (--no-phase2).")
        return

    np, pd, mannwhitneyu = _import_phase2_deps()

    if not replicas_csv.exists():
        raise FileNotFoundError(f"Phase 2 requires {replicas_csv}")

    if fdr_method.lower() != "bh":
        raise ValueError("Only fdr_method='bh' is implemented.")

    df = pd.read_csv(replicas_csv)
    if df.empty:
        vlog(1, verbose, "[WARN] replicas_parsed.csv is empty; generating empty Phase 2 outputs.")

    score_col = _resolve_col(df, SCORE_SYNONYMS)
    cluster_col = _resolve_col(df, CLUSTER_SYNONYMS)
    qa_col = _resolve_col(df, QA_SYNONYMS)

    for c in [score_col, cluster_col]:
        df[c] = pd.to_numeric(df[c], errors="coerce")

    df["_excluded_reason"] = df.apply(
        lambda r: _excluded_reason(r, score_col, cluster_col, qa_col, exclude_peptide_suffix, sanity_cfg),
        axis=1,
    )
    df["_used"] = df["_excluded_reason"] == "used"

    skipped_legacy_groups = sorted(df.loc[df["_excluded_reason"] == "legacy", "peptide_id"].dropna().unique().tolist())

    group_rows = []
    set_warnings: List[str] = []
    peptides_csv, peptide_meta, peptide_to_set, poi_metadata_ids = _load_peptide_metadata(exp_root, poi_id)
    grouped = df.groupby(["protein_id", "peptide_id"], sort=True, dropna=False)
    for (protein_id, peptide_id), g in grouped:
        n_detected = int(len(g))
        g_used = g[g["_used"]]
        n_used = int(len(g_used))

        excluded_counts = g.loc[~g["_used"], "_excluded_reason"].value_counts().to_dict()
        e = g_used[score_col].dropna().astype(float).to_numpy()
        cl = g_used[cluster_col].dropna().astype(int).to_numpy()

        if n_used > 0 and len(e) > 0:
            eint_median = float(np.median(e))
            eint_iqr = float(np.percentile(e, 75) - np.percentile(e, 25))
            eint_mean = float(np.mean(e))
            eint_std = float(np.std(e, ddof=1)) if len(e) > 1 else float("nan")
            eint_p10 = float(np.percentile(e, 10))
            eint_p90 = float(np.percentile(e, 90))
        else:
            eint_median = eint_iqr = eint_mean = eint_std = eint_p10 = eint_p90 = float("nan")

        if len(cl) > 0:
            vals, cnts = np.unique(cl, return_counts=True)
            top_idx = int(np.argmax(cnts))
            cluster_top = int(vals[top_idx])
            f1 = float(cnts[top_idx] / len(cl))
            probs = cnts / cnts.sum()
            h = float(-np.sum(probs * np.log(probs)))
            neff = float(np.exp(h))
        else:
            cluster_top = np.nan
            f1 = h = neff = float("nan")

        group_rows.append({
            "protein_id": protein_id,
            "peptide_id": peptide_id,
            "n_detected": n_detected,
            "n_used": n_used,
            "fail_rate": float(1.0 - (n_used / n_detected)) if n_detected else float("nan"),
            "n_excluded_openmm_missing": int(excluded_counts.get("openmm_missing", 0)),
            "n_excluded_legacy": int(excluded_counts.get("legacy", 0)),
            "n_excluded_other": int(excluded_counts.get("other", 0)),
            "Eint_median": eint_median,
            "Eint_IQR": eint_iqr,
            "Eint_mean": eint_mean,
            "Eint_std": eint_std,
            "Eint_p10": eint_p10,
            "Eint_p90": eint_p90,
            "cluster_mode_top": cluster_top,
            "f1_cluster": f1,
            "H_cluster": h,
            "Neff_cluster": neff,
            "flag_low_n": bool(n_used < min_replicas),
            "flag_no_data": bool(n_used == 0),
        })

    group_summary = pd.DataFrame(group_rows, columns=PHASE2_GROUP_SUMMARY_FIELDS)
    group_summary.to_csv(outdir / "group_summary.csv", index=False) #, float_format="%.2f")

    # Set-level analysis uses metadata from peptides.csv and keeps the existing
    # peptide-level outputs unchanged. The primary analysis is sequence-balanced:
    # each control peptide contributes one summary statistic so large control sets
    # do not dominate solely because they contain more raw replicas.
    df_nonlegacy = df[~df["peptide_id"].astype(str).str.endswith(exclude_peptide_suffix)].copy()
    df_nonlegacy["_peptide_set"] = df_nonlegacy["peptide_id"].map(peptide_to_set)
    missing_meta_ids = sorted(
        df_nonlegacy.loc[df_nonlegacy["_peptide_set"].isna(), "peptide_id"].dropna().astype(str).unique().tolist()
    )
    if missing_meta_ids:
        set_warnings.append(
            "Excluded from set-level analysis due to missing peptides.csv metadata: "
            + ", ".join(missing_meta_ids[:20])
            + (" ..." if len(missing_meta_ids) > 20 else "")
        )
    df_set = df_nonlegacy[df_nonlegacy["_peptide_set"].notna()].copy()
    df_set["_peptide_set"] = df_set["_peptide_set"].astype(str)
    df_set_used = df_set[df_set["_used"]].copy()

    peptide_summary = group_summary.copy()
    peptide_summary = peptide_summary[~peptide_summary["peptide_id"].astype(str).str.endswith(exclude_peptide_suffix)].copy()
    peptide_summary["_peptide_set"] = peptide_summary["peptide_id"].map(peptide_to_set)
    peptide_summary = peptide_summary[peptide_summary["_peptide_set"].notna()].copy()

    set_rows = []
    if not df_set.empty:
        detected = (
            df_set.groupby(["protein_id", "_peptide_set"])["peptide_id"]
            .agg(lambda s: sorted(set(map(str, s.dropna().tolist()))))
            .to_dict()
        )
        used_members = (
            df_set_used.groupby(["protein_id", "_peptide_set"])["peptide_id"]
            .agg(lambda s: sorted(set(map(str, s.dropna().tolist()))))
            .to_dict()
        )
        for key_tuple in sorted(detected.keys()):
            protein_id, peptide_set = key_tuple
            detected_peptides = detected.get(key_tuple, [])
            used_peptides = used_members.get(key_tuple, [])
            pep_sub = peptide_summary[
                (peptide_summary["protein_id"] == protein_id) & (peptide_summary["_peptide_set"] == peptide_set)
            ].copy()
            pep_sub = pep_sub.sort_values("peptide_id")
            raw_sub = df_set_used[(df_set_used["protein_id"] == protein_id) & (df_set_used["_peptide_set"] == peptide_set)].copy()
            seq_e = pep_sub["Eint_median"].dropna().astype(float).to_numpy()
            raw_e = raw_sub[score_col].dropna().astype(float).to_numpy()
            seq_f1 = pep_sub["f1_cluster"].dropna().astype(float).to_numpy()
            seq_neff = pep_sub["Neff_cluster"].dropna().astype(float).to_numpy()
            pool_f1, pool_neff = _cluster_metrics_from_values(raw_sub[cluster_col].dropna().astype(float).to_numpy())
            low_n_set = (len(used_peptides) < min_set_peptides) or (len(raw_sub) < min_set_replicas)
            set_rows.append({
                "protein_id": protein_id,
                "peptide_set": peptide_set,
                "n_peptides_detected": len(detected_peptides),
                "n_peptides_used": len(used_peptides),
                "n_replicas_used": int(len(raw_sub)),
                "member_peptides": ";".join(used_peptides),
                "flag_low_n_set": bool(low_n_set),
                "Eint_seq_median": float(np.median(seq_e)) if len(seq_e) else np.nan,
                "Eint_seq_IQR": float(np.percentile(seq_e, 75) - np.percentile(seq_e, 25)) if len(seq_e) else np.nan,
                "Eint_seq_mean": float(np.mean(seq_e)) if len(seq_e) else np.nan,
                "Eint_seq_std": float(np.std(seq_e, ddof=1)) if len(seq_e) > 1 else np.nan,
                "Eint_raw_median": float(np.median(raw_e)) if len(raw_e) else np.nan,
                "Eint_raw_IQR": float(np.percentile(raw_e, 75) - np.percentile(raw_e, 25)) if len(raw_e) else np.nan,
                "Eint_raw_mean": float(np.mean(raw_e)) if len(raw_e) else np.nan,
                "Eint_raw_std": float(np.std(raw_e, ddof=1)) if len(raw_e) > 1 else np.nan,
                "f1_seq_median": float(np.median(seq_f1)) if len(seq_f1) else np.nan,
                "f1_seq_IQR": float(np.percentile(seq_f1, 75) - np.percentile(seq_f1, 25)) if len(seq_f1) else np.nan,
                "Neff_seq_median": float(np.median(seq_neff)) if len(seq_neff) else np.nan,
                "Neff_seq_IQR": float(np.percentile(seq_neff, 75) - np.percentile(seq_neff, 25)) if len(seq_neff) else np.nan,
                "f1_raw_pool": pool_f1,
                "Neff_raw_pool": pool_neff,
            })
    set_summary = pd.DataFrame(set_rows, columns=SET_SUMMARY_FIELDS)
    set_summary.to_csv(outdir / "set_summary.csv", index=False)

    gidx = {(r["protein_id"], r["peptide_id"]): r for r in group_rows}
    discr_rows = []
    rng = np.random.default_rng(seed)

    for protein_id, gp in df.groupby("protein_id", sort=True):
        gp_nonlegacy = gp[~gp["peptide_id"].astype(str).str.endswith(exclude_peptide_suffix)]
        poi_used = gp_nonlegacy[(gp_nonlegacy["peptide_id"] == poi_id) & (gp_nonlegacy["_used"])][score_col].dropna().astype(float).to_numpy()
        if len(poi_used) == 0:
            continue

        for control_id in sorted(gp_nonlegacy["peptide_id"].dropna().astype(str).unique().tolist()):
            if control_id == poi_id:
                continue
            ctrl_used = gp_nonlegacy[(gp_nonlegacy["peptide_id"] == control_id) & (gp_nonlegacy["_used"])][score_col].dropna().astype(float).to_numpy()
            n_poi = int(len(poi_used))
            n_ctrl = int(len(ctrl_used))
            low_n = (n_poi < min_replicas) or (n_ctrl < min_replicas)
            skipped = bool(low_n and not allow_low_n_comparisons)

            poi_summary = gidx.get((protein_id, poi_id), {})
            ctrl_summary = gidx.get((protein_id, control_id), {})
            base = {
                "protein_id": protein_id, "poi_id": poi_id, "control_id": control_id,
                "n_poi": n_poi, "n_control": n_ctrl,
                "median_poi": float(np.median(poi_used)) if n_poi else np.nan,
                "median_control": float(np.median(ctrl_used)) if n_ctrl else np.nan,
                "delta_median": (float(np.median(poi_used) - np.median(ctrl_used)) if (n_poi and n_ctrl) else np.nan),
                "auc": np.nan, "cliffs_delta": np.nan, "p_mwu": np.nan, "q_fdr": np.nan, "p_perm": np.nan,
                "f1_cluster_poi": poi_summary.get("f1_cluster", np.nan),
                "neff_cluster_poi": poi_summary.get("Neff_cluster", np.nan),
                "f1_cluster_control": ctrl_summary.get("f1_cluster", np.nan),
                "neff_cluster_control": ctrl_summary.get("Neff_cluster", np.nan),
                "delta_f1": (poi_summary.get("f1_cluster", np.nan) - ctrl_summary.get("f1_cluster", np.nan)
                             if (poi_summary and ctrl_summary) else np.nan),
                "delta_neff": (poi_summary.get("Neff_cluster", np.nan) - ctrl_summary.get("Neff_cluster", np.nan)
                               if (poi_summary and ctrl_summary) else np.nan),
                "flag_low_n": bool(low_n),
                "flag_skipped": bool(skipped),
            }

            if skipped or n_poi == 0 or n_ctrl == 0:
                discr_rows.append(base)
                continue

            mwu = mannwhitneyu(poi_used, ctrl_used, alternative="two-sided", method="auto")
            _, auc, cliffs = _compute_auc_and_cliffs(poi_used, ctrl_used)
            base["auc"] = float(auc)
            base["cliffs_delta"] = float(cliffs)
            base["p_mwu"] = float(mwu.pvalue)
            if permutation_tests:
                base["p_perm"] = float(_perm_pvalue_delta_median(poi_used, ctrl_used, n_perm=n_perm, rng=rng))
            discr_rows.append(base)

    discr_df = pd.DataFrame(discr_rows, columns=PHASE2_DISCRIMINATION_FIELDS)
    if not discr_df.empty:
        out_pieces = []
        for protein_id, sub in discr_df.groupby("protein_id", sort=True):
            sub = sub.copy()
            mask = sub["p_mwu"].notna()
            if mask.any():
                qvals = _bh_fdr(sub.loc[mask, "p_mwu"].astype(float).tolist())
                sub.loc[mask, "q_fdr"] = qvals
            out_pieces.append(sub)
        discr_df = pd.concat(out_pieces, ignore_index=True)
    discr_df.to_csv(outdir / "protein_discrimination.csv", index=False) #, float_format="%.2f")

    set_discr_rows = []
    set_lookup = {}
    if not set_summary.empty:
        set_lookup = {
            (str(r["protein_id"]), str(r["peptide_set"])): r
            for _, r in set_summary.iterrows()
        }

    for protein_id, gp in df_set.groupby("protein_id", sort=True):
        gp_used = gp[gp["_used"]].copy()
        poi_candidates = gp_used[gp_used["peptide_id"] == poi_id]["peptide_id"].dropna().astype(str).unique().tolist()
        if not poi_candidates:
            if poi_metadata_ids:
                poi_candidates = [p for p in poi_metadata_ids if p in gp_used["peptide_id"].astype(str).unique().tolist()]
        if not poi_candidates:
            continue
        poi_peptide_id = sorted(poi_candidates)[0]
        poi_set = str(peptide_to_set.get(poi_peptide_id, ""))
        poi_raw = gp_used[gp_used["peptide_id"] == poi_peptide_id][score_col].dropna().astype(float).to_numpy()
        poi_clusters = gp_used[gp_used["peptide_id"] == poi_peptide_id][cluster_col].dropna().astype(float).to_numpy()
        if len(poi_raw) == 0:
            continue
        poi_summary_row = gidx.get((protein_id, poi_peptide_id), {})
        poi_seq_summary = float(np.median(poi_raw))
        poi_f1 = float(poi_summary_row.get("f1_cluster", np.nan))
        poi_neff = float(poi_summary_row.get("Neff_cluster", np.nan))

        sets_present = sorted(gp["_peptide_set"].dropna().astype(str).unique().tolist())
        for control_set in sets_present:
            if poi_set and control_set == poi_set:
                continue
            control_members = gp_used[(gp_used["_peptide_set"] == control_set) & (gp_used["peptide_id"] != poi_peptide_id)].copy()
            control_peptides = sorted(control_members["peptide_id"].dropna().astype(str).unique().tolist())
            set_info = set_lookup.get((protein_id, control_set), {})
            low_n_set = bool(
                (len(control_peptides) < min_set_peptides) or (len(control_members) < min_set_replicas)
            )
            skipped = bool(len(control_peptides) == 0 or len(control_members) == 0)
            base = {
                "protein_id": protein_id,
                "poi_id": poi_peptide_id,
                "control_set": control_set,
                "n_poi_replicas": int(len(poi_raw)),
                "n_control_peptides": int(len(control_peptides)),
                "n_control_replicas_total": int(len(control_members)),
                "poi_Eint_seq_summary": poi_seq_summary,
                "control_Eint_seq_median": np.nan,
                "control_Eint_seq_IQR": np.nan,
                "delta_seq_median": np.nan,
                "poi_rank_among_control_seq": np.nan,
                "poi_percentile_in_control_seq": np.nan,
                "auc_seq": np.nan,
                "p_empirical_seq": np.nan,
                "poi_f1": poi_f1,
                "control_f1_seq_median": np.nan,
                "delta_f1_seq": np.nan,
                "poi_Neff": poi_neff,
                "control_Neff_seq_median": np.nan,
                "delta_Neff_seq": np.nan,
                "n_poi_replicas_raw": int(len(poi_raw)),
                "n_control_replicas_raw": int(len(control_members)),
                "poi_raw_median": float(np.median(poi_raw)) if len(poi_raw) else np.nan,
                "control_raw_median": float(np.median(control_members[score_col].dropna().astype(float))) if len(control_members) else np.nan,
                "delta_raw_median": np.nan,
                "auc_raw": np.nan,
                "cliffs_delta_raw": np.nan,
                "p_mwu_raw": np.nan,
                "q_fdr_raw": np.nan,
                "p_perm_raw": np.nan,
                "delta_f1_rawpool": np.nan,
                "delta_Neff_rawpool": np.nan,
                "flag_low_n_set": low_n_set,
                "flag_skipped": skipped,
            }
            if skipped:
                set_discr_rows.append(base)
                continue
            control_seq = (
                peptide_summary[
                    (peptide_summary["protein_id"] == protein_id)
                    & (peptide_summary["_peptide_set"] == control_set)
                    & (peptide_summary["peptide_id"] != poi_peptide_id)
                ]
                .copy()
                .sort_values("peptide_id")
            )
            control_seq_vals = control_seq["Eint_median"].dropna().astype(float).to_numpy()
            control_seq_f1 = control_seq["f1_cluster"].dropna().astype(float).to_numpy()
            control_seq_neff = control_seq["Neff_cluster"].dropna().astype(float).to_numpy()
            control_raw = control_members[score_col].dropna().astype(float).to_numpy()
            if len(control_seq_vals):
                n_better_or_equal = int(np.sum(control_seq_vals >= poi_seq_summary))
                n_strict_better = int(np.sum(control_seq_vals < poi_seq_summary))
                percentile_seq = float(n_better_or_equal / len(control_seq_vals))
                base["control_Eint_seq_median"] = float(np.median(control_seq_vals))
                base["control_Eint_seq_IQR"] = float(np.percentile(control_seq_vals, 75) - np.percentile(control_seq_vals, 25))
                base["delta_seq_median"] = float(poi_seq_summary - np.median(control_seq_vals))
                base["poi_rank_among_control_seq"] = int(n_strict_better + 1)
                base["poi_percentile_in_control_seq"] = percentile_seq
                base["auc_seq"] = percentile_seq
                base["p_empirical_seq"] = float((np.sum(control_seq_vals <= poi_seq_summary) + 1.0) / (len(control_seq_vals) + 1.0))
                base["control_f1_seq_median"] = float(np.median(control_seq_f1)) if len(control_seq_f1) else np.nan
                base["delta_f1_seq"] = float(poi_f1 - np.median(control_seq_f1)) if len(control_seq_f1) and math.isfinite(poi_f1) else np.nan
                base["control_Neff_seq_median"] = float(np.median(control_seq_neff)) if len(control_seq_neff) else np.nan
                base["delta_Neff_seq"] = float(poi_neff - np.median(control_seq_neff)) if len(control_seq_neff) and math.isfinite(poi_neff) else np.nan
            if len(control_raw):
                base["delta_raw_median"] = float(np.median(poi_raw) - np.median(control_raw))
                mwu = mannwhitneyu(poi_raw, control_raw, alternative="two-sided", method="auto")
                _, auc_raw, cliffs_raw = _compute_auc_and_cliffs(poi_raw, control_raw)
                base["auc_raw"] = float(auc_raw)
                base["cliffs_delta_raw"] = float(cliffs_raw)
                base["p_mwu_raw"] = float(mwu.pvalue)
                if permutation_tests:
                    base["p_perm_raw"] = float(_perm_pvalue_delta_median(poi_raw, control_raw, n_perm=n_perm, rng=rng))
                control_pool_f1 = _as_float(set_info.get("f1_raw_pool", np.nan))
                control_pool_neff = _as_float(set_info.get("Neff_raw_pool", np.nan))
                if math.isfinite(control_pool_f1) and math.isfinite(poi_f1):
                    base["delta_f1_rawpool"] = float(poi_f1 - control_pool_f1)
                if math.isfinite(control_pool_neff) and math.isfinite(poi_neff):
                    base["delta_Neff_rawpool"] = float(poi_neff - control_pool_neff)
            set_discr_rows.append(base)

    set_discr_df = pd.DataFrame(set_discr_rows, columns=PROTEIN_SET_DISCRIMINATION_FIELDS)
    if not set_discr_df.empty:
        pieces = []
        for protein_id, sub in set_discr_df.groupby("protein_id", sort=True):
            sub = sub.copy()
            mask = sub["p_mwu_raw"].notna()
            if mask.any():
                sub.loc[mask, "q_fdr_raw"] = _bh_fdr(sub.loc[mask, "p_mwu_raw"].astype(float).tolist())
            pieces.append(sub)
        set_discr_df = pd.concat(pieces, ignore_index=True)
    set_discr_df.to_csv(outdir / "protein_set_discrimination.csv", index=False)

    set_rank_rows = []
    for protein_id, sub in set_discr_df.groupby("protein_id", sort=True):
        sub2 = sub[sub["flag_skipped"] == False].copy()
        if sub2.empty:
            continue
        sub2["_priority"] = sub2["control_set"].map(lambda s: _control_set_priority(str(s))[0])
        sub2["_primary_p"] = pd.to_numeric(sub2["p_empirical_seq"], errors="coerce").fillna(1.0)
        sub2["_secondary_q"] = pd.to_numeric(sub2["q_fdr_raw"], errors="coerce").fillna(1.0)
        sub2["_percentile_sort"] = pd.to_numeric(sub2["poi_percentile_in_control_seq"], errors="coerce").fillna(-1.0)
        best = sub2.sort_values(
            ["_priority", "_primary_p", "_secondary_q", "_percentile_sort"],
            ascending=[True, True, True, False],
        ).iloc[0]
        primary_percentile = float(best["poi_percentile_in_control_seq"]) if pd.notna(best["poi_percentile_in_control_seq"]) else np.nan
        primary_empirical_p = float(best["p_empirical_seq"]) if pd.notna(best["p_empirical_seq"]) else np.nan
        secondary_auc_raw = float(best["auc_raw"]) if pd.notna(best["auc_raw"]) else np.nan
        secondary_q_fdr_raw = float(best["q_fdr_raw"]) if pd.notna(best["q_fdr_raw"]) else np.nan
        delta_f1_seq = float(best["delta_f1_seq"]) if pd.notna(best["delta_f1_seq"]) else np.nan
        delta_f1_rawpool = float(best["delta_f1_rawpool"]) if pd.notna(best["delta_f1_rawpool"]) else np.nan
        verdict = "WEAK_SET_SIGNAL"
        if bool(best.get("flag_low_n_set", False)):
            verdict = "CONFOUNDED_SET_SIGNAL"
        elif (
            math.isfinite(primary_percentile) and primary_percentile >= SET_VERDICT_THRESHOLDS["strong_percentile"]
            and math.isfinite(primary_empirical_p) and primary_empirical_p <= SET_VERDICT_THRESHOLDS["strong_empirical_p"]
            and math.isfinite(secondary_auc_raw) and secondary_auc_raw >= SET_VERDICT_THRESHOLDS["strong_auc_raw"]
            and math.isfinite(delta_f1_seq) and delta_f1_seq >= SET_VERDICT_THRESHOLDS["strong_delta_f1"]
        ):
            verdict = "STRONG_SET_SIGNAL"
        elif (
            math.isfinite(primary_percentile) and primary_percentile >= SET_VERDICT_THRESHOLDS["moderate_percentile"]
            and math.isfinite(primary_empirical_p) and primary_empirical_p <= SET_VERDICT_THRESHOLDS["moderate_empirical_p"]
            and (
                (math.isfinite(secondary_auc_raw) and secondary_auc_raw >= SET_VERDICT_THRESHOLDS["moderate_auc_raw"])
                or (math.isfinite(delta_f1_seq) and delta_f1_seq >= SET_VERDICT_THRESHOLDS["moderate_delta_f1"])
            )
        ):
            verdict = "MODERATE_SET_SIGNAL"
        elif (
            (math.isfinite(secondary_q_fdr_raw) and secondary_q_fdr_raw < alpha and math.isfinite(primary_percentile) and primary_percentile < 0.5)
            or (math.isfinite(delta_f1_seq) and delta_f1_seq < 0)
        ):
            verdict = "CONFOUNDED_SET_SIGNAL"
        set_rank_rows.append({
            "protein_id": protein_id,
            "best_control_set": best["control_set"],
            "primary_percentile": primary_percentile,
            "primary_empirical_p": primary_empirical_p,
            "secondary_auc_raw": secondary_auc_raw,
            "secondary_q_fdr_raw": secondary_q_fdr_raw,
            "delta_f1_seq": delta_f1_seq,
            "delta_f1_rawpool": delta_f1_rawpool,
            "verdict_set_level": verdict,
        })
    set_rank_df = pd.DataFrame(set_rank_rows, columns=PROTEIN_SET_RANK_FIELDS)
    set_rank_df.to_csv(outdir / "protein_set_rank.csv", index=False)

    global_rows = []
    for protein_id, sub in discr_df.groupby("protein_id", sort=True):
        sub2 = sub[sub["flag_skipped"] == False].copy()
        if sub2.empty:
            continue
        sub2["_q_sort"] = sub2["q_fdr"].fillna(1.0)
        sub2["_auc_sort"] = sub2["auc"].fillna(-1.0)
        best = sub2.sort_values(["_q_sort", "_auc_sort"], ascending=[True, False]).iloc[0]

        poi_g = group_summary[(group_summary["protein_id"] == protein_id) & (group_summary["peptide_id"] == poi_id)]
        fail_rate_poi = float(poi_g.iloc[0]["fail_rate"]) if not poi_g.empty else np.nan
        ctrl_fail = group_summary[(group_summary["protein_id"] == protein_id) & (group_summary["peptide_id"] != poi_id)]["fail_rate"].dropna().astype(float)
        fail_rate_controls_median = float(ctrl_fail.median()) if not ctrl_fail.empty else np.nan

        auc_best = float(best["auc"]) if pd.notna(best["auc"]) else np.nan
        q_best = float(best["q_fdr"]) if pd.notna(best["q_fdr"]) else np.nan
        delta_f1_best = float(best["delta_f1"]) if pd.notna(best["delta_f1"]) else np.nan

        verdict = "WEAK"
        polya = sub2[sub2["control_id"].astype(str).str.lower() == "polya"]
        if not polya.empty:
            p = polya.iloc[0]
            if (pd.notna(p["auc"]) and float(p["auc"]) < VERDICT_THRESHOLDS["polya_auc_confounded"]) or                (pd.notna(p["delta_median"]) and float(p["delta_median"]) > VERDICT_THRESHOLDS["confounded_delta_median"]):
                verdict = "CONFOUNDED"

        if verdict != "CONFOUNDED":
            if (pd.notna(q_best) and q_best < alpha and pd.notna(auc_best) and auc_best >= VERDICT_THRESHOLDS["strong_auc"] and
                pd.notna(delta_f1_best) and delta_f1_best >= VERDICT_THRESHOLDS["strong_delta_f1"] and
                pd.notna(fail_rate_poi) and fail_rate_poi <= VERDICT_THRESHOLDS["poi_fail_rate_max"]):
                verdict = "STRONG"
            elif (pd.notna(q_best) and q_best < alpha and
                  ((pd.notna(auc_best) and auc_best >= VERDICT_THRESHOLDS["moderate_auc"]) or
                   (pd.notna(delta_f1_best) and delta_f1_best >= VERDICT_THRESHOLDS["moderate_delta_f1"]))):
                verdict = "MODERATE"
            elif (pd.isna(q_best) or q_best >= alpha or
                  (pd.notna(auc_best) and VERDICT_THRESHOLDS["weak_auc_low"] <= auc_best <= VERDICT_THRESHOLDS["weak_auc_high"])):
                verdict = "WEAK"

        global_rows.append({
            "protein_id": protein_id,
            "best_control_type": best["control_id"],
            "auc_best": auc_best,
            "q_best": q_best,
            "delta_f1_best": delta_f1_best,
            "fail_rate_poi": fail_rate_poi,
            "fail_rate_controls_median": fail_rate_controls_median,
            "verdict": verdict,
        })

    global_df = pd.DataFrame(global_rows, columns=PHASE2_GLOBAL_RANK_FIELDS)
    global_df.to_csv(outdir / "global_rank_by_protein.csv", index=False)

    parse_report = {}
    if parse_report_json.exists():
        try:
            parse_report = json.loads(parse_report_json.read_text(encoding="utf-8"))
        except Exception:
            parse_report = {}

    lines = []
    lines.append("# Phase 2 Metrics Report")
    lines.append("")
    lines.append("## Run metadata")
    lines.append(f"- exp_root: `{exp_root}`")
    lines.append(f"- outdir: `{outdir}`")
    lines.append(f"- timestamp: `{now_iso()}`")
    lines.append(f"- poi_id: `{poi_id}`")
    lines.append(f"- min_replicas: `{min_replicas}`")
    lines.append(f"- min_set_peptides: `{min_set_peptides}`")
    lines.append(f"- min_set_replicas: `{min_set_replicas}`")
    lines.append(f"- allow_low_n_comparisons: `{allow_low_n_comparisons}`")
    lines.append(f"- alpha: `{alpha}`")
    lines.append(f"- fdr_method: `{fdr_method}`")
    lines.append(f"- permutation_tests: `{permutation_tests}`")
    if permutation_tests:
        lines.append(f"- n_perm: `{n_perm}`")
        lines.append(f"- seed: `{seed}`")
    lines.append("")
    lines.append("## Set-level analysis summary")
    lines.append("")
    lines.append(f"- peptides.csv: `{peptides_csv}`")
    lines.append(f"- POI metadata candidates from peptide_set=={poi_id!r}: {', '.join(poi_metadata_ids) if poi_metadata_ids else 'none'}")
    lines.append(f"- set rows: `{len(set_summary)}`; protein-set comparisons: `{len(set_discr_df)}`; ranked proteins: `{len(set_rank_df)}`")
    if missing_meta_ids:
        lines.append(f"- Missing peptide_set metadata excluded from set-level analysis: {', '.join(missing_meta_ids)}")
    if set_warnings:
        for warning in set_warnings:
            lines.append(f"- WARNING: {warning}")
    lines.append("")
    if not set_summary.empty:
        lines.append("| protein_id | peptide_set | n_peptides_detected | n_peptides_used | n_replicas_used | flag_low_n_set | member_peptides |")
        lines.append("|---|---|---:|---:|---:|---:|---|")
        for _, r in set_summary.sort_values(["protein_id", "peptide_set"]).iterrows():
            lines.append(
                f"| {r['protein_id']} | {r['peptide_set']} | {int(r['n_peptides_detected'])} | {int(r['n_peptides_used'])} | {int(r['n_replicas_used'])} | {bool(r['flag_low_n_set'])} | {r['member_peptides']} |"
            )
        lines.append("")

    lines.append("## Data coverage")
    lines.append("")
    lines.append("| protein_id | peptide_id | n_detected | n_used | fail_rate | flag_low_n | flag_no_data |")
    lines.append("|---|---:|---:|---:|---:|---:|---:|")
    for _, r in group_summary.sort_values(["protein_id", "peptide_id"]).iterrows():
        lines.append(f"| {r['protein_id']} | {r['peptide_id']} | {int(r['n_detected'])} | {int(r['n_used'])} | {float(r['fail_rate']):.3f} | {bool(r['flag_low_n'])} | {bool(r['flag_no_data'])} |")
    lines.append("")
    lines.append(f"Skipped legacy groups (suffix `{exclude_peptide_suffix}`): {', '.join(skipped_legacy_groups) if skipped_legacy_groups else 'none'}")
    lines.append("")

    verdict_map = {r['protein_id']: r for r in global_rows}
    lines.append("## Per-protein discrimination")
    for protein_id in sorted(df["protein_id"].dropna().astype(str).unique().tolist()):
        lines.append("")
        lines.append(f"### {protein_id}")
        poi = group_summary[(group_summary["protein_id"] == protein_id) & (group_summary["peptide_id"] == poi_id)]
        if poi.empty:
            lines.append(f"- WARNING: POI group `{poi_id}` not present or has no usable replicas.")
            continue
        pr = poi.iloc[0]
        lines.append(f"- POI summary: Eint_median={pr['Eint_median']:.3f}, IQR={pr['Eint_IQR']:.3f}, f1_cluster={pr['f1_cluster']:.3f}, Neff={pr['Neff_cluster']:.3f}, fail_rate={pr['fail_rate']:.3f}")
        sub = discr_df[discr_df["protein_id"] == protein_id].copy()
        if sub.empty:
            lines.append("- No control comparisons available.")
            continue
        sub["_qsort"] = sub["q_fdr"].fillna(1.0)
        sub["_asort"] = sub["auc"].fillna(-1.0)
        sub = sub.sort_values(["_qsort", "_asort"], ascending=[True, False])
        lines.append("| control_id | n_poi | n_control | median_poi | median_control | delta_median | auc | q_fdr | delta_f1 | delta_neff | low_n | skipped |")
        lines.append("|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|")
        for _, r in sub.iterrows():
            lines.append(
                f"| {r['control_id']} | {int(r['n_poi'])} | {int(r['n_control'])} | {r['median_poi']:.3f} | {r['median_control']:.3f} | {r['delta_median']:.3f} | {r['auc']:.3f} | {r['q_fdr']:.3g} | {r['delta_f1']:.3f} | {r['delta_neff']:.3f} | {bool(r['flag_low_n'])} | {bool(r['flag_skipped'])} |"
            )
        gv = verdict_map.get(protein_id)
        if gv:
            lines.append(f"- Best candidate control: `{gv['best_control_type']}`; verdict: **{gv['verdict']}**")

    lines.append("")
    lines.append("## Set-level analysis")
    lines.append("")
    lines.append("Primary set-level analysis is sequence-balanced: each control peptide contributes one median winner_omm_dE_interaction value, which prevents large pooled control sets from dominating by raw replica count alone.")
    lines.append("Secondary set-level analysis pools raw replicas across the control set for a higher-power but potentially control-heavy sensitivity view.")
    for protein_id in sorted(set_discr_df["protein_id"].dropna().astype(str).unique().tolist()) if not set_discr_df.empty else []:
        lines.append("")
        lines.append(f"### {protein_id}")
        sub = set_discr_df[set_discr_df["protein_id"] == protein_id].copy()
        sub["_priority"] = sub["control_set"].map(lambda s: _control_set_priority(str(s))[0])
        sub = sub.sort_values(["_priority", "p_empirical_seq", "q_fdr_raw"], ascending=[True, True, True])
        lines.append("| control_set | n_control_peptides | n_control_replicas_total | poi_Eint_seq_summary | control_Eint_seq_median | control_Eint_seq_IQR | poi_percentile_in_control_seq | p_empirical_seq | delta_f1_seq | auc_raw | q_fdr_raw | delta_f1_rawpool | skipped |")
        lines.append("|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|")
        for _, r in sub.iterrows():
            lines.append(
                f"| {r['control_set']} | {int(r['n_control_peptides'])} | {int(r['n_control_replicas_total'])} | {r['poi_Eint_seq_summary']:.3f} | {r['control_Eint_seq_median']:.3f} | {r['control_Eint_seq_IQR']:.3f} | {r['poi_percentile_in_control_seq']:.3f} | {r['p_empirical_seq']:.3g} | {r['delta_f1_seq']:.3f} | {r['auc_raw']:.3f} | {r['q_fdr_raw']:.3g} | {r['delta_f1_rawpool']:.3f} | {bool(r['flag_skipped'])} |"
            )
        rv = set_rank_df[set_rank_df["protein_id"] == protein_id]
        if not rv.empty:
            rr = rv.iloc[0]
            lines.append(
                f"- Best set-level comparison: `{rr['best_control_set']}`; verdict: **{rr['verdict_set_level']}** "
                f"(primary_percentile={rr['primary_percentile']:.3f}, primary_empirical_p={rr['primary_empirical_p']:.3g}, secondary_auc_raw={rr['secondary_auc_raw']:.3f}, secondary_q_fdr_raw={rr['secondary_q_fdr_raw']:.3g})."
            )

    lines.append("")
    lines.append("## QA exclusions")
    excl = df["_excluded_reason"].value_counts().to_dict()
    lines.append(f"- openmm_missing: {int(excl.get('openmm_missing', 0))}")
    lines.append(f"- legacy: {int(excl.get('legacy', 0))}")
    lines.append(f"- other: {int(excl.get('other', 0))}")

    failures = parse_report.get("replica_failures", []) if isinstance(parse_report, dict) else []
    if failures:
        lines.append("- Top 10 replica failures from parse_report:")
        for item in failures[:10]:
            key = item.get("key", {})
            lines.append(f"  - {key.get('protein_id','?')}/{key.get('peptide_id','?')}/replica_{key.get('replica_id','?')}: {item.get('reason','')} ")

    lines.append("")
    lines.append("## Interpretation note")
    lines.append("These metrics are protocol-discrimination diagnostics and are **not** binding free energies (not ΔG estimates).")
    lines.append("Set-level primary percentile values are reported on a 0-1 scale, where larger values mean the POI is more favorable (more negative/better) than a larger fraction of control-sequence medians.")
    lines.append("The primary set-level p_empirical_seq is an empirical percentile-based tail probability surrogate for the sequence-balanced comparison, not a standard two-sample p-value.")

    write_text(outdir / "metrics_report.md", "\n".join(lines) + "\n")

    run_info = {}
    if run_info_json.exists():
        try:
            run_info = json.loads(run_info_json.read_text(encoding="utf-8"))
        except Exception:
            run_info = {}
    run_info["phase2"] = {
        "timestamp": now_iso(),
        "enabled": True,
        "params": {
            "poi_id": poi_id,
            "exclude_peptide_suffix": exclude_peptide_suffix,
            "min_replicas": min_replicas,
            "allow_low_n_comparisons": allow_low_n_comparisons,
            "alpha": alpha,
            "fdr_method": fdr_method,
            "permutation_tests": permutation_tests,
            "n_perm": n_perm,
            "seed": seed,
            "min_set_peptides": min_set_peptides,
            "min_set_replicas": min_set_replicas,
            "sanity": asdict(sanity_cfg),
        },
        "outputs": {
            "group_summary_csv": str(outdir / "group_summary.csv"),
            "protein_discrimination_csv": str(outdir / "protein_discrimination.csv"),
            "global_rank_by_protein_csv": str(outdir / "global_rank_by_protein.csv"),
            "set_summary_csv": str(outdir / "set_summary.csv"),
            "protein_set_discrimination_csv": str(outdir / "protein_set_discrimination.csv"),
            "protein_set_rank_csv": str(outdir / "protein_set_rank.csv"),
            "metrics_report_md": str(outdir / "metrics_report.md"),
        },
    }
    with run_info_json.open("w", encoding="utf-8") as f:
        json.dump(run_info, f, indent=2, sort_keys=True)

    vlog(1, verbose, "[INFO] Phase 2 outputs written: group_summary.csv, protein_discrimination.csv, global_rank_by_protein.csv, set_summary.csv, protein_set_discrimination.csv, protein_set_rank.csv, metrics_report.md")


STATUS_ORDER = {"PASS": 0, "WARN": 1, "FAIL": 2}


def _combine_status(*statuses: str) -> str:
    return max(statuses, key=lambda s: STATUS_ORDER.get(s, 2)) if statuses else "PASS"


def _as_float(v) -> float:
    try:
        return float(v)
    except Exception:
        return float("nan")


def _count_model_lines_streaming(path: Path) -> int:
    n = 0
    with path.open("r", encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith("MODEL"):
                n += 1
    return n


def _phase2_checklist_report(
    outdir: Path,
    exp_root: Path,
    parse_report_json: Path,
    parse_report_txt: Path,
    metrics_report_md: Path,
    poi_id: str,
    exclude_peptide_suffix: str,
    sanity_cfg: SanityConfig,
    topk_k: int,
    min_replicas: int,
    alpha: float,
    checklist_max_detail: int,
    min_set_peptides: int,
    min_set_replicas: int,
    verbose: int = 0,
) -> Tuple[str, str]:
    np, pd, _ = _import_phase2_deps()
    eps = 1e-6
    ts = now_iso()

    def _fmt_replica_label(replica_id: object) -> str:
        rid = str(replica_id)
        try:
            rid = f"{int(rid):03d}"
        except Exception:
            pass
        return f"replica_{rid}"

    group_csv = outdir / "group_summary.csv"
    discr_csv = outdir / "protein_discrimination.csv"
    set_summary_csv = outdir / "set_summary.csv"
    set_discr_csv = outdir / "protein_set_discrimination.csv"
    topk_index_csv = outdir / "topk_concat_index.csv"
    topk_parsed_csv = outdir / "topk_parsed.csv"
    replicas_csv = outdir / "replicas_parsed.csv"

    for req in [group_csv, discr_csv, set_summary_csv, set_discr_csv, topk_index_csv, topk_parsed_csv, replicas_csv]:
        if not req.exists():
            raise FileNotFoundError(f"Checklist requires {req}")

    group = pd.read_csv(group_csv)
    discr = pd.read_csv(discr_csv)
    set_summary = pd.read_csv(set_summary_csv)
    set_discr = pd.read_csv(set_discr_csv)
    topk_idx = pd.read_csv(topk_index_csv)
    topk_parsed = pd.read_csv(topk_parsed_csv)
    replicas = pd.read_csv(replicas_csv)

    parse_groups = {}
    if parse_report_json.exists():
        try:
            parse_obj = json.loads(parse_report_json.read_text(encoding="utf-8"))
            parse_groups = parse_obj.get("groups_discovered", {}) if isinstance(parse_obj, dict) else {}
        except Exception as ex:
            vlog(1, verbose, f"[WARN] Checklist could not read parse_report.json: {ex}")
    elif parse_report_txt.exists():
        pat = re.compile(r"^\s*(\S+/\S+):\s+detected=(\d+)\s+summary=(\d+)\s+rescored=(\d+)\s+out=(\d+)\s*$")
        in_block = False
        for line in parse_report_txt.read_text(encoding="utf-8", errors="replace").splitlines():
            if line.strip() == "Groups discovered (protein/peptide):":
                in_block = True
                continue
            if in_block:
                m = pat.match(line)
                if not m:
                    if line.strip() == "":
                        continue
                    break
                parse_groups[m.group(1)] = {
                    "replicas_detected": int(m.group(2)),
                    "replicas_with_summary": int(m.group(3)),
                    "replicas_with_rescored": int(m.group(4)),
                    "replicas_with_out": int(m.group(5)),
                }

    group["protein_id"] = group["protein_id"].astype(str)
    group["peptide_id"] = group["peptide_id"].astype(str)
    discr["protein_id"] = discr.get("protein_id", pd.Series(dtype=str)).astype(str)
    if "poi_id" in discr.columns:
        discr["poi_id"] = discr["poi_id"].astype(str)
    if "control_id" in discr.columns:
        discr["control_id"] = discr["control_id"].astype(str)
    for df in [topk_idx, topk_parsed, replicas]:
        if "protein_id" in df.columns:
            df["protein_id"] = df["protein_id"].astype(str)
        if "peptide_id" in df.columns:
            df["peptide_id"] = df["peptide_id"].astype(str)
        if "replica_id" in df.columns:
            df["replica_id"] = df["replica_id"].astype(str)

    def _truncate_path(path_txt: str, max_len: int = 120) -> str:
        ptxt = str(path_txt or "").strip()
        if len(ptxt) <= max_len:
            return ptxt
        head_len = min(30, max_len - 3)
        tail_len = max_len - 3 - head_len
        if tail_len <= 0:
            return ptxt[:max_len]
        return f"{ptxt[:head_len]}...{ptxt[-tail_len:]}"

    group_idx = {(r["protein_id"], r["peptide_id"]): r for _, r in group.iterrows()}

    details_1 = []
    details_2 = []
    details_3 = []
    details_4 = []
    details_5 = []
    details_6 = []

    # 1A
    s1a = "PASS"
    fail_1a = []
    warn_1a = []
    for _, r in group.iterrows():
        pid, pep = str(r["protein_id"]), str(r["peptide_id"])
        n_det = int(_as_float(r.get("n_detected", float("nan"))) or 0)
        n_used = int(_as_float(r.get("n_used", float("nan"))) or 0)
        fr = _as_float(r.get("fail_rate", float("nan")))
        if n_used < 0 or n_det < 0 or n_used > n_det:
            fail_1a.append(f"{pid}/{pep}: n_used={n_used}, n_detected={n_det}")
        if math.isfinite(fr) and (fr < -eps or fr > 1.0 + eps):
            fail_1a.append(f"{pid}/{pep}: fail_rate={fr:.6f} outside [0,1]")
        if n_det == 0 and n_used != 0:
            fail_1a.append(f"{pid}/{pep}: n_detected=0 but n_used={n_used}")
        if pep.endswith(exclude_peptide_suffix) and n_used != 0:
            fail_1a.append(f"{pid}/{pep}: legacy group has n_used={n_used} (expected 0)")

    gs_keys = {f"{p}/{q}" for p, q in group[["protein_id", "peptide_id"]].itertuples(index=False, name=None)}
    parse_keys = set(parse_groups.keys())
    missing_in_group = sorted(parse_keys - gs_keys)
    if missing_in_group:
        warn_1a.append(f"Groups in parse_report but missing in group_summary: {', '.join(missing_in_group[:checklist_max_detail])}")

    if fail_1a:
        s1a = "FAIL"
        details_1.append("1A FAIL: group invariants violated.")
        details_1.extend([f"- {x}" for x in fail_1a[:checklist_max_detail]])
    elif warn_1a:
        s1a = "WARN"
        details_1.append("1A WARN: minor coverage differences.")
        details_1.extend([f"- {x}" for x in warn_1a[:checklist_max_detail]])
    else:
        details_1.append("1A PASS: group_summary invariants and parse coverage look consistent.")

    # 1B
    s1b = "PASS"
    fail_1b = []
    warn_1b = []
    fail_1b_offenders = []
    proteins_with_poi = sorted(group[(group["peptide_id"] == poi_id) & (pd.to_numeric(group["n_used"], errors="coerce").fillna(0) > 0)]["protein_id"].unique().tolist())
    for prot in proteins_with_poi:
        sub = discr[(discr["protein_id"] == prot) & (discr.get("poi_id", "") == poi_id)]
        if sub.empty:
            fail_1b.append(f"{prot}: POI present but no discrimination rows found")
            continue
        if bool(sub["flag_skipped"].fillna(False).astype(bool).all()):
            warn_1b.append(f"{prot}: all POI comparisons skipped by low-n constraints")

    for _, r in discr.iterrows():
        pid, poi, ctrl = str(r.get("protein_id", "")), str(r.get("poi_id", "")), str(r.get("control_id", ""))
        if poi.endswith(exclude_peptide_suffix) or ctrl.endswith(exclude_peptide_suffix):
            fail_1b.append(f"{pid}: legacy group appears in discrimination ({poi} vs {ctrl})")
        g_poi = group_idx.get((pid, poi))
        g_ctrl = group_idx.get((pid, ctrl))
        if g_poi is None or g_ctrl is None:
            fail_1b.append(f"{pid}: missing group_summary row for comparison {poi} vs {ctrl}")
            continue
        n_poi = int(_as_float(r.get("n_poi", 0)) or 0)
        n_ctrl = int(_as_float(r.get("n_control", 0)) or 0)
        exp_poi = int(_as_float(g_poi.get("n_used", 0)) or 0)
        exp_ctrl = int(_as_float(g_ctrl.get("n_used", 0)) or 0)
        if n_poi != exp_poi or n_ctrl != exp_ctrl:
            fail_1b.append(f"{pid}:{poi} vs {ctrl} n_poi/n_control=({n_poi},{n_ctrl}) expected ({exp_poi},{exp_ctrl})")
            fail_1b_offenders.append({
                "protein_id": pid,
                "poi_id": poi,
                "control_id": ctrl,
                "n_poi": n_poi,
                "n_used_poi": exp_poi,
                "n_control": n_ctrl,
                "n_used_control": exp_ctrl,
            })

    if fail_1b:
        s1b = "FAIL"
        details_1.append("1B FAIL: join inconsistencies between group_summary and discrimination.")
        details_1.extend([f"- {x}" for x in fail_1b[:checklist_max_detail]])
    elif warn_1b:
        s1b = "WARN"
        details_1.append("1B WARN: proteins with POI but all comparisons skipped (low-n).")
        details_1.extend([f"- {x}" for x in warn_1b[:checklist_max_detail]])
    else:
        details_1.append("1B PASS: discrimination joins and legacy exclusions are coherent.")

    s1 = _combine_status(s1a, s1b)

    # 2A + 2B
    s2a = "PASS"
    mild_2a = []
    fail_2a = []
    fail_2a_offenders = []
    active = discr[~discr.get("flag_skipped", False).fillna(False).astype(bool)].copy() if not discr.empty else discr.copy()
    for _, r in active.iterrows():
        pid = str(r["protein_id"])
        ctrl = str(r["control_id"])
        mp, mc = _as_float(r.get("median_poi")), _as_float(r.get("median_control"))
        delta, auc = _as_float(r.get("delta_median")), _as_float(r.get("auc"))
        if not (math.isfinite(mp) and math.isfinite(mc) and math.isfinite(delta) and math.isfinite(auc)):
            continue
        if abs(mp - mc) <= eps:
            if abs(delta) > 1e-3 or abs(auc - 0.5) > 0.05:
                mild_2a.append(f"{pid}/{ctrl}: tie medians but delta={delta:.3g}, auc={auc:.3f}")
            continue
        poi_better = mp < mc
        if poi_better and delta > eps:
            fail_2a.append(f"{pid}/{ctrl}: median_poi < median_control but delta_median={delta:.3g} > 0")
        if (not poi_better) and delta < -eps:
            fail_2a.append(f"{pid}/{ctrl}: median_poi > median_control but delta_median={delta:.3g} < 0")
        if poi_better and auc < 0.45:
            fail_2a.append(f"{pid}/{ctrl}: POI better by medians but auc={auc:.3f} < 0.45")
            fail_2a_offenders.append({"protein_id": pid, "control_id": ctrl, "median_poi": mp, "median_control": mc, "delta_median": delta, "auc": auc})
        elif poi_better and auc < 0.5:
            mild_2a.append(f"{pid}/{ctrl}: POI better by medians but auc={auc:.3f} mildly contradictory")
        if (not poi_better) and auc > 0.55:
            fail_2a.append(f"{pid}/{ctrl}: POI worse by medians but auc={auc:.3f} > 0.55")
            fail_2a_offenders.append({"protein_id": pid, "control_id": ctrl, "median_poi": mp, "median_control": mc, "delta_median": delta, "auc": auc})
        elif (not poi_better) and auc > 0.5:
            mild_2a.append(f"{pid}/{ctrl}: POI worse by medians but auc={auc:.3f} mildly contradictory")

    if fail_2a:
        s2a = "FAIL"
        details_2.append("2A FAIL: effect-direction contradictions detected.")
        details_2.extend([f"- {x}" for x in fail_2a[:checklist_max_detail]])
    elif mild_2a:
        s2a = "WARN"
        details_2.append("2A WARN: mild AUC/median orientation contradictions.")
        details_2.extend([f"- {x}" for x in mild_2a[:checklist_max_detail]])
    else:
        details_2.append("2A PASS: effect direction is coherent (medians, delta_median, AUC).")

    s2b = "PASS"
    fail_2b = []
    for _, r in discr.iterrows():
        auc, cd = _as_float(r.get("auc")), _as_float(r.get("cliffs_delta"))
        if math.isfinite(auc) and math.isfinite(cd):
            if abs(cd - (2.0 * auc - 1.0)) > 1e-3:
                fail_2b.append(f"{r.get('protein_id','?')}/{r.get('control_id','?')}: cliffs={cd:.6f}, 2*auc-1={(2.0*auc-1.0):.6f}")
    if fail_2b:
        s2b = "FAIL"
        details_2.append("2B FAIL: Cliff's delta != (2*auc-1).")
        details_2.extend([f"- {x}" for x in fail_2b[:checklist_max_detail]])
    else:
        details_2.append("2B PASS: Cliff's delta consistency checks passed.")
    s2 = _combine_status(s2a, s2b)

    # 3A + 3B
    s3a = "PASS"
    fail_3a = []
    non_skipped = discr[~discr.get("flag_skipped", False).fillna(False).astype(bool)].copy() if not discr.empty else discr.copy()
    for _, r in non_skipped.iterrows():
        p, q = _as_float(r.get("p_mwu")), _as_float(r.get("q_fdr"))
        pid, ctrl = str(r.get("protein_id", "?")), str(r.get("control_id", "?"))
        if not math.isfinite(p) or not (0.0 - eps <= p <= 1.0 + eps):
            fail_3a.append(f"{pid}/{ctrl}: invalid p_mwu={p}")
        if not math.isfinite(q) or not (0.0 - eps <= q <= 1.0 + eps):
            fail_3a.append(f"{pid}/{ctrl}: invalid q_fdr={q}")
        if math.isfinite(p) and math.isfinite(q) and q + 1e-9 < p:
            fail_3a.append(f"{pid}/{ctrl}: q_fdr={q:.6g} < p_mwu={p:.6g}")
    if fail_3a:
        s3a = "FAIL"
        details_3.append("3A FAIL: p/q bounds or BH monotonic property violated.")
        details_3.extend([f"- {x}" for x in fail_3a[:checklist_max_detail]])
    else:
        details_3.append("3A PASS: p-values and q-values are in range and q>=p.")

    s3b = "PASS"
    fail_3b = []
    fail_3b_offenders = []
    for pid, sub in non_skipped.groupby("protein_id", sort=True):
        sub = sub.copy()
        mask = sub["p_mwu"].notna()
        if not bool(mask.any()):
            continue
        pvals = sub.loc[mask, "p_mwu"].astype(float).tolist()
        q_expected = _bh_fdr(pvals)
        q_observed = sub.loc[mask, "q_fdr"].astype(float).tolist()
        for idx, (qo, qe) in enumerate(zip(q_observed, q_expected), start=1):
            if abs(qo - qe) > 1e-6:
                fail_3b.append(f"{pid}: row#{idx} q_fdr={qo:.6g} expected_BH={qe:.6g}")
                fail_3b_offenders.append({"protein_id": pid, "row": idx, "p_mwu": pvals[idx - 1], "q_reported": qo, "q_expected": qe})
    if fail_3b:
        s3b = "FAIL"
        details_3.append("3B FAIL: FDR appears not applied per protein with BH.")
        details_3.extend([f"- {x}" for x in fail_3b[:checklist_max_detail]])
    else:
        details_3.append("3B PASS: per-protein BH recomputation matches q_fdr.")
    s3 = _combine_status(s3a, s3b)

    # 4
    s4 = "PASS"
    fail_4 = []
    for _, r in group.iterrows():
        n_used = int(_as_float(r.get("n_used", 0)) or 0)
        pid, pep = str(r["protein_id"]), str(r["peptide_id"])
        f1, neff = _as_float(r.get("f1_cluster")), _as_float(r.get("Neff_cluster"))
        if n_used > 0:
            if not (math.isfinite(f1) and -eps <= f1 <= 1.0 + eps):
                fail_4.append(f"{pid}/{pep}: f1_cluster={f1} outside [0,1]")
            if not (math.isfinite(neff) and 1.0 - eps <= neff <= n_used + eps):
                fail_4.append(f"{pid}/{pep}: Neff_cluster={neff} outside [1,n_used={n_used}]")
    for _, r in discr.iterrows():
        f1p, f1c, d1 = _as_float(r.get("f1_cluster_poi")), _as_float(r.get("f1_cluster_control")), _as_float(r.get("delta_f1"))
        np_p, np_c, dn = _as_float(r.get("neff_cluster_poi")), _as_float(r.get("neff_cluster_control")), _as_float(r.get("delta_neff"))
        pid, ctrl = str(r.get("protein_id", "?")), str(r.get("control_id", "?"))
        if math.isfinite(f1p) and math.isfinite(f1c) and math.isfinite(d1) and abs(d1 - (f1p - f1c)) > 1e-6:
            fail_4.append(f"{pid}/{ctrl}: delta_f1={d1:.6f} != f1_poi-f1_ctrl={(f1p-f1c):.6f}")
        if math.isfinite(np_p) and math.isfinite(np_c) and math.isfinite(dn) and abs(dn - (np_p - np_c)) > 1e-6:
            fail_4.append(f"{pid}/{ctrl}: delta_neff={dn:.6f} != neff_poi-neff_ctrl={(np_p-np_c):.6f}")
    if fail_4:
        s4 = "FAIL"
        details_4.append("4 FAIL: cluster convergence metrics are inconsistent.")
        details_4.extend([f"- {x}" for x in fail_4[:checklist_max_detail]])
    else:
        details_4.append("4 PASS: f1/Neff ranges and deltas are consistent.")

    # 5A + 5B
    s5a = "PASS"
    warn_5a, fail_5a = [], []
    has_underfilled_group_5a = False

    replica_topk_counts = {}
    if not topk_parsed.empty:
        replica_topk_counts = {
            (pid, pep, rid): int(cnt)
            for (pid, pep, rid), cnt in topk_parsed.groupby(["protein_id", "peptide_id", "replica_id"]).size().items()
        }

    replicas_by_group = {}
    replica_path_by_key = {}
    expected_winners = {}
    if not replicas.empty:
        score_col = _resolve_col(replicas, SCORE_SYNONYMS)
        cluster_col = _resolve_col(replicas, CLUSTER_SYNONYMS)
        qa_col = _resolve_col(replicas, QA_SYNONYMS)
        if score_col and cluster_col and qa_col:
            replicas[score_col] = pd.to_numeric(replicas[score_col], errors="coerce")
            replicas[cluster_col] = pd.to_numeric(replicas[cluster_col], errors="coerce")
            replicas["_excluded_reason"] = replicas.apply(
                lambda rr: _excluded_reason(rr, score_col, cluster_col, qa_col, exclude_peptide_suffix, sanity_cfg), axis=1
            )
            for (pid, pep), sub in replicas.groupby(["protein_id", "peptide_id"], sort=False):
                if str(pep).endswith(exclude_peptide_suffix):
                    continue
                used_sub = sub[sub["_excluded_reason"] == "used"].copy()
                replicas_by_group[(pid, pep)] = used_sub
                exp = int(used_sub["rescored_pdb_path"].astype(str).str.strip().ne("").sum())
                expected_winners[(pid, pep)] = exp
                for _, rr in used_sub.iterrows():
                    replica_path_by_key[(pid, pep, str(rr.get("replica_id", "")))] = str(rr.get("rescored_pdb_path", ""))

    details_5.append("5A policy: PDB concat mismatch policy: WARN if observed<expected; FAIL if observed>expected or file missing.")
    topk_group_counts = {}
    if not topk_parsed.empty:
        topk_group_counts = {
            (pid, pep): int(cnt)
            for (pid, pep), cnt in topk_parsed.groupby(["protein_id", "peptide_id"]).size().items()
        }

    concat_model_count_by_group = {}
    legacy_index_rows = []
    for _, r in topk_idx.iterrows():
        pid, pep = str(r.get("protein_id", "")), str(r.get("peptide_id", ""))
        if pep.endswith(exclude_peptide_suffix):
            legacy_index_rows.append(f"{pid}/{pep}")
            continue
        g = group_idx.get((pid, pep))
        if g is None:
            fail_5a.append(f"{pid}/{pep}: missing group_summary row for topk index")
            continue
        n_used = int(_as_float(g.get("n_used", 0)) or 0)
        kpr = int(_as_float(r.get("k_per_replica", 0)) or 0)
        expected = kpr * n_used
        indexed_obs = int(_as_float(r.get("n_models_total", 0)) or 0)
        path = Path(str(r.get("concat_pdb_path", "")))
        if not path.exists():
            fail_5a.append(f"{pid}/{pep}: concat PDB missing ({path})")
            continue
        streamed = _count_model_lines_streaming(path)
        observed = streamed
        concat_model_count_by_group[(pid, pep)] = observed
        if indexed_obs != streamed:
            warn_5a.append(f"{pid}/{pep}: index n_models_total={indexed_obs}, MODEL count={streamed}; using streamed count")
        if observed < expected:
            group_missing = expected - observed
            has_underfilled_group_5a = True
            details_5.append(f"- {pid}/{pep}: expected={expected}, observed={observed}, missing={group_missing} (WARN)")
            if kpr > 0 and (pid, pep) in replicas_by_group:
                incomplete = []
                for _, rr in replicas_by_group[(pid, pep)].iterrows():
                    rid = str(rr.get("replica_id", ""))
                    c = int(replica_topk_counts.get((pid, pep, rid), 0))
                    if c < kpr:
                        ptxt = replica_path_by_key.get((pid, pep, rid), "")
                        pp = Path(str(ptxt)) if str(ptxt).strip() else None
                        incomplete.append((rid, c, ptxt, pp.exists() if pp else False, bool(pp)))
                incomplete.sort(key=lambda x: (x[1], x[0]))
                listed = incomplete[:checklist_max_detail]
                sum_missing = 0
                has_missing_file = any((not has_path) or (has_path and not exists) for _, _, _, exists, has_path in incomplete)
                for rid, c, ptxt, exists, has_path in listed:
                    miss = kpr - c
                    sum_missing += miss
                    extra_path = f"; rescored_pdb={_truncate_path(ptxt)}" if str(ptxt).strip() else ""
                    details_5.append(f"  - {_fmt_replica_label(rid)}: {c}/{kpr} (missing {miss}){extra_path}")
                details_5.append(f"  - Missing models accounted for by listed replicas: {sum_missing}/{group_missing}")
                remaining_missing = group_missing - sum_missing
                if remaining_missing > 0:
                    details_5.append(f"  - ... plus {remaining_missing} missing models not shown (increase --checklist-max-detail)")
                if has_missing_file:
                    details_5.append("  - Likely cause: missing rescored PDB for some replicas (OpenMM stage missing).")
                else:
                    details_5.append("  - Likely cause: rescored PDB contains fewer than k MODEL blocks for one or more replicas (OpenMM stage produced <k outputs or output truncated).")
            else:
                details_5.append("  - Likely cause: rescored PDB contains fewer than k MODEL blocks for one or more replicas (OpenMM stage produced <k outputs or output truncated).")
        elif observed > expected:
            fail_5a.append(f"{pid}/{pep}: expected={expected}, observed={observed} (unexpected extra models)")

    concat_topk_crosscheck_warnings = []
    for _, r in topk_idx.iterrows():
        pid, pep = str(r.get("protein_id", "")), str(r.get("peptide_id", ""))
        if pep.endswith(exclude_peptide_suffix):
            continue
        observed_models_topk = int(topk_group_counts.get((pid, pep), 0))
        observed_models_concat = concat_model_count_by_group.get((pid, pep))
        if observed_models_concat is None:
            continue
        if observed_models_concat != observed_models_topk:
            msg = (
                f"{pid}/{pep}: concat MODEL count ({observed_models_concat}) differs from topk_parsed rows ({observed_models_topk}). Possible export bug; inspect concat generation."
            )
            warn_5a.append(msg)
            concat_topk_crosscheck_warnings.append(msg)

    if legacy_index_rows:
        details_5.append("5A WARN: Legacy export files detected from prior runs; ignored in current-run validation.")
        details_5.extend([f"- ignored legacy concat index row: {x}" for x in legacy_index_rows[:checklist_max_detail]])
    if not concat_topk_crosscheck_warnings:
        details_5.append("5A PASS: concat MODEL counts match topk_parsed row counts for all concat groups.")
    if fail_5a:
        s5a = "FAIL"
        details_5.append("5A FAIL: topk_concat export structural mismatches.")
        details_5.extend([f"- {x}" for x in fail_5a[:checklist_max_detail]])
    elif has_underfilled_group_5a or warn_5a:
        s5a = "WARN"
        if has_underfilled_group_5a:
            details_5.append("5A WARN: topk_concat has fewer models than expected in some groups.")
        else:
            details_5.append("5A WARN: topk_concat cross-check warnings detected.")
        details_5.extend([f"- {x}" for x in warn_5a[:checklist_max_detail]])
        if has_underfilled_group_5a:
            details_5.append("- Action: if needed, rerun replicas with incomplete rescored MODEL blocks.")
    else:
        details_5.append("5A PASS: topk_concat MODEL counts match expected k_per_replica*n_used.")

    s5b = "PASS"
    warn_5b, fail_5b = [], []
    for _, g in group.iterrows():
        pid, pep = str(g["protein_id"]), str(g["peptide_id"])
        if pep.endswith(exclude_peptide_suffix):
            continue
        n_used = int(_as_float(g.get("n_used", 0)) or 0)
        if n_used <= 0:
            continue
        winner_dir = outdir / "poses_winners" / pid / pep
        exp = expected_winners.get((pid, pep), n_used)
        if not winner_dir.exists():
            fail_5b.append(f"{pid}/{pep}: poses_winners directory missing but n_used={n_used}")
            continue
        obs = len(list(winner_dir.glob("replica_*_model1.pdb")))
        if obs < exp:
            warn_5b.append(f"{pid}/{pep}: winners observed={obs}, expected~={exp} (possible missing exports)")
    if fail_5b:
        s5b = "FAIL"
        details_5.append("5B FAIL: winner export directories missing where expected.")
        details_5.extend([f"- {x}" for x in fail_5b[:checklist_max_detail]])
    elif warn_5b:
        s5b = "WARN"
        details_5.append("5B WARN: fewer winner exports than expected in some groups.")
        details_5.extend([f"- {x}" for x in warn_5b[:checklist_max_detail]])
        details_5.append("- Action: inspect missing replicas and rerun winner export if required.")
    else:
        details_5.append("5B PASS: winner export presence looks consistent for n_used>0 groups.")
    s5 = _combine_status(s5a, s5b)

    s6 = "PASS"
    fail_6 = []
    warn_6 = []
    peptide_to_set = {}
    if exp_root.exists():
        _, _, peptide_to_set, _ = _load_peptide_metadata(exp_root, poi_id)
    nonlegacy = replicas[~replicas["peptide_id"].astype(str).str.endswith(exclude_peptide_suffix)].copy()
    missing_map = sorted(nonlegacy.loc[~nonlegacy["peptide_id"].astype(str).isin(set(peptide_to_set.keys())), "peptide_id"].dropna().astype(str).unique().tolist())
    if missing_map:
        warn_6.append(f"Non-legacy peptides missing peptide_set metadata: {', '.join(missing_map[:checklist_max_detail])}")
    self_comparison_found = False
    if not set_discr.empty and not set_summary.empty:
        set_idx = {(str(r["protein_id"]), str(r["peptide_set"])): r for _, r in set_summary.iterrows()}
        poi_set_by_protein = {}
        poi_rows = set_summary[set_summary["peptide_set"].astype(str).str.casefold() == str(poi_id).casefold()].copy()
        for _, poi_row in poi_rows.iterrows():
            poi_set_by_protein[str(poi_row["protein_id"])] = str(poi_row["peptide_set"])
        for _, r in set_discr.iterrows():
            pid = str(r.get("protein_id", ""))
            control_set = str(r.get("control_set", ""))
            poi_set = poi_set_by_protein.get(pid, "")
            if poi_set and control_set == poi_set:
                self_comparison_found = True
                fail_6.append(f"{pid}/{control_set}: self-comparison row should not be present in protein_set_discrimination")
                continue
            key = (pid, control_set)
            if key not in set_idx and not bool(r.get("flag_skipped", False)):
                fail_6.append(f"{pid}/{control_set}: missing set_summary row for set-level comparison")
                continue
            if key in set_idx:
                expected = int(_as_float(set_idx[key].get("n_peptides_used", 0)) or 0)
                observed = int(_as_float(r.get("n_control_peptides", 0)) or 0)
                # The control set count may exclude the POI peptide if the POI itself
                # is annotated inside the same peptide_set, so accept either exact
                # match or expected-1.
                if observed not in {expected, max(0, expected - 1)}:
                    fail_6.append(f"{pid}/{control_set}: n_control_peptides={observed} expected={expected} (or {max(0, expected - 1)} if POI is inside the set)")
            pct = _as_float(r.get("poi_percentile_in_control_seq"))
            if math.isfinite(pct) and (pct < -eps or pct > 1.0 + eps):
                fail_6.append(f"{pid}/{control_set}: percentile={pct:.6f} outside [0,1]")
            p_raw = _as_float(r.get("p_mwu_raw"))
            q_raw = _as_float(r.get("q_fdr_raw"))
            if math.isfinite(p_raw) and math.isfinite(q_raw) and q_raw + 1e-9 < p_raw:
                fail_6.append(f"{pid}/{control_set}: q_fdr_raw={q_raw:.6g} < p_mwu_raw={p_raw:.6g}")
        for pid, sub in set_discr.groupby("protein_id", sort=True):
            mask = sub["p_mwu_raw"].notna()
            if mask.any():
                expected_q = _bh_fdr(sub.loc[mask, "p_mwu_raw"].astype(float).tolist())
                observed_q = sub.loc[mask, "q_fdr_raw"].astype(float).tolist()
                for qo, qe in zip(observed_q, expected_q):
                    if abs(qo - qe) > 1e-6:
                        fail_6.append(f"{pid}: raw pooled BH mismatch q_reported={qo:.6g} expected={qe:.6g}")
                        break
        matched_present = set_summary[
            set_summary["peptide_set"].astype(str).str.contains("matched", case=False, na=False)
            & (pd.to_numeric(set_summary["n_peptides_used"], errors="coerce").fillna(0) >= min_set_peptides)
        ]
        for _, row in matched_present.iterrows():
            pid = str(row["protein_id"])
            control_set = str(row["peptide_set"])
            found = not set_discr[(set_discr["protein_id"].astype(str) == pid) & (set_discr["control_set"].astype(str) == control_set)].empty
            if not found:
                fail_6.append(f"{pid}/{control_set}: matched-like set missing from protein_set_discrimination")
    if fail_6:
        s6 = "FAIL"
        details_6.append("6 FAIL: set-level analysis integrity checks failed.")
        details_6.extend([f"- {x}" for x in fail_6[:checklist_max_detail]])
    elif warn_6:
        if not self_comparison_found:
            details_6.append("6 PASS note: self-comparisons correctly excluded from set-level discrimination.")
        s6 = "WARN"
        details_6.append("6 WARN: set-level analysis integrity has warnings.")
        details_6.extend([f"- {x}" for x in warn_6[:checklist_max_detail]])
    else:
        details_6.append("6 PASS: set-level mapping, counts, percentiles, and BH q-values look consistent.")
        details_6.append("6 PASS note: self-comparisons correctly excluded from set-level discrimination.")

    def _quick_control_label(name: object) -> str:
        ctrl = str(name or "?")
        lname = ctrl.lower()
        if "matched" in lname:
            return "matched"
        if "decoy" in lname:
            return "decoy"
        if "polya" in lname:
            return "polya"
        return ctrl

    def _fmt_set_quick_chunk(rr) -> str:
        ctrl = _quick_control_label(rr.get("control_set", "?"))
        pct = _as_float(rr.get("poi_percentile_in_control_seq"))
        p_emp = _as_float(rr.get("p_empirical_seq"))
        auc_raw = _as_float(rr.get("auc_raw"))
        q_raw = _as_float(rr.get("q_fdr_raw"))
        d1_seq = _as_float(rr.get("delta_f1_seq"))
        parts = [f"POI vs {ctrl}"]
        if math.isfinite(pct):
            parts.append(f"percentile={pct:.3f}")
        if math.isfinite(p_emp):
            parts.append(f"p_emp={p_emp:.3g}")
        if math.isfinite(auc_raw):
            parts.append(f"AUC_raw={auc_raw:.3f}")
        if math.isfinite(q_raw):
            parts.append(f"q_raw={q_raw:.3g}")
        if math.isfinite(d1_seq):
            parts.append(f"delta_f1_seq={d1_seq:+.2f}")
        return " ".join([parts[0], *parts[1:]])

    def _fmt_individual_quick_chunk(rr) -> str:
        ctrl = _quick_control_label(rr.get("control_id", "?"))
        auc = _as_float(rr.get("auc"))
        qv = _as_float(rr.get("q_fdr"))
        d1 = _as_float(rr.get("delta_f1"))
        parts = [f"POI vs {ctrl}"]
        if math.isfinite(auc):
            parts.append(f"AUC={auc:.3f}")
        if math.isfinite(qv):
            parts.append(f"q={qv:.3g}")
        if math.isfinite(d1):
            parts.append(f"delta_f1={d1:+.2f}")
        return " ".join([parts[0], *parts[1:]])

    quick_lines = []
    for prot in proteins_with_poi:
        set_sub = set_discr[(set_discr.get("protein_id", "") == prot) & (set_discr.get("poi_id", "") == poi_id)].copy()
        if not set_sub.empty:
            set_sub = set_sub[~set_sub.get("flag_skipped", False).fillna(False).astype(bool)].copy()
        indiv_sub = discr[(discr.get("protein_id", "") == prot) & (discr.get("poi_id", "") == poi_id)].copy()
        if not indiv_sub.empty:
            indiv_sub = indiv_sub[~indiv_sub.get("flag_skipped", False).fillna(False).astype(bool)].copy()

        if set_sub.empty and indiv_sub.empty:
            continue

        if not set_sub.empty:
            set_sub["control_set"] = set_sub["control_set"].astype(str)
            set_sub["p_empirical_seq_sort"] = pd.to_numeric(set_sub.get("p_empirical_seq"), errors="coerce").fillna(float("inf"))
            set_sub["q_fdr_raw_sort"] = pd.to_numeric(set_sub.get("q_fdr_raw"), errors="coerce").fillna(float("inf"))
            set_sub["percentile_sort"] = pd.to_numeric(set_sub.get("poi_percentile_in_control_seq"), errors="coerce").fillna(-1.0)
            set_sub["_priority"] = set_sub["control_set"].map(lambda s: _control_set_priority(str(s))[0])
            set_sub = set_sub.sort_values(["_priority", "p_empirical_seq_sort", "q_fdr_raw_sort", "percentile_sort", "control_set"], ascending=[True, True, True, False, True], kind="stable")

        if not indiv_sub.empty:
            indiv_sub["control_id"] = indiv_sub["control_id"].astype(str)
            indiv_sub["q_fdr_sort"] = pd.to_numeric(indiv_sub.get("q_fdr"), errors="coerce").fillna(float("inf"))
            indiv_sub["auc_sort"] = pd.to_numeric(indiv_sub.get("auc"), errors="coerce").fillna(-1.0)
            indiv_sub = indiv_sub.sort_values(["q_fdr_sort", "auc_sort", "control_id"], ascending=[True, False, True], kind="stable")

        chunks = []
        used_set_controls = set()

        def _first_set_match(pred):
            if set_sub.empty:
                return None
            matched = set_sub[set_sub["control_set"].map(pred)]
            if matched.empty:
                return None
            rr = matched.iloc[0]
            used_set_controls.add(str(rr["control_set"]))
            return _fmt_set_quick_chunk(rr)

        def _first_individual_match(pred):
            if indiv_sub.empty:
                return None
            matched = indiv_sub[indiv_sub["control_id"].map(pred)]
            if matched.empty:
                return None
            return _fmt_individual_quick_chunk(matched.iloc[0])

        primary_specs = [
            (lambda s: "matched" in s.lower(), None),
            (lambda s: "decoy" in s.lower(), lambda s: "decoy" in s.lower()),
            (lambda s: "polya" in s.lower(), lambda s: s.lower() == "polya"),
        ]
        for set_pred, indiv_pred in primary_specs:
            chunk = _first_set_match(set_pred)
            if chunk is None and indiv_pred is not None:
                chunk = _first_individual_match(indiv_pred)
            if chunk is not None:
                chunks.append(chunk)

        if not set_sub.empty:
            for _, rr in set_sub.iterrows():
                ctrl = str(rr.get("control_set", ""))
                if ctrl in used_set_controls:
                    continue
                chunks.append(_fmt_set_quick_chunk(rr))

        if chunks:
            quick_lines.append(f"{prot}: " + "; ".join(chunks))

    groups_in_group_summary = int(len(group))
    groups_analyzed_nonlegacy = int(((~group["peptide_id"].astype(str).str.endswith(exclude_peptide_suffix)) & (pd.to_numeric(group["n_used"], errors="coerce").fillna(0) > 0)).sum())
    comparisons_total = int(len(discr))
    comparisons_skipped = int(discr.get("flag_skipped", False).fillna(False).astype(bool).sum()) if not discr.empty else 0
    comparisons_run = int(max(0, comparisons_total - comparisons_skipped))
    excluded_openmm_missing = int(pd.to_numeric(group.get("n_excluded_openmm_missing", 0), errors="coerce").fillna(0).sum())
    excluded_legacy = int(pd.to_numeric(group.get("n_excluded_legacy", 0), errors="coerce").fillna(0).sum())
    excluded_other = int(pd.to_numeric(group.get("n_excluded_other", 0), errors="coerce").fillna(0).sum())
    replicas_total_detected = int(pd.to_numeric(group.get("n_detected", 0), errors="coerce").fillna(0).sum())
    replicas_total_used = int(pd.to_numeric(group.get("n_used", 0), errors="coerce").fillna(0).sum())

    offender_lines = []
    if _combine_status(s1, s2, s3, s4, s5, s6) == "FAIL":
        offender_lines.append("Top offenders (first N):")
        if s1b == "FAIL" and fail_1b_offenders:
            offender_lines.append("- 1B join inconsistencies:")
            for it in fail_1b_offenders[:checklist_max_detail]:
                offender_lines.append(f"  - {it['protein_id']}, {it['poi_id']} vs {it['control_id']}, n_poi={it['n_poi']}, n_used_poi={it['n_used_poi']}, n_control={it['n_control']}, n_used_control={it['n_used_control']}")
            offender_lines.append("  Action: revisar filtros de inclusión o merge keys.")
        if s2a == "FAIL" and fail_2a_offenders:
            offender_lines.append("- 2A orientation contradictions:")
            for it in fail_2a_offenders[:checklist_max_detail]:
                offender_lines.append(f"  - {it['protein_id']}/{it['control_id']}: median_poi={it['median_poi']:.3f}, median_control={it['median_control']:.3f}, delta_median={it['delta_median']:.3f}, auc={it['auc']:.3f}")
            offender_lines.append("  Action: revisar bug de orientación en AUC.")
        if s3b == "FAIL" and fail_3b_offenders:
            offender_lines.append("- 3B FDR per-protein BH mismatches:")
            for it in fail_3b_offenders[:checklist_max_detail]:
                offender_lines.append(f"  - {it['protein_id']} row#{it['row']}: p={it['p_mwu']:.3g}, q_reported={it['q_reported']:.3g}, q_expected={it['q_expected']:.3g}")
            offender_lines.append("  Action: recomputar BH por proteína y sobrescribir q_fdr por bloque protein_id.")

    overall = _combine_status(s1, s2, s3, s4, s5, s6)

    stdout_lines = [
        f"Script: {SCRIPT_NAME} {SCRIPT_VERSION}",
        "===== Phase 2 Checklist (dockanalrep) =====",
        f"Timestamp: {ts}",
        f"Outdir: {outdir}",
        f"Config: k={topk_k}, min_replicas={min_replicas}, min_set_peptides={min_set_peptides}, min_set_replicas={min_set_replicas}, alpha={alpha}, FDR=BH-per-protein, legacy_suffix=\"{exclude_peptide_suffix}\"",
        f"1) Cobertura y coherencia básica: {s1}",
        *[f"   {x}" for x in details_1],
        f"2) Sentido de métricas: {s2}",
        *[f"   {x}" for x in details_2],
        f"3) Significancia y FDR: {s3}",
        *[f"   {x}" for x in details_3],
        f"4) Convergencia por clusters: {s4}",
        *[f"   {x}" for x in details_4],
        f"5) Verificación de exports PDB: {s5}",
        *[f"   {x}" for x in details_5],
        f"6) Set-level analysis integrity: {s6}",
        *[f"   {x}" for x in details_6],
        f"Coverage: groups_in_group_summary={groups_in_group_summary}, groups_analyzed_nonlegacy={groups_analyzed_nonlegacy}, comparisons_total={comparisons_total}, comparisons_run={comparisons_run}, comparisons_skipped={comparisons_skipped}",
        f"Replicas: total_detected={replicas_total_detected}, total_used={replicas_total_used}, excluded_legacy={excluded_legacy}, excluded_openmm_missing={excluded_openmm_missing}, excluded_other={excluded_other}",
        "Protein summaries (quick):",
        *[f"   - {x}" for x in quick_lines],
        *[f"   {x}" for x in offender_lines],
        f"Overall: {overall}",
        "==========================================",
    ]

    md_lines = [
        f"Script: `{SCRIPT_NAME} {SCRIPT_VERSION}`",
        "## Phase 2 Checklist",
        f"_Run: {ts}_",
        f"(outdir: `{outdir}`)",
        f"Config: `k={topk_k}, min_replicas={min_replicas}, min_set_peptides={min_set_peptides}, min_set_replicas={min_set_replicas}, alpha={alpha}, FDR=BH-per-protein, legacy_suffix=\"{exclude_peptide_suffix}\"`",
        "",
        f"### 1) Cobertura y coherencia básica — {s1}",
        *details_1,
        "",
        f"### 2) Sentido de métricas — {s2}",
        *details_2,
        "",
        f"### 3) Significancia y FDR — {s3}",
        *details_3,
        "",
        f"### 4) Convergencia por clusters — {s4}",
        *details_4,
        "",
        f"### 5) Verificación de exports PDB — {s5}",
        *details_5,
        "",
        f"### 6) Set-level analysis integrity — {s6}",
        *details_6,
        "",
        f"Coverage: groups_in_group_summary={groups_in_group_summary}, groups_analyzed_nonlegacy={groups_analyzed_nonlegacy}, comparisons_total={comparisons_total}, comparisons_run={comparisons_run}, comparisons_skipped={comparisons_skipped}",
        f"Replicas: total_detected={replicas_total_detected}, total_used={replicas_total_used}, excluded_legacy={excluded_legacy}, excluded_openmm_missing={excluded_openmm_missing}, excluded_other={excluded_other}",
        "",
        "Protein summaries (quick):",
        *[f"- {x}" for x in quick_lines],
        "",
        *offender_lines,
        "",
        f"**Overall: {overall}**",
        "",
    ]

    checklist_stdout = "\n".join(stdout_lines)
    checklist_md = "\n".join(md_lines)
    print(checklist_stdout)

    previous = metrics_report_md.read_text(encoding="utf-8") if metrics_report_md.exists() else ""
    with metrics_report_md.open("w", encoding="utf-8") as f:
        if previous:
            f.write(previous.rstrip() + "\n\n")
        f.write(checklist_md)

    return overall, checklist_stdout


# ----------------------------
# Pipeline principal: parseo + consolidación + export PDB + QA
# ----------------------------

def run_parse_and_consolidate(exp_root: Path, outdir: Path, topk_k: int, sanity_cfg: SanityConfig, exclude_peptide_suffix: str, verbose: int = 0) -> None:
    safe_mkdir(outdir)
    vlog(1, verbose, f"[INFO] Discovering replicas under: {exp_root / 'runs'}")

    discovered, group_counts = discover_replicas(exp_root)
    vlog(1, verbose, f"[INFO] Discovered {len(discovered)} replica directories across {len(group_counts)} groups.")

    # Report skeleton
    report = {
        "timestamp": now_iso(),
        "exp_root": str(exp_root.resolve()),
        "outdir": str(outdir.resolve()),
        "topk_k": topk_k,
        "sanity": asdict(sanity_cfg),
        "groups_discovered": group_counts,
        "replicas_total_discovered": len(discovered),
        "replicas_parsed_ok": 0,
        "replicas_parsed_fail": 0,
        "replica_failures": [],  # list of {key, reason}
        "warnings": [],
    }

    replicas_rows: List[Dict[str, object]] = []
    topk_rows: List[Dict[str, object]] = []
    clusters_rows: List[Dict[str, object]] = []

    # For exports
    winners_dir = outdir / "poses_winners"
    concat_dir = outdir / "poses_topk_concat"
    safe_mkdir(winners_dir)
    safe_mkdir(concat_dir)

    # To build concatenated topk pdb per group
    group_to_concat_blocks: Dict[Tuple[str, str, str], List[str]] = {}
    group_to_pose_only_blocks: Dict[Tuple[str, str, str], List[str]] = {}
    group_to_first_target_lines: Dict[Tuple[str, str, str], List[str]] = {}
    group_to_ppi_unrestrained_residues: Dict[Tuple[str, str, str], List[Tuple[str, str, int]]] = {}
    group_to_deint_values: Dict[Tuple[str, str, str], List[float]] = {}
    split_strategy_counts: Dict[str, int] = {
        "first_ter": 0,
        "chain_fallback": 0,
        "unsplit": 0,
        "empty": 0,
    }

    # For summary group-level later (placeholder)
    parsed_replicas: List[ParsedReplica] = []

    for idx, rp in enumerate(discovered, start=1):
        key = rp.key
        qa_status = "OK"
        qa_msgs: List[str] = []
        replica_tag = f"{key.protein_id}/{key.peptide_id}/replica_{key.replica_id}"
        vlog(1, verbose, f"[INFO] [{idx}/{len(discovered)}] Processing {replica_tag}")

        if rp.summary_path is None or not rp.summary_path.exists():
            qa_status = "MISSING_SUMMARY"
            qa_msgs.append("No *_summary.dlg found for replica.")
            report["replicas_parsed_fail"] += 1
            report["replica_failures"].append({"key": asdict(key), "reason": "; ".join(qa_msgs)})
            # record minimal row for traceability (optional)
            replicas_rows.append({
                "experiment_id": key.experiment_id,
                "protein_id": key.protein_id,
                "peptide_id": key.peptide_id,
                "replica_id": key.replica_id,
                "summary_path": "",
                "rescored_pdb_path": str(rp.rescored_pdb_path) if rp.rescored_pdb_path else "",
                "out_pdb_path": str(rp.out_pdb_path) if rp.out_pdb_path else "",
                "qa_status": qa_status,
                "qa_message": "; ".join(qa_msgs),
            })
            vlog(1, verbose, f"[WARN] {replica_tag}: missing *_summary.dlg; skipping parse.")
            continue

        try:
            text = read_text(rp.summary_path)
            meta, clusters, topk, warnings = parse_summary_dlg(text)
            if warnings:
                qa_msgs.extend([f"PARSE_WARN: {w}" for w in warnings])

            if not topk:
                qa_status = "INCOMPLETE"
                qa_msgs.append("No top-k (OMM Ranking) rows found.")
        except Exception as ex:
            qa_status = "PARSE_ERROR"
            qa_msgs.append(f"Exception parsing summary: {repr(ex)}")
            report["replicas_parsed_fail"] += 1
            report["replica_failures"].append({"key": asdict(key), "reason": "; ".join(qa_msgs)})
            replicas_rows.append({
                "experiment_id": key.experiment_id,
                "protein_id": key.protein_id,
                "peptide_id": key.peptide_id,
                "replica_id": key.replica_id,
                "summary_path": str(rp.summary_path),
                "rescored_pdb_path": str(rp.rescored_pdb_path) if rp.rescored_pdb_path else "",
                "out_pdb_path": str(rp.out_pdb_path) if rp.out_pdb_path else "",
                "qa_status": qa_status,
                "qa_message": "; ".join(qa_msgs),
            })
            vlog(1, verbose, f"[WARN] {replica_tag}: parse error in summary ({repr(ex)}).")
            continue

        # Winner selection: prefer rank_openmm==1, else smallest rank_openmm
        winner_row: Optional[ParsedTopKRow] = None
        if topk:
            # rank_openmm should be 1..k; pick min rank_openmm
            winner_row = sorted(topk, key=lambda r: r.rank_openmm)[0]
            # In well-formed output, winner_row.rank_openmm == 1 and model_id == 1
            if winner_row.rank_openmm != 1:
                qa_msgs.append(f"Unexpected: winner rank_openmm is {winner_row.rank_openmm} (expected 1).")

        # Sanity filter
        if winner_row is not None:
            ok, msg = sanity_check_eint(winner_row.omm_dE_interaction, sanity_cfg)
            if not ok:
                qa_status = "SANITY_FAIL" if qa_status == "OK" else qa_status
                qa_msgs.append(f"SANITY: {msg}")

        # Export winner pose
        winner_export_path = ""
        if rp.rescored_pdb_path and rp.rescored_pdb_path.exists():
            try:
                block1 = extract_model1(rp.rescored_pdb_path)
                if block1:
                    dst = winners_dir / key.protein_id / key.peptide_id / f"replica_{key.replica_id}_model1.pdb"
                    safe_mkdir(dst.parent)
                    write_text(dst, block1)
                    winner_export_path = str(dst)
                else:
                    qa_msgs.append("Could not extract MODEL 1 from rescored PDB.")
            except Exception as ex:
                qa_msgs.append(f"Winner export failed: {repr(ex)}")
        else:
                qa_msgs.append("No rescored PDB found; skipping winner pose export.")

        # Accumulate concat top-k blocks per group (if rescored exists)
        if rp.rescored_pdb_path and rp.rescored_pdb_path.exists():
            try:
                model_entries: List[Tuple[int, str]] = []
                for idx, (mid, block) in enumerate(iter_pdb_models(read_text(rp.rescored_pdb_path)), start=1):
                    model_entries.append((mid, block))
                    if idx >= topk_k:
                        break

                if model_entries:
                    blocks = [block for _, block in model_entries]
                    gk = key.group_key()
                    group_to_concat_blocks.setdefault(gk, []).extend(blocks)

                    model_to_deint = {r.model_id: r.omm_dE_interaction for r in topk}
                    for mid, _ in model_entries:
                        if mid in model_to_deint:
                            group_to_deint_values.setdefault(gk, []).append(model_to_deint[mid])
                        else:
                            qa_msgs.append(f"No dE_Interaction found for MODEL {mid} in top-k table.")

                    pose_blocks: List[str] = []
                    for i, block in enumerate(blocks):
                        target_lines, pose_lines, split_strategy = split_model_target_and_pose(block)
                        split_strategy_counts[split_strategy] = split_strategy_counts.get(split_strategy, 0) + 1
                        if i == 0 and gk not in group_to_first_target_lines and target_lines:
                            group_to_first_target_lines[gk] = target_lines
                            residues = extract_unrestrained_receptor_residues(target_lines)
                            if residues:
                                group_to_ppi_unrestrained_residues[gk] = residues

                        if split_strategy == "chain_fallback":
                            qa_msgs.append("Target/pose split used chain-based fallback (no TER found).")
                        elif split_strategy == "unsplit":
                            qa_msgs.append("Could not split target/pose (no TER and no chain transition).")

                        if pose_lines:
                            pose_block = f"MODEL        1\n{''.join(pose_lines)}ENDMDL\n"
                            pose_blocks.append(pose_block)

                    if pose_blocks:
                        group_to_pose_only_blocks.setdefault(gk, []).extend(pose_blocks)
                    else:
                        qa_msgs.append("Could not split target/pose for pose-only concat export.")
                else:
                    qa_msgs.append("No MODEL blocks found for concat export.")
            except Exception as ex:
                qa_msgs.append(f"Concat export read failed: {repr(ex)}")

        # Map winner cluster_mode_id already computed (Definition B)
        winner_cluster_mode_id = winner_row.cluster_mode_id if winner_row else None
        if winner_row and winner_cluster_mode_id is None:
            qa_status = "CLUSTER_MAP_FAIL" if qa_status == "OK" else qa_status
            qa_msgs.append("Could not map winner_best_run_pose_id to cluster_mode_id (Definition B).")

        # Construct ParsedReplica object (for later group-level summary)
        parsed = ParsedReplica(
            experiment_id=key.experiment_id,
            protein_id=key.protein_id,
            peptide_id=key.peptide_id,
            replica_id=key.replica_id,
            summary_path=str(rp.summary_path),
            rescored_pdb_path=str(rp.rescored_pdb_path) if rp.rescored_pdb_path else "",
            out_pdb_path=str(rp.out_pdb_path) if rp.out_pdb_path else "",
            sequence=str(meta.get("sequence", "")),
            N_runs=int(meta.get("N_runs", 0) or 0),
            n_evals=int(meta.get("n_evals", 0) or 0),
            nmin=int(meta.get("nmin", 0) or 0),
            env=str(meta.get("env", "")),
            nitr=int(meta.get("nitr", 0) or 0),
            adcp_bestEnergy=float(meta.get("adcp_bestEnergy", float("nan"))),
            adcp_bestEnergy_run=int(meta.get("adcp_bestEnergy_run", 0) or 0),
            adcp_bestEnergies_list=str(meta.get("adcp_bestEnergies_list", "")),
            winner_rank_openmm=int(winner_row.rank_openmm if winner_row else 0),
            winner_model_id=int(winner_row.model_id if winner_row else 0),
            winner_omm_dE_interaction=float(winner_row.omm_dE_interaction if winner_row else float("nan")),
            winner_omm_dE_complex_minus_receptor=float(winner_row.omm_score_complex_minus_receptor if winner_row else float("nan")),
            winner_adcp_affinity=float(winner_row.adcp_affinity_kcalmol if winner_row else float("nan")),
            winner_rank_adcp=int(winner_row.rank_adcp if winner_row else 0),
            winner_cluster_size=int(winner_row.cluster_size if winner_row else 0),
            winner_best_run_pose_id=str(winner_row.best_run_pose_id if winner_row else ""),
            winner_cluster_mode_id=winner_cluster_mode_id,
            qa_status=qa_status,
            qa_message="; ".join(qa_msgs).strip(),
        )
        parsed_replicas.append(parsed)
        vlog(2, verbose, f"[DEBUG] {replica_tag}: qa_status={qa_status}, winner_model={parsed.winner_model_id}, winner_dEint={parsed.winner_omm_dE_interaction}")

        # replicas_parsed.csv row
        replicas_rows.append({
            **asdict(parsed),
            "winner_pose_export_path": winner_export_path,
        })

        # topk_parsed.csv rows
        for r in topk:
            topk_rows.append({
                "experiment_id": key.experiment_id,
                "protein_id": key.protein_id,
                "peptide_id": key.peptide_id,
                "replica_id": key.replica_id,
                "summary_path": str(rp.summary_path),
                "rescored_pdb_path": str(rp.rescored_pdb_path) if rp.rescored_pdb_path else "",
                "model_id": r.model_id,
                "rank_openmm": r.rank_openmm,
                "rank_adcp": r.rank_adcp,
                "omm_score_complex_minus_receptor": r.omm_score_complex_minus_receptor,
                "omm_dE_interaction": r.omm_dE_interaction,
                "adcp_affinity_kcalmol": r.adcp_affinity_kcalmol,
                "cluster_size": r.cluster_size,
                "best_run_pose_id": r.best_run_pose_id,
                "cluster_mode_id": r.cluster_mode_id if r.cluster_mode_id is not None else "",
                "qa_status": qa_status,
            })

        # clusters_parsed.csv rows
        for c in clusters:
            clusters_rows.append({
                "experiment_id": key.experiment_id,
                "protein_id": key.protein_id,
                "peptide_id": key.peptide_id,
                "replica_id": key.replica_id,
                "summary_path": str(rp.summary_path),
                "mode": c.mode,
                "adcp_affinity_kcalmol": c.adcp_affinity_kcalmol,
                "cluster_size": c.cluster_size,
                "best_run_pose_id": c.best_run_pose_id,
                "qa_status": qa_status,
            })

        if qa_status == "OK":
            report["replicas_parsed_ok"] += 1
            vlog(1, verbose, f"[OK]   {replica_tag}: parsed successfully.")
        else:
            report["replicas_parsed_fail"] += 1
            report["replica_failures"].append({"key": asdict(key), "reason": parsed.qa_message})
            vlog(1, verbose, f"[WARN] {replica_tag}: completed with QA status {qa_status}.")

    # Write CSVs
    replicas_fields = [
        "experiment_id", "protein_id", "peptide_id", "replica_id",
        "summary_path", "rescored_pdb_path", "out_pdb_path",
        "sequence", "N_runs", "n_evals", "nmin", "env", "nitr",
        "adcp_bestEnergy", "adcp_bestEnergy_run", "adcp_bestEnergies_list",
        "winner_rank_openmm", "winner_model_id",
        "winner_omm_dE_interaction", "winner_omm_dE_complex_minus_receptor",
        "winner_adcp_affinity", "winner_rank_adcp", "winner_cluster_size",
        "winner_best_run_pose_id", "winner_cluster_mode_id",
        "qa_status", "qa_message",
        "winner_pose_export_path",
    ]
    write_csv(outdir / "replicas_parsed.csv", replicas_rows, replicas_fields)
    vlog(1, verbose, f"[INFO] Wrote {outdir / 'replicas_parsed.csv'} ({len(replicas_rows)} rows)")

    topk_fields = [
        "experiment_id", "protein_id", "peptide_id", "replica_id",
        "summary_path", "rescored_pdb_path",
        "model_id", "rank_openmm", "rank_adcp",
        "omm_score_complex_minus_receptor", "omm_dE_interaction",
        "adcp_affinity_kcalmol", "cluster_size", "best_run_pose_id", "cluster_mode_id",
        "qa_status",
    ]
    write_csv(outdir / "topk_parsed.csv", topk_rows, topk_fields)
    vlog(1, verbose, f"[INFO] Wrote {outdir / 'topk_parsed.csv'} ({len(topk_rows)} rows)")

    clusters_fields = [
        "experiment_id", "protein_id", "peptide_id", "replica_id",
        "summary_path",
        "mode", "adcp_affinity_kcalmol", "cluster_size", "best_run_pose_id",
        "qa_status",
    ]
    write_csv(outdir / "clusters_parsed.csv", clusters_rows, clusters_fields)
    vlog(1, verbose, f"[INFO] Wrote {outdir / 'clusters_parsed.csv'} ({len(clusters_rows)} rows)")

    # Export concatenated top-k PDB per group
    concat_index_rows: List[Dict[str, object]] = []
    for gk, blocks in group_to_concat_blocks.items():
        exp_id, protein_id, peptide_id = gk
        if peptide_id.endswith(exclude_peptide_suffix):
            continue
        # renumber models from 1..N
        concat_text = renumber_models(blocks, start_at=1)
        out_path = concat_dir / protein_id / peptide_id / f"topk_concat_{protein_id}_{peptide_id}.pdb"
        safe_mkdir(out_path.parent)
        write_text(out_path, concat_text)

        pose_only_path = ""
        pose_blocks = group_to_pose_only_blocks.get(gk, [])
        if pose_blocks:
            pose_text = renumber_models(pose_blocks, start_at=1)
            pose_out_path = concat_dir / protein_id / peptide_id / f"topk_concat_poses_{protein_id}_{peptide_id}.pdb"
            write_text(pose_out_path, pose_text)
            pose_only_path = str(pose_out_path)

        first_target_path = ""
        first_target_lines = group_to_first_target_lines.get(gk, [])
        if first_target_lines:
            target_out_path = concat_dir / protein_id / peptide_id / f"topk_first_target_{protein_id}_{peptide_id}.pdb"
            write_text(target_out_path, "".join(first_target_lines))
            first_target_path = str(target_out_path)

        ppi_unrestrained_path = ""
        ppi_residues = group_to_ppi_unrestrained_residues.get(gk, [])
        if ppi_residues:
            ppi_out_path = concat_dir / protein_id / peptide_id / f"topk_first_target_PPI_{protein_id}_{peptide_id}.dat"
            write_text(ppi_out_path, format_unrestrained_receptor_residues_dat(ppi_residues))
            ppi_unrestrained_path = str(ppi_out_path)

        deint_path = ""
        deint_values = group_to_deint_values.get(gk, [])
        if deint_values:
            deint_out_path = concat_dir / protein_id / peptide_id / f"topk_omm_dEinter_{protein_id}_{peptide_id}.dat"
            write_text(deint_out_path, "\n".join(str(v) for v in deint_values) + "\n")
            deint_path = str(deint_out_path)

        concat_index_rows.append({
            "experiment_id": exp_id,
            "protein_id": protein_id,
            "peptide_id": peptide_id,
            "k_per_replica": topk_k,
            "n_models_total": len(blocks),
            "concat_pdb_path": str(out_path),
            "concat_pose_only_pdb_path": pose_only_path,
            "first_target_pdb_path": first_target_path,
            "first_target_ppi_dat_path": ppi_unrestrained_path,
            "deinter_dat_path": deint_path,
        })
    write_csv(outdir / "topk_concat_index.csv", concat_index_rows,
              [
                  "experiment_id", "protein_id", "peptide_id", "k_per_replica", "n_models_total",
                  "concat_pdb_path", "concat_pose_only_pdb_path", "first_target_pdb_path", "first_target_ppi_dat_path", "deinter_dat_path",
              ])
    vlog(1, verbose, f"[INFO] Wrote {outdir / 'topk_concat_index.csv'} ({len(concat_index_rows)} groups)")
    print(
        "[INFO] Target/pose split strategies used: "
        f"first_ter={split_strategy_counts.get('first_ter', 0)}, "
        f"chain_fallback={split_strategy_counts.get('chain_fallback', 0)}, "
        f"unsplit={split_strategy_counts.get('unsplit', 0)}, "
        f"empty={split_strategy_counts.get('empty', 0)}"
    )


    # Write parse_report.json + parse_report.txt
    report_path = outdir / "parse_report.json"
    safe_mkdir(report_path.parent)
    with report_path.open("w", encoding="utf-8") as f:
        json.dump(report, f, indent=2, sort_keys=True)

    # Human-readable short report
    lines: List[str] = []
    lines.append(f"anal_adcp_agfr_replicas_campaign parse report @ {report['timestamp']}")
    lines.append(f"exp_root: {report['exp_root']}")
    lines.append(f"outdir:   {report['outdir']}")
    lines.append(f"replicas discovered: {report['replicas_total_discovered']}")
    lines.append(f"replicas OK:         {report['replicas_parsed_ok']}")
    lines.append(f"replicas FAIL:       {report['replicas_parsed_fail']}")
    lines.append("")
    lines.append("Groups discovered (protein/peptide):")
    for gk, c in sorted(group_counts.items()):
        lines.append(f"  {gk}: detected={c['replicas_detected']} summary={c['replicas_with_summary']} rescored={c['replicas_with_rescored']} out={c['replicas_with_out']}")
    lines.append("")
    if report["replica_failures"]:
        lines.append("Replica failures (first 25):")
        for item in report["replica_failures"][:25]:
            k = item["key"]
            lines.append(f"  {k['protein_id']}/{k['peptide_id']}/replica_{k['replica_id']}: {item['reason']}")
    else:
        lines.append("No replica failures recorded.")

    write_text(outdir / "parse_report.txt", "\n".join(lines) + "\n")
    vlog(1, verbose, f"[INFO] Wrote reports: {outdir / 'parse_report.json'} and {outdir / 'parse_report.txt'}")


# ----------------------------
# CLI
# ----------------------------

def build_argparser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="anal_adcp_agfr_replicas_campaign",
        description="Parse + consolidate ADCP/OpenMM replica outputs; export PDB; QA.",
    )
    p.add_argument("--exp-root", required=True, help="Experiment root directory (contains runs/).")
    p.add_argument("--outdir", required=True, help="Output directory (e.g., analysis).")
    p.add_argument("--topk", type=int, default=5, help="How many top OpenMM-ranked models per replica to export/concat (default: 5).")
    p.add_argument("-v", "--verbose", action="count", default=0,
                   help="Increase progress verbosity (-v: info, -vv: debug).")

    # Sanity filter options
    p.add_argument("--sanity-disable", action="store_true", help="Disable sanity checks on energies.")
    p.add_argument("--sanity-eint-min", type=float, default=-1.0e6, help="Min allowed winner dE_interaction (default: -1e6).")
    p.add_argument("--sanity-eint-max", type=float, default= 1.0e6, help="Max allowed winner dE_interaction (default: 1e6).")

    # Phase 2 switches
    p.add_argument("--run-phase2", dest="run_phase2", action="store_true", default=True,
                   help="Run Phase 2 metrics/significance after Phase 1 parsing (default: enabled).")
    p.add_argument("--no-phase2", dest="run_phase2", action="store_false",
                   help="Disable Phase 2 metrics/significance.")
    p.add_argument("--poi-id", default="POI", help="Peptide ID treated as POI within each protein (default: POI).")
    p.add_argument("--exclude-peptide-suffix", default="_old",
                   help="Exclude peptide IDs ending with this suffix from Phase 2 (default: _old).")
    p.add_argument("--min-replicas", type=int, default=10,
                   help="Minimum n_used per group for unflagged comparisons (default: 10).")
    p.add_argument("--min-set-peptides", type=int, default=1,
                   help="Minimum number of control peptides required for unflagged set-level comparisons (default: 1).")
    p.add_argument("--min-set-replicas", type=int, default=1,
                   help="Minimum number of pooled control replicas required for unflagged set-level comparisons (default: 1).")
    p.add_argument("--allow-low-n-comparisons", action="store_true",
                   help="Allow POI-vs-control comparisons even when n_used < min-replicas.")
    p.add_argument("--alpha", type=float, default=0.05, help="Alpha for reporting thresholds (default: 0.05).")
    p.add_argument("--fdr-method", default="bh", choices=["bh"], help="FDR correction method (default: bh).")
    p.add_argument("--permutation-tests", action="store_true",
                   help="Enable optional permutation p-value robustness check.")
    p.add_argument("--n-perm", type=int, default=10000,
                   help="Number of permutations when --permutation-tests is enabled (default: 10000).")
    p.add_argument("--seed", type=int, default=None,
                   help="Optional RNG seed for permutation tests.")
    p.add_argument("--run-checklist", dest="run_checklist", action="store_true", default=True,
                   help="Run fast Phase 2 checklist validation at end of run (default: enabled).")
    p.add_argument("--no-checklist", dest="run_checklist", action="store_false",
                   help="Disable Phase 2 checklist validation.")
    p.add_argument("--checklist-max-detail", type=int, default=10,
                   help="Maximum failing/warning examples per checklist subsection (default: 10).")
    return p


def main(argv: Optional[List[str]] = None) -> int:
    args = build_argparser().parse_args(argv)

    exp_root = Path(args.exp_root).expanduser().resolve()
    outdir = Path(args.outdir).expanduser().resolve()

    if not exp_root.exists():
        eprint(f"ERROR: exp-root does not exist: {exp_root}")
        return 2

    sanity_cfg = SanityConfig(
        enable=not args.sanity_disable,
        eint_min=float(args.sanity_eint_min),
        eint_max=float(args.sanity_eint_max),
    )

    # Ensure outdir is inside exp_root if user passed "analysis" relative
    safe_mkdir(outdir)

    # Write run_info.json
    run_info = {
        "timestamp": now_iso(),
        "argv": sys.argv,
        "exp_root": str(exp_root),
        "outdir": str(outdir),
        "topk": int(args.topk),
        "sanity": asdict(sanity_cfg),
        "phase2": {
            "timestamp": now_iso(),
            "enabled": bool(args.run_phase2),
            "params": {
                "poi_id": str(args.poi_id),
                "exclude_peptide_suffix": str(args.exclude_peptide_suffix),
                "min_replicas": int(args.min_replicas),
                "min_set_peptides": int(args.min_set_peptides),
                "min_set_replicas": int(args.min_set_replicas),
                "allow_low_n_comparisons": bool(args.allow_low_n_comparisons),
                "alpha": float(args.alpha),
                "fdr_method": str(args.fdr_method),
                "permutation_tests": bool(args.permutation_tests),
                "n_perm": int(args.n_perm),
                "seed": args.seed,
            },
        },
        "notes": {
            "phases_implemented": [
                "parse_summary_dlg",
                "consolidate_csvs(replicas_parsed, topk_parsed, clusters_parsed)",
                "export_pose_winners",
                "export_topk_concat_pdb",
                "basic_QA_and_parse_report",
                "phase2_group_summary_discrimination_global_rank_metrics_report",
            ],
        },
    }
    with (outdir / "run_info.json").open("w", encoding="utf-8") as f:
        json.dump(run_info, f, indent=2, sort_keys=True)

    run_parse_and_consolidate(
        exp_root,
        outdir,
        topk_k=int(args.topk),
        sanity_cfg=sanity_cfg,
        exclude_peptide_suffix=str(args.exclude_peptide_suffix),
        verbose=int(args.verbose),
    )
    if args.verbose:
        eprint(f"[INFO] Completed parse+consolidation for {exp_root} -> {outdir}")

    run_phase2_from_replicas_csv(
        replicas_csv=outdir / "replicas_parsed.csv",
        outdir=outdir,
        parse_report_json=outdir / "parse_report.json",
        run_info_json=outdir / "run_info.json",
        sanity_cfg=sanity_cfg,
        run_phase2=bool(args.run_phase2),
        poi_id=str(args.poi_id),
        exclude_peptide_suffix=str(args.exclude_peptide_suffix),
        min_replicas=int(args.min_replicas),
        allow_low_n_comparisons=bool(args.allow_low_n_comparisons),
        alpha=float(args.alpha),
        fdr_method=str(args.fdr_method),
        permutation_tests=bool(args.permutation_tests),
        n_perm=int(args.n_perm),
        seed=args.seed,
        exp_root=exp_root,
        min_set_peptides=int(args.min_set_peptides),
        min_set_replicas=int(args.min_set_replicas),
        verbose=int(args.verbose),
    )

    if args.run_checklist:
        _phase2_checklist_report(
            outdir=outdir,
            exp_root=exp_root,
            parse_report_json=outdir / "parse_report.json",
            parse_report_txt=outdir / "parse_report.txt",
            metrics_report_md=outdir / "metrics_report.md",
            poi_id=str(args.poi_id),
            exclude_peptide_suffix=str(args.exclude_peptide_suffix),
            sanity_cfg=sanity_cfg,
            topk_k=int(args.topk),
            min_replicas=int(args.min_replicas),
            min_set_peptides=int(args.min_set_peptides),
            min_set_replicas=int(args.min_set_replicas),
            alpha=float(args.alpha),
            checklist_max_detail=max(1, int(args.checklist_max_detail)),
            verbose=int(args.verbose),
        )
    else:
        vlog(1, int(args.verbose), "[INFO] Phase 2 checklist disabled by CLI (--no-checklist).")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
