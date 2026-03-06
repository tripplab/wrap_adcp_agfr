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


# ----------------------------
# Pipeline principal: parseo + consolidación + export PDB + QA
# ----------------------------

def run_parse_and_consolidate(exp_root: Path, outdir: Path, topk_k: int, sanity_cfg: SanityConfig, verbose: int = 0) -> None:
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

    # ----------------------------
    # PHASE 2 (placeholder): Esquema de tablas -> métricas -> tests -> reporte
    # ----------------------------
    #
    # En esta implementación entregamos:
    #   - Parsing (replicas_parsed.csv, topk_parsed.csv, clusters_parsed.csv)
    #   - Exports de poses (winners + concatenado top-k)
    #   - QA básico y parse_report
    #
    # A continuación dejamos EN EL CÓDIGO los puntos exactos donde se añadirán
    # las fases siguientes:
    #
    #   (A) group_summary.csv
    #       Propósito: resumir por (protein, peptide) las métricas robustas:
    #         - n_ok/fail_rate
    #         - Eint_median/IQR (winner_omm_dE_interaction)
    #         - convergencia: f1_cluster y Neff_cluster usando winner_cluster_mode_id
    #
    #   (B) discrimination.csv
    #       Propósito: comparar dentro de cada proteína POI vs controles:
    #         - wins/AUC (Mann–Whitney U normalizado)
    #         - Cliff's delta
    #         - p-values (Mann–Whitney + permutación opcional)
    #         - ajuste por múltiples comparaciones (FDR)
    #         - diagnóstico automático (verdict)
    #
    #   (C) Reporte/QA ampliado
    #       Propósito: un reporte legible (txt/json) con:
    #         - grupos vacíos
    #         - réplicas fallidas y causas
    #         - sanity outliers
    #         - resumen por proteína
    #
    # Notas:
    #   - La información necesaria YA está en replicas_parsed.csv
    #     (score winner + cluster_mode_id + QA).
    #   - Solo falta codificar los agregados y tests.
    #
    # Implementaremos esas fases en funciones separadas para mantener limpio:
    #   compute_group_summary(parsed_replicas, outdir)
    #   compute_discrimination(group_summary, parsed_replicas, outdir)
    #   write_human_report(report, outdir)
    #

    # ---- PHASE 2A: group_summary.csv (no implementado aún) ----
    # group_summary_rows = compute_group_summary(parsed_replicas)
    # write_csv(outdir / "group_summary.csv", group_summary_rows, group_summary_fields)

    # ---- PHASE 2B: discrimination.csv (no implementado aún) ----
    # discrimination_rows = compute_discrimination(parsed_replicas, group_summary_rows)
    # write_csv(outdir / "discrimination.csv", discrimination_rows, discrimination_fields)

    # ---- PHASE 2C: reporte/QA ampliado (no implementado aún) ----
    # write_text(outdir / "analysis_report.txt", render_human_report(...))

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

    # Future-phase switches (placeholders)
    p.add_argument("--compute-metrics", action="store_true",
                   help="(Future) Compute group_summary/discrimination metrics after parsing (not implemented yet).")
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
        "notes": {
            "phases_implemented": [
                "parse_summary_dlg",
                "consolidate_csvs(replicas_parsed, topk_parsed, clusters_parsed)",
                "export_pose_winners",
                "export_topk_concat_pdb",
                "basic_QA_and_parse_report",
            ],
            "phases_pending": [
                "group_summary.csv (aggregations per protein-peptide)",
                "discrimination.csv (wins/AUC, MWU, permutation, FDR, verdict)",
                "expanded human report / QA diagnostics summary",
            ],
        },
    }
    with (outdir / "run_info.json").open("w", encoding="utf-8") as f:
        json.dump(run_info, f, indent=2, sort_keys=True)

    run_parse_and_consolidate(exp_root, outdir, topk_k=int(args.topk), sanity_cfg=sanity_cfg, verbose=int(args.verbose))
    if args.verbose:
        eprint(f"[INFO] Completed parse+consolidation for {exp_root} -> {outdir}")

    if args.compute_metrics:
        # Placeholder hook for future phases (intentionally not implemented).
        eprint("NOTE: --compute-metrics is not implemented yet. Parsing outputs are ready in outdir.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
