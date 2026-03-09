#!/usr/bin/env python3
"""
pepctrl.py — Peptide controls: matched decoys, scrambled, and random sequences.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import random
import re
import sys
import time
from collections import Counter
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Sequence, Set, Tuple

AA20 = set("ACDEFGHIKLMNPQRSTVWY")
AA_LIST = sorted(AA20)

KD_HYDRO = {
    "I": 4.5, "V": 4.2, "L": 3.8, "F": 2.8, "C": 2.5,
    "M": 1.9, "A": 1.8, "G": -0.4, "T": -0.7, "S": -0.8,
    "W": -0.9, "Y": -1.3, "P": -1.6, "H": -3.2, "E": -3.5,
    "Q": -3.5, "D": -3.5, "N": -3.5, "K": -3.9, "R": -4.5,
}
AROM = set("FYW")
POS = set("KR")
NEG = set("DE")

CLASS_MAP = {
    "A": "hydrophobic", "V": "hydrophobic", "L": "hydrophobic", "I": "hydrophobic", "M": "hydrophobic",
    "F": "aromatic", "Y": "aromatic", "W": "aromatic",
    "K": "positive", "R": "positive", "H": "polar",
    "D": "negative", "E": "negative",
    "S": "polar", "T": "polar", "N": "polar", "Q": "polar",
    "C": "special", "P": "special", "G": "special",
}
CLASS_TO_AA: Dict[str, List[str]] = {}
for aa, cls in CLASS_MAP.items():
    CLASS_TO_AA.setdefault(cls, []).append(aa)

COARSE_GROUPS = {
    "hydrophobic": list("AVLIM"),
    "aromatic": list("FYW"),
    "positive": list("KRH"),
    "negative": list("DE"),
    "polar": list("STNQ"),
    "special": list("CPG"),
}
AA_TO_COARSE_GROUP = {aa: group for group, residues in COARSE_GROUPS.items() for aa in residues}
SOFT_KEYS = [
    "charge_at_ph",
    "pI",
    "hydro_mean",
    "hydro_var",
    "aromatic_burden",
    "cpg_burden",
    "complexity",
    "charge_patterning",
    "hydro_blockiness",
]

PKA = {"Nterm": 9.69, "Cterm": 2.34, "K": 10.5, "R": 12.5, "H": 6.0, "D": 3.9, "E": 4.1, "C": 8.3, "Y": 10.1}


@dataclass
class AnnealSettings:
    steps: int
    t_start: float
    t_end: float


@dataclass
class ObjectiveWeights:
    charge_at_ph: float
    pI: float
    hydro_mean: float
    hydro_var: float
    aromatic_burden: float
    cpg_burden: float
    complexity: float
    charge_patterning: float
    hydro_blockiness: float


@dataclass
class Descriptor:
    length: int
    net_charge: float
    charge_at_ph: float
    pI: float
    hydro_mean: float
    hydro_var: float
    aromatic_burden: float
    cpg_burden: float
    complexity: float
    charge_patterning: float
    hydro_blockiness: float

    def as_dict(self) -> Dict[str, float]:
        return {
            "length": self.length,
            "net_charge": self.net_charge,
            "charge_at_ph": self.charge_at_ph,
            "pI": self.pI,
            "hydro_mean": self.hydro_mean,
            "hydro_var": self.hydro_var,
            "aromatic_burden": self.aromatic_burden,
            "cpg_burden": self.cpg_burden,
            "complexity": self.complexity,
            "charge_patterning": self.charge_patterning,
            "hydro_blockiness": self.hydro_blockiness,
        }


@dataclass
class DecoyResult:
    sequence: str
    decoy_class: str
    seed: int
    objective_score: float
    descriptor: Descriptor
    descriptor_delta: Dict[str, float]
    identity_to_poi: float
    matches_to_poi: int
    relaxation_stage_used: int
    flags: List[str]


@dataclass
class Config:
    poi: str
    seq_type: str
    decoy_class: str
    n: int
    out: str
    fmt: str
    no_header: bool
    include_poi_row: bool
    extra_cols: bool
    allow_partial: bool
    seed: Optional[int]
    max_identity: float
    max_matches: int
    max_identity_between_decoys: float
    min_complexity: float
    max_run: int
    no_simple_repeats: bool
    exclude_literals: List[str]
    exclude_regexes: List[re.Pattern]
    ph: float
    include_termini: bool
    anneal: AnnealSettings
    oversample_factor: int
    progress_every: int
    relaxation: bool
    max_relax_stage: int
    property_tolerances: Dict[str, float]
    weights: ObjectiveWeights
    run_report_json: Optional[str]


@dataclass
class RunStats:
    generated_attempts: int = 0
    valid_candidates: int = 0
    selected_final: int = 0
    rejections: Dict[str, int] = field(default_factory=dict)
    relax_log: List[Dict[str, Any]] = field(default_factory=list)
    timings: Dict[str, float] = field(default_factory=dict)
    anneal_runs_started: int = 0
    anneal_runs_completed: int = 0
    total_proposals_generated: int = 0
    total_sequences_described: int = 0
    total_hard_valid_proposals: int = 0
    total_hard_failures: int = 0
    total_soft_failures: int = 0
    total_soft_valid_candidates: int = 0
    total_best_of_run_candidates: int = 0
    total_pool_admissions: int = 0
    total_duplicate_rejections: int = 0
    target_pool_expected: int = 0
    pool_valid_achieved: int = 0
    written_files: List[Dict[str, str]] = field(default_factory=list)


def safe_rate(num: float, den: float) -> float:
    return num / den if den else 0.0


def build_breakdown(counts: Dict[str, int], denom: int) -> Dict[str, Dict[str, float]]:
    return {
        key: {
            "count": val,
            "pct_of_failures": safe_rate(val, sum(counts.values())),
            "pct_of_total_proposals": safe_rate(val, denom),
        }
        for key, val in sorted(counts.items())
    }


def build_sampling_summary(stats: RunStats) -> Dict[str, float]:
    return {
        "anneal_runs_started": stats.anneal_runs_started,
        "anneal_runs_completed": stats.anneal_runs_completed,
        "total_proposals_generated": stats.total_proposals_generated,
        "total_sequences_described": stats.total_sequences_described,
        "total_hard_valid_proposals": stats.total_hard_valid_proposals,
        "total_soft_valid_candidates": stats.total_soft_valid_candidates,
        "total_best_of_run_candidates": stats.total_best_of_run_candidates,
        "total_pool_admissions": stats.total_pool_admissions,
        "total_final_selected": stats.selected_final,
        "target_pool": stats.target_pool_expected,
        "achieved_pool": stats.pool_valid_achieved,
        "acceptance_rate_hard_valid": safe_rate(stats.total_hard_valid_proposals, stats.total_proposals_generated),
        "acceptance_rate_soft_valid": safe_rate(stats.total_soft_valid_candidates, stats.total_hard_valid_proposals),
        "acceptance_rate_pool": safe_rate(stats.total_pool_admissions, stats.total_best_of_run_candidates),
        "final_selection_rate": safe_rate(stats.selected_final, stats.total_pool_admissions),
    }


def parse_args(argv: Optional[List[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("--poi", type=str)
    p.add_argument("--poi-file", type=str)
    p.add_argument("--type", dest="seq_type", choices=["decoy", "scramble", "random"], required=True)
    p.add_argument("--decoy-class", choices=["exact-composition", "class-matched", "property-matched"], default="exact-composition")
    p.add_argument("--n", type=int, required=True)
    p.add_argument("--out", type=str, default="-")
    p.add_argument("--format", choices=["csv", "tsv", "jsonl", "fasta"], default="csv")
    p.add_argument("--no-header", action="store_true")
    p.add_argument("--extra-cols", action="store_true")
    p.add_argument("--no-poi-row", action="store_true")
    p.add_argument("--allow-partial", action="store_true")
    p.add_argument("--seed", type=int)

    p.add_argument("--ph", type=float, default=7.0)
    p.add_argument("--include-termini", action="store_true")

    p.add_argument("--max-identity", type=float, default=0.33)
    p.add_argument("--max-matches", type=int)
    p.add_argument("--max-identity-between-decoys", type=float, default=0.5)
    p.add_argument("--min-complexity", type=float, default=0.55)
    p.add_argument("--max-run", type=int, default=3)
    p.add_argument("--no-simple-repeats", action="store_true", default=True)
    p.add_argument("--allow-simple-repeats", action="store_true")

    p.add_argument("--exclude-motif", action="append", default=[])
    p.add_argument("--exclude-motif-file", type=str)

    p.add_argument("--anneal-steps", type=int, default=3000)
    p.add_argument("--anneal-temp-start", type=float, default=2.0)
    p.add_argument("--anneal-temp-end", type=float, default=0.02)
    p.add_argument("--oversample-factor", type=int, default=6)
    p.add_argument("--progress-every", type=int, default=1)

    p.add_argument("--relax-soft", action="store_true")
    p.add_argument("--max-relax-stage", type=int, default=3)

    p.add_argument("--tol-charge-at-ph", type=float, default=1.0)
    p.add_argument("--tol-pi", type=float, default=1.0)
    p.add_argument("--tol-hydro-mean", type=float, default=0.8)
    p.add_argument("--tol-hydro-var", type=float, default=1.0)
    p.add_argument("--tol-aromatic-burden", type=float, default=0.2)
    p.add_argument("--tol-cpg-burden", type=float, default=0.2)
    p.add_argument("--tol-complexity", type=float, default=0.3)
    p.add_argument("--tol-charge-patterning", type=float, default=0.25)
    p.add_argument("--tol-hydro-blockiness", type=float, default=0.25)

    p.add_argument("--w-charge-at-ph", type=float, default=1.0)
    p.add_argument("--w-pi", type=float, default=1.0)
    p.add_argument("--w-hydro-mean", type=float, default=1.0)
    p.add_argument("--w-hydro-var", type=float, default=0.8)
    p.add_argument("--w-aromatic-burden", type=float, default=0.8)
    p.add_argument("--w-cpg-burden", type=float, default=0.6)
    p.add_argument("--w-complexity", type=float, default=0.8)
    p.add_argument("--w-charge-patterning", type=float, default=1.0)
    p.add_argument("--w-hydro-blockiness", type=float, default=0.8)

    p.add_argument("--run-report-json", type=str, default=None)

    # deprecated legacy options: accepted for compatibility, no-op
    p.add_argument("--method", type=str, default=None)
    p.add_argument("--refine-steps", type=int, default=None)
    p.add_argument("--relax", action="store_true")

    return p.parse_args(argv)


def validate_sequence(seq: str) -> None:
    bad = sorted(set(seq) - AA20)
    if bad:
        raise ValueError(f"Invalid residue(s): {bad}")


def load_poi(ns: argparse.Namespace) -> str:
    seq = ns.poi
    if not seq and ns.poi_file:
        with open(ns.poi_file, "r", encoding="utf-8") as f:
            for line in f:
                if line.strip():
                    seq = line.strip()
                    break
    if not seq:
        raise ValueError("Use --poi or --poi-file")
    seq = seq.upper().strip()
    validate_sequence(seq)
    return seq


def compile_excludes(values: List[str], file_path: Optional[str]) -> Tuple[List[str], List[re.Pattern]]:
    motifs = list(values)
    if file_path:
        with open(file_path, "r", encoding="utf-8") as f:
            motifs.extend([x.strip() for x in f if x.strip() and not x.startswith("#")])
    literals, regexes = [], []
    for m in motifs:
        if m.startswith("re:"):
            regexes.append(re.compile(m[3:]))
        else:
            literals.append(m.upper())
    return literals, regexes


def resolve_config(ns: argparse.Namespace, poi: str) -> Config:
    if ns.method is not None or ns.refine_steps is not None or ns.relax:
        print("[pepctrl] DEPRECATED: --method/--refine-steps/--relax are deprecated under annealing decoy framework.")
    max_matches = ns.max_matches if ns.max_matches is not None else int(math.floor(ns.max_identity * len(poi)))
    literals, regexes = compile_excludes(ns.exclude_motif or [], ns.exclude_motif_file)
    return Config(
        poi=poi,
        seq_type=ns.seq_type,
        decoy_class=ns.decoy_class,
        n=ns.n,
        out=ns.out,
        fmt=ns.format,
        no_header=ns.no_header,
        include_poi_row=(not ns.no_poi_row),
        extra_cols=ns.extra_cols,
        allow_partial=ns.allow_partial,
        seed=ns.seed,
        max_identity=ns.max_identity,
        max_matches=max_matches,
        max_identity_between_decoys=ns.max_identity_between_decoys,
        min_complexity=ns.min_complexity,
        max_run=ns.max_run,
        no_simple_repeats=(ns.no_simple_repeats and not ns.allow_simple_repeats),
        exclude_literals=literals,
        exclude_regexes=regexes,
        ph=ns.ph,
        include_termini=ns.include_termini,
        anneal=AnnealSettings(ns.anneal_steps, ns.anneal_temp_start, ns.anneal_temp_end),
        oversample_factor=max(1, ns.oversample_factor),
        progress_every=max(1, ns.progress_every),
        relaxation=ns.relax_soft,
        max_relax_stage=max(0, ns.max_relax_stage),
        property_tolerances={
            "charge_at_ph": ns.tol_charge_at_ph,
            "pI": ns.tol_pi,
            "hydro_mean": ns.tol_hydro_mean,
            "hydro_var": ns.tol_hydro_var,
            "aromatic_burden": ns.tol_aromatic_burden,
            "cpg_burden": ns.tol_cpg_burden,
            "complexity": ns.tol_complexity,
            "charge_patterning": ns.tol_charge_patterning,
            "hydro_blockiness": ns.tol_hydro_blockiness,
        },
        weights=ObjectiveWeights(
            ns.w_charge_at_ph, ns.w_pi, ns.w_hydro_mean, ns.w_hydro_var,
            ns.w_aromatic_burden, ns.w_cpg_burden, ns.w_complexity,
            ns.w_charge_patterning, ns.w_hydro_blockiness,
        ),
        run_report_json=ns.run_report_json,
    )


def approx_charge_at_ph(seq: str, ph: float, include_termini: bool) -> float:
    counts = Counter(seq)
    pos = counts["K"] / (1 + 10 ** (ph - PKA["K"]))
    pos += counts["R"] / (1 + 10 ** (ph - PKA["R"]))
    pos += counts["H"] / (1 + 10 ** (ph - PKA["H"]))
    neg = counts["D"] / (1 + 10 ** (PKA["D"] - ph))
    neg += counts["E"] / (1 + 10 ** (PKA["E"] - ph))
    neg += counts["C"] / (1 + 10 ** (PKA["C"] - ph))
    neg += counts["Y"] / (1 + 10 ** (PKA["Y"] - ph))
    if include_termini and seq:
        pos += 1 / (1 + 10 ** (ph - PKA["Nterm"]))
        neg += 1 / (1 + 10 ** (PKA["Cterm"] - ph))
    return pos - neg


def approx_pI(seq: str, include_termini: bool) -> float:
    lo, hi = 0.0, 14.0
    for _ in range(40):
        mid = (lo + hi) / 2
        q = approx_charge_at_ph(seq, mid, include_termini)
        if q > 0:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


def shannon_complexity(seq: str) -> float:
    L = len(seq)
    if L == 0:
        return 0.0
    cnt = Counter(seq)
    ent = 0.0
    for c in cnt.values():
        p = c / L
        ent -= p * math.log2(p)
    return ent / math.log2(min(20, L))


def charge_patterning(seq: str, ph: float, include_termini: bool) -> float:
    q = [approx_charge_at_ph(aa, ph, False) for aa in seq]
    if include_termini and q:
        q[0] += 1.0
        q[-1] -= 1.0
    if len(q) < 2:
        return 0.0
    return sum(abs(q[i] - q[i - 1]) for i in range(1, len(q))) / (len(q) - 1)


def hydro_blockiness(seq: str) -> float:
    vals = [KD_HYDRO[a] for a in seq]
    if not vals:
        return 0.0
    m = sum(vals) / len(vals)
    mask = [1 if v >= m else 0 for v in vals]
    transitions = sum(1 for i in range(1, len(mask)) if mask[i] != mask[i - 1])
    return 1.0 - transitions / max(1, len(mask) - 1)


def describe(seq: str, cfg: Config) -> Descriptor:
    vals = [KD_HYDRO[a] for a in seq]
    m = sum(vals) / len(vals)
    var = sum((x - m) ** 2 for x in vals) / len(vals)
    ch = approx_charge_at_ph(seq, cfg.ph, cfg.include_termini)
    return Descriptor(
        length=len(seq),
        net_charge=round(ch, 6),
        charge_at_ph=ch,
        pI=approx_pI(seq, cfg.include_termini),
        hydro_mean=m,
        hydro_var=var,
        aromatic_burden=sum(1 for a in seq if a in AROM) / len(seq),
        cpg_burden=sum(1 for a in seq if a in {"C", "P", "G"}) / len(seq),
        complexity=shannon_complexity(seq),
        charge_patterning=charge_patterning(seq, cfg.ph, cfg.include_termini),
        hydro_blockiness=hydro_blockiness(seq),
    )


def identity_stats(a: str, b: str) -> Tuple[int, float]:
    matches = sum(1 for x, y in zip(a, b) if x == y)
    return matches, matches / len(a)


def has_simple_repeats(seq: str) -> bool:
    for k in (2, 3):
        for i in range(len(seq) - 2 * k + 1):
            if seq[i:i + k] == seq[i + k:i + 2 * k]:
                return True
    return False


def max_run_length(seq: str) -> int:
    best = cur = 1
    for i in range(1, len(seq)):
        cur = cur + 1 if seq[i] == seq[i - 1] else 1
        best = max(best, cur)
    return best


def hard_check(seq: str, cfg: Config, poi: str, seen: Set[str], mode_meta: Dict[str, Any]) -> Tuple[bool, str]:
    if set(seq) - AA20:
        return False, "invalid_aa"
    if len(seq) != len(poi):
        return False, "length"
    if seq == poi:
        return False, "poi_exclusion"
    if seq in seen:
        return False, "duplicate"
    m, ident = identity_stats(seq, poi)
    if m > cfg.max_matches or ident > cfg.max_identity:
        return False, "identity_to_poi"
    if any(lit in seq for lit in cfg.exclude_literals):
        return False, "forbidden_motif"
    if any(rx.search(seq) for rx in cfg.exclude_regexes):
        return False, "forbidden_motif"
    if max_run_length(seq) > cfg.max_run:
        return False, "low_complexity"
    if shannon_complexity(seq) < cfg.min_complexity:
        return False, "low_complexity"
    if cfg.no_simple_repeats and has_simple_repeats(seq):
        return False, "low_complexity"
    if cfg.decoy_class == "exact-composition":
        if Counter(seq) != mode_meta["poi_counts"]:
            return False, "exact_composition"
    if cfg.decoy_class == "class-matched":
        cls = Counter(CLASS_MAP[a] for a in seq)
        if cls != mode_meta["poi_class_counts"]:
            return False, "exact_class_counts"
    return True, "ok"


def objective(desc: Descriptor, poi_desc: Descriptor, w: ObjectiveWeights) -> float:
    d = {
        "charge_at_ph": abs(desc.charge_at_ph - poi_desc.charge_at_ph),
        "pI": abs(desc.pI - poi_desc.pI),
        "hydro_mean": abs(desc.hydro_mean - poi_desc.hydro_mean),
        "hydro_var": abs(desc.hydro_var - poi_desc.hydro_var),
        "aromatic_burden": abs(desc.aromatic_burden - poi_desc.aromatic_burden),
        "cpg_burden": abs(desc.cpg_burden - poi_desc.cpg_burden),
        "complexity": abs(desc.complexity - poi_desc.complexity),
        "charge_patterning": abs(desc.charge_patterning - poi_desc.charge_patterning),
        "hydro_blockiness": abs(desc.hydro_blockiness - poi_desc.hydro_blockiness),
    }
    return (
        w.charge_at_ph * d["charge_at_ph"] + w.pI * d["pI"] + w.hydro_mean * d["hydro_mean"] +
        w.hydro_var * d["hydro_var"] + w.aromatic_burden * d["aromatic_burden"] + w.cpg_burden * d["cpg_burden"] +
        w.complexity * d["complexity"] + w.charge_patterning * d["charge_patterning"] + w.hydro_blockiness * d["hydro_blockiness"]
    )


def descriptor_deltas(desc: Descriptor, poi: Descriptor) -> Dict[str, float]:
    out = {}
    for k, v in desc.as_dict().items():
        out[k] = v - poi.as_dict()[k]
    return out


def current_soft_scale(stage: int) -> float:
    return 1.0 + 0.25 * stage


def soft_tolerance_failures(desc: Descriptor, poi: Descriptor, cfg: Config, stage: int) -> List[str]:
    if cfg.decoy_class != "property-matched":
        return []
    scale = current_soft_scale(stage)
    desc_map = desc.as_dict()
    poi_map = poi.as_dict()
    failed = []
    for key in SOFT_KEYS:
        tol = cfg.property_tolerances[key]
        delta = abs(desc_map[key] - poi_map[key])
        if delta > tol * scale + 1e-12:
            failed.append(key)
    return failed


def init_candidate_property_matched(poi: str, rng: random.Random) -> str:
    L = len(poi)
    poi_counts = Counter(poi)
    target_arom = poi_counts["F"] + poi_counts["Y"] + poi_counts["W"]
    target_cpg = poi_counts["C"] + poi_counts["P"] + poi_counts["G"]
    target_pos = poi_counts["K"] + poi_counts["R"]
    target_neg = poi_counts["D"] + poi_counts["E"]
    pos_bias = max(0, target_pos - target_neg)
    neg_bias = max(0, target_neg - target_pos)
    neutral_target = max(0, L - pos_bias - neg_bias)

    pool: List[str] = []
    pool.extend(rng.choices(["F", "Y", "W"], k=target_arom))
    pool.extend(rng.choices(["C", "P", "G"], k=target_cpg))
    pool.extend(rng.choices(["K", "R", "H"], k=pos_bias))
    pool.extend(rng.choices(["D", "E"], k=neg_bias))

    neutral_pool = ["A", "V", "L", "I", "M", "T", "S", "N", "Q", "G", "P"]
    pool.extend(rng.choices(neutral_pool, k=neutral_target))

    if len(pool) > L:
        rng.shuffle(pool)
        pool = pool[:L]
    elif len(pool) < L:
        # Prefer filling with moderate hydrophobicity residues to stay near typical POI manifold.
        pool.extend(rng.choices(["A", "V", "L", "T", "S", "Q"], k=L - len(pool)))

    rng.shuffle(pool)
    return "".join(pool)


def init_candidate(poi: str, cfg: Config, rng: random.Random) -> str:
    if cfg.decoy_class == "exact-composition":
        s = list(poi)
        rng.shuffle(s)
        return "".join(s)
    if cfg.decoy_class == "class-matched":
        labels = [CLASS_MAP[a] for a in poi]
        rng.shuffle(labels)
        return "".join(rng.choice(CLASS_TO_AA[c]) for c in labels)
    return init_candidate_property_matched(poi, rng)


def move(seq: str, cfg: Config, rng: random.Random) -> str:
    s = list(seq)
    L = len(s)
    if cfg.decoy_class == "exact-composition":
        i, j = rng.randrange(L), rng.randrange(L)
        s[i], s[j] = s[j], s[i]
    elif cfg.decoy_class == "class-matched":
        if rng.random() < 0.5:
            i, j = rng.randrange(L), rng.randrange(L)
            s[i], s[j] = s[j], s[i]
        else:
            i = rng.randrange(L)
            cls = CLASS_MAP[s[i]]
            s[i] = rng.choice(CLASS_TO_AA[cls])
    else:
        r = rng.random()
        if r < 0.45:
            i = rng.randrange(L)
            s[i] = rng.choice(AA_LIST)
        elif r < 0.7:
            i, j = rng.randrange(L), rng.randrange(L)
            s[i], s[j] = s[j], s[i]
        else:
            i = rng.randrange(L)
            group = AA_TO_COARSE_GROUP[s[i]]
            s[i] = rng.choice(COARSE_GROUPS[group])
    return "".join(s)


def anneal_one(poi: str, poi_desc: Descriptor, cfg: Config, seed_i: int,
               mode_meta: Dict[str, Any], stage: int, run: RunStats, seen: Set[str]) -> Optional[DecoyResult]:
    run.anneal_runs_started += 1
    local_rng = random.Random(seed_i)
    start = init_candidate(poi, cfg, local_rng)
    run.total_sequences_described += 1
    current = start
    current_desc = describe(current, cfg)
    current_score = objective(current_desc, poi_desc, cfg.weights)
    best = (current, current_desc, current_score)
    best_score = current_score

    for step in range(cfg.anneal.steps):
        t = cfg.anneal.t_start * ((cfg.anneal.t_end / cfg.anneal.t_start) ** (step / max(1, cfg.anneal.steps - 1)))
        prop = move(current, cfg, local_rng)
        run.generated_attempts += 1
        run.total_proposals_generated += 1

        ok, reason = hard_check(prop, cfg, poi, set(), mode_meta)
        if not ok:
            run.total_hard_failures += 1
            run.rejections[reason] = run.rejections.get(reason, 0) + 1
            continue

        run.total_hard_valid_proposals += 1
        run.total_sequences_described += 1
        d = describe(prop, cfg)
        sc = objective(d, poi_desc, cfg.weights)

        accept = sc < current_score or local_rng.random() < math.exp(-(sc - current_score) / max(t, 1e-9))
        if accept:
            current, current_desc, current_score = prop, d, sc
        if sc < best_score:
            best = (prop, d, sc)
            best_score = sc

    seq, d, sc = best
    run.total_best_of_run_candidates += 1
    ok, reason = hard_check(seq, cfg, poi, seen, mode_meta)
    if not ok:
        run.total_hard_failures += 1
        if reason == "duplicate":
            run.total_duplicate_rejections += 1
        run.rejections[reason] = run.rejections.get(reason, 0) + 1
        run.anneal_runs_completed += 1
        return None

    soft_failed = soft_tolerance_failures(d, poi_desc, cfg, stage)
    if soft_failed:
        run.total_soft_failures += 1
        run.rejections["soft_tolerance"] = run.rejections.get("soft_tolerance", 0) + 1
        for key in soft_failed:
            bucket = f"soft_tolerance_{key.lower()}"
            run.rejections[bucket] = run.rejections.get(bucket, 0) + 1
        run.anneal_runs_completed += 1
        return None

    run.total_soft_valid_candidates += 1
    matches, ident = identity_stats(seq, poi)
    run.valid_candidates += 1
    run.anneal_runs_completed += 1
    return DecoyResult(
        sequence=seq,
        decoy_class=cfg.decoy_class,
        seed=seed_i,
        objective_score=sc,
        descriptor=d,
        descriptor_delta=descriptor_deltas(d, poi_desc),
        identity_to_poi=ident,
        matches_to_poi=matches,
        relaxation_stage_used=stage,
        flags=[],
    )


def select_diverse(pool: List[DecoyResult], n: int, max_ident: float) -> Tuple[List[DecoyResult], Dict[str, Any]]:
    chosen: List[DecoyResult] = []
    diversity_rejections = 0
    for cand in sorted(pool, key=lambda x: x.objective_score):
        if any(identity_stats(cand.sequence, x.sequence)[1] > max_ident + 1e-12 for x in chosen):
            diversity_rejections += 1
            continue
        chosen.append(cand)
        if len(chosen) == n:
            break
    min_pair = 1.0
    for i in range(len(chosen)):
        for j in range(i + 1, len(chosen)):
            min_pair = min(min_pair, identity_stats(chosen[i].sequence, chosen[j].sequence)[1])
    return chosen, {
        "strict_max_pairwise_identity": max_ident,
        "observed_min_pairwise_identity": min_pair if len(chosen) > 1 else None,
        "rejected_by_diversity_filter": diversity_rejections,
    }


def write_decoys(decoys: List[DecoyResult], poi: str, poi_desc: Descriptor, cfg: Config, stats: RunStats) -> None:
    cols = ["decoy_id", "decoy_class", "sequence", "seed", "objective_score", "descriptors", "descriptor_deltas_vs_poi",
            "identity_to_poi", "matches_to_poi", "relaxation_stage_used", "flags"]
    h = sys.stdout if cfg.out == "-" else open(cfg.out, "w", encoding="utf-8", newline="")
    try:
        if cfg.fmt in ("csv", "tsv"):
            delim = "," if cfg.fmt == "csv" else "\t"
            writer = csv.DictWriter(h, fieldnames=cols, delimiter=delim)
            if not cfg.no_header:
                writer.writeheader()
            if cfg.include_poi_row:
                writer.writerow({"decoy_id": "POI", "decoy_class": "poi", "sequence": poi, "seed": "", "objective_score": "",
                                 "descriptors": json.dumps(poi_desc.as_dict(), sort_keys=True), "descriptor_deltas_vs_poi": json.dumps({}, sort_keys=True),
                                 "identity_to_poi": 1.0, "matches_to_poi": len(poi), "relaxation_stage_used": "", "flags": ""})
            for i, d in enumerate(decoys, 1):
                writer.writerow({"decoy_id": f"D{i:04d}", "decoy_class": d.decoy_class, "sequence": d.sequence, "seed": d.seed,
                                 "objective_score": f"{d.objective_score:.6g}", "descriptors": json.dumps(d.descriptor.as_dict(), sort_keys=True),
                                 "descriptor_deltas_vs_poi": json.dumps(d.descriptor_delta, sort_keys=True), "identity_to_poi": f"{d.identity_to_poi:.6g}",
                                 "matches_to_poi": d.matches_to_poi, "relaxation_stage_used": d.relaxation_stage_used, "flags": ";".join(d.flags)})
        else:
            for i, d in enumerate(decoys, 1):
                h.write(json.dumps({"decoy_id": f"D{i:04d}", **d.__dict__}, default=lambda x: x.as_dict() if hasattr(x, "as_dict") else str(x)) + "\n")
    finally:
        if h is not sys.stdout:
            h.close()
    stats.written_files.append({
        "path": cfg.out,
        "kind": "accepted_sequences_table",
        "format": cfg.fmt,
        "description": "tabla de secuencias aceptadas con descriptores y deltas vs POI",
    })




def accepted_sequences_summary(decoys: List[DecoyResult]) -> List[Dict[str, Any]]:
    out = []
    for i, d in enumerate(decoys, 1):
        out.append({
            "decoy_id": f"D{i:04d}",
            "sequence": d.sequence,
            "objective_score": d.objective_score,
            "identity_to_poi": d.identity_to_poi,
            "matches_to_poi": d.matches_to_poi,
            "charge_at_ph": d.descriptor.charge_at_ph,
            "pI": d.descriptor.pI,
            "hydro_mean": d.descriptor.hydro_mean,
            "hydro_var": d.descriptor.hydro_var,
            "aromatic_burden": d.descriptor.aromatic_burden,
            "cpg_burden": d.descriptor.cpg_burden,
            "complexity": d.descriptor.complexity,
            "charge_patterning": d.descriptor.charge_patterning,
            "hydro_blockiness": d.descriptor.hydro_blockiness,
            "relaxation_stage_used": d.relaxation_stage_used,
        })
    return out


def print_stdout_summary(cfg: Config, stats: RunStats, poi_desc: Descriptor, decoys: List[DecoyResult], partial: bool, failed: bool, seed: int) -> None:
    sampling = build_sampling_summary(stats)
    hard_counts = {k: v for k, v in stats.rejections.items() if not k.startswith("soft_tolerance")}
    soft_counts = {k: v for k, v in stats.rejections.items() if k.startswith("soft_tolerance_")}
    hard_breakdown = build_breakdown(hard_counts, stats.total_proposals_generated)
    soft_breakdown = build_breakdown(soft_counts, stats.total_proposals_generated)

    print("[pepctrl] sampling_summary:")
    for k, v in sampling.items():
        print(f"[pepctrl]   {k}: {v}")

    print(f"[pepctrl] rejection_totals: total_hard_failures={stats.total_hard_failures} total_soft_failures={stats.total_soft_failures}")
    print("[pepctrl] hard_failure_breakdown:")
    for k, v in hard_breakdown.items():
        print(f"[pepctrl]   {k}: count={v['count']} pct_of_failures={v['pct_of_failures']:.4f} pct_of_total_proposals={v['pct_of_total_proposals']:.4f}")
    print("[pepctrl] soft_failure_breakdown (counts are failures per descriptor; soft_tolerance counts failed sequences):")
    for k, v in soft_breakdown.items():
        print(f"[pepctrl]   {k}: count={v['count']} pct_of_failures={v['pct_of_failures']:.4f} pct_of_total_proposals={v['pct_of_total_proposals']:.4f}")

    print("[pepctrl] written_files:")
    for wf in stats.written_files:
        print(f"[pepctrl]   path={wf['path']} kind={wf['kind']} format={wf['format']} description={wf['description']}")

    print(f"[pepctrl] accepted_sequences: {len(decoys)}")
    if decoys:
        for i, d in enumerate(decoys, 1):
            print(
                f"[pepctrl] D{i:04d} seq={d.sequence} obj={d.objective_score:.6g} id={d.identity_to_poi:.4f} "
                f"matches={d.matches_to_poi} charge_at_ph={d.descriptor.charge_at_ph:.4f} pI={d.descriptor.pI:.4f} "
                f"hydro_mean={d.descriptor.hydro_mean:.4f} hydro_var={d.descriptor.hydro_var:.4f} "
                f"aromatic_burden={d.descriptor.aromatic_burden:.4f} cpg_burden={d.descriptor.cpg_burden:.4f} "
                f"complexity={d.descriptor.complexity:.4f} charge_patterning={d.descriptor.charge_patterning:.4f} "
                f"hydro_blockiness={d.descriptor.hydro_blockiness:.4f} relaxation_stage_used={d.relaxation_stage_used}"
            )
    else:
        print("[pepctrl] No se aceptaron secuencias en la selección final.")

    print(
        f"[pepctrl] final_poi_reference sequence={cfg.poi} charge_at_ph={poi_desc.charge_at_ph:.4f} pI={poi_desc.pI:.4f} "
        f"hydro_mean={poi_desc.hydro_mean:.4f}"
    )
    print(
        f"[pepctrl] final_status requested={cfg.n} target_pool={stats.target_pool_expected} pool_valid={stats.pool_valid_achieved} "
        f"selected={stats.selected_final} partial={partial} failed={failed} total_s={stats.timings.get('total_s', 0.0):.4f} seed={seed}"
    )

def write_report(cfg: Config, poi_desc: Descriptor, stats: RunStats, pool_size: int, selected: int,
                 diversity: Dict[str, Any], partial: bool, failed: bool, seed: int, decoys: List[DecoyResult]) -> None:
    if not cfg.run_report_json:
        return
    soft_breakdown_legacy = {f"soft_tolerance_{k.lower()}": stats.rejections.get(f"soft_tolerance_{k.lower()}", 0) for k in SOFT_KEYS}
    hard_counts = {k: v for k, v in stats.rejections.items() if not k.startswith("soft_tolerance")}
    soft_counts = {k: v for k, v in stats.rejections.items() if k.startswith("soft_tolerance_")}
    sampling_summary = build_sampling_summary(stats)
    report = {
        "poi_summary": {"sequence": cfg.poi, "descriptor": poi_desc.as_dict()},
        "decoy_class": cfg.decoy_class,
        "config_cli": sys.argv[1:],
        "seed": seed,
        "hard_constraints": {
            "aa20_validity": True, "length": len(cfg.poi), "poi_exclusion": True, "duplicate_exclusion": True,
            "max_identity_to_poi": cfg.max_identity, "strict_inter_decoy_identity": cfg.max_identity_between_decoys,
            "forbidden_motifs": cfg.exclude_literals + [f"re:{r.pattern}" for r in cfg.exclude_regexes],
            "minimum_complexity": cfg.min_complexity,
            "exact_composition": cfg.decoy_class == "exact-composition",
            "exact_class_counts": cfg.decoy_class == "class-matched",
        },
        "soft_constraints": cfg.property_tolerances,
        "objective_weights": cfg.weights.__dict__,
        "annealing_settings": cfg.anneal.__dict__,
        "requested_generated_selected": {
            "requested": cfg.n,
            "target_pool": stats.target_pool_expected,
            "pool_valid": pool_size,
            "selected": selected,
        },
        "rejection_reasons": stats.rejections,
        "total_hard_failures": stats.total_hard_failures,
        "total_soft_failures": stats.total_soft_failures,
        "hard_failure_breakdown": build_breakdown(hard_counts, stats.total_proposals_generated),
        "soft_failure_breakdown": build_breakdown(soft_counts, stats.total_proposals_generated),
        "soft_tolerance_breakdown": soft_breakdown_legacy,
        "soft_tolerance_count_definition": {
            "soft_tolerance": "counts sequences rejected by soft tolerance (sequence-level)",
            "soft_tolerance_<descriptor>": "counts descriptor-level tolerance failures; one sequence can increment multiple descriptors",
        },
        "sampling_summary": sampling_summary,
        "partial": partial,
        "failed": failed,
        "relaxation_used": stats.relax_log,
        "final_diversity_summary": diversity,
        "timings": stats.timings,
        "written_files": stats.written_files,
        "accepted_sequences_summary": accepted_sequences_summary(decoys),
    }
    with open(cfg.run_report_json, "w", encoding="utf-8") as f:
        json.dump(report, f, indent=2)


def generate_decoys(cfg: Config, seed: int, stats: RunStats) -> Tuple[List[DecoyResult], Dict[str, Any], bool, bool, int]:
    poi_desc = describe(cfg.poi, cfg)
    print(f"[pepctrl] POI descriptor: {json.dumps(poi_desc.as_dict(), sort_keys=True)}")
    mode_meta = {
        "poi_counts": Counter(cfg.poi),
        "poi_class_counts": Counter(CLASS_MAP[a] for a in cfg.poi),
    }

    pool: List[DecoyResult] = []
    seen: Set[str] = set()
    target_pool = cfg.n * cfg.oversample_factor
    stats.target_pool_expected = target_pool
    stage_max = cfg.max_relax_stage if cfg.relaxation else 0

    for stage in range(stage_max + 1):
        if stage > 0 and cfg.decoy_class == "property-matched":
            stats.relax_log.append({
                "stage": stage,
                "scale": current_soft_scale(stage),
                "relaxed_parameters": sorted(cfg.property_tolerances.keys()),
            })
        while len(pool) < target_pool:
            idx = len(pool) + 1
            one_seed = seed + stage * 100000 + idx
            cand = anneal_one(cfg.poi, poi_desc, cfg, one_seed, mode_meta, stage, stats, seen)
            if cand:
                pool.append(cand)
                stats.total_pool_admissions += 1
                seen.add(cand.sequence)
                if len(pool) % cfg.progress_every == 0:
                    print(f"[pepctrl] progress valid_pool={len(pool)}/{target_pool} stage={stage}")
            if stats.generated_attempts > cfg.anneal.steps * target_pool * (stage + 1) * 4:
                break
        if len(pool) >= target_pool:
            break

    chosen, diversity = select_diverse(pool, cfg.n, cfg.max_identity_between_decoys)
    partial = len(chosen) < cfg.n
    failed = partial and not cfg.allow_partial
    stats.pool_valid_achieved = len(pool)
    stats.selected_final = len(chosen)
    stats.total_duplicate_rejections += diversity.get("rejected_by_diversity_filter", 0)
    return chosen, diversity, partial, failed, len(pool)


def generate_scramble_set(poi: str, n: int, rng: random.Random) -> List[str]:
    out, seen = [], set()
    while len(out) < n:
        s = list(poi)
        rng.shuffle(s)
        c = "".join(s)
        if c != poi and c not in seen:
            out.append(c)
            seen.add(c)
    return out


def generate_random_set(L: int, n: int, rng: random.Random) -> List[str]:
    out, seen = [], set()
    while len(out) < n:
        c = "".join(rng.choice(AA_LIST) for _ in range(L))
        if c not in seen:
            out.append(c)
            seen.add(c)
    return out


def main(argv: Optional[List[str]] = None) -> int:
    t0 = time.perf_counter()
    try:
        ns = parse_args(argv)
        poi = load_poi(ns)
        cfg = resolve_config(ns, poi)
    except Exception as e:
        print(f"[pepctrl] ERROR: {e}", file=sys.stderr)
        return 3

    seed = cfg.seed if cfg.seed is not None else random.SystemRandom().randint(1, 2**31 - 1)
    stats = RunStats()

    if cfg.seq_type == "decoy":
        t1 = time.perf_counter()
        decoys, diversity, partial, failed, pool_size = generate_decoys(cfg, seed, stats)
        stats.timings["generate_s"] = time.perf_counter() - t1
        poi_desc = describe(cfg.poi, cfg)
        t2 = time.perf_counter()
        write_decoys(decoys, cfg.poi, poi_desc, cfg, stats)
        stats.timings["write_s"] = time.perf_counter() - t2
        stats.timings["total_s"] = time.perf_counter() - t0
        if cfg.run_report_json:
            stats.written_files.append({
                "path": cfg.run_report_json,
                "kind": "run_report_json",
                "format": "json",
                "description": "reporte de ejecución con métricas de muestreo, rechazos y resumen final",
            })
        write_report(cfg, poi_desc, stats, pool_size, len(decoys), diversity, partial, failed, seed, decoys)
        print_stdout_summary(cfg, stats, poi_desc, decoys, partial, failed, seed)
        if failed:
            print(f"[pepctrl] partial failure: requested={cfg.n} selected={len(decoys)} allow with --allow-partial")
            return 2
        if partial:
            print(f"[pepctrl] partial output: requested={cfg.n} selected={len(decoys)}")
            return 2
        return 0

    seqs = generate_scramble_set(cfg.poi, cfg.n, rng) if cfg.seq_type == "scramble" else generate_random_set(len(cfg.poi), cfg.n, rng)
    for i, s in enumerate(seqs, 1):
        print(f">{cfg.seq_type[0].upper()}{i:04d}\n{s}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
