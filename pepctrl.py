#!/usr/bin/env python3
"""
pepctrl.py — Peptide controls: matched decoys, scrambled, and random sequences

Modes:
- decoy: multi-criteria matched decoys (fair controls)
- scramble: permutation of POI (same multiset; same length)
- random: random AA20 sequence (same length)

NEW: when output format is CSV/TSV, the script prepends a POI row so you can see
the POI metrics alongside generated sequences.

POI row conventions:
- For decoy output: decoy_id = "POI", type="poi", score="", refine_steps_used="", relax_stage="", flags="".
- For scramble/random output: seq_id="POI", type="poi", matches_to_poi=length, identity_to_poi=1.0

Default output: CSV (TSV optional). Manifest/logs go to stderr.

Author: trippm@tripplab.com [Feb 2026]
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
from dataclasses import dataclass, field, replace
from typing import Dict, List, Optional, Sequence, Set, Tuple


# ----------------------------
# Constants / scales (v1: KD)
# ----------------------------
AA20 = set("ACDEFGHIKLMNPQRSTVWY")

KD_HYDRO = {
    "I": 4.5, "V": 4.2, "L": 3.8, "F": 2.8, "C": 2.5,
    "M": 1.9, "A": 1.8, "G": -0.4, "T": -0.7, "S": -0.8,
    "W": -0.9, "Y": -1.3, "P": -1.6, "H": -3.2, "E": -3.5,
    "Q": -3.5, "D": -3.5, "N": -3.5, "K": -3.9, "R": -4.5,
}

AROM = set("FYW")
SPECIAL = {"C", "P", "G"}
POS = set("KR")
NEG = set("DE")


# ----------------------------
# Data classes
# ----------------------------
@dataclass
class Config:
    poi_seq: str
    n: int
    seq_type: str  # decoy|scramble|random
    out_path: str
    fmt: str
    no_header: bool

    # charge model (used in decoy and also for reporting in other modes)
    ph: float
    his_charge: float
    include_termini: bool

    # hydrophobicity (decoy)
    hydro_scale: str
    tol_hydro_mean: float

    # matching (decoy)
    match_arom: str  # exact|frac
    tol_arom_frac: float
    match_special: str  # exact|set|none

    # quality (decoy and optionally for others)
    max_identity: float
    max_matches: int
    max_identity_between_decoys: float
    max_run: int
    min_shannon: float
    no_simple_repeats: bool

    # motifs
    exclude_motifs_raw: List[str]
    exclude_regexes: List[re.Pattern] = field(default_factory=list, init=False)
    exclude_literals: List[str] = field(default_factory=list, init=False)

    # charge distribution (decoy)
    charge_shape_mode: str  # off|l1|auto
    tol_charge_shape: float
    tol_charge_polarization: float
    pol_margin: float

    # generation (decoy)
    method: str  # mutate|compose|hybrid
    batch_size: int
    max_attempts: int
    refine_steps: int
    min_mutations: int
    max_mutations: int
    seed: Optional[int]

    # convergence
    relax: bool
    relax_max_hydro: float
    relax_max_shape: float
    relax_arom: bool
    allow_partial: bool

    # logging
    verbosity: int
    timing: bool
    progress_every: int
    extra_cols: bool

    # optional quality filtering for scramble/random
    filter_scramble_random: bool

    # POI row in CSV/TSV
    include_poi_row: bool


@dataclass
class POIStats:
    L: int
    charge_vec: List[float]
    net_charge: float
    hydro_vec: List[float]
    hydro_mean: float
    arom_count: int
    c_count: int
    p_count: int
    g_count: int
    charge_abs_sum: float
    charge_pol: float
    charge_shape_norm: List[float]


@dataclass
class Candidate:
    seq: str
    score: float
    matches_to_poi: int
    identity_to_poi: float
    net_charge: float
    hydro_mean: float
    arom_count: int
    c_count: int
    p_count: int
    g_count: int
    charge_shape_l1: float
    charge_pol: float
    refine_steps_used: int
    relax_stage: int
    flags: Set[str] = field(default_factory=set)


@dataclass
class RunStats:
    generated_total: int = 0
    rejected_motif: int = 0
    rejected_low_complexity: int = 0
    rejected_identity: int = 0
    rejected_charge: int = 0
    rejected_arom: int = 0
    rejected_special: int = 0
    rejected_hydro: int = 0
    rejected_charge_shape: int = 0
    accepted_pool: int = 0

    # for scramble/random bookkeeping
    rejected_duplicate: int = 0

    t_parse: float = 0.0
    t_poi_stats: float = 0.0
    t_generate: float = 0.0
    t_select: float = 0.0
    t_write: float = 0.0


# ----------------------------
# CLI / config
# ----------------------------
def parse_args(argv: Optional[List[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        prog="pepctrl.py",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Generate peptide controls for a POI: matched decoys, scrambles, or random sequences."
    )

    # input/output
    p.add_argument("--poi", type=str, default=None, help="POI sequence (AA20).")
    p.add_argument("--poi-file", type=str, default=None, help="File containing POI sequence (first non-empty line).")
    p.add_argument("--type", dest="seq_type", choices=["decoy", "scramble", "random"], required=True,
                   help="Type of sequences to generate.")
    p.add_argument("--n", type=int, required=True, help="Number of sequences desired.")
    p.add_argument("--out", type=str, default="-", help="Output path, or '-' for stdout.")
    p.add_argument("--format", type=str, choices=["csv", "tsv", "fasta", "jsonl"], default="csv", help="Output format.")
    p.add_argument("--no-header", action="store_true", help="Disable header row (csv/tsv).")
    p.add_argument("--extra-cols", action="store_true", help="Include extra columns (seed, type).")
    p.add_argument("--no-poi-row", action="store_true",
                   help="Do not prepend a POI row (CSV/TSV only). Default: include POI row.")

    # charge model
    p.add_argument("--ph", type=float, default=7.0, help="pH (used only to report; model is simplified in v1).")
    p.add_argument("--his-charge", type=float, default=0.1, help="Charge assigned to histidine at this pH (simplified).")
    p.add_argument("--include-termini", action="store_true", help="Include N(+1) and C(-1) termini charges.")

    # hydrophobicity (decoy)
    p.add_argument("--hydro-scale", type=str, choices=["kd"], default="kd", help="Hydrophobicity scale (v1: kd only).")
    p.add_argument("--tol-hydro-mean", type=float, default=None, help="Tolerance for mean hydrophobicity match (auto if omitted).")

    # matching (decoy)
    p.add_argument("--match-arom", choices=["exact", "frac"], default="exact", help="Match aromatics by exact count or fraction.")
    p.add_argument("--tol-arom-frac", type=float, default=0.05, help="Tolerance on aromatic fraction (only for match-arom=frac).")
    p.add_argument("--match-special", choices=["exact", "set", "none"], default="exact",
                   help="Match special residues: exact (C,P,G individually), set (total C+P+G), none.")

    # quality (decoy; optionally used to filter scramble/random if enabled)
    p.add_argument("--max-identity", type=float, default=0.33, help="Max position-wise identity to POI (decoy).")
    p.add_argument("--max-matches", type=int, default=None, help="Max identical positions to POI (decoy; overrides max-identity).")
    p.add_argument("--max-identity-between-decoys", type=float, default=0.50, help="Max identity allowed between output sequences.")
    p.add_argument("--max-run", type=int, default=3, help="Max allowed run length of identical residue.")
    p.add_argument("--min-shannon", type=float, default=None, help="Min Shannon entropy (auto if omitted).")
    p.add_argument("--no-simple-repeats", action="store_true", default=True, help="Disallow simple k-mer repeats (k=2,3).")
    p.add_argument("--allow-simple-repeats", action="store_true", help="Allow simple repeats (overrides --no-simple-repeats).")

    # motifs
    p.add_argument("--exclude-motif", action="append", default=[],
                   help="Motif to exclude. Literal substring by default; regex if prefixed with 're:'. Repeatable.")
    p.add_argument("--exclude-motif-file", type=str, default=None,
                   help="File with one motif per line (literal; use 're:' prefix for regex).")

    # charge distribution (decoy)
    p.add_argument("--charge-shape", choices=["auto", "off", "l1"], default="auto",
                   help="Charge shape matching (auto=off if few charges; else l1).")
    p.add_argument("--tol-charge-shape", type=float, default=None, help="Tolerance for charge-shape L1 distance (auto if omitted).")
    p.add_argument("--tol-charge-polarization", type=float, default=0.15,
                   help="Allowed extra polarization beyond POI (units: charge/pos).")
    p.add_argument("--pol-margin", type=float, default=0.15, help="Alias of tol-charge-polarization.")

    # generation (decoy)
    p.add_argument("--method", choices=["mutate", "compose", "hybrid"], default="hybrid", help="Generation method (decoy).")
    p.add_argument("--batch-size", type=int, default=2000, help="Batch size for candidate generation (decoy).")
    p.add_argument("--max-attempts", type=int, default=None, help="Max generated candidates (decoy; auto if omitted).")
    p.add_argument("--refine-steps", type=int, default=50, help="Local refinement steps per candidate (decoy).")
    p.add_argument("--min-mutations", type=int, default=None, help="Min mutations for mutate/hybrid (decoy; auto if omitted).")
    p.add_argument("--max-mutations", type=int, default=None, help="Max mutations for mutate/hybrid (decoy; auto if omitted).")
    p.add_argument("--seed", type=int, default=None, help="RNG seed for reproducibility.")

    # convergence (decoy)
    p.add_argument("--relax", action="store_true", help="Enable controlled tolerance relaxation if insufficient decoys.")
    p.add_argument("--relax-max-hydro", type=float, default=0.30, help="Max hydro tolerance under relaxation.")
    p.add_argument("--relax-max-shape", type=float, default=1.20, help="Max shape tolerance under relaxation.")
    p.add_argument("--relax-arom", action="store_true", help="Allow relaxing aromatic matching (exact->frac) under relax.")
    p.add_argument("--allow-partial", action="store_true", help="Write partial results if insufficient sequences (exit code 2).")

    # scramble/random filtering (optional)
    p.add_argument("--filter-sr", action="store_true",
                   help="Apply motif/low-complexity/identity filters also to scramble/random outputs.")

    # logging/timing
    p.add_argument("-v", "--verbose", action="count", default=0, help="Verbosity (-v, -vv, -vvv).")
    p.add_argument("--timing", action="store_true", help="Print timing per stage to stderr.")
    p.add_argument("--progress-every", type=int, default=20000, help="Progress report cadence (generated candidates).")

    return p.parse_args(argv)


def load_poi_sequence(ns: argparse.Namespace) -> str:
    seq = None
    if ns.poi:
        seq = ns.poi.strip()
    elif ns.poi_file:
        with open(ns.poi_file, "r", encoding="utf-8") as f:
            for line in f:
                line = line.strip()
                if line:
                    seq = line
                    break
    if not seq:
        raise ValueError("POI sequence not provided. Use --poi or --poi-file.")
    seq = seq.strip().upper()
    validate_sequence(seq)
    return seq


def validate_sequence(seq: str) -> None:
    bad = sorted(set(seq) - AA20)
    if bad:
        raise ValueError(f"Invalid residue(s) in sequence: {bad}. Allowed: {''.join(sorted(AA20))}")


def auto_tol_hydro_mean(L: int) -> float:
    if L <= 10:
        return 0.18
    if L >= 25:
        return 0.10
    return 0.12


def auto_min_shannon(L: int) -> float:
    if L <= 12:
        return 2.2
    if L <= 25:
        return 2.5
    return 2.7


def auto_tol_charge_shape(abs_sum: float) -> float:
    if abs_sum <= 1.0:
        return 999.0
    if abs_sum <= 3.0:
        return 0.8
    return 0.6


def compile_excludes(excludes: List[str]) -> Tuple[List[str], List[re.Pattern]]:
    literals: List[str] = []
    regexes: List[re.Pattern] = []
    for m in excludes:
        m = m.strip()
        if not m:
            continue
        if m.startswith("re:"):
            regexes.append(re.compile(m[3:]))
        else:
            literals.append(m.upper())
    return literals, regexes


def resolve_config(ns: argparse.Namespace, poi_seq: str) -> Config:
    L = len(poi_seq)

    tol_h = ns.tol_hydro_mean if ns.tol_hydro_mean is not None else auto_tol_hydro_mean(L)
    min_sh = ns.min_shannon if ns.min_shannon is not None else auto_min_shannon(L)
    max_attempts = ns.max_attempts if ns.max_attempts is not None else max(200000, 5000 * ns.n)

    if ns.max_matches is not None:
        max_matches = ns.max_matches
    else:
        max_matches = int(math.floor(ns.max_identity * L + 1e-9))

    if ns.min_mutations is not None:
        min_mut = ns.min_mutations
    else:
        min_mut = max(1, int(math.ceil(0.5 * L)))

    if ns.max_mutations is not None:
        max_mut = ns.max_mutations
    else:
        max_mut = max(min_mut, min(L, int(math.ceil(0.8 * L))))

    motifs = list(ns.exclude_motif or [])
    if ns.exclude_motif_file:
        with open(ns.exclude_motif_file, "r", encoding="utf-8") as f:
            for line in f:
                line = line.strip()
                if line and not line.startswith("#"):
                    motifs.append(line)

    literals, regexes = compile_excludes(motifs)

    no_simple_repeats = ns.no_simple_repeats and (not ns.allow_simple_repeats)
    pol_margin = ns.pol_margin if ns.pol_margin is not None else ns.tol_charge_polarization
    tol_shape = ns.tol_charge_shape if ns.tol_charge_shape is not None else -1.0

    cfg = Config(
        poi_seq=poi_seq,
        n=ns.n,
        seq_type=ns.seq_type,
        out_path=ns.out,
        fmt=ns.format,
        no_header=ns.no_header,

        ph=ns.ph,
        his_charge=ns.his_charge,
        include_termini=ns.include_termini,

        hydro_scale=ns.hydro_scale,
        tol_hydro_mean=tol_h,

        match_arom=ns.match_arom,
        tol_arom_frac=ns.tol_arom_frac,
        match_special=ns.match_special,

        max_identity=ns.max_identity,
        max_matches=max_matches,
        max_identity_between_decoys=ns.max_identity_between_decoys,
        max_run=ns.max_run,
        min_shannon=min_sh,
        no_simple_repeats=no_simple_repeats,

        exclude_motifs_raw=motifs,

        charge_shape_mode=ns.charge_shape,
        tol_charge_shape=tol_shape,
        tol_charge_polarization=ns.tol_charge_polarization,
        pol_margin=pol_margin,

        method=ns.method,
        batch_size=ns.batch_size,
        max_attempts=max_attempts,
        refine_steps=ns.refine_steps,
        min_mutations=min_mut,
        max_mutations=max_mut,
        seed=ns.seed,

        relax=ns.relax,
        relax_max_hydro=ns.relax_max_hydro,
        relax_max_shape=ns.relax_max_shape,
        relax_arom=ns.relax_arom,
        allow_partial=ns.allow_partial,

        verbosity=ns.verbose or 0,
        timing=ns.timing,
        progress_every=ns.progress_every,
        extra_cols=ns.extra_cols,

        filter_scramble_random=ns.filter_sr,
        include_poi_row=(not ns.no_poi_row),
    )

    cfg.exclude_regexes = regexes
    cfg.exclude_literals = literals
    return cfg


# ----------------------------
# Metrics
# ----------------------------
def charge_vector(seq: str, his_charge: float, include_termini: bool) -> List[float]:
    q = []
    for aa in seq:
        if aa in POS:
            q.append(1.0)
        elif aa in NEG:
            q.append(-1.0)
        elif aa == "H":
            q.append(float(his_charge))
        else:
            q.append(0.0)
    if include_termini and q:
        q[0] += 1.0
        q[-1] -= 1.0
    return q


def net_charge(q: Sequence[float]) -> float:
    return float(sum(q))


def hydro_profile(seq: str, scale: str) -> List[float]:
    if scale != "kd":
        raise ValueError(f"Unsupported hydro scale: {scale}")
    return [KD_HYDRO[aa] for aa in seq]


def mean(xs: Sequence[float]) -> float:
    return float(sum(xs) / len(xs)) if xs else 0.0


def aromatic_count(seq: str) -> int:
    return sum(1 for aa in seq if aa in AROM)


def special_counts(seq: str) -> Tuple[int, int, int]:
    return (seq.count("C"), seq.count("P"), seq.count("G"))


def identity_stats(seq: str, poi: str) -> Tuple[int, float]:
    matches = sum(1 for a, b in zip(seq, poi) if a == b)
    return matches, matches / len(poi)


def max_run_length(seq: str) -> int:
    if not seq:
        return 0
    best = 1
    cur = 1
    for i in range(1, len(seq)):
        if seq[i] == seq[i - 1]:
            cur += 1
            best = max(best, cur)
        else:
            cur = 1
    return best


def has_simple_repeats(seq: str) -> bool:
    L = len(seq)
    if L >= 6:
        for i in range(L - 5):
            k2 = seq[i:i + 2]
            if seq[i:i + 6] == k2 * 3:
                return True
    if L >= 6:
        for i in range(L - 5):
            k3 = seq[i:i + 3]
            if seq[i:i + 6] == k3 * 2:
                return True
    return False


def shannon_entropy(seq: str) -> float:
    L = len(seq)
    if L == 0:
        return 0.0
    counts: Dict[str, int] = {}
    for aa in seq:
        counts[aa] = counts.get(aa, 0) + 1
    ent = 0.0
    for c in counts.values():
        p = c / L
        ent -= p * math.log2(p)
    return ent


def charge_shape_norm(q: Sequence[float]) -> List[float]:
    denom = sum(abs(x) for x in q) + 1e-12
    return [x / denom for x in q]


def charge_shape_l1(q_norm_decoy: Sequence[float], q_norm_poi: Sequence[float]) -> float:
    return float(sum(abs(a - b) for a, b in zip(q_norm_decoy, q_norm_poi)))


def charge_polarization(q: Sequence[float]) -> float:
    L = len(q)
    if L == 0:
        return 0.0
    mid = L // 2
    a = mean(q[:mid]) if mid > 0 else 0.0
    b = mean(q[mid:]) if mid < L else 0.0
    return abs(a - b)


def compute_poi_stats(poi: str, cfg: Config) -> POIStats:
    q = charge_vector(poi, cfg.his_charge, cfg.include_termini)
    netq = net_charge(q)
    hvec = hydro_profile(poi, cfg.hydro_scale)
    hmean = mean(hvec)
    aromc = aromatic_count(poi)
    cc, pc, gc = special_counts(poi)
    abs_sum = sum(abs(x) for x in q)
    pol = charge_polarization(q)
    qn = charge_shape_norm(q)
    return POIStats(
        L=len(poi),
        charge_vec=q,
        net_charge=netq,
        hydro_vec=hvec,
        hydro_mean=hmean,
        arom_count=aromc,
        c_count=cc,
        p_count=pc,
        g_count=gc,
        charge_abs_sum=abs_sum,
        charge_pol=pol,
        charge_shape_norm=qn,
    )


# ----------------------------
# Filters / scoring (decoy + optional for S/R)
# ----------------------------
def contains_excluded_motif(seq: str, cfg: Config) -> bool:
    s = seq.upper()
    for lit in cfg.exclude_literals:
        if lit and lit in s:
            return True
    for rx in cfg.exclude_regexes:
        if rx.search(seq):
            return True
    return False


def passes_low_complexity(seq: str, cfg: Config) -> bool:
    if max_run_length(seq) > cfg.max_run:
        return False
    if shannon_entropy(seq) < cfg.min_shannon:
        return False
    if cfg.no_simple_repeats and has_simple_repeats(seq):
        return False
    return True


def passes_identity(seq: str, poi: str, cfg: Config) -> Tuple[bool, int, float]:
    matches, ident = identity_stats(seq, poi)
    if matches > cfg.max_matches:
        return False, matches, ident
    if ident > cfg.max_identity + 1e-12:
        return False, matches, ident
    return True, matches, ident


def arom_ok(seq: str, poi_stats: POIStats, cfg: Config) -> bool:
    a = aromatic_count(seq)
    if cfg.match_arom == "exact":
        return a == poi_stats.arom_count
    frac = a / poi_stats.L
    poi_frac = poi_stats.arom_count / poi_stats.L
    return abs(frac - poi_frac) <= cfg.tol_arom_frac + 1e-12


def special_ok(seq: str, poi_stats: POIStats, cfg: Config) -> bool:
    c, p, g = special_counts(seq)
    if cfg.match_special == "none":
        return True
    if cfg.match_special == "set":
        return (c + p + g) == (poi_stats.c_count + poi_stats.p_count + poi_stats.g_count)
    return (c == poi_stats.c_count) and (p == poi_stats.p_count) and (g == poi_stats.g_count)


def score_candidate(
    hydro_mean: float,
    charge_shape: float,
    pol: float,
    ident: float,
    poi_stats: POIStats,
    cfg: Config,
) -> float:
    w_hydro = 1.0
    w_shape = 0.8 if poi_stats.charge_abs_sum > 1.0 else 0.0
    w_pol = 0.5 if poi_stats.charge_abs_sum > 1.0 else 0.0
    w_ident = 0.3

    d_h = abs(hydro_mean - poi_stats.hydro_mean)
    pol_excess = max(0.0, pol - (poi_stats.charge_pol + cfg.pol_margin))
    return w_hydro * d_h + w_shape * charge_shape + w_pol * pol_excess + w_ident * ident


def evaluate_decoy_candidate(
    seq: str,
    poi: str,
    poi_stats: POIStats,
    cfg: Config,
    run: RunStats,
    refine_steps_used: int,
    relax_stage: int,
    flags: Optional[Set[str]] = None
) -> Optional[Candidate]:
    run.generated_total += 1

    if contains_excluded_motif(seq, cfg):
        run.rejected_motif += 1
        return None

    if not passes_low_complexity(seq, cfg):
        run.rejected_low_complexity += 1
        return None

    ok_id, matches, ident = passes_identity(seq, poi, cfg)
    if not ok_id:
        run.rejected_identity += 1
        return None

    q = charge_vector(seq, cfg.his_charge, cfg.include_termini)
    netq = net_charge(q)
    if abs(netq - poi_stats.net_charge) > 1e-9:
        run.rejected_charge += 1
        return None

    if not arom_ok(seq, poi_stats, cfg):
        run.rejected_arom += 1
        return None

    if not special_ok(seq, poi_stats, cfg):
        run.rejected_special += 1
        return None

    hmean = mean(hydro_profile(seq, cfg.hydro_scale))
    if abs(hmean - poi_stats.hydro_mean) > cfg.tol_hydro_mean + 1e-12:
        run.rejected_hydro += 1
        return None

    pol = charge_polarization(q)

    if cfg.charge_shape_mode == "off":
        shape = 0.0
    else:
        if poi_stats.charge_abs_sum <= 1.0 and cfg.charge_shape_mode in ("auto", "l1"):
            shape = 0.0
        else:
            qn = charge_shape_norm(q)
            shape = charge_shape_l1(qn, poi_stats.charge_shape_norm)
            if shape > cfg.tol_charge_shape + 1e-12:
                run.rejected_charge_shape += 1
                return None
            if pol > (poi_stats.charge_pol + cfg.tol_charge_polarization + 1e-12):
                run.rejected_charge_shape += 1
                return None

    sc = score_candidate(hmean, shape, pol, ident, poi_stats, cfg)
    return Candidate(
        seq=seq,
        score=sc,
        matches_to_poi=matches,
        identity_to_poi=ident,
        net_charge=netq,
        hydro_mean=hmean,
        arom_count=aromatic_count(seq),
        c_count=seq.count("C"),
        p_count=seq.count("P"),
        g_count=seq.count("G"),
        charge_shape_l1=shape,
        charge_pol=pol,
        refine_steps_used=refine_steps_used,
        relax_stage=relax_stage,
        flags=set(flags or []),
    )


# ----------------------------
# Scramble (S) and Random (R)
# ----------------------------
def generate_scramble(poi: str, rng: random.Random, max_tries: int = 200) -> str:
    s = list(poi)
    for _ in range(max_tries):
        rng.shuffle(s)
        cand = "".join(s)
        if cand != poi:
            return cand
    return "".join(s)


def generate_random(length: int, rng: random.Random) -> str:
    aa_list = sorted(AA20)
    return "".join(rng.choice(aa_list) for _ in range(length))


# ----------------------------
# Decoy generation helpers
# ----------------------------
def target_counts(poi: str) -> Dict[str, int]:
    return {
        "L": len(poi),
        "arom": aromatic_count(poi),
        "C": poi.count("C"),
        "P": poi.count("P"),
        "G": poi.count("G"),
        "pos": sum(1 for aa in poi if aa in POS),
        "neg": sum(1 for aa in poi if aa in NEG),
        "H": poi.count("H"),
    }


def build_candidate_from_composition(tgt: Dict[str, int], rng: random.Random) -> str:
    L = tgt["L"]
    seq_list: List[str] = []

    seq_list += ["C"] * tgt["C"]
    seq_list += ["P"] * tgt["P"]
    seq_list += ["G"] * tgt["G"]

    for _ in range(tgt["arom"]):
        seq_list.append(rng.choice(["F", "Y", "W"]))

    for _ in range(tgt["pos"]):
        seq_list.append(rng.choice(["K", "R"]))
    for _ in range(tgt["neg"]):
        seq_list.append(rng.choice(["D", "E"]))

    seq_list += ["H"] * tgt["H"]

    remaining = L - len(seq_list)
    if remaining < 0:
        raise RuntimeError("Internal error: negative remaining length during composition build.")

    forbidden = set("CPGFYWKRDEH")
    neutral_pool = [aa for aa in sorted(AA20 - forbidden)]
    for _ in range(remaining):
        seq_list.append(rng.choice(neutral_pool))

    rng.shuffle(seq_list)
    return "".join(seq_list)


def mutate_within_class(seq: str, rng: random.Random) -> str:
    L = len(seq)
    if L < 2:
        return seq
    op = rng.random()
    s = list(seq)

    if op < 0.45:
        i, j = rng.randrange(L), rng.randrange(L)
        if i != j:
            s[i], s[j] = s[j], s[i]
        return "".join(s)

    i = rng.randrange(L)
    aa = s[i]
    if aa in POS:
        s[i] = "K" if aa == "R" else "R"
    elif aa in NEG:
        s[i] = "D" if aa == "E" else "E"
    elif aa in AROM:
        choices = ["F", "Y", "W"]
        choices.remove(aa)
        s[i] = rng.choice(choices)
    else:
        if aa in SPECIAL:
            return seq
        forbidden = set("CPGFYWKRDEH")
        neutral_pool = [x for x in sorted(AA20 - forbidden)]
        if aa in neutral_pool and len(neutral_pool) > 1:
            cand = rng.choice(neutral_pool)
            while cand == aa:
                cand = rng.choice(neutral_pool)
            s[i] = cand
    return "".join(s)


def provisional_metrics(seq: str, poi: str, poi_stats: POIStats, cfg: Config) -> Dict[str, float]:
    hmean = mean(hydro_profile(seq, cfg.hydro_scale))
    q = charge_vector(seq, cfg.his_charge, cfg.include_termini)
    pol = charge_polarization(q)
    shape = 0.0
    if poi_stats.charge_abs_sum > 1.0:
        shape = charge_shape_l1(charge_shape_norm(q), poi_stats.charge_shape_norm)
    _, ident = identity_stats(seq, poi)
    prov = score_candidate(hmean, shape, pol, ident, poi_stats, cfg)
    return {"hydro_mean": hmean, "shape": shape, "pol": pol, "ident": ident, "prov_score": prov}


def refine_candidate(seq: str, poi: str, poi_stats: POIStats, cfg: Config, rng: random.Random) -> Tuple[str, int, Set[str]]:
    flags: Set[str] = set()
    best = seq
    bestm = provisional_metrics(best, poi, poi_stats, cfg)

    steps_used = 0
    for step in range(cfg.refine_steps):
        prop = mutate_within_class(best, rng)

        if contains_excluded_motif(prop, cfg):
            continue
        if not passes_low_complexity(prop, cfg):
            continue
        ok_id, _, _ = passes_identity(prop, poi, cfg)
        if not ok_id:
            continue
        q = charge_vector(prop, cfg.his_charge, cfg.include_termini)
        if abs(net_charge(q) - poi_stats.net_charge) > 1e-9:
            continue

        m = provisional_metrics(prop, poi, poi_stats, cfg)
        if m["prov_score"] + 1e-12 < bestm["prov_score"]:
            best = prop
            bestm = m
            steps_used = step + 1
            flags.add("refined")
            if abs(m["hydro_mean"] - poi_stats.hydro_mean) <= min(0.02, cfg.tol_hydro_mean * 0.25):
                flags.add("early_stop")
                break

    return best, steps_used, flags


def relax_config(cfg: Config, stage: int) -> Config:
    if stage == 0:
        return cfg
    new = replace(cfg)
    new.tol_hydro_mean = min(cfg.relax_max_hydro, cfg.tol_hydro_mean + 0.03 * stage)
    if cfg.charge_shape_mode != "off" and cfg.tol_charge_shape < 900:
        new.tol_charge_shape = min(cfg.relax_max_shape, cfg.tol_charge_shape + 0.1 * stage)
    if cfg.relax_arom and stage >= 3 and cfg.match_arom == "exact":
        new.match_arom = "frac"
    return new


def generate_decoy_pool(
    poi: str,
    poi_stats: POIStats,
    base_cfg: Config,
    rng: random.Random,
    run: RunStats
) -> Tuple[List[Candidate], int]:
    tgt = target_counts(poi)
    pool: List[Candidate] = []
    seen: Set[str] = set()

    cfg0 = replace(base_cfg)
    if cfg0.charge_shape_mode == "auto":
        cfg0.charge_shape_mode = "off" if poi_stats.charge_abs_sum <= 1.0 else "l1"
    if cfg0.tol_charge_shape < 0:
        cfg0.tol_charge_shape = auto_tol_charge_shape(poi_stats.charge_abs_sum)

    max_stage = 0 if not cfg0.relax else 6
    relax_stage_used = 0

    attempts = 0
    stage = 0
    while stage <= max_stage:
        cfg = relax_config(cfg0, stage)
        relax_stage_used = stage

        while attempts < cfg.max_attempts:
            if cfg.verbosity >= 1 and attempts > 0 and (attempts % cfg.progress_every == 0):
                print(f"[pepctrl] generated={attempts} pool={len(pool)} stage={stage}", file=sys.stderr)

            for _ in range(cfg.batch_size):
                attempts += 1
                if attempts > cfg.max_attempts:
                    break

                if cfg.method in ("compose", "hybrid"):
                    seq = build_candidate_from_composition(tgt, rng)
                else:
                    seq = poi

                if cfg.method in ("mutate", "hybrid"):
                    k = rng.randint(cfg.min_mutations, cfg.max_mutations)
                    for _m in range(k):
                        seq = mutate_within_class(seq, rng)

                refine_steps_used = 0
                flags: Set[str] = set()
                if cfg.method in ("hybrid", "mutate") and cfg.refine_steps > 0:
                    seq, refine_steps_used, flags = refine_candidate(seq, poi, poi_stats, cfg, rng)

                if seq in seen:
                    continue
                seen.add(seq)

                cand = evaluate_decoy_candidate(
                    seq, poi, poi_stats, cfg, run,
                    refine_steps_used=refine_steps_used,
                    relax_stage=stage,
                    flags=flags
                )
                if cand is not None:
                    pool.append(cand)

            if pool:
                pool.sort(key=lambda c: c.score)
                pool_max = max(5 * cfg.n, 200)
                if len(pool) > pool_max:
                    pool = pool[:pool_max]

            if len(pool) >= max(cfg.n, 20):
                break

        if len(pool) >= max(cfg.n, 20):
            break

        stage += 1

    run.accepted_pool = len(pool)
    return pool, relax_stage_used


# ----------------------------
# Diversity helpers
# ----------------------------
def hamming_distance(a: str, b: str) -> int:
    return sum(1 for x, y in zip(a, b) if x != y)


def identity_between(a: str, b: str) -> float:
    L = len(a)
    if L == 0:
        return 0.0
    matches = sum(1 for x, y in zip(a, b) if x == y)
    return matches / L


def select_diverse_sequences(seqs: List[str], n: int, max_ident_between: float) -> List[str]:
    if not seqs:
        return []
    chosen = [seqs[0]]
    remaining = seqs[1:]

    while len(chosen) < n and remaining:
        best_idx = None
        best_min_dist = -1
        for i, s in enumerate(remaining):
            if any(identity_between(s, c) > max_ident_between + 1e-12 for c in chosen):
                continue
            min_dist = min(hamming_distance(s, c) for c in chosen)
            if min_dist > best_min_dist:
                best_min_dist = min_dist
                best_idx = i
        if best_idx is None:
            chosen.append(remaining.pop(0))
        else:
            chosen.append(remaining.pop(best_idx))

    return chosen[:n]


def select_diverse_decoys(pool: List[Candidate], n: int, max_ident_between: float) -> List[Candidate]:
    if not pool:
        return []
    pool_sorted = sorted(pool, key=lambda c: c.score)
    chosen = [pool_sorted[0]]
    remaining = pool_sorted[1:]

    while len(chosen) < n and remaining:
        best_idx = None
        best_min_dist = -1
        for i, cand in enumerate(remaining):
            if any(identity_between(cand.seq, c.seq) > max_ident_between + 1e-12 for c in chosen):
                continue
            min_dist = min(hamming_distance(cand.seq, c.seq) for c in chosen)
            if min_dist > best_min_dist:
                best_min_dist = min_dist
                best_idx = i
        if best_idx is None:
            chosen.append(remaining.pop(0))
        else:
            chosen.append(remaining.pop(best_idx))

    return chosen[:n]


# ----------------------------
# S/R set generation
# ----------------------------
def generate_scramble_set(poi: str, cfg: Config, rng: random.Random, run: RunStats) -> List[str]:
    seen: Set[str] = set()
    out: List[str] = []
    attempts = 0
    max_attempts = max(cfg.max_attempts, 20000)

    while len(out) < cfg.n and attempts < max_attempts:
        attempts += 1
        run.generated_total += 1
        s = generate_scramble(poi, rng)

        if s in seen:
            run.rejected_duplicate += 1
            continue
        seen.add(s)

        if cfg.filter_scramble_random:
            if contains_excluded_motif(s, cfg):
                run.rejected_motif += 1
                continue
            if not passes_low_complexity(s, cfg):
                run.rejected_low_complexity += 1
                continue
            ok_id, _, _ = passes_identity(s, poi, cfg)
            if not ok_id:
                run.rejected_identity += 1
                continue

        out.append(s)

        if cfg.verbosity >= 1 and attempts % cfg.progress_every == 0:
            print(f"[pepctrl] (scramble) attempts={attempts} collected={len(out)}", file=sys.stderr)

    return out


def generate_random_set(poi: str, cfg: Config, rng: random.Random, run: RunStats) -> List[str]:
    seen: Set[str] = set()
    out: List[str] = []
    attempts = 0
    max_attempts = max(cfg.max_attempts, 20000)

    while len(out) < cfg.n and attempts < max_attempts:
        attempts += 1
        run.generated_total += 1
        s = generate_random(len(poi), rng)

        if s in seen:
            run.rejected_duplicate += 1
            continue
        seen.add(s)

        if cfg.filter_scramble_random:
            if contains_excluded_motif(s, cfg):
                run.rejected_motif += 1
                continue
            if not passes_low_complexity(s, cfg):
                run.rejected_low_complexity += 1
                continue
            ok_id, _, _ = passes_identity(s, poi, cfg)
            if not ok_id:
                run.rejected_identity += 1
                continue

        out.append(s)

        if cfg.verbosity >= 1 and attempts % cfg.progress_every == 0:
            print(f"[pepctrl] (random) attempts={attempts} collected={len(out)}", file=sys.stderr)

    return out


# ----------------------------
# Output helpers (with POI row)
# ----------------------------
def open_out(path: str):
    if path == "-" or path == "":
        return sys.stdout
    return open(path, "w", encoding="utf-8", newline="")


def write_output_sequences(seqs: List[str], poi: str, cfg: Config, seed_used: int) -> None:
    cols = ["seq_id", "type", "sequence", "length", "matches_to_poi", "identity_to_poi"]
    if cfg.extra_cols:
        cols += ["seed"]

    out_handle = open_out(cfg.out_path)
    close_after = out_handle is not sys.stdout

    try:
        if cfg.fmt in ("csv", "tsv"):
            delimiter = "," if cfg.fmt == "csv" else "\t"
            quoting = csv.QUOTE_MINIMAL if cfg.fmt == "csv" else csv.QUOTE_NONE
            writer = csv.DictWriter(
                out_handle,
                fieldnames=cols,
                delimiter=delimiter,
                quotechar='"',
                quoting=quoting,
                escapechar="\\" if cfg.fmt == "tsv" else None,
                lineterminator="\n",
            )
            if not cfg.no_header:
                writer.writeheader()

            # POI row
            if cfg.include_poi_row:
                row = {
                    "seq_id": "POI",
                    "type": "poi",
                    "sequence": poi,
                    "length": len(poi),
                    "matches_to_poi": len(poi),
                    "identity_to_poi": "1",
                }
                if cfg.extra_cols:
                    row["seed"] = seed_used
                writer.writerow(row)

            for i, s in enumerate(seqs, start=1):
                matches, ident = identity_stats(s, poi)
                row = {
                    "seq_id": f"{cfg.seq_type[0].upper()}{i:04d}",
                    "type": cfg.seq_type,
                    "sequence": s,
                    "length": len(poi),
                    "matches_to_poi": matches,
                    "identity_to_poi": f"{ident:.6g}",
                }
                if cfg.extra_cols:
                    row["seed"] = seed_used
                writer.writerow(row)

        elif cfg.fmt == "fasta":
            for i, s in enumerate(seqs, start=1):
                matches, ident = identity_stats(s, poi)
                header = f">{cfg.seq_type[0].upper()}{i:04d} type={cfg.seq_type} ident={ident:.3f} matches={matches}"
                out_handle.write(header + "\n")
                out_handle.write(s + "\n")

        elif cfg.fmt == "jsonl":
            # POI row
            if cfg.include_poi_row:
                obj = {
                    "seq_id": "POI",
                    "type": "poi",
                    "sequence": poi,
                    "length": len(poi),
                    "matches_to_poi": len(poi),
                    "identity_to_poi": 1.0,
                }
                if cfg.extra_cols:
                    obj["seed"] = seed_used
                out_handle.write(json.dumps(obj) + "\n")

            for i, s in enumerate(seqs, start=1):
                matches, ident = identity_stats(s, poi)
                obj = {
                    "seq_id": f"{cfg.seq_type[0].upper()}{i:04d}",
                    "type": cfg.seq_type,
                    "sequence": s,
                    "length": len(poi),
                    "matches_to_poi": matches,
                    "identity_to_poi": ident,
                }
                if cfg.extra_cols:
                    obj["seed"] = seed_used
                out_handle.write(json.dumps(obj) + "\n")
        else:
            raise ValueError(f"Unsupported format: {cfg.fmt}")
    finally:
        if close_after:
            out_handle.close()


def write_output_decoys(decoys: List[Candidate], poi_stats: POIStats, cfg: Config, seed_used: int) -> None:
    cols = [
        "decoy_id", "type", "sequence", "score", "length",
        "matches_to_poi", "identity_to_poi",
        "net_charge", "hydro_mean",
        "arom_count", "c_count", "p_count", "g_count",
        "charge_shape_l1", "charge_pol",
        "refine_steps_used", "relax_stage", "flags",
    ]
    if cfg.extra_cols:
        cols += ["seed"]

    out_handle = open_out(cfg.out_path)
    close_after = out_handle is not sys.stdout

    try:
        if cfg.fmt in ("csv", "tsv"):
            delimiter = "," if cfg.fmt == "csv" else "\t"
            quoting = csv.QUOTE_MINIMAL if cfg.fmt == "csv" else csv.QUOTE_NONE
            writer = csv.DictWriter(
                out_handle,
                fieldnames=cols,
                delimiter=delimiter,
                quotechar='"',
                quoting=quoting,
                escapechar="\\" if cfg.fmt == "tsv" else None,
                lineterminator="\n",
            )
            if not cfg.no_header:
                writer.writeheader()

            # POI row (metrics)
            if cfg.include_poi_row:
                matches, ident = identity_stats(cfg.poi_seq, cfg.poi_seq)
                row = {
                    "decoy_id": "POI",
                    "type": "poi",
                    "sequence": cfg.poi_seq,
                    "score": "",
                    "length": poi_stats.L,
                    "matches_to_poi": matches,
                    "identity_to_poi": f"{ident:.6g}",
                    "net_charge": f"{poi_stats.net_charge:.6g}",
                    "hydro_mean": f"{poi_stats.hydro_mean:.6g}",
                    "arom_count": poi_stats.arom_count,
                    "c_count": poi_stats.c_count,
                    "p_count": poi_stats.p_count,
                    "g_count": poi_stats.g_count,
                    "charge_shape_l1": "0",
                    "charge_pol": f"{poi_stats.charge_pol:.6g}",
                    "refine_steps_used": "",
                    "relax_stage": "",
                    "flags": "",
                }
                if cfg.extra_cols:
                    row["seed"] = seed_used
                writer.writerow(row)

            for i, d in enumerate(decoys, start=1):
                row = {
                    "decoy_id": f"D{i:04d}",
                    "type": "decoy",
                    "sequence": d.seq,
                    "score": f"{d.score:.6g}",
                    "length": poi_stats.L,
                    "matches_to_poi": d.matches_to_poi,
                    "identity_to_poi": f"{d.identity_to_poi:.6g}",
                    "net_charge": f"{d.net_charge:.6g}",
                    "hydro_mean": f"{d.hydro_mean:.6g}",
                    "arom_count": d.arom_count,
                    "c_count": d.c_count,
                    "p_count": d.p_count,
                    "g_count": d.g_count,
                    "charge_shape_l1": f"{d.charge_shape_l1:.6g}",
                    "charge_pol": f"{d.charge_pol:.6g}",
                    "refine_steps_used": d.refine_steps_used,
                    "relax_stage": d.relax_stage,
                    "flags": ";".join(sorted(d.flags)) if d.flags else "",
                }
                if cfg.extra_cols:
                    row["seed"] = seed_used
                writer.writerow(row)

        elif cfg.fmt == "fasta":
            for i, d in enumerate(decoys, start=1):
                header = (
                    f">D{i:04d} type=decoy score={d.score:.4g} ident={d.identity_to_poi:.3f} "
                    f"netq={d.net_charge:.3g}"
                )
                out_handle.write(header + "\n")
                out_handle.write(d.seq + "\n")

        elif cfg.fmt == "jsonl":
            if cfg.include_poi_row:
                obj = {
                    "decoy_id": "POI",
                    "type": "poi",
                    "sequence": cfg.poi_seq,
                    "length": poi_stats.L,
                    "matches_to_poi": poi_stats.L,
                    "identity_to_poi": 1.0,
                    "net_charge": poi_stats.net_charge,
                    "hydro_mean": poi_stats.hydro_mean,
                    "arom_count": poi_stats.arom_count,
                    "c_count": poi_stats.c_count,
                    "p_count": poi_stats.p_count,
                    "g_count": poi_stats.g_count,
                    "charge_shape_l1": 0.0,
                    "charge_pol": poi_stats.charge_pol,
                }
                if cfg.extra_cols:
                    obj["seed"] = seed_used
                out_handle.write(json.dumps(obj) + "\n")

            for i, d in enumerate(decoys, start=1):
                obj = {
                    "decoy_id": f"D{i:04d}",
                    "type": "decoy",
                    "sequence": d.seq,
                    "score": d.score,
                    "length": poi_stats.L,
                    "matches_to_poi": d.matches_to_poi,
                    "identity_to_poi": d.identity_to_poi,
                    "net_charge": d.net_charge,
                    "hydro_mean": d.hydro_mean,
                    "arom_count": d.arom_count,
                    "c_count": d.c_count,
                    "p_count": d.p_count,
                    "g_count": d.g_count,
                    "charge_shape_l1": d.charge_shape_l1,
                    "charge_pol": d.charge_pol,
                    "refine_steps_used": d.refine_steps_used,
                    "relax_stage": d.relax_stage,
                    "flags": sorted(d.flags),
                }
                if cfg.extra_cols:
                    obj["seed"] = seed_used
                out_handle.write(json.dumps(obj) + "\n")
        else:
            raise ValueError(f"Unsupported format: {cfg.fmt}")
    finally:
        if close_after:
            out_handle.close()


# ----------------------------
# Manifest
# ----------------------------
def print_manifest(run: RunStats, poi_stats: POIStats, cfg: Config, seed_used: int, relax_stage_used: int) -> None:
    if cfg.verbosity <= 0:
        return

    def kv(k, v):
        print(f"[pepctrl] {k}: {v}", file=sys.stderr)

    kv("type", cfg.seq_type)
    kv("poi_length", poi_stats.L)
    kv("seed", seed_used)
    kv("format", cfg.fmt)
    kv("include_poi_row", cfg.include_poi_row)
    kv("filter_sr", cfg.filter_scramble_random)

    if cfg.seq_type == "decoy":
        kv("tol_hydro_mean", cfg.tol_hydro_mean)
        kv("max_identity", cfg.max_identity)
        kv("max_matches", cfg.max_matches)
        kv("match_arom", cfg.match_arom)
        kv("match_special", cfg.match_special)
        kv("include_termini", cfg.include_termini)
        kv("his_charge", cfg.his_charge)
        kv("charge_shape_mode", cfg.charge_shape_mode)
        kv("tol_charge_shape", cfg.tol_charge_shape)
        kv("tol_charge_polarization", cfg.tol_charge_polarization)
        kv("relax_enabled", cfg.relax)
        kv("relax_stage_used", relax_stage_used)

    if cfg.exclude_motifs_raw:
        kv("exclude_motifs", cfg.exclude_motifs_raw)

    kv("generated_total", run.generated_total)
    kv("rejected_duplicate", run.rejected_duplicate)

    if cfg.seq_type == "decoy":
        kv("accepted_pool", run.accepted_pool)
        kv("rejected_motif", run.rejected_motif)
        kv("rejected_low_complexity", run.rejected_low_complexity)
        kv("rejected_identity", run.rejected_identity)
        kv("rejected_charge", run.rejected_charge)
        kv("rejected_arom", run.rejected_arom)
        kv("rejected_special", run.rejected_special)
        kv("rejected_hydro", run.rejected_hydro)
        kv("rejected_charge_shape", run.rejected_charge_shape)
    else:
        if cfg.filter_scramble_random:
            kv("rejected_motif", run.rejected_motif)
            kv("rejected_low_complexity", run.rejected_low_complexity)
            kv("rejected_identity", run.rejected_identity)

    if cfg.timing:
        kv("t_parse_s", f"{run.t_parse:.4f}")
        kv("t_poi_stats_s", f"{run.t_poi_stats:.4f}")
        kv("t_generate_s", f"{run.t_generate:.4f}")
        kv("t_select_s", f"{run.t_select:.4f}")
        kv("t_write_s", f"{run.t_write:.4f}")


# ----------------------------
# Main
# ----------------------------
def main(argv: Optional[List[str]] = None) -> int:
    run = RunStats()
    t0 = time.perf_counter()
    try:
        ns = parse_args(argv)
        poi = load_poi_sequence(ns)
        cfg = resolve_config(ns, poi)
    except Exception as e:
        print(f"[pepctrl] ERROR: {e}", file=sys.stderr)
        return 3
    run.t_parse = time.perf_counter() - t0

    # RNG seed
    if cfg.seed is None:
        seed_used = random.SystemRandom().randint(1, 2**31 - 1)
    else:
        seed_used = cfg.seed
    rng = random.Random(seed_used)

    # POI stats
    t1 = time.perf_counter()
    poi_stats = compute_poi_stats(poi, cfg)
    run.t_poi_stats = time.perf_counter() - t1

    # finalize charge-shape default tol if needed (decoy mode)
    if cfg.seq_type == "decoy":
        if cfg.charge_shape_mode == "auto":
            cfg.charge_shape_mode = "off" if poi_stats.charge_abs_sum <= 1.0 else "l1"
        if cfg.tol_charge_shape < 0:
            cfg.tol_charge_shape = auto_tol_charge_shape(poi_stats.charge_abs_sum)

    # generate + write
    t2 = time.perf_counter()
    relax_stage_used = 0
    exit_code = 0

    if cfg.seq_type == "decoy":
        pool, relax_stage_used = generate_decoy_pool(poi, poi_stats, cfg, rng, run)
        chosen = select_diverse_decoys(pool, cfg.n, cfg.max_identity_between_decoys)
        run.t_generate = time.perf_counter() - t2

        if len(chosen) < cfg.n and not cfg.allow_partial:
            if cfg.fmt in ("csv", "tsv"):
                t4 = time.perf_counter()
                write_output_decoys([], poi_stats, replace(cfg), seed_used)
                run.t_write = time.perf_counter() - t4
            print_manifest(run, poi_stats, cfg, seed_used, relax_stage_used)
            if cfg.verbosity >= 1:
                print(f"[pepctrl] Insufficient decoys: requested={cfg.n} found={len(chosen)} (use --allow-partial or --relax).",
                      file=sys.stderr)
            return 2

        t4 = time.perf_counter()
        write_output_decoys(chosen, poi_stats, cfg, seed_used)
        run.t_write = time.perf_counter() - t4
        if len(chosen) < cfg.n:
            exit_code = 2

    elif cfg.seq_type == "scramble":
        seqs = generate_scramble_set(poi, cfg, rng, run)
        run.t_generate = time.perf_counter() - t2
        seqs = select_diverse_sequences(seqs, min(cfg.n, len(seqs)), cfg.max_identity_between_decoys)

        if len(seqs) < cfg.n and not cfg.allow_partial:
            if cfg.fmt in ("csv", "tsv"):
                t4 = time.perf_counter()
                write_output_sequences([], poi, replace(cfg), seed_used)
                run.t_write = time.perf_counter() - t4
            print_manifest(run, poi_stats, cfg, seed_used, relax_stage_used)
            if cfg.verbosity >= 1:
                print(f"[pepctrl] Insufficient scrambles: requested={cfg.n} found={len(seqs)} (use --allow-partial).",
                      file=sys.stderr)
            return 2

        t4 = time.perf_counter()
        write_output_sequences(seqs, poi, cfg, seed_used)
        run.t_write = time.perf_counter() - t4
        if len(seqs) < cfg.n:
            exit_code = 2

    elif cfg.seq_type == "random":
        seqs = generate_random_set(poi, cfg, rng, run)
        run.t_generate = time.perf_counter() - t2
        seqs = select_diverse_sequences(seqs, min(cfg.n, len(seqs)), cfg.max_identity_between_decoys)

        if len(seqs) < cfg.n and not cfg.allow_partial:
            if cfg.fmt in ("csv", "tsv"):
                t4 = time.perf_counter()
                write_output_sequences([], poi, replace(cfg), seed_used)
                run.t_write = time.perf_counter() - t4
            print_manifest(run, poi_stats, cfg, seed_used, relax_stage_used)
            if cfg.verbosity >= 1:
                print(f"[pepctrl] Insufficient random sequences: requested={cfg.n} found={len(seqs)} (use --allow-partial).",
                      file=sys.stderr)
            return 2

        t4 = time.perf_counter()
        write_output_sequences(seqs, poi, cfg, seed_used)
        run.t_write = time.perf_counter() - t4
        if len(seqs) < cfg.n:
            exit_code = 2

    else:
        print(f"[pepctrl] ERROR: Unknown type {cfg.seq_type}", file=sys.stderr)
        return 3

    run.t_select = 0.0  # kept for compatibility in manifest timing fields
    print_manifest(run, poi_stats, cfg, seed_used, relax_stage_used)
    return exit_code


if __name__ == "__main__":
    raise SystemExit(main())

