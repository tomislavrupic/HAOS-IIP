#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import json
import random
from dataclasses import dataclass
from datetime import datetime
from itertools import combinations
from pathlib import Path
from typing import Callable


SCRIPT_PATH = Path(__file__).resolve()
REPO_ROOT = SCRIPT_PATH.parents[2]
DATA_DIR = REPO_ROOT / "data"
EXPERIMENT_DIR = REPO_ROOT / "experiments" / "toy_sequence_controls"
RESULTS_TABLE_PATH = EXPERIMENT_DIR / "results_table.md"

SEQUENCE_LENGTH = 9
PROTECTED_SIZE = 6
SWAP_COUNT = 3
CORRIDOR_WIDTH = 4
BURST_WIDTH = 4
TOTAL_PROTECTED_PAIRS = PROTECTED_SIZE * (PROTECTED_SIZE - 1) // 2

# No canonical index rule for the Paper 25.2 toy subset was located in the workspace.
# This fixed candidate remains the default only because, under the original global-swap
# protocol, its baseline mean reproduces the externally reported ~0.688 toy score closely.
DEFAULT_REPORTED_SET = (0, 1, 2, 5, 6, 8)
DEFAULT_STORY_FAMILY = "anisotropic_local_disruptions"

ALL_PROTECTED_SETS = list(combinations(range(SEQUENCE_LENGTH), PROTECTED_SIZE))
CONTIGUOUS_PROTECTED_SETS = [
    tuple(range(start, start + PROTECTED_SIZE))
    for start in range(SEQUENCE_LENGTH - PROTECTED_SIZE + 1)
]


@dataclass(frozen=True)
class PerturbationFamily:
    name: str
    description: str
    generator: Callable[[random.Random], list[int]]


@dataclass
class RunningStats:
    count: int = 0
    total: float = 0.0
    total_sq: float = 0.0
    exact_recoveries: int = 0

    def add(self, score: float) -> None:
        self.count += 1
        self.total += score
        self.total_sq += score * score
        if score >= 1.0 - 1.0e-12:
            self.exact_recoveries += 1

    def snapshot(self) -> dict[str, float]:
        mean = self.total / self.count if self.count else 0.0
        variance = (self.total_sq / self.count - mean * mean) if self.count else 0.0
        return {
            "mean_protected_pairwise_coherence": mean,
            "variance": max(variance, 0.0),
            "exact_recoverability_rate": self.exact_recoveries / self.count if self.count else 0.0,
            "trial_count": float(self.count),
        }


def swap_positions(sequence: list[int], left: int, right: int) -> None:
    sequence[left], sequence[right] = sequence[right], sequence[left]


def random_swaps(rng: random.Random) -> list[int]:
    sequence = list(range(SEQUENCE_LENGTH))
    for _ in range(SWAP_COUNT):
        left, right = rng.sample(range(SEQUENCE_LENGTH), 2)
        swap_positions(sequence, left, right)
    return sequence


def adjacent_swaps(rng: random.Random) -> list[int]:
    sequence = list(range(SEQUENCE_LENGTH))
    for _ in range(SWAP_COUNT):
        left = rng.randrange(SEQUENCE_LENGTH - 1)
        swap_positions(sequence, left, left + 1)
    return sequence


def corridor_limited_swaps(rng: random.Random) -> list[int]:
    sequence = list(range(SEQUENCE_LENGTH))
    start = rng.randrange(SEQUENCE_LENGTH - CORRIDOR_WIDTH + 1)
    corridor = list(range(start, start + CORRIDOR_WIDTH))
    for _ in range(SWAP_COUNT):
        left, right = rng.sample(corridor, 2)
        swap_positions(sequence, left, right)
    return sequence


def burst_errors(rng: random.Random) -> list[int]:
    sequence = list(range(SEQUENCE_LENGTH))
    start = rng.randrange(SEQUENCE_LENGTH - BURST_WIDTH + 1)
    burst_indices = list(range(start, start + BURST_WIDTH))
    burst_values = [sequence[index] for index in burst_indices]
    while True:
        shuffled = burst_values[:]
        rng.shuffle(shuffled)
        if shuffled != burst_values:
            break
    for index, value in zip(burst_indices, shuffled):
        sequence[index] = value
    return sequence


def anisotropic_local_disruptions(rng: random.Random) -> list[int]:
    sequence = list(range(SEQUENCE_LENGTH))
    anchor = rng.randrange(1, SEQUENCE_LENGTH - 1)
    direction = rng.choice((-1, 1))
    position = anchor
    for _ in range(SWAP_COUNT):
        nxt = position + direction
        if nxt < 0 or nxt >= SEQUENCE_LENGTH:
            direction *= -1
            nxt = position + direction
        swap_positions(sequence, position, nxt)
        position = nxt
    return sequence


PERTURBATION_FAMILIES = {
    family.name: family
    for family in [
        PerturbationFamily(
            name="random_swaps",
            description="three unconstrained swaps anywhere in the sequence",
            generator=random_swaps,
        ),
        PerturbationFamily(
            name="adjacent_swaps",
            description="three adjacent swaps at arbitrary positions",
            generator=adjacent_swaps,
        ),
        PerturbationFamily(
            name="corridor_limited_swaps",
            description=f"three swaps restricted to one random contiguous corridor of width {CORRIDOR_WIDTH}",
            generator=corridor_limited_swaps,
        ),
        PerturbationFamily(
            name="burst_errors",
            description=f"one random contiguous burst of width {BURST_WIDTH}, internally permuted",
            generator=burst_errors,
        ),
        PerturbationFamily(
            name="anisotropic_local_disruptions",
            description="one directional local displacement built from three consecutive adjacent swaps",
            generator=anisotropic_local_disruptions,
        ),
    ]
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run the multi-family 9-event toy recoverability control suite."
    )
    parser.add_argument("--trials", type=int, default=50000)
    parser.add_argument("--ranking-trials", type=int, default=20000)
    parser.add_argument("--seed", type=int, default=20260316)
    parser.add_argument(
        "--reported-set",
        type=str,
        default=",".join(str(item) for item in DEFAULT_REPORTED_SET),
        help="Comma-separated protected subset used as the fixed reported candidate.",
    )
    parser.add_argument(
        "--haos-set",
        type=str,
        default="",
        help="Optional comma-separated protected subset for a separately supplied HAOS selection rule.",
    )
    parser.add_argument(
        "--story-family",
        type=str,
        default=DEFAULT_STORY_FAMILY,
        help="Perturbation family treated as the default HAOS-story proxy.",
    )
    return parser.parse_args()


def parse_protected_set(raw: str, label: str) -> tuple[int, ...]:
    values = tuple(int(part.strip()) for part in raw.split(",") if part.strip())
    if len(values) != PROTECTED_SIZE:
        raise ValueError(f"{label} must contain exactly {PROTECTED_SIZE} positions")
    if len(set(values)) != PROTECTED_SIZE:
        raise ValueError(f"{label} must not repeat positions")
    if any(value < 0 or value >= SEQUENCE_LENGTH for value in values):
        raise ValueError(f"{label} entries must lie in [0, {SEQUENCE_LENGTH - 1}]")
    return tuple(sorted(values))


def protected_subsequence(sequence: list[int], protected_lookup: frozenset[int]) -> list[int]:
    return [value for value in sequence if value in protected_lookup]


def protected_pairwise_coherence(ordering: list[int]) -> float:
    preserved = 0
    for left_index, left_value in enumerate(ordering):
        for right_value in ordering[left_index + 1 :]:
            if left_value < right_value:
                preserved += 1
    return preserved / TOTAL_PROTECTED_PAIRS


def full_repair(ordering: list[int]) -> list[int]:
    return sorted(ordering)


def local_repair(ordering: list[int]) -> list[int]:
    repaired = list(ordering)
    for index in range(len(repaired) - 1):
        if repaired[index] > repaired[index + 1]:
            repaired[index], repaired[index + 1] = repaired[index + 1], repaired[index]
            break
    return repaired


def evaluate_conditions_for_family(
    trials: int,
    seed: int,
    family: PerturbationFamily,
    reported_set: tuple[int, ...],
    haos_set: tuple[int, ...] | None,
) -> dict[str, RunningStats]:
    repairs: dict[str, Callable[[list[int]], list[int]]] = {
        "baseline": lambda ordering: list(ordering),
        "full_repair": full_repair,
        "local_repair": local_repair,
    }
    stats = {
        f"{prefix}_{repair_name}": RunningStats()
        for prefix in ("reported", "random", "contiguous")
        for repair_name in repairs
    }
    if haos_set is not None:
        for repair_name in repairs:
            stats[f"haos_{repair_name}"] = RunningStats()

    excluded_random_sets = {reported_set}
    if haos_set is not None:
        excluded_random_sets.add(haos_set)
    random_pool = [
        protected_set
        for protected_set in ALL_PROTECTED_SETS
        if protected_set not in excluded_random_sets and protected_set not in CONTIGUOUS_PROTECTED_SETS
    ]
    contiguous_pool = list(CONTIGUOUS_PROTECTED_SETS)
    reported_lookup = frozenset(reported_set)
    haos_lookup = frozenset(haos_set) if haos_set is not None else None

    rng = random.Random(seed)
    for _ in range(trials):
        sequence = family.generator(rng)
        random_set = random_pool[rng.randrange(len(random_pool))]
        contiguous_set = contiguous_pool[rng.randrange(len(contiguous_pool))]

        lookups = {
            "reported": reported_lookup,
            "random": frozenset(random_set),
            "contiguous": frozenset(contiguous_set),
        }
        if haos_lookup is not None:
            lookups["haos"] = haos_lookup

        for prefix, lookup in lookups.items():
            ordering = protected_subsequence(sequence, lookup)
            for repair_name, repair_fn in repairs.items():
                score = protected_pairwise_coherence(repair_fn(ordering))
                stats[f"{prefix}_{repair_name}"].add(score)

    return stats


def exhaustive_rankings_for_family(
    trials: int,
    seed: int,
    family: PerturbationFamily,
    candidate_sets: dict[str, tuple[int, ...]],
) -> dict[str, dict[str, float | int | str]]:
    baseline_totals = {protected_set: 0.0 for protected_set in ALL_PROTECTED_SETS}
    local_totals = {protected_set: 0.0 for protected_set in ALL_PROTECTED_SETS}
    protected_lookups = {
        protected_set: frozenset(protected_set) for protected_set in ALL_PROTECTED_SETS
    }

    rng = random.Random(seed + 1)
    for _ in range(trials):
        sequence = family.generator(rng)
        for protected_set in ALL_PROTECTED_SETS:
            ordering = protected_subsequence(sequence, protected_lookups[protected_set])
            baseline_totals[protected_set] += protected_pairwise_coherence(ordering)
            local_totals[protected_set] += protected_pairwise_coherence(local_repair(ordering))

    baseline_means = {
        protected_set: total / trials for protected_set, total in baseline_totals.items()
    }
    local_means = {
        protected_set: total / trials for protected_set, total in local_totals.items()
    }
    baseline_ranked = sorted(
        baseline_means.items(),
        key=lambda item: (item[1], item[0]),
        reverse=True,
    )
    local_ranked = sorted(
        local_means.items(),
        key=lambda item: (item[1], item[0]),
        reverse=True,
    )

    results: dict[str, dict[str, float | int | str]] = {}
    for label, protected_set in candidate_sets.items():
        baseline_rank = next(
            index + 1 for index, item in enumerate(baseline_ranked) if item[0] == protected_set
        )
        local_rank = next(
            index + 1 for index, item in enumerate(local_ranked) if item[0] == protected_set
        )
        results[label] = {
            "protected_set": str(protected_set),
            "baseline_mean": baseline_means[protected_set],
            "baseline_rank_of_84": baseline_rank,
            "local_mean": local_means[protected_set],
            "local_rank_of_84": local_rank,
            "baseline_top_set": str(baseline_ranked[0][0]),
            "baseline_top_mean": baseline_ranked[0][1],
            "local_top_set": str(local_ranked[0][0]),
            "local_top_mean": local_ranked[0][1],
        }
    return results


def long_form_rows_for_family(
    family_name: str,
    family_description: str,
    stats: dict[str, RunningStats],
    reported_set: tuple[int, ...],
    haos_set: tuple[int, ...] | None,
) -> list[dict[str, object]]:
    snapshots = {name: value.snapshot() for name, value in stats.items()}
    rows: list[dict[str, object]] = []
    subset_families = [
        ("reported", reported_set, "fixed candidate subset"),
        ("random", "", "uniform random non-contiguous control"),
        ("contiguous", "", "uniform contiguous 6-position control"),
    ]
    if haos_set is not None:
        subset_families.append(("haos", haos_set, "externally supplied HAOS subset"))

    for subset_family, protected_set, note in subset_families:
        for repair_mode in ("baseline", "full_repair", "local_repair"):
            snapshot = snapshots[f"{subset_family}_{repair_mode}"]
            row: dict[str, object] = {
                "perturbation_family": family_name,
                "perturbation_description": family_description,
                "subset_family": subset_family,
                "repair_mode": repair_mode,
                "protected_set": protected_set,
                "mean_protected_pairwise_coherence": snapshot["mean_protected_pairwise_coherence"],
                "variance": snapshot["variance"],
                "exact_recoverability_rate": snapshot["exact_recoverability_rate"],
                "trial_count": int(snapshot["trial_count"]),
                "delta_vs_random_mean": "",
                "delta_vs_contiguous_mean": "",
                "beats_random_mean": "",
                "beats_contiguous_mean": "",
                "notes": note,
            }
            if subset_family in {"reported", "haos"}:
                random_mean = snapshots[f"random_{repair_mode}"]["mean_protected_pairwise_coherence"]
                contiguous_mean = snapshots[f"contiguous_{repair_mode}"]["mean_protected_pairwise_coherence"]
                delta_vs_random = snapshot["mean_protected_pairwise_coherence"] - random_mean
                delta_vs_contiguous = (
                    snapshot["mean_protected_pairwise_coherence"] - contiguous_mean
                )
                row["delta_vs_random_mean"] = delta_vs_random
                row["delta_vs_contiguous_mean"] = delta_vs_contiguous
                row["beats_random_mean"] = delta_vs_random > 0.0
                row["beats_contiguous_mean"] = delta_vs_contiguous > 0.0
            rows.append(row)
    return rows


def family_summary_row(
    family_name: str,
    family_description: str,
    story_family: str,
    stats: dict[str, RunningStats],
    rankings: dict[str, dict[str, float | int | str]],
) -> dict[str, object]:
    snapshots = {name: value.snapshot() for name, value in stats.items()}
    ranking = rankings["reported"]
    return {
        "perturbation_family": family_name,
        "description": family_description,
        "story_family": family_name == story_family,
        "reported_baseline_mean": snapshots["reported_baseline"]["mean_protected_pairwise_coherence"],
        "reported_baseline_rank_of_84": int(ranking["baseline_rank_of_84"]),
        "reported_local_mean": snapshots["reported_local_repair"]["mean_protected_pairwise_coherence"],
        "reported_local_rank_of_84": int(ranking["local_rank_of_84"]),
        "baseline_delta_vs_random_mean": snapshots["reported_baseline"]["mean_protected_pairwise_coherence"]
        - snapshots["random_baseline"]["mean_protected_pairwise_coherence"],
        "baseline_delta_vs_contiguous_mean": snapshots["reported_baseline"]["mean_protected_pairwise_coherence"]
        - snapshots["contiguous_baseline"]["mean_protected_pairwise_coherence"],
        "local_delta_vs_random_mean": snapshots["reported_local_repair"]["mean_protected_pairwise_coherence"]
        - snapshots["random_local_repair"]["mean_protected_pairwise_coherence"],
        "local_delta_vs_contiguous_mean": snapshots["reported_local_repair"]["mean_protected_pairwise_coherence"]
        - snapshots["contiguous_local_repair"]["mean_protected_pairwise_coherence"],
        "baseline_beats_random": snapshots["reported_baseline"]["mean_protected_pairwise_coherence"]
        > snapshots["random_baseline"]["mean_protected_pairwise_coherence"],
        "baseline_beats_contiguous": snapshots["reported_baseline"]["mean_protected_pairwise_coherence"]
        > snapshots["contiguous_baseline"]["mean_protected_pairwise_coherence"],
        "local_beats_random": snapshots["reported_local_repair"]["mean_protected_pairwise_coherence"]
        > snapshots["random_local_repair"]["mean_protected_pairwise_coherence"],
        "local_beats_contiguous": snapshots["reported_local_repair"]["mean_protected_pairwise_coherence"]
        > snapshots["contiguous_local_repair"]["mean_protected_pairwise_coherence"],
        "baseline_top_set": ranking["baseline_top_set"],
        "local_top_set": ranking["local_top_set"],
    }


def decide_question(summary_rows: list[dict[str, object]], story_family: str) -> dict[str, object]:
    story_row = next(row for row in summary_rows if row["perturbation_family"] == story_family)
    other_rows = [row for row in summary_rows if row["perturbation_family"] != story_family]
    best_baseline_rank = min(int(row["reported_baseline_rank_of_84"]) for row in summary_rows)
    best_local_rank = min(int(row["reported_local_rank_of_84"]) for row in summary_rows)
    answer = (
        int(story_row["reported_baseline_rank_of_84"]) <= best_baseline_rank
        and int(story_row["reported_local_rank_of_84"]) <= best_local_rank
        and bool(story_row["baseline_beats_random"])
        and bool(story_row["baseline_beats_contiguous"])
        and bool(story_row["local_beats_random"])
        and bool(story_row["local_beats_contiguous"])
    )
    return {
        "question": (
            f"Does the fixed candidate become consistently more special under `{story_family}`, "
            "the chosen HAOS-story perturbation proxy?"
        ),
        "answer": "yes" if answer else "no",
        "decision_rule": (
            "Answer yes only if the story family gives the fixed candidate the best or tied-best "
            "exhaustive rank among all 84 subsets in both baseline and bounded local repair, "
            "while also beating matched random and contiguous controls in both modes."
        ),
        "story_family_baseline_rank_of_84": int(story_row["reported_baseline_rank_of_84"]),
        "story_family_local_rank_of_84": int(story_row["reported_local_rank_of_84"]),
        "next_best_baseline_rank_of_84": min(
            int(row["reported_baseline_rank_of_84"]) for row in other_rows
        ) if other_rows else int(story_row["reported_baseline_rank_of_84"]),
        "next_best_local_rank_of_84": min(
            int(row["reported_local_rank_of_84"]) for row in other_rows
        ) if other_rows else int(story_row["reported_local_rank_of_84"]),
    }


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    fieldnames = [
        "perturbation_family",
        "perturbation_description",
        "subset_family",
        "repair_mode",
        "protected_set",
        "mean_protected_pairwise_coherence",
        "variance",
        "exact_recoverability_rate",
        "trial_count",
        "delta_vs_random_mean",
        "delta_vs_contiguous_mean",
        "beats_random_mean",
        "beats_contiguous_mean",
        "notes",
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def format_float(value: object) -> str:
    if value == "":
        return ""
    return f"{float(value):.6f}"


def write_results_table(
    path: Path,
    json_rel: str,
    csv_rel: str,
    summary_rows: list[dict[str, object]],
    question_result: dict[str, object],
    story_family: str,
    reported_set: tuple[int, ...],
    selection_rule_available: bool,
) -> None:
    lines = [
        "# Toy Sequence Perturbation Family Comparison",
        "",
        f"- Timestamped JSON: `{json_rel}`",
        f"- Timestamped CSV: `{csv_rel}`",
        f"- Sequence length: `{SEQUENCE_LENGTH}`",
        f"- Protected positions: `{PROTECTED_SIZE}`",
        f"- Fixed candidate subset: `{reported_set}`",
        f"- Explicit HAOS protected-set selection rule available in workspace: `{selection_rule_available}`",
        f"- Default HAOS-story perturbation proxy: `{story_family}`",
        "",
        "## Family definitions",
    ]
    for family_name in PERTURBATION_FAMILIES:
        marker = " (story proxy)" if family_name == story_family else ""
        lines.append(
            f"- `{family_name}`{marker}: {PERTURBATION_FAMILIES[family_name].description}"
        )

    lines.extend(
        [
            "",
            "## Question",
            "",
            f"- {question_result['question']}",
            f"- Answer: `{question_result['answer']}`",
            f"- Rule: {question_result['decision_rule']}",
            "",
            "## Family comparison",
            "",
            "| family | story proxy | baseline rank / 84 | local rank / 84 | baseline delta vs random | baseline delta vs contiguous | local delta vs random | local delta vs contiguous |",
            "| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |",
        ]
    )
    ordered_rows = sorted(
        summary_rows,
        key=lambda row: (
            0 if row["perturbation_family"] == story_family else 1,
            int(row["reported_baseline_rank_of_84"]),
            int(row["reported_local_rank_of_84"]),
            str(row["perturbation_family"]),
        ),
    )
    for row in ordered_rows:
        lines.append(
            "| "
            + f"{row['perturbation_family']} | "
            + f"{row['story_family']} | "
            + f"{int(row['reported_baseline_rank_of_84'])} | "
            + f"{int(row['reported_local_rank_of_84'])} | "
            + f"{format_float(row['baseline_delta_vs_random_mean'])} | "
            + f"{format_float(row['baseline_delta_vs_contiguous_mean'])} | "
            + f"{format_float(row['local_delta_vs_random_mean'])} | "
            + f"{format_float(row['local_delta_vs_contiguous_mean'])} |"
        )

    lines.extend(
        [
            "",
            "## Note",
            "",
            "- Full repair is omitted from the comparison table because it restores exact protected order by construction and therefore does not distinguish perturbation families.",
            "- The exhaustive ranking remains the main nontrivial test: each family ranks the fixed candidate against all 84 possible 6-position subsets.",
        ]
    )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    args = parse_args()
    if args.story_family not in PERTURBATION_FAMILIES:
        raise ValueError(
            f"story-family must be one of {sorted(PERTURBATION_FAMILIES)}"
        )

    DATA_DIR.mkdir(parents=True, exist_ok=True)
    EXPERIMENT_DIR.mkdir(parents=True, exist_ok=True)

    reported_set = parse_protected_set(args.reported_set, "reported-set")
    haos_set = parse_protected_set(args.haos_set, "haos-set") if args.haos_set else None
    selection_rule_available = haos_set is not None

    long_form_rows: list[dict[str, object]] = []
    family_rankings: dict[str, dict[str, dict[str, float | int | str]]] = {}
    family_summaries: list[dict[str, object]] = []
    ranking_candidates = {"reported": reported_set}
    if haos_set is not None:
        ranking_candidates["haos"] = haos_set

    for family_name, family in PERTURBATION_FAMILIES.items():
        family_seed = args.seed + 1000 * list(PERTURBATION_FAMILIES).index(family_name)
        stats = evaluate_conditions_for_family(
            trials=args.trials,
            seed=family_seed,
            family=family,
            reported_set=reported_set,
            haos_set=haos_set,
        )
        rankings = exhaustive_rankings_for_family(
            trials=args.ranking_trials,
            seed=family_seed,
            family=family,
            candidate_sets=ranking_candidates,
        )
        family_rankings[family_name] = rankings
        long_form_rows.extend(
            long_form_rows_for_family(
                family_name=family_name,
                family_description=family.description,
                stats=stats,
                reported_set=reported_set,
                haos_set=haos_set,
            )
        )
        family_summaries.append(
            family_summary_row(
                family_name=family_name,
                family_description=family.description,
                story_family=args.story_family,
                stats=stats,
                rankings=rankings,
            )
        )

    question_result = decide_question(
        summary_rows=family_summaries,
        story_family=args.story_family,
    )

    stamp = datetime.utcnow().strftime("%Y%m%d_%H%M%S")
    json_path = DATA_DIR / f"{stamp}_toy_sequence_recoverability_controls.json"
    csv_path = DATA_DIR / f"{stamp}_toy_sequence_recoverability_controls.csv"

    summary = {
        "experiment": "toy_sequence_recoverability_controls",
        "sequence_length": SEQUENCE_LENGTH,
        "protected_size": PROTECTED_SIZE,
        "swap_count": SWAP_COUNT,
        "corridor_width": CORRIDOR_WIDTH,
        "burst_width": BURST_WIDTH,
        "trials_per_family": args.trials,
        "ranking_trials_per_family": args.ranking_trials,
        "seed": args.seed,
        "selection_rule_available": selection_rule_available,
        "reported_set": reported_set,
        "haos_set": haos_set,
        "story_family": args.story_family,
        "story_family_assumption": (
            "The default story-family proxy is anisotropic_local_disruptions because it is the most structured, "
            "directional, and local perturbation family in the requested set. No canonical HAOS perturbation family "
            "definition was located in the workspace."
        ),
        "perturbation_families": {
            family_name: family.description for family_name, family in PERTURBATION_FAMILIES.items()
        },
        "question_result": question_result,
        "family_summaries": family_summaries,
        "rows": long_form_rows,
        "rankings": family_rankings,
        "reported_set_note": (
            "Default fixed candidate chosen to reproduce the externally reported baseline under the original "
            "global-swap toy protocol. No explicit HAOS selection rule for the subset was located in the workspace."
        ),
        "local_repair_note": (
            "Local repair is a single left-to-right adjacent inversion swap on the protected subsequence. "
            "It is bounded and does not guarantee full restoration."
        ),
        "full_repair_note": (
            "Full repair sorts the protected subsequence into reference order. "
            "It guarantees exact protected-order recovery by construction."
        ),
    }
    json_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    write_csv(csv_path, long_form_rows)
    write_results_table(
        path=RESULTS_TABLE_PATH,
        json_rel=json_path.relative_to(REPO_ROOT).as_posix(),
        csv_rel=csv_path.relative_to(REPO_ROOT).as_posix(),
        summary_rows=family_summaries,
        question_result=question_result,
        story_family=args.story_family,
        reported_set=reported_set,
        selection_rule_available=selection_rule_available,
    )

    print(f"Wrote {json_path}")
    print(f"Wrote {csv_path}")
    print(f"Wrote {RESULTS_TABLE_PATH}")


if __name__ == "__main__":
    main()
