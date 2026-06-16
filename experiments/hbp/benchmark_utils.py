from __future__ import annotations

import hashlib
from pathlib import Path
from typing import Any, Iterable

import numpy as np


def stable_hash(payload: Any, prefix: str) -> str:
    import json

    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":"), default=str).encode("utf-8")
    return prefix + hashlib.sha256(encoded).hexdigest()[:24]


def summarize_array(values: np.ndarray | Iterable[float], prefix: str) -> dict[str, float]:
    arr = np.asarray(values, dtype=float)
    arr = np.nan_to_num(arr, nan=0.0, posinf=0.0, neginf=0.0).reshape(-1)
    if arr.size == 0:
        return {
            f"{prefix}_mean": 0.0,
            f"{prefix}_std": 0.0,
            f"{prefix}_min": 0.0,
            f"{prefix}_max": 0.0,
            f"{prefix}_l1": 0.0,
            f"{prefix}_l2": 0.0,
            f"{prefix}_zero_fraction": 1.0,
        }
    return {
        f"{prefix}_mean": float(np.mean(arr)),
        f"{prefix}_std": float(np.std(arr)),
        f"{prefix}_min": float(np.min(arr)),
        f"{prefix}_max": float(np.max(arr)),
        f"{prefix}_l1": float(np.sum(np.abs(arr))),
        f"{prefix}_l2": float(np.linalg.norm(arr)),
        f"{prefix}_zero_fraction": float(np.mean(np.isclose(arr, 0.0))),
    }


def standardize(train: np.ndarray, values: np.ndarray) -> tuple[np.ndarray, dict[str, list[float]]]:
    mean = np.mean(train, axis=0)
    std = np.std(train, axis=0)
    std = np.where(std < 1.0e-12, 1.0, std)
    return (values - mean) / std, {"mean": mean.tolist(), "std": std.tolist()}


def fit_ridge(train_x: np.ndarray, train_y: np.ndarray, alpha: float = 1.0e-4) -> np.ndarray:
    x = np.asarray(train_x, dtype=float)
    y = np.asarray(train_y, dtype=float).reshape(-1, 1)
    x_aug = np.concatenate([np.ones((x.shape[0], 1), dtype=float), x], axis=1)
    gram = x_aug.T @ x_aug
    gram += alpha * np.eye(gram.shape[0], dtype=float)
    gram[0, 0] -= alpha
    beta = np.linalg.solve(gram, x_aug.T @ y)
    return beta.reshape(-1)


def predict_ridge(beta: np.ndarray, x: np.ndarray) -> np.ndarray:
    x_aug = np.concatenate([np.ones((x.shape[0], 1), dtype=float), x], axis=1)
    return (x_aug @ beta.reshape(-1, 1)).reshape(-1)


def rankdata(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    order = np.argsort(values, kind="mergesort")
    ranks = np.empty_like(order, dtype=float)
    ranks[order] = np.arange(len(values), dtype=float)
    _, counts = np.unique(values[order], return_counts=True)
    start = 0
    for count in counts:
        end = start + count
        avg = float(np.mean(np.arange(start, end, dtype=float)))
        ranks[order[start:end]] = avg
        start = end
    return ranks


def spearman(x: np.ndarray, y: np.ndarray) -> float:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if x.size < 2 or y.size < 2:
        return 0.0
    rx = rankdata(x)
    ry = rankdata(y)
    if np.std(rx) < 1.0e-12 or np.std(ry) < 1.0e-12:
        return 0.0
    return float(np.corrcoef(rx, ry)[0, 1])


def mae(x: np.ndarray, y: np.ndarray) -> float:
    return float(np.mean(np.abs(np.asarray(x) - np.asarray(y))))


def rmse(x: np.ndarray, y: np.ndarray) -> float:
    return float(np.sqrt(np.mean((np.asarray(x) - np.asarray(y)) ** 2)))


def top_k_precision(pred: np.ndarray, truth: np.ndarray, k: int = 100) -> float:
    pred = np.asarray(pred, dtype=float)
    truth = np.asarray(truth, dtype=float)
    k = max(1, min(int(k), pred.size))
    pred_top = set(np.argsort(pred)[-k:])
    truth_top = set(np.argsort(truth)[-k:])
    return float(len(pred_top & truth_top) / k)


def build_feature_vector(arrays: Iterable[np.ndarray | Iterable[float]], prefixes: Iterable[str]) -> np.ndarray:
    features: list[float] = []
    for arr, prefix in zip(arrays, prefixes, strict=True):
        summary = summarize_array(arr, prefix)
        features.extend(summary.values())
    return np.asarray(features, dtype=float)

