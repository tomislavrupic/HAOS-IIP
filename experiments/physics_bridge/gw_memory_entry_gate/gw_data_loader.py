#!/usr/bin/env python3
"""Phase 64 strain-data loader for the GW memory entry gate.

The loader has two roles:

1. provide a deterministic GW150914-like surrogate so the phase can be run
   offline and reproduced without external downloads;
2. accept local GWOSC-style HDF5, CSV/TXT, or NPY strain files when the user
   wants to test real files.

This module does not download data and does not claim that any loaded strain
contains a real gravitational-memory observable. It only prepares strain-like
time series for the Phase 64 proxy telemetry.
"""

from __future__ import annotations

import csv
import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np


@dataclass(frozen=True)
class StrainRecord:
    """A single detector-like strain segment.

    `time` and `event_time` are seconds in the same coordinate. For generated
    surrogates they are relative seconds. For HDF5 files they may start as GPS
    seconds, but `prepare_record` returns an analysis window with time reset to
    zero and `event_time` shifted into that window.
    """

    time: np.ndarray
    strain: np.ndarray
    sample_rate: float
    source_label: str
    source_kind: str
    detector: str
    event_time: float
    metadata: dict[str, Any]


def generate_gw150914_like_surrogate(
    *,
    duration: float = 2.0,
    sample_rate: float = 1024.0,
    event_time: float = 1.15,
    seed: int = 6401,
    signal_to_noise: float = 14.0,
    detector: str = "SYNTH",
) -> StrainRecord:
    """Generate a deterministic compact-binary strain surrogate.

    The signal is a chirp with a short damped ringdown in colored Gaussian
    background noise. It is intentionally only GW150914-like: useful for testing
    the entry-gate telemetry without network access, but not a physical waveform
    model or NRSurrogate replacement.
    """

    rng = np.random.default_rng(seed)
    sample_count = int(round(duration * sample_rate))
    time = np.arange(sample_count, dtype=float) / float(sample_rate)
    centered = time - float(event_time)

    # A light colored-noise background, deterministic under the seed.
    white = rng.normal(size=sample_count)
    low = np.cumsum(rng.normal(scale=0.12, size=sample_count))
    low = low - float(np.mean(low))
    low = low / max(float(np.std(low)), 1.0e-12)
    noise = 0.82 * white + 0.18 * low
    noise = noise / max(float(np.std(noise)), 1.0e-12)

    # Chirp before merger. Frequency rises from tens to a few hundred Hz.
    pre = np.clip((centered + 0.42) / 0.42, 0.0, 1.0)
    frequency = 34.0 + 210.0 * pre**2.35
    phase = 2.0 * math.pi * np.cumsum(frequency) / float(sample_rate)
    envelope = np.exp(-((centered + 0.035) / 0.155) ** 2) * (0.20 + 1.10 * pre**1.7)
    chirp = envelope * np.sin(phase)

    # Damped post-merger ringdown.
    post = centered > 0.0
    ringdown = np.zeros_like(time)
    ringdown[post] = (
        np.exp(-centered[post] / 0.045)
        * np.sin(2.0 * math.pi * 235.0 * centered[post] + 0.7)
    )
    waveform = chirp + 0.42 * ringdown
    waveform = waveform / max(float(np.std(waveform)), 1.0e-12)

    # SNR here is a controlled analysis-scale ratio, not detector SNR.
    injected = noise + waveform * (signal_to_noise / 18.0)
    injected = injected - float(np.mean(injected))

    return StrainRecord(
        time=time,
        strain=injected.astype(float),
        sample_rate=float(sample_rate),
        source_label="gw150914_like_surrogate",
        source_kind="synthetic_surrogate",
        detector=detector,
        event_time=float(event_time),
        metadata={
            "duration": float(duration),
            "seed": int(seed),
            "signal_to_noise": float(signal_to_noise),
            "claim_label": "strain-derived proxy input; synthetic surrogate",
        },
    )


def load_strain_record(
    path: Path | None,
    *,
    sample_rate: float = 1024.0,
    duration: float = 2.0,
    event_time: float | None = None,
    seed: int = 6401,
    signal_to_noise: float = 14.0,
    detector: str = "SYNTH",
) -> StrainRecord:
    """Load a local strain record or return the deterministic surrogate."""

    if path is None:
        default_event = duration * 0.575 if event_time is None else float(event_time)
        return generate_gw150914_like_surrogate(
            duration=duration,
            sample_rate=sample_rate,
            event_time=default_event,
            seed=seed,
            signal_to_noise=signal_to_noise,
            detector=detector,
        )

    path = path.expanduser().resolve()
    suffix = path.suffix.lower()
    if suffix in {".hdf5", ".h5"}:
        return load_hdf5_strain(path, event_time=event_time)
    if suffix in {".csv", ".txt", ".dat"}:
        return load_text_strain(path, sample_rate=sample_rate, event_time=event_time, detector=detector)
    if suffix == ".npy":
        return load_npy_strain(path, sample_rate=sample_rate, event_time=event_time, detector=detector)
    raise ValueError(f"Unsupported strain file type: {path}")


def load_hdf5_strain(path: Path, *, event_time: float | None = None) -> StrainRecord:
    """Load a GWOSC-style HDF5 strain file.

    The common GWOSC layout contains `strain/Strain` plus metadata such as
    `meta/GPSstart` and `meta/Duration`. The function also searches for the
    longest one-dimensional numeric dataset as a fallback for local variants.
    """

    try:
        import h5py  # type: ignore[import-not-found]
    except ImportError as exc:  # pragma: no cover - depends on optional h5py
        raise ImportError("Loading HDF5 strain files requires h5py.") from exc

    with h5py.File(path, "r") as handle:
        dataset = _find_hdf5_strain_dataset(handle)
        strain = np.asarray(dataset, dtype=float)
        strain = np.ravel(strain)
        if strain.size < 16:
            raise ValueError(f"HDF5 strain dataset is too short: {path}")

        dt = _hdf5_float_attr(dataset, "Xspacing")
        if dt is None and "strain" in handle:
            dt = _hdf5_float_attr(handle["strain"], "Xspacing")
        sample_rate = 1.0 / dt if dt and dt > 0 else _infer_sample_rate_from_meta(handle, strain.size)

        gps_start = _read_hdf5_scalar(handle, "meta/GPSstart")
        duration = _read_hdf5_scalar(handle, "meta/Duration")
        if gps_start is None:
            gps_start = 0.0
        if duration is not None and duration > 0:
            sample_rate = float(strain.size) / float(duration)

        time = float(gps_start) + np.arange(strain.size, dtype=float) / float(sample_rate)
        detector = _read_hdf5_string(handle, "meta/Detector") or path.stem.split("-")[0]

    default_event = float(time[len(time) // 2]) if event_time is None else float(event_time)
    return StrainRecord(
        time=time,
        strain=strain,
        sample_rate=float(sample_rate),
        source_label=path.name,
        source_kind="gwosc_hdf5",
        detector=str(detector),
        event_time=default_event,
        metadata={
            "path": str(path),
            "claim_label": "strain-derived proxy input; local HDF5 file",
        },
    )


def load_text_strain(
    path: Path,
    *,
    sample_rate: float,
    event_time: float | None,
    detector: str,
) -> StrainRecord:
    """Load CSV/TXT strain.

    A two-column file is interpreted as `time,strain`; a one-column file is
    interpreted as strain sampled at `sample_rate`.
    """

    delimiter = "," if path.suffix.lower() == ".csv" else None
    try:
        data = np.loadtxt(path, delimiter=delimiter)
    except ValueError:
        data = np.loadtxt(path, delimiter=delimiter, skiprows=1)
    data = np.asarray(data, dtype=float)
    if data.ndim == 1:
        strain = data
        time = np.arange(strain.size, dtype=float) / float(sample_rate)
    elif data.shape[1] >= 2:
        time = data[:, 0]
        strain = data[:, 1]
        sample_rate = infer_sample_rate(time)
    else:
        raise ValueError(f"Cannot parse text strain file: {path}")

    default_event = float(time[len(time) // 2]) if event_time is None else float(event_time)
    return StrainRecord(
        time=time,
        strain=strain,
        sample_rate=float(sample_rate),
        source_label=path.name,
        source_kind="local_text",
        detector=detector,
        event_time=default_event,
        metadata={
            "path": str(path),
            "claim_label": "strain-derived proxy input; local text file",
        },
    )


def load_npy_strain(
    path: Path,
    *,
    sample_rate: float,
    event_time: float | None,
    detector: str,
) -> StrainRecord:
    """Load NPY strain.

    A shape `(N,)` array is strain only. A shape `(N, 2)` array is interpreted
    as `time,strain`.
    """

    data = np.asarray(np.load(path), dtype=float)
    if data.ndim == 1:
        strain = data
        time = np.arange(strain.size, dtype=float) / float(sample_rate)
    elif data.ndim == 2 and data.shape[1] >= 2:
        time = data[:, 0]
        strain = data[:, 1]
        sample_rate = infer_sample_rate(time)
    else:
        raise ValueError(f"Cannot parse NPY strain file: {path}")

    default_event = float(time[len(time) // 2]) if event_time is None else float(event_time)
    return StrainRecord(
        time=time,
        strain=strain,
        sample_rate=float(sample_rate),
        source_label=path.name,
        source_kind="local_npy",
        detector=detector,
        event_time=default_event,
        metadata={
            "path": str(path),
            "claim_label": "strain-derived proxy input; local NPY file",
        },
    )


def prepare_record(
    record: StrainRecord,
    *,
    duration: float,
    analysis_sample_rate: float | None = None,
    event_time: float | None = None,
) -> StrainRecord:
    """Crop around the event and optionally resample to an analysis rate."""

    raw_event = record.event_time if event_time is None else float(event_time)
    start = raw_event - duration / 2.0
    stop = raw_event + duration / 2.0
    record_start = float(record.time[0])
    record_stop = float(record.time[-1])
    mask = (record.time >= start) & (record.time < stop)
    if start >= record_start and stop <= record_stop and np.count_nonzero(mask) >= 16:
        time = record.time[mask]
        strain = record.strain[mask]
    else:
        time = record.time
        strain = record.strain
        duration = float(time[-1] - time[0]) + 1.0 / max(infer_sample_rate(time), 1.0e-12)

    if analysis_sample_rate is not None and analysis_sample_rate > 0:
        new_count = max(32, int(round(duration * analysis_sample_rate)))
        new_time_abs = float(time[0]) + np.arange(new_count, dtype=float) / float(analysis_sample_rate)
        strain = np.interp(new_time_abs, time, strain)
        time = new_time_abs
        sample_rate = float(analysis_sample_rate)
    else:
        sample_rate = infer_sample_rate(time)

    relative_time = time - float(time[0])
    relative_event = float(raw_event - float(time[0]))
    return StrainRecord(
        time=relative_time.astype(float),
        strain=np.asarray(strain, dtype=float),
        sample_rate=float(sample_rate),
        source_label=record.source_label,
        source_kind=record.source_kind,
        detector=record.detector,
        event_time=relative_event,
        metadata={
            **record.metadata,
            "analysis_duration": float(relative_time[-1] - relative_time[0]),
            "analysis_sample_rate": float(sample_rate),
            "analysis_event_time": float(relative_event),
        },
    )


def write_record_preview(path: Path, record: StrainRecord) -> None:
    """Write a compact CSV preview of the prepared strain segment."""

    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(["time_seconds", "strain"])
        for t_value, h_value in zip(record.time, record.strain, strict=True):
            writer.writerow([f"{float(t_value):.9f}", f"{float(h_value):.12e}"])


def write_record_metadata(path: Path, record: StrainRecord) -> None:
    payload = {
        "source_label": record.source_label,
        "source_kind": record.source_kind,
        "detector": record.detector,
        "sample_rate": float(record.sample_rate),
        "sample_count": int(record.strain.size),
        "duration": float(record.time[-1] - record.time[0]) if record.time.size else 0.0,
        "event_time": float(record.event_time),
        "metadata": record.metadata,
    }
    path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")


def infer_sample_rate(time: np.ndarray) -> float:
    if len(time) < 2:
        return 1.0
    spacing = np.diff(np.asarray(time, dtype=float))
    spacing = spacing[np.isfinite(spacing) & (spacing > 0)]
    if spacing.size == 0:
        return 1.0
    return float(1.0 / np.median(spacing))


def _find_hdf5_strain_dataset(handle: Any) -> Any:
    for candidate in ("strain/Strain", "Strain", "strain"):
        if candidate in handle:
            obj = handle[candidate]
            if hasattr(obj, "shape") and len(obj.shape) == 1:
                return obj

    datasets: list[Any] = []

    def visitor(_name: str, obj: Any) -> None:
        if not hasattr(obj, "shape") or not hasattr(obj, "dtype"):
            return
        if len(obj.shape) != 1:
            return
        if not np.issubdtype(obj.dtype, np.number):
            return
        datasets.append(obj)

    handle.visititems(visitor)
    if not datasets:
        raise ValueError("No one-dimensional numeric strain dataset found in HDF5 file.")
    return max(datasets, key=lambda ds: int(ds.shape[0]))


def _hdf5_float_attr(obj: Any, key: str) -> float | None:
    if key not in getattr(obj, "attrs", {}):
        return None
    value = obj.attrs[key]
    try:
        return float(np.asarray(value).ravel()[0])
    except (TypeError, ValueError, IndexError):
        return None


def _read_hdf5_scalar(handle: Any, path: str) -> float | None:
    if path not in handle:
        return None
    try:
        return float(np.asarray(handle[path]).ravel()[0])
    except (TypeError, ValueError, IndexError):
        return None


def _read_hdf5_string(handle: Any, path: str) -> str | None:
    if path not in handle:
        return None
    value = np.asarray(handle[path]).ravel()[0]
    if isinstance(value, bytes):
        return value.decode("utf-8", errors="replace")
    return str(value)


def _infer_sample_rate_from_meta(handle: Any, sample_count: int) -> float:
    duration = _read_hdf5_scalar(handle, "meta/Duration")
    if duration is not None and duration > 0:
        return float(sample_count) / float(duration)
    return 4096.0
