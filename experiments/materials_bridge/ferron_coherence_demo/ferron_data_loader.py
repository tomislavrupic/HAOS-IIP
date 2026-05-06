#!/usr/bin/env python3
"""Data loader for the ferron coherence materials bridge.

This module is intentionally sidecar-only. It tries to fetch the public
Figshare source-data record for the Nature Materials ferron paper, stores raw
downloads untouched, records provenance and hashes, and exposes conservative
parsers for common numeric data formats.

If public data cannot be downloaded, the module can build a clearly labeled
synthetic smoke-test dataset for software verification only.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
import re
import tarfile
import urllib.error
import urllib.request
import zipfile
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np


ROOT = Path(__file__).resolve().parent
OUTPUTS_DIR = ROOT / "outputs"
RAW_DIR = OUTPUTS_DIR / "raw"
EXTRACTED_DIR = OUTPUTS_DIR / "derived_extracted"
SMOKE_DIR = OUTPUTS_DIR / "smoke_test_dataset"
MANIFEST_PATH = ROOT / "data_manifest.json"

PRIMARY_PAPER_DOI = "10.1038/s41563-026-02597-4"
PRIMARY_DATA_DOI = "10.6084/m9.figshare.31895293"
PRIMARY_DATA_URL = "https://doi.org/10.6084/m9.figshare.31895293"
FIGSHARE_ARTICLE_ID = "31895293"
FIGSHARE_API_URL = f"https://api.figshare.com/v2/articles/{FIGSHARE_ARTICLE_ID}"
FIGSHARE_SEARCH_URL = "https://api.figshare.com/v2/articles/search"

SUPPORTED_EXTENSIONS = {
    ".csv",
    ".txt",
    ".tsv",
    ".dat",
    ".json",
    ".npy",
    ".npz",
    ".xlsx",
}
ARCHIVE_EXTENSIONS = {".zip", ".tar", ".gz", ".tgz", ".bz2", ".xz"}


def utc_now_iso() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def compute_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _request_json(url: str, *, data: dict[str, Any] | None = None) -> dict[str, Any] | list[Any]:
    payload = None
    headers = {
        "Accept": "application/json",
        "User-Agent": "haos-iip-ferron-materials-bridge/1.0",
    }
    if data is not None:
        payload = json.dumps(data).encode("utf-8")
        headers["Content-Type"] = "application/json"
    request = urllib.request.Request(url, data=payload, headers=headers)
    with urllib.request.urlopen(request, timeout=45) as response:
        return json.loads(response.read().decode("utf-8"))


def _safe_name(name: str, fallback: str) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9._ -]+", "_", name).strip()
    cleaned = cleaned.replace("/", "_").replace("\\", "_")
    return cleaned or fallback


def _download_file(url: str, destination: Path) -> None:
    request = urllib.request.Request(
        url,
        headers={"User-Agent": "haos-iip-ferron-materials-bridge/1.0"},
    )
    with urllib.request.urlopen(request, timeout=120) as response:
        with destination.open("wb") as handle:
            for chunk in iter(lambda: response.read(1024 * 1024), b""):
                handle.write(chunk)


def _resolve_figshare_metadata() -> tuple[dict[str, Any] | None, list[dict[str, Any]]]:
    attempts: list[dict[str, Any]] = []
    try:
        metadata = _request_json(FIGSHARE_API_URL)
        if isinstance(metadata, dict):
            attempts.append({"url": FIGSHARE_API_URL, "status": "ok"})
            return metadata, attempts
    except (urllib.error.URLError, urllib.error.HTTPError, TimeoutError, json.JSONDecodeError) as exc:
        attempts.append({"url": FIGSHARE_API_URL, "status": "failed", "error": str(exc)})

    try:
        search = _request_json(
            FIGSHARE_SEARCH_URL,
            data={"search_for": PRIMARY_DATA_DOI, "page": 1, "page_size": 5},
        )
        attempts.append({"url": FIGSHARE_SEARCH_URL, "status": "ok"})
        if isinstance(search, list) and search:
            article_id = str(search[0].get("id", "")).strip()
            if article_id:
                article_url = f"https://api.figshare.com/v2/articles/{article_id}"
                metadata = _request_json(article_url)
                if isinstance(metadata, dict):
                    attempts.append({"url": article_url, "status": "ok"})
                    return metadata, attempts
    except (urllib.error.URLError, urllib.error.HTTPError, TimeoutError, json.JSONDecodeError) as exc:
        attempts.append({"url": FIGSHARE_SEARCH_URL, "status": "failed", "error": str(exc)})

    return None, attempts


def download_or_locate_dataset(root: Path = ROOT, *, allow_download: bool = True) -> dict[str, Any]:
    """Download the Figshare dataset or locate existing raw files.

    Raw files are stored under ``outputs/raw`` and never modified. The returned
    dictionary is also written through ``write_data_manifest``.
    """

    raw_dir = root / "outputs" / "raw"
    raw_dir.mkdir(parents=True, exist_ok=True)
    manifest_path = root / "data_manifest.json"
    existing = sorted(path for path in raw_dir.iterdir() if path.is_file())
    downloaded_at = utc_now_iso()

    manifest: dict[str, Any] = {
        "experiment": "Materials Bridge Line A - Ferron Coherence Recoverability Probe",
        "primary_paper_doi": PRIMARY_PAPER_DOI,
        "source_doi": PRIMARY_DATA_DOI,
        "source_url": PRIMARY_DATA_URL,
        "figshare_api_url": FIGSHARE_API_URL,
        "download_timestamp_utc": downloaded_at,
        "status": "NOT_ATTEMPTED",
        "files": [],
        "download_attempts": [],
    }

    if existing:
        prior_records = _prior_manifest_records(manifest_path)
        manifest["status"] = "RAW_FILES_ALREADY_PRESENT"
        for path in existing:
            prior = prior_records.get(path.name, {})
            source_url = str(prior.get("source_url") or PRIMARY_DATA_URL)
            file_timestamp = str(prior.get("download_timestamp_utc") or downloaded_at)
            record = _manifest_file_record(path, source_url, file_timestamp)
            for key in (
                "figshare_file_id",
                "figshare_reported_size_bytes",
                "figshare_reported_md5",
            ):
                if key in prior:
                    record[key] = prior[key]
            manifest["files"].append(record)
        write_data_manifest(manifest, manifest_path)
        return manifest

    if not allow_download:
        manifest["status"] = "DOWNLOAD_SKIPPED"
        write_data_manifest(manifest, manifest_path)
        return manifest

    metadata, attempts = _resolve_figshare_metadata()
    manifest["download_attempts"].extend(attempts)
    if not metadata:
        manifest["status"] = "DOWNLOAD_FAILED"
        write_data_manifest(manifest, manifest_path)
        return manifest

    manifest["figshare_title"] = metadata.get("title")
    manifest["figshare_id"] = metadata.get("id")
    manifest["figshare_published_date"] = metadata.get("published_date")
    files = metadata.get("files") or []
    if not isinstance(files, list) or not files:
        manifest["status"] = "NO_FILES_LISTED"
        write_data_manifest(manifest, manifest_path)
        return manifest

    for index, file_info in enumerate(files):
        if not isinstance(file_info, dict):
            continue
        source_url = file_info.get("download_url")
        if not source_url and file_info.get("id"):
            source_url = f"https://ndownloader.figshare.com/files/{file_info['id']}"
        if not source_url:
            manifest["download_attempts"].append(
                {"file": file_info.get("name", f"file_{index}"), "status": "no_download_url"}
            )
            continue

        name = _safe_name(str(file_info.get("name") or f"figshare_file_{index}"), f"figshare_file_{index}")
        destination = raw_dir / name
        if destination.exists():
            destination = raw_dir / f"{index}_{name}"
        try:
            _download_file(str(source_url), destination)
            record = _manifest_file_record(destination, str(source_url), downloaded_at)
            record["figshare_file_id"] = file_info.get("id")
            record["figshare_reported_size_bytes"] = file_info.get("size")
            record["figshare_reported_md5"] = file_info.get("computed_md5") or file_info.get("md5")
            manifest["files"].append(record)
            manifest["download_attempts"].append(
                {"file": name, "source_url": source_url, "status": "downloaded"}
            )
        except (urllib.error.URLError, urllib.error.HTTPError, TimeoutError, OSError) as exc:
            manifest["download_attempts"].append(
                {"file": name, "source_url": source_url, "status": "failed", "error": str(exc)}
            )

    manifest["status"] = "DOWNLOADED" if manifest["files"] else "DOWNLOAD_FAILED"
    write_data_manifest(manifest, manifest_path)
    return manifest


def _manifest_file_record(path: Path, source_url: str, downloaded_at: str) -> dict[str, Any]:
    return {
        "file_name": path.name,
        "relative_path": str(path.relative_to(ROOT)),
        "source_url": source_url,
        "source_doi": PRIMARY_DATA_DOI,
        "download_timestamp_utc": downloaded_at,
        "size_bytes": path.stat().st_size,
        "sha256": compute_sha256(path),
    }


def _prior_manifest_records(path: Path) -> dict[str, dict[str, Any]]:
    if not path.exists():
        return {}
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return {}
    records: dict[str, dict[str, Any]] = {}
    for file_record in payload.get("files", []):
        if isinstance(file_record, dict) and file_record.get("file_name"):
            records[str(file_record["file_name"])] = file_record
    return records


def write_data_manifest(manifest: dict[str, Any], path: Path = MANIFEST_PATH) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def inspect_dataset_files(root: Path = ROOT) -> dict[str, Any]:
    """Inspect raw and derived files, extracting archives into derived storage."""

    raw_dir = root / "outputs" / "raw"
    extracted_dir = root / "outputs" / "derived_extracted"
    extracted_dir.mkdir(parents=True, exist_ok=True)
    raw_dir.mkdir(parents=True, exist_ok=True)

    files: list[dict[str, Any]] = []
    extraction_events: list[dict[str, Any]] = []
    for path in sorted(raw_dir.iterdir()):
        if not path.is_file():
            continue
        info = _file_info(path, "raw")
        files.append(info)
        if _is_archive(path):
            event = _extract_archive(path, extracted_dir / path.stem)
            extraction_events.append(event)

    for path in sorted(extracted_dir.rglob("*")):
        if path.is_file():
            files.append(_file_info(path, "derived_extracted"))

    inspection = {
        "source_doi": PRIMARY_DATA_DOI,
        "inspected_at_utc": utc_now_iso(),
        "files": files,
        "extraction_events": extraction_events,
        "supported_extensions": sorted(SUPPORTED_EXTENSIONS),
    }
    (root / "outputs").mkdir(parents=True, exist_ok=True)
    (root / "outputs" / "data_inspection.json").write_text(
        json.dumps(inspection, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return inspection


def _file_info(path: Path, role: str) -> dict[str, Any]:
    suffix = path.suffix.lower()
    return {
        "relative_path": str(path.relative_to(ROOT)),
        "file_name": path.name,
        "role": role,
        "extension": suffix,
        "size_bytes": path.stat().st_size,
        "sha256": compute_sha256(path),
        "parse_status": "candidate" if suffix in SUPPORTED_EXTENSIONS else "unsupported",
    }


def _is_archive(path: Path) -> bool:
    name = path.name.lower()
    return path.suffix.lower() in ARCHIVE_EXTENSIONS or name.endswith((".tar.gz", ".tar.bz2", ".tar.xz"))


def _extract_archive(path: Path, destination: Path) -> dict[str, Any]:
    destination.mkdir(parents=True, exist_ok=True)
    try:
        if zipfile.is_zipfile(path):
            with zipfile.ZipFile(path) as archive:
                extracted = _safe_extract_zip(archive, destination)
            return {"archive": str(path.relative_to(ROOT)), "status": "extracted", "files": extracted}
        if tarfile.is_tarfile(path):
            with tarfile.open(path) as archive:
                extracted = _safe_extract_tar(archive, destination)
            return {"archive": str(path.relative_to(ROOT)), "status": "extracted", "files": extracted}
        return {"archive": str(path.relative_to(ROOT)), "status": "unsupported_archive"}
    except (OSError, zipfile.BadZipFile, tarfile.TarError) as exc:
        return {"archive": str(path.relative_to(ROOT)), "status": "extract_failed", "error": str(exc)}


def _safe_extract_zip(archive: zipfile.ZipFile, destination: Path) -> list[str]:
    extracted: list[str] = []
    base = destination.resolve()
    for member in archive.infolist():
        if member.is_dir():
            continue
        target = (destination / member.filename).resolve()
        if not str(target).startswith(str(base)):
            continue
        target.parent.mkdir(parents=True, exist_ok=True)
        with archive.open(member) as source, target.open("wb") as output:
            output.write(source.read())
        extracted.append(str(target.relative_to(ROOT)))
    return extracted


def _safe_extract_tar(archive: tarfile.TarFile, destination: Path) -> list[str]:
    extracted: list[str] = []
    base = destination.resolve()
    for member in archive.getmembers():
        if not member.isfile():
            continue
        target = (destination / member.name).resolve()
        if not str(target).startswith(str(base)):
            continue
        source = archive.extractfile(member)
        if source is None:
            continue
        target.parent.mkdir(parents=True, exist_ok=True)
        with source, target.open("wb") as output:
            output.write(source.read())
        extracted.append(str(target.relative_to(ROOT)))
    return extracted


def load_time_traces(root: Path = ROOT, inspection: dict[str, Any] | None = None) -> list[dict[str, Any]]:
    traces: list[dict[str, Any]] = []
    for table in _iter_tables(root, inspection):
        columns = table["columns"]
        time_col = _find_column(columns, TIME_ALIASES)
        signal_col = _find_column(columns, SIGNAL_ALIASES)
        freq_col = _find_column(columns, FREQUENCY_ALIASES)
        if time_col is None or signal_col is None or freq_col is not None:
            continue
        rows = _numeric_rows(table["rows"])
        grouped = _group_rows(rows, columns, excluded={time_col, signal_col})
        for key, group_rows in grouped.items():
            time = _column_array(group_rows, time_col)
            signal = _column_array(group_rows, signal_col)
            if time.size < 8 or signal.size != time.size:
                continue
            order = np.argsort(time)
            metadata = _metadata_from_key(key)
            traces.append(
                {
                    "kind": "time_trace",
                    "source_file": table["relative_path"],
                    "sample_id": metadata.get("sample_id", "unknown"),
                    "condition_type": _condition_type(metadata),
                    "condition_value": _condition_value(metadata),
                    "distance_um": _float_or_none(metadata.get("distance_um")),
                    "time_ps": time[order],
                    "signal": signal[order],
                    "metadata": metadata,
                }
            )
    return traces


def load_frequency_spectra(root: Path = ROOT, inspection: dict[str, Any] | None = None) -> list[dict[str, Any]]:
    spectra: list[dict[str, Any]] = []
    for table in _iter_tables(root, inspection):
        columns = table["columns"]
        freq_col = _find_column(columns, FREQUENCY_ALIASES)
        amp_col = _find_column(columns, AMPLITUDE_ALIASES)
        peak_col = _find_column(columns, PEAK_AMPLITUDE_ALIASES)
        if freq_col is not None and amp_col is not None:
            rows = _numeric_rows(table["rows"])
            grouped = _group_rows(rows, columns, excluded={freq_col, amp_col})
            for key, group_rows in grouped.items():
                frequency = _column_array(group_rows, freq_col)
                amplitude = _column_array(group_rows, amp_col)
                if frequency.size < 4 or amplitude.size != frequency.size:
                    continue
                order = np.argsort(frequency)
                metadata = _metadata_from_key(key)
                spectra.append(
                    {
                        "kind": "frequency_spectrum",
                        "source_file": table["relative_path"],
                        "sample_id": metadata.get("sample_id", "unknown"),
                        "condition_type": _condition_type(metadata),
                        "condition_value": _condition_value(metadata),
                        "distance_um": _float_or_none(metadata.get("distance_um")),
                        "frequency_THz": frequency[order],
                        "amplitude": amplitude[order],
                        "metadata": metadata,
                    }
                )
            continue

        if peak_col is not None:
            rows = _numeric_rows(table["rows"])
            for row_index, row in enumerate(rows):
                metadata = _metadata_from_row(row, columns)
                peak_value = _float_or_none(row.get(peak_col))
                if peak_value is None:
                    continue
                spectra.append(
                    {
                        "kind": "peak_record",
                        "source_file": table["relative_path"],
                        "sample_id": str(metadata.get("sample_id", "unknown")),
                        "condition_type": _condition_type(metadata, row_index=row_index),
                        "condition_value": _condition_value(metadata, row_index=row_index),
                        "distance_um": _float_or_none(metadata.get("distance_um")),
                        "peak_amplitude": float(peak_value),
                        "frequency_THz": _float_or_none(metadata.get("frequency_THz")),
                        "coherence_time_ps": _float_or_none(metadata.get("coherence_time_ps")),
                        "metadata": metadata,
                    }
                )
    return spectra


def load_stft_amplitudes(root: Path = ROOT, inspection: dict[str, Any] | None = None) -> list[dict[str, Any]]:
    stft_records: list[dict[str, Any]] = []
    for table in _iter_tables(root, inspection):
        columns = table["columns"]
        time_col = _find_column(columns, TIME_ALIASES)
        freq_col = _find_column(columns, FREQUENCY_ALIASES)
        amp_col = _find_column(columns, AMPLITUDE_ALIASES)
        if time_col is None or freq_col is None or amp_col is None:
            continue
        rows = _numeric_rows(table["rows"])
        grouped = _group_rows(rows, columns, excluded={time_col, freq_col, amp_col})
        for key, group_rows in grouped.items():
            time = _column_array(group_rows, time_col)
            frequency = _column_array(group_rows, freq_col)
            amplitude = _column_array(group_rows, amp_col)
            if time.size < 8 or frequency.size != time.size or amplitude.size != time.size:
                continue
            metadata = _metadata_from_key(key)
            stft_records.append(
                {
                    "kind": "stft_amplitude",
                    "source_file": table["relative_path"],
                    "sample_id": metadata.get("sample_id", "unknown"),
                    "condition_type": _condition_type(metadata),
                    "condition_value": _condition_value(metadata),
                    "distance_um": _float_or_none(metadata.get("distance_um")),
                    "time_ps": time,
                    "frequency_THz": frequency,
                    "amplitude": amplitude,
                    "metadata": metadata,
                }
            )
    return stft_records


def build_smoke_test_dataset_if_needed(root: Path = ROOT) -> dict[str, Any]:
    """Create a deterministic synthetic dataset for pipeline verification only."""

    smoke_dir = root / "outputs" / "smoke_test_dataset"
    smoke_dir.mkdir(parents=True, exist_ok=True)
    frequency = np.linspace(2.45, 3.80, 240)
    distances = np.array([0.0, 2.0, 4.0, 6.0, 8.0], dtype=float)
    amplitudes = np.array([1.00, 0.92, 0.78, 0.55, 0.34], dtype=float)
    linewidths = np.array([0.055, 0.062, 0.082, 0.155, 0.235], dtype=float)
    drifts = np.array([0.000, 0.006, -0.011, 0.022, 0.052], dtype=float)

    spectra_rows: list[dict[str, Any]] = []
    for distance, amp, width, drift in zip(distances, amplitudes, linewidths, drifts):
        center = 3.13 + drift
        background = 0.045 + 0.010 * distance
        profile = background + amp * np.exp(-0.5 * ((frequency - center) / width) ** 2)
        shoulder = 0.030 * np.exp(-0.5 * ((frequency - 2.83) / 0.09) ** 2)
        for freq, value in zip(frequency, profile + shoulder):
            spectra_rows.append(
                {
                    "sample_id": "synthetic_NbOI2_smoke",
                    "distance_um": f"{distance:.3f}",
                    "frequency_THz": f"{freq:.6f}",
                    "amplitude": f"{value:.9f}",
                }
            )

    time = np.linspace(0.0, 14.0, 700)
    trace_rows: list[dict[str, Any]] = []
    for distance, amp in zip(distances, amplitudes):
        tau = max(2.2, 9.0 - 0.7 * distance)
        phase = 0.12 * distance
        signal = amp * np.sin(2.0 * math.pi * 3.13 * time + phase) * np.exp(-time / tau)
        signal += 0.018 * np.sin(2.0 * math.pi * 1.20 * time + 0.4)
        for t_value, s_value in zip(time, signal):
            trace_rows.append(
                {
                    "sample_id": "synthetic_NbOI2_smoke",
                    "distance_um": f"{distance:.3f}",
                    "time_ps": f"{t_value:.6f}",
                    "signal": f"{s_value:.9f}",
                }
            )

    _write_csv(smoke_dir / "smoke_spectra.csv", spectra_rows)
    _write_csv(smoke_dir / "smoke_time_traces.csv", trace_rows)

    inspection = {
        "source_doi": None,
        "inspected_at_utc": utc_now_iso(),
        "files": [
            _file_info(smoke_dir / "smoke_spectra.csv", "synthetic_smoke_test"),
            _file_info(smoke_dir / "smoke_time_traces.csv", "synthetic_smoke_test"),
        ],
        "extraction_events": [],
        "supported_extensions": sorted(SUPPORTED_EXTENSIONS),
        "smoke_test_only": True,
    }
    return {
        "inspection": inspection,
        "time_traces": load_time_traces(root, inspection),
        "frequency_spectra": load_frequency_spectra(root, inspection),
        "stft_amplitudes": load_stft_amplitudes(root, inspection),
        "files": [str(smoke_dir / "smoke_spectra.csv"), str(smoke_dir / "smoke_time_traces.csv")],
    }


def _write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        return
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def _iter_tables(root: Path, inspection: dict[str, Any] | None) -> list[dict[str, Any]]:
    if inspection is None:
        inspection = inspect_dataset_files(root)
    tables: list[dict[str, Any]] = []
    for file_info in inspection.get("files", []):
        if file_info.get("parse_status") != "candidate":
            continue
        path = ROOT / str(file_info["relative_path"])
        try:
            tables.extend(_load_table(path))
        except (OSError, ValueError, json.JSONDecodeError, ImportError):
            continue
    return tables


def _load_table(path: Path) -> list[dict[str, Any]]:
    suffix = path.suffix.lower()
    if suffix in {".csv", ".tsv", ".txt", ".dat"}:
        return [_read_text_table(path)]
    if suffix == ".json":
        return _read_json_tables(path)
    if suffix == ".npy":
        return [_read_npy_table(path)]
    if suffix == ".npz":
        return _read_npz_tables(path)
    if suffix == ".xlsx":
        return _read_xlsx_tables(path)
    raise ValueError(f"Unsupported table extension: {path}")


def _read_text_table(path: Path) -> dict[str, Any]:
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    lines = [line for line in lines if line.strip() and not line.lstrip().startswith("#")]
    if not lines:
        raise ValueError(f"Empty text table: {path}")
    delimiter = _guess_delimiter(path, lines[0])
    rows: list[dict[str, Any]] = []
    if delimiter == "whitespace":
        first = re.split(r"\s+", lines[0].strip())
        has_header = not all(_looks_float(token) for token in first)
        columns = [_normalize_header(value, index) for index, value in enumerate(first if has_header else [])]
        if not columns:
            columns = [f"col_{index}" for index in range(len(first))]
        data_lines = lines[1:] if has_header else lines
        for line in data_lines:
            values = re.split(r"\s+", line.strip())
            if len(values) >= len(columns):
                rows.append(dict(zip(columns, values[: len(columns)])))
    else:
        reader = csv.DictReader(lines, delimiter=delimiter)
        if reader.fieldnames is None:
            raise ValueError(f"No header in text table: {path}")
        rename = {name: _normalize_header(name, index) for index, name in enumerate(reader.fieldnames)}
        columns = [rename[name] for name in reader.fieldnames]
        for row in reader:
            rows.append({rename[key]: value for key, value in row.items() if key in rename})

    return {"relative_path": str(path.relative_to(ROOT)), "columns": columns, "rows": rows}


def _read_json_tables(path: Path) -> list[dict[str, Any]]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    candidates: list[Any]
    if isinstance(payload, list):
        candidates = [payload]
    elif isinstance(payload, dict):
        candidates = []
        for key in ("data", "rows", "values"):
            if isinstance(payload.get(key), list):
                candidates.append(payload[key])
        if _dict_of_arrays(payload):
            candidates.append(payload)
    else:
        raise ValueError(f"Unsupported JSON payload: {path}")

    tables: list[dict[str, Any]] = []
    for index, candidate in enumerate(candidates):
        rows = _rows_from_json_candidate(candidate)
        if not rows:
            continue
        columns = list(rows[0].keys())
        rel = str(path.relative_to(ROOT))
        if len(candidates) > 1:
            rel = f"{rel}#table_{index}"
        tables.append({"relative_path": rel, "columns": columns, "rows": rows})
    return tables


def _read_npy_table(path: Path) -> dict[str, Any]:
    array = np.asarray(np.load(path, allow_pickle=False))
    return _array_to_table(path, array)


def _read_npz_tables(path: Path) -> list[dict[str, Any]]:
    loaded = np.load(path, allow_pickle=False)
    keys = list(loaded.keys())
    keyed = {key: np.asarray(loaded[key]) for key in keys}
    if {"time_ps", "frequency_THz", "amplitude"}.issubset(keyed):
        return [_dict_arrays_to_table(path, keyed, keys=["time_ps", "frequency_THz", "amplitude"])]
    if {"frequency_THz", "amplitude"}.issubset(keyed):
        return [_dict_arrays_to_table(path, keyed, keys=["frequency_THz", "amplitude"])]
    if {"time_ps", "signal"}.issubset(keyed):
        return [_dict_arrays_to_table(path, keyed, keys=["time_ps", "signal"])]
    tables = []
    for key in keys:
        array = keyed[key]
        if array.ndim in {1, 2}:
            table = _array_to_table(path, array)
            table["relative_path"] = f"{path.relative_to(ROOT)}#{key}"
            tables.append(table)
    return tables


def _read_xlsx_tables(path: Path) -> list[dict[str, Any]]:
    try:
        import openpyxl  # type: ignore[import-not-found]
    except ImportError as exc:
        raise ImportError("XLSX parsing requires openpyxl.") from exc

    workbook = openpyxl.load_workbook(path, read_only=True, data_only=True)
    tables: list[dict[str, Any]] = []
    for sheet in workbook.worksheets:
        raw_rows = list(sheet.iter_rows(values_only=True))
        if not raw_rows:
            continue
        expanded = _expand_xy_sheet(path, sheet.title, raw_rows)
        if expanded:
            tables.extend(expanded)
            continue

        header = raw_rows[0]
        columns = _unique_columns(
            [_normalize_header(str(value), index) for index, value in enumerate(header) if value is not None]
        )
        rows: list[dict[str, Any]] = []
        for raw_row in raw_rows[1:]:
            values = list(raw_row)[: len(columns)]
            if not any(value is not None for value in values):
                continue
            rows.append({column: value for column, value in zip(columns, values)})
        if rows:
            tables.append(
                {
                    "relative_path": f"{path.relative_to(ROOT)}#{sheet.title}",
                    "columns": columns,
                    "rows": rows,
                }
            )
    return tables


def _expand_xy_sheet(path: Path, sheet_title: str, raw_rows: list[tuple[Any, ...]]) -> list[dict[str, Any]]:
    """Unpivot common source-data XLSX layouts into tidy two-column tables."""

    header = ["" if value is None else str(value) for value in raw_rows[0]]
    data_rows = raw_rows[1:]
    tables: list[dict[str, Any]] = []
    used_pairs: set[tuple[int, int]] = set()

    # Pattern A: repeated x/y pairs, for example:
    # dt (ps), dR/R 0um, dt (ps), dR/R 2um, ...
    index = 0
    while index < len(header) - 1:
        x_column = _classify_x_header(header[index])
        y_column = _classify_y_header(header[index + 1])
        if x_column and y_column:
            table = _xy_pair_to_table(
                path,
                sheet_title,
                data_rows,
                x_index=index,
                y_index=index + 1,
                x_column=x_column,
                y_column=y_column,
                y_label=header[index + 1],
            )
            if table is not None:
                tables.append(table)
                used_pairs.add((index, index + 1))
            index += 2
            continue
        index += 1

    if tables:
        return tables

    # Pattern B: one x column followed by multiple y columns, for example:
    # Frequency, b-axis, c-axis.
    if header:
        x_column = _classify_x_header(header[0])
        if x_column:
            for y_index in range(1, len(header)):
                y_column = _classify_y_header(header[y_index])
                if not y_column:
                    continue
                table = _xy_pair_to_table(
                    path,
                    sheet_title,
                    data_rows,
                    x_index=0,
                    y_index=y_index,
                    x_column=x_column,
                    y_column=y_column,
                    y_label=header[y_index],
                )
                if table is not None:
                    tables.append(table)

    return tables


def _xy_pair_to_table(
    path: Path,
    sheet_title: str,
    data_rows: list[tuple[Any, ...]],
    *,
    x_index: int,
    y_index: int,
    x_column: str,
    y_column: str,
    y_label: str,
) -> dict[str, Any] | None:
    metadata = _metadata_from_label(y_label)
    metadata["sample_id"] = metadata.get("sample_id") or _safe_name(f"{path.stem}_{sheet_title}", "sample")
    columns = [x_column, y_column]
    for metadata_column in ("sample_id", "distance_um", "temperature_K", "fluence_mJ_cm2", "thickness_nm"):
        if metadata.get(metadata_column) is not None:
            columns.append(metadata_column)

    rows: list[dict[str, Any]] = []
    for raw_row in data_rows:
        if x_index >= len(raw_row) or y_index >= len(raw_row):
            continue
        x_value = _float_or_none(raw_row[x_index])
        y_value = _float_or_none(raw_row[y_index])
        if x_value is None or y_value is None:
            continue
        row: dict[str, Any] = {x_column: x_value, y_column: y_value}
        for metadata_column in columns[2:]:
            row[metadata_column] = metadata.get(metadata_column)
        rows.append(row)

    if len(rows) < 4:
        return None
    suffix = _safe_name(y_label, f"series_{y_index}")
    return {
        "relative_path": f"{path.relative_to(ROOT)}#{sheet_title}:{suffix}",
        "columns": columns,
        "rows": rows,
    }


def _classify_x_header(label: str) -> str | None:
    lower = label.strip().lower()
    if not lower:
        return None
    if "time" in lower or "delay" in lower or lower.startswith("dt"):
        return "time_ps"
    if "frequency" in lower or lower in {"freq", "thz"}:
        return "frequency_THz"
    return None


def _classify_y_header(label: str) -> str | None:
    lower = label.strip().lower()
    if not lower:
        return None
    if "d" in lower and "r" in lower and "/" in lower:
        return "signal"
    if "signal" in lower or "reflectance" in lower or "deflection" in lower:
        return "signal"
    if "intensity" in lower or "amplitude" in lower or "power" in lower or "im(" in lower:
        return "amplitude"
    if lower in {"b-axis", "c-axis", "b axis", "c axis"}:
        return "amplitude"
    return None


def _metadata_from_label(label: str) -> dict[str, Any]:
    metadata: dict[str, Any] = {}
    distance_match = re.search(r"([0-9]+(?:\.[0-9]+)?)\s*u?m\b", label, flags=re.IGNORECASE)
    if distance_match:
        metadata["distance_um"] = float(distance_match.group(1))
    thickness_match = re.search(r"([0-9]+(?:\.[0-9]+)?)\s*nm\b", label, flags=re.IGNORECASE)
    if thickness_match:
        metadata["thickness_nm"] = float(thickness_match.group(1))
    temperature_match = re.search(r"([0-9]+(?:\.[0-9]+)?)\s*k\b", label, flags=re.IGNORECASE)
    if temperature_match:
        metadata["temperature_K"] = float(temperature_match.group(1))
    fluence_match = re.search(r"([0-9]+(?:\.[0-9]+)?)\s*mj", label, flags=re.IGNORECASE)
    if fluence_match:
        metadata["fluence_mJ_cm2"] = float(fluence_match.group(1))
    return metadata


def _unique_columns(columns: list[str]) -> list[str]:
    counts: dict[str, int] = {}
    unique: list[str] = []
    for column in columns:
        count = counts.get(column, 0)
        counts[column] = count + 1
        unique.append(column if count == 0 else f"{column}_{count + 1}")
    return unique


def _array_to_table(path: Path, array: np.ndarray) -> dict[str, Any]:
    if array.ndim == 1:
        columns = ["time_ps", "signal"]
        rows = [{"time_ps": float(index), "signal": float(value)} for index, value in enumerate(array)]
    elif array.ndim == 2 and array.shape[1] >= 3:
        lower_name = path.name.lower()
        if "stft" in lower_name or "timefreq" in lower_name or "time_frequency" in lower_name:
            columns = ["time_ps", "frequency_THz", "amplitude"]
        else:
            columns = ["distance_um", "frequency_THz", "amplitude"]
        rows = [
            {columns[0]: float(row[0]), columns[1]: float(row[1]), columns[2]: float(row[2])}
            for row in array[:, :3]
        ]
    elif array.ndim == 2 and array.shape[1] >= 2:
        lower_name = path.name.lower()
        if "freq" in lower_name or "spectrum" in lower_name or "fft" in lower_name or "thz" in lower_name:
            columns = ["frequency_THz", "amplitude"]
        else:
            columns = ["time_ps", "signal"]
        rows = [{columns[0]: float(row[0]), columns[1]: float(row[1])} for row in array[:, :2]]
    else:
        raise ValueError(f"Unsupported array shape {array.shape}: {path}")
    return {"relative_path": str(path.relative_to(ROOT)), "columns": columns, "rows": rows}


def _dict_arrays_to_table(path: Path, arrays: dict[str, np.ndarray], *, keys: list[str]) -> dict[str, Any]:
    size = min(arrays[key].size for key in keys)
    rows = []
    for index in range(size):
        rows.append({key: float(np.ravel(arrays[key])[index]) for key in keys})
    return {"relative_path": str(path.relative_to(ROOT)), "columns": keys, "rows": rows}


def _rows_from_json_candidate(candidate: Any) -> list[dict[str, Any]]:
    if isinstance(candidate, list) and candidate and isinstance(candidate[0], dict):
        return [
            {_normalize_header(str(key), index): value for index, (key, value) in enumerate(row.items())}
            for row in candidate
        ]
    if isinstance(candidate, dict) and _dict_of_arrays(candidate):
        keys = list(candidate.keys())
        size = min(len(candidate[key]) for key in keys)
        rows = []
        for row_index in range(size):
            rows.append(
                {
                    _normalize_header(str(key), index): candidate[key][row_index]
                    for index, key in enumerate(keys)
                }
            )
        return rows
    return []


def _dict_of_arrays(payload: dict[str, Any]) -> bool:
    values = list(payload.values())
    return bool(values) and all(isinstance(value, list) for value in values)


def _guess_delimiter(path: Path, first_line: str) -> str:
    suffix = path.suffix.lower()
    if suffix == ".csv":
        return ","
    if suffix == ".tsv":
        return "\t"
    if "," in first_line:
        return ","
    if "\t" in first_line:
        return "\t"
    return "whitespace"


def _normalize_header(name: str, index: int) -> str:
    if not name or name.lower() == "none":
        return f"col_{index}"
    normalized = name.strip().lower()
    normalized = normalized.replace("delta", "delta")
    normalized = normalized.replace("dr/r0", "delta_r_over_r")
    normalized = normalized.replace("dr/r_0", "delta_r_over_r")
    normalized = normalized.replace("frequency (thz)", "frequency_thz")
    normalized = normalized.replace("time (ps)", "time_ps")
    normalized = normalized.replace("delay (ps)", "delay_ps")
    normalized = normalized.replace("distance (um)", "distance_um")
    normalized = re.sub(r"[^a-z0-9]+", "_", normalized).strip("_")
    aliases = {
        "frequency_thz": "frequency_THz",
        "freq_thz": "frequency_THz",
        "time_ps": "time_ps",
        "delay_ps": "delay_ps",
        "distance_um": "distance_um",
        "temperature_k": "temperature_K",
        "fluence_mj_cm2": "fluence_mJ_cm2",
        "thickness_nm": "thickness_nm",
        "coherence_time_ps": "coherence_time_ps",
    }
    return aliases.get(normalized, normalized or f"col_{index}")


TIME_ALIASES = {
    "time_ps",
    "delay_ps",
    "pump_probe_delay_ps",
    "pump_probe_delay",
    "delay",
    "time",
}
FREQUENCY_ALIASES = {"frequency_THz", "freq_THz", "frequency", "freq", "thz"}
SIGNAL_ALIASES = {
    "delta_r_over_r",
    "delta_r_r0",
    "dr_over_r0",
    "dr_r0",
    "signal",
    "amplitude",
    "intensity",
    "reflectance",
    "deflection",
    "polarization",
}
AMPLITUDE_ALIASES = {"amplitude", "power", "intensity", "signal", "stft_amplitude"}
PEAK_AMPLITUDE_ALIASES = {"peak_amplitude", "peak_amp", "normalized_amplitude", "amplitude_3_13_thz"}
METADATA_COLUMNS = {
    "sample_id",
    "distance_um",
    "temperature_K",
    "fluence_mJ_cm2",
    "thickness_nm",
    "frequency_THz",
    "coherence_time_ps",
}


def _find_column(columns: list[str], aliases: set[str]) -> str | None:
    normalized_aliases = {_alias_key(alias) for alias in aliases}
    for column in columns:
        if _alias_key(column) in normalized_aliases:
            return column
    return None


def _alias_key(value: str) -> str:
    return re.sub(r"[^a-z0-9]+", "", value.lower())


def _numeric_rows(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    return rows


def _group_rows(
    rows: list[dict[str, Any]],
    columns: list[str],
    *,
    excluded: set[str],
) -> dict[tuple[tuple[str, Any], ...], list[dict[str, Any]]]:
    metadata_columns = [column for column in columns if column not in excluded and column in METADATA_COLUMNS]
    if not metadata_columns:
        metadata_columns = ["sample_id"] if "sample_id" in columns else []
    grouped: dict[tuple[tuple[str, Any], ...], list[dict[str, Any]]] = {}
    for row in rows:
        key_items = []
        for column in metadata_columns:
            value = row.get(column)
            if value not in (None, ""):
                key_items.append((column, _stable_value(value)))
        key = tuple(key_items) if key_items else (("sample_id", "unknown"),)
        grouped.setdefault(key, []).append(row)
    return grouped


def _metadata_from_key(key: tuple[tuple[str, Any], ...]) -> dict[str, Any]:
    return {column: value for column, value in key}


def _metadata_from_row(row: dict[str, Any], columns: list[str]) -> dict[str, Any]:
    metadata: dict[str, Any] = {}
    for column in columns:
        if column in METADATA_COLUMNS and row.get(column) not in (None, ""):
            metadata[column] = _stable_value(row[column])
    return metadata


def _condition_type(metadata: dict[str, Any], *, row_index: int | None = None) -> str:
    for column in ("distance_um", "fluence_mJ_cm2", "temperature_K", "thickness_nm", "delay_ps", "time_ps"):
        if metadata.get(column) not in (None, ""):
            return column
    return "condition_index" if row_index is not None else "dataset"


def _condition_value(metadata: dict[str, Any], *, row_index: int | None = None) -> Any:
    condition = _condition_type(metadata, row_index=row_index)
    if condition == "condition_index":
        return row_index
    if condition == "dataset":
        return 0.0
    return metadata.get(condition)


def _stable_value(value: Any) -> Any:
    numeric = _float_or_none(value)
    if numeric is not None and math.isfinite(numeric):
        return float(numeric)
    return str(value)


def _float_or_none(value: Any) -> float | None:
    if value in (None, ""):
        return None
    try:
        numeric = float(value)
    except (TypeError, ValueError):
        return None
    if not math.isfinite(numeric):
        return None
    return numeric


def _column_array(rows: list[dict[str, Any]], column: str) -> np.ndarray:
    values = [_float_or_none(row.get(column)) for row in rows]
    return np.asarray([value for value in values if value is not None], dtype=float)


def _looks_float(value: str) -> bool:
    try:
        float(value)
        return True
    except ValueError:
        return False
