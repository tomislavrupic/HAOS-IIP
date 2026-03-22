#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import re
import shutil
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PAPERS_DIR = ROOT / "papers"
PDF_RELEASES_DIR = PAPERS_DIR / "pdf_releases"
INDEX_PATH = PDF_RELEASES_DIR / "INDEX.md"
RELEASE_TARGET_PATH = ROOT / "RELEASE_TARGET.json"


def run(cmd: list[str], *, cwd: Path | None = None, capture_output: bool = False) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        cmd,
        cwd=str(cwd) if cwd else None,
        check=True,
        text=True,
        capture_output=capture_output,
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Compile a phase paper, create the next numbered PDF release, commit, and push.",
    )
    parser.add_argument("--phase-dir", help="Phase directory to stage, relative to repo root.")
    parser.add_argument("--paper", help="Paper .tex path, relative to repo root.")
    parser.add_argument("--title", help="Release title. Defaults to \\title{...} from the TeX source.")
    parser.add_argument("--commit-message", help="Commit message. Defaults to a freeze-style message.")
    parser.add_argument(
        "--extra",
        action="append",
        default=[],
        help="Additional repo-relative path to stage. Can be repeated.",
    )
    parser.add_argument("--write-target", action="store_true", help="Write RELEASE_TARGET.json from the resolved inputs and exit.")
    parser.add_argument("--dry-run", action="store_true", help="Show the computed release plan without writing, committing, or pushing.")
    return parser.parse_args()


def repo_rel(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def require_path(path: Path, kind: str) -> Path:
    if not path.exists():
        raise FileNotFoundError(f"{kind} not found: {repo_rel(path) if path.is_absolute() and path.exists() else path}")
    return path


def extract_title(tex_path: Path) -> str:
    text = tex_path.read_text(encoding="utf-8")
    match = re.search(r"\\title\{(.+?)\}", text, flags=re.DOTALL)
    if not match:
        raise RuntimeError(f"Could not extract \\title{{...}} from {repo_rel(tex_path)}")
    title = re.sub(r"\s+", " ", match.group(1)).strip()
    return title


def load_release_target() -> dict:
    if not RELEASE_TARGET_PATH.exists():
        return {}
    with RELEASE_TARGET_PATH.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def write_release_target(payload: dict) -> None:
    RELEASE_TARGET_PATH.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")


def git_dirty_paths() -> list[str]:
    result = run(
        ["git", "-c", "core.fsmonitor=false", "status", "--short", "--untracked-files=all"],
        cwd=ROOT,
        capture_output=True,
    )
    paths: list[str] = []
    for line in result.stdout.splitlines():
        if not line.strip():
            continue
        path = line[3:]
        if " -> " in path:
            path = path.split(" -> ", 1)[1]
        paths.append(path.strip())
    return paths


def infer_phase_dir(dirty_paths: list[str]) -> str | None:
    candidates = set()
    for path in dirty_paths:
        parts = Path(path).parts
        if not parts:
            continue
        top = parts[0]
        if top.startswith("phase") and top not in {"phase.yaml"}:
            candidates.add(top)
    if len(candidates) == 1:
        return next(iter(candidates))
    return None


def infer_paper(dirty_paths: list[str]) -> str | None:
    candidates = sorted(
        path
        for path in dirty_paths
        if path.startswith("papers/") and path.endswith(".tex") and "/drafts/" not in path
    )
    if len(candidates) == 1:
        return candidates[0]
    return None


def resolve_release_inputs(args: argparse.Namespace) -> tuple[str, str, str | None, str | None, list[str]]:
    target = load_release_target()
    dirty_paths = git_dirty_paths()

    phase_dir = args.phase_dir or target.get("phase_dir") or infer_phase_dir(dirty_paths)
    paper = args.paper or target.get("paper") or infer_paper(dirty_paths)
    title = args.title or target.get("title")
    commit_message = args.commit_message or target.get("commit_message")
    extra = list(target.get("extra", [])) + list(args.extra or [])

    if not phase_dir or not paper:
        hints = []
        inferred_phase = infer_phase_dir(dirty_paths)
        inferred_paper = infer_paper(dirty_paths)
        if inferred_phase:
            hints.append(f"inferred phase_dir={inferred_phase}")
        if inferred_paper:
            hints.append(f"inferred paper={inferred_paper}")
        hint_text = "; ".join(hints) if hints else "no unique dirty phase dir or paper .tex detected"
        raise RuntimeError(
            "Release target is ambiguous. Set RELEASE_TARGET.json or pass --phase-dir and --paper. "
            f"Current detection: {hint_text}."
        )

    deduped_extra: list[str] = []
    seen = set()
    for item in extra:
        if item not in seen:
            deduped_extra.append(item)
            seen.add(item)
    return phase_dir, paper, title, commit_message, deduped_extra


def sanitize_filename(text: str) -> str:
    cleaned = text.replace("/", "-").replace(":", " -")
    cleaned = re.sub(r"\s+", " ", cleaned).strip()
    cleaned = re.sub(r'[<>:"/\\\\|?*]+', "-", cleaned)
    return cleaned


def next_release_label() -> str:
    pattern = re.compile(r"^(\d+)\.(\d+)\s+.+\.pdf$")
    max_major = 0
    for path in PDF_RELEASES_DIR.iterdir():
        if not path.is_file() or path.name == "INDEX.md":
            continue
        match = pattern.match(path.name)
        if match:
            max_major = max(max_major, int(match.group(1)))
    return f"{max_major + 1}.1"


def compile_tex(tex_path: Path) -> Path:
    cmd = ["pdflatex", "-interaction=nonstopmode", "-halt-on-error", tex_path.name]
    run(cmd, cwd=tex_path.parent)
    run(cmd, cwd=tex_path.parent)
    pdf_path = tex_path.with_suffix(".pdf")
    require_path(pdf_path, "compiled PDF")
    return pdf_path


def update_index(release_label: str, title: str) -> None:
    line = f"- `{release_label}` {title}"
    text = INDEX_PATH.read_text(encoding="utf-8")
    if line in text:
        return

    marker = "\nRule for future papers:\n"
    if marker in text:
        text = text.replace(marker, f"{line}\n{marker}")
    else:
        text = text.rstrip() + "\n" + line + "\n"
    INDEX_PATH.write_text(text, encoding="utf-8")


def current_branch() -> str:
    result = run(["git", "rev-parse", "--abbrev-ref", "HEAD"], cwd=ROOT, capture_output=True)
    branch = result.stdout.strip()
    if branch == "HEAD":
        raise RuntimeError("Release command requires a named branch, not detached HEAD.")
    return branch


def stage_paths(paths: list[Path]) -> None:
    cmd = ["git", "-c", "core.fsmonitor=false", "add", "--"]
    cmd.extend(str(path) for path in paths)
    run(cmd, cwd=ROOT)


def commit(commit_message: str) -> str:
    run(["git", "-c", "core.fsmonitor=false", "commit", "-m", commit_message], cwd=ROOT)
    result = run(["git", "rev-parse", "--short", "HEAD"], cwd=ROOT, capture_output=True)
    return result.stdout.strip()


def push(branch: str) -> None:
    run(["git", "-c", "core.fsmonitor=false", "push", "origin", branch], cwd=ROOT)


def main() -> int:
    args = parse_args()
    phase_dir_value, paper_value, title_override, commit_override, extra_values = resolve_release_inputs(args)

    phase_dir = require_path(ROOT / phase_dir_value, "phase directory")
    tex_path = require_path(ROOT / paper_value, "paper source")
    if tex_path.suffix != ".tex":
        raise RuntimeError(f"--paper must point to a .tex file, got {repo_rel(tex_path)}")

    title = title_override or extract_title(tex_path)
    release_label = next_release_label()
    release_filename = f"{release_label} {sanitize_filename(title)}.pdf"
    release_pdf_path = PDF_RELEASES_DIR / release_filename
    source_pdf_path = tex_path.with_suffix(".pdf")
    commit_message = commit_override or f"Freeze {phase_dir.name} release {release_label}"
    extra_paths = [ROOT / extra for extra in extra_values]

    if args.write_target:
        payload = {
            "phase_dir": phase_dir_value,
            "paper": paper_value,
            "title": title,
            "commit_message": commit_message,
            "extra": extra_values,
        }
        write_release_target(payload)
        print(f"wrote_target={repo_rel(RELEASE_TARGET_PATH)}")
        return 0

    if args.dry_run:
        print(f"phase_dir={repo_rel(phase_dir)}")
        print(f"paper={repo_rel(tex_path)}")
        print(f"title={title}")
        print(f"release_label={release_label}")
        print(f"release_pdf={repo_rel(release_pdf_path)}")
        print(f"commit_message={commit_message}")
        if extra_paths:
            print("extra_paths=" + ",".join(repo_rel(path) for path in extra_paths))
        return 0

    compiled_pdf_path = compile_tex(tex_path)
    if compiled_pdf_path != source_pdf_path:
        raise RuntimeError("Compiled PDF path mismatch.")

    shutil.copy2(source_pdf_path, release_pdf_path)
    update_index(release_label, title)

    paths_to_stage = [phase_dir, tex_path, source_pdf_path, release_pdf_path, INDEX_PATH]
    for path in extra_paths:
        require_path(path, "extra path")
        paths_to_stage.append(path)

    stage_paths(paths_to_stage)
    branch = current_branch()
    commit_hash = commit(commit_message)
    push(branch)

    print(f"release_label={release_label}")
    print(f"commit={commit_hash}")
    print(f"branch={branch}")
    print(f"release_pdf={repo_rel(release_pdf_path)}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except subprocess.CalledProcessError as exc:
        print(f"release_phase.py failed while running: {' '.join(exc.cmd)}", file=sys.stderr)
        raise
    except Exception as exc:  # noqa: BLE001
        print(f"release_phase.py error: {exc}", file=sys.stderr)
        raise SystemExit(1)
