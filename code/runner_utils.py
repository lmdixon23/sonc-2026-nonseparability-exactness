#!/usr/bin/env python3
"""Shared deterministic subprocess runner for release verification."""
from __future__ import annotations

from pathlib import Path
from typing import Iterable
import json
import subprocess
import sys


def run_scripts(
    base: Path,
    scripts: Iterable[str],
    *,
    timeout_seconds: int,
    output_json: Path,
    echo: bool = True,
) -> list[dict[str, object]]:
    records: list[dict[str, object]] = []
    for name in scripts:
        path = base / name
        if not path.is_file():
            raise FileNotFoundError(path)
        completed = subprocess.run(
            [sys.executable, str(path)],
            cwd=str(base.parent),
            capture_output=True,
            text=True,
            timeout=timeout_seconds,
            check=False,
        )
        record = {
            "script": name,
            "returncode": completed.returncode,
            "stdout": completed.stdout,
            "stderr": completed.stderr,
        }
        records.append(record)
        if echo:
            print(f"=== {name} ===")
            print(completed.stdout, end="")
            if completed.stderr:
                print(completed.stderr, end="", file=sys.stderr)
        if completed.returncode != 0:
            output_json.parent.mkdir(parents=True, exist_ok=True)
            output_json.write_text(json.dumps(records, indent=2) + "\n", encoding="utf-8")
            raise SystemExit(completed.returncode)
    output_json.parent.mkdir(parents=True, exist_ok=True)
    output_json.write_text(json.dumps(records, indent=2) + "\n", encoding="utf-8")
    return records
