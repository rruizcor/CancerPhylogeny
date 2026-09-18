from __future__ import annotations

import logging
from pathlib import Path
from typing import Iterable


def configure_logging() -> None:
    logging.basicConfig(
        level=logging.INFO,
        format="[%(asctime)s] %(levelname)-5s %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )


def find_project_root(start: Path | None = None) -> Path:
    current = (start or Path.cwd()).resolve()
    for candidate in [current, *current.parents]:
        if (candidate / "config" / "tcga_projects.yaml").exists():
            return candidate
    raise FileNotFoundError("Could not locate project root containing config/tcga_projects.yaml")


def project_path(*parts: str, root: Path | None = None) -> Path:
    return (root or find_project_root()).joinpath(*parts)


def require_file(path: Path, description: str = "required file") -> Path:
    if not path.exists():
        raise FileNotFoundError(f"Missing {description}: {path}")
    return path


def ensure_dir(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    return path


def not_implemented(step: str, next_action: str) -> None:
    raise NotImplementedError(f"{step} is a starter placeholder. Next implementation step: {next_action}")


def report_missing(paths: Iterable[Path]) -> list[Path]:
    return [path for path in paths if not path.exists()]
