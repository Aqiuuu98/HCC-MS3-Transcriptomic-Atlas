"""Path adapter for preserved analysis sources; performs no analysis or writes.

HCC_PROJECT_ROOT must identify an isolated restored workspace. Frozen inputs
should be mounted read-only there and derived outputs directed to its writable
copy. This module deliberately refuses the original project and this release.
"""
from pathlib import Path
import os


def project_path(relative=""):
    value = os.environ.get("HCC_PROJECT_ROOT", "").strip()
    if not value:
        raise RuntimeError("Set HCC_PROJECT_ROOT to an isolated restored analysis workspace.")
    normalized = value.replace("\\", "/").rstrip("/").casefold()
    if any(normalized == item or normalized.startswith(item + "/")
           for item in ("e:/tcga", "/mnt/e/tcga")):
        raise RuntimeError("The original frozen project cannot be used as the execution workspace.")
    root = Path(value).expanduser().resolve()
    release = Path(__file__).resolve().parent
    if root == release or root.is_relative_to(release):
        raise RuntimeError("Choose an execution workspace outside the release archive.")
    if not root.is_dir():
        raise RuntimeError("HCC_PROJECT_ROOT must be an existing restored workspace.")
    suffix = relative.replace("\\", "/").lstrip("/")
    result = (root / suffix).resolve()
    if not result.is_relative_to(root):
        raise RuntimeError("A project path must remain inside HCC_PROJECT_ROOT.")
    return str(result)
