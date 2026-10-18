"""
Sync shared, repo-level files into the R and Python packages.

The repo root holds the single source of truth for:
  - VERSION      -> r/DESCRIPTION (Version:), python/pyproject.toml (version =)
  - fixtures/    -> r/inst/extdata/fixtures/ (R CMD check cannot see ../fixtures)

Usage:
  python tools/sync.py           # write the copies
  python tools/sync.py --check   # exit 1 if any copy is out of date; changes nothing
"""

import argparse
import filecmp
import re
import shutil
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent

VERSION_FILE = ROOT / "VERSION"
R_DESCRIPTION = ROOT / "r" / "DESCRIPTION"
PY_PYPROJECT = ROOT / "python" / "pyproject.toml"

FIXTURES_SRC = ROOT / "fixtures"
FIXTURES_DST = ROOT / "r" / "inst" / "extdata" / "fixtures"
FIXTURE_ITEMS = ["test_manifest.tsv", "inputs", "expected"]


def read_version():
    """
    Read VERSION, which must be plain X.Y.Z (valid for both R and PEP 440).
    """
    version = VERSION_FILE.read_text().strip()
    if not re.fullmatch(r"\d+\.\d+\.\d+", version):
        sys.exit(f"[sync] VERSION must be X.Y.Z (e.g. 0.1.0); got {version!r}")
    return version


def sync_version_field(path, pattern, version, check):
    """
    Rewrite the version field in a metadata file.

    Returns True if it was out of date.
    """
    if not path.exists():
        return False
    text = path.read_text()
    match = re.search(pattern, text, flags=re.MULTILINE)
    if match is None:
        sys.exit(f"[sync] no version field found in {path.relative_to(ROOT)}")
    if match.group(1) == version:
        return False
    if not check:
        path.write_text(text[: match.start(1)] + version + text[match.end(1) :])
    print(f"[sync] {path.relative_to(ROOT)}: version {match.group(1)} -> {version}")
    return True


def visible_files(root):
    """
    Relative paths of all files under root, skipping hidden files and directories.
    """
    return {
        p.relative_to(root)
        for p in root.rglob("*")
        if p.is_file()
        and not any(part.startswith(".") for part in p.relative_to(root).parts)
    }


def dirs_differ(src, dst):
    """
    True if the directory trees src and dst differ in any (non-hidden) file name or
    content.
    """
    if not dst.exists():
        return True
    src_files, dst_files = visible_files(src), visible_files(dst)
    if src_files != dst_files:
        return True
    return any(not filecmp.cmp(src / f, dst / f, shallow=False) for f in src_files)


def sync_fixtures(check):
    """
    Mirror the fixture items into the R package.

    Returns True if they were out of date.
    """
    stale = False
    for item in FIXTURE_ITEMS:
        src, dst = FIXTURES_SRC / item, FIXTURES_DST / item
        differs = (
            dirs_differ(src, dst)
            if src.is_dir()
            else (not dst.exists() or not filecmp.cmp(src, dst, shallow=False))
        )
        if not differs:
            continue
        stale = True
        print(f"[sync] {dst.relative_to(ROOT)} is out of date")
        if check:
            continue
        if dst.is_dir():
            shutil.rmtree(dst)
        dst.parent.mkdir(parents=True, exist_ok=True)
        if src.is_dir():
            shutil.copytree(src, dst, ignore=shutil.ignore_patterns(".*"))
        else:
            shutil.copy2(src, dst)
    return stale


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "--check",
        action="store_true",
        help="report out-of-date copies without changing them",
    )
    args = parser.parse_args()

    version = read_version()
    stale = [
        sync_version_field(
            R_DESCRIPTION, r"^Version:\s*(\S+)\s*$", version, args.check
        ),
        sync_version_field(
            PY_PYPROJECT, r'^version\s*=\s*"([^"]+)"', version, args.check
        ),
        sync_fixtures(args.check),
    ]

    if args.check and any(stale):
        sys.exit("[sync] out of sync; run `pixi run -e dev sync`")
    print(
        f"[sync] {'all copies match' if not any(stale) else 'synced'} (VERSION {version})"
    )


if __name__ == "__main__":
    main()
