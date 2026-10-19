# Where the shared fixtures are: tests/fixtures/ inside an unpacked sdist (the
# build copies them there), otherwise the repo-level fixtures/ directory.

from pathlib import Path

_TESTS_DIR = Path(__file__).resolve().parent
FIXTURES_DIR = next(
    path for path in [_TESTS_DIR / "fixtures", _TESTS_DIR.parents[1] / "fixtures"]
    if (path / "test_manifest.tsv").exists()
)
