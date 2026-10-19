"""
SARDINE: epigenetic clock and phenotype associations from EWAS summary statistics.

Estimates the association between an epigenetic clock and a phenotype using only EWAS
summary statistics, the clock's CpG weights, and a reference DNA methylation covariance
matrix. No individual-level data is required.
"""

from ._version import __version__
from .result import SardineResult
from .sardine import example, sardine
from .utils import SardineError, SardineWarning

__all__ = [
    "SardineError",
    "SardineResult",
    "SardineWarning",
    "__version__",
    "example",
    "sardine",
]
