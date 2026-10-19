# =============================================================================
# python/src/sardine/utils.py — logging and small numeric helpers shared across
# the package (mirrors r/R/utils.R)
# =============================================================================

import math
import numbers
import re
import sys
import warnings

import numpy as np
import pandas as pd
from scipy import stats

from ._version import __version__


class SardineError(ValueError):
    """
    Raised when SARDINE's inputs are invalid (the counterpart of R's stop()).
    """


class SardineWarning(UserWarning):
    """
    Issued for inputs that are allowed but need the user's attention.
    """


# Matches one sprintf() conversion specification, e.g. "%d", "%.6g", "%s" or "%%"
_FORMAT_SPEC = re.compile(r"%(?:%|[-+ 0#]*\d*(?:\.\d+)?[diouxXeEfFgGs])")


def format_number(x):
    """
    Format a number as R's as.character() does: integers without a decimal point, other
    values with up to 15 significant digits, NA/Inf spelled as in R.

    The formatting needs to match so that we can test that the behavior of the R and
    Python packages is identical.
    """
    if isinstance(x, numbers.Integral):
        return str(int(x))
    x = float(x)
    if math.isnan(x):
        return "NA"
    if math.isinf(x):
        return "Inf" if x > 0 else "-Inf"
    return f"{x:.15g}"


def sprintf(fmt, *args):
    """
    Format like R's sprintf(): numbers passed to %s print as in R, and non-finite
    numbers print as NA / Inf / -Inf in any numeric conversion.

    The formatting needs to match so that we can test that the behavior of the R and
    Python packages is identical.
    """
    converted = list(args)
    pieces = []
    last = 0
    i = 0
    for match in _FORMAT_SPEC.finditer(fmt):
        spec = match.group(0)
        if spec != "%%":
            arg = converted[i]
            is_number = isinstance(arg, numbers.Real) and not isinstance(arg, bool)
            if spec.endswith("s") and is_number:
                converted[i] = format_number(arg)
            elif is_number and not math.isfinite(float(arg)):
                spec, converted[i] = "%s", format_number(arg)
            i += 1
        pieces.append(fmt[last : match.start()] + spec)
        last = match.end()
    pieces.append(fmt[last:])
    return "".join(pieces) % tuple(converted)


def log_message(verbose, level, fmt, *args):
    """
    Emit a tagged log message, honouring the verbose flag.

    Every diagnostic message in the package goes through this function, so it all goes
    to stderr and `verbose=False` reliably silences all except for warnings.

    level is one of "INFO", "NOTE" or "WARN", or None to print the message without a tag
    (e.g. the banner or a blank line). fmt is a sprintf() format string and args are the
    values interpolated into it.
    """
    if level is not None and level.upper() == "WARN":
        warnings.warn(
            sprintf("[%s] %s", level, sprintf(fmt, *args)), SardineWarning, stacklevel=3
        )
        return

    if not verbose:
        return
    text = sprintf(fmt, *args)
    print(text if level is None else sprintf("[%s] %s", level, text), file=sys.stderr)


def log_table(verbose, label, x, digits=None):
    """
    Log a preview of a table: a blank line, an "[INFO] <label>:" line, then the table.
    """
    if not verbose:
        return
    float_format = (lambda v: f"{v:.{digits}g}") if digits is not None else None
    table_text = x.to_string(float_format=float_format)
    log_message(verbose, None, "")
    log_message(verbose, "INFO", "%s:\n%s", label, table_text)


def is_missing(x):
    """
    True for None or NaN (R's NA).
    """
    return x is None or (isinstance(x, float) and math.isnan(x))


def append_note(x, msg):
    """
    Append a note to an existing note string without producing "NA; message".
    """
    if is_missing(msg):
        return x
    return msg if is_missing(x) else f"{x}; {msg}"


def sqrt_or_na(x):
    """
    Square root that returns NaN rather than warning on negative input.

    Note: Negative variance estimates are possible if Var(Y) is too large
    relative to the clock variance.
    """
    with np.errstate(invalid="ignore"):
        return np.sqrt(np.float64(x))


def two_sided_normal_pvalue(z):
    """
    Two-sided normal p-value from a Z statistic.
    """
    return 2 * stats.norm.sf(np.abs(z))


def two_sided_neg_log10_pvalue(z):
    """
    Two-sided -log10 p-value, computed in log space to avoid underflow to zero.
    """
    return -(stats.norm.logsf(np.abs(z)) + np.log(2)) / np.log(10)


def format_elapsed_time(elapsed_seconds):
    """
    Format an elapsed duration in seconds as "1 hour, 2 minutes, 3 seconds".
    """
    hours = math.floor(elapsed_seconds / 3600)
    minutes = math.floor((elapsed_seconds - hours * 3600) / 60)
    seconds = round(elapsed_seconds - hours * 3600 - minutes * 60)

    parts = []
    if hours > 0:
        parts.append(f"{hours} hour{'' if hours == 1 else 's'}")
    if minutes > 0:
        parts.append(f"{minutes} minute{'' if minutes == 1 else 's'}")
    if seconds > 0 or len(parts) == 0:
        parts.append(f"{seconds} second{'' if seconds == 1 else 's'}")
    return ", ".join(parts)


def match_option(value, choices, arg_name):
    """
    Match an option argument to its canonical value.

    Matching ignores case and whitespace, so "y ~ m" matches "Y~M".
    """

    def normalize(x):
        return re.sub(r"\s", "", x).upper()

    if isinstance(value, str):
        matched = [
            choice for choice in choices if normalize(choice) == normalize(value)
        ]
        if len(matched) == 1:
            return matched[0]

    got = f'"{value}"' if isinstance(value, str) else repr(value)
    raise SardineError(
        sprintf(
            "[ERROR] '%s' must be one of %s; got %s.",
            arg_name,
            ", ".join(f'"{c}"' for c in choices),
            got,
        )
    )


def sardine_version():
    """
    SARDINE package version.
    """
    return __version__


def match(x, table):
    """
    Positions of the first match of each element of x in table, -1 where there is none
    (akin to R's match()).
    """
    first = {}
    for i, value in enumerate(table):
        first.setdefault(value, i)
    return np.array([first.get(value, -1) for value in x], dtype=int)


def intersect(x, y):
    """
    Unique elements of x that are also in y, in the order they appear in x (equivalent
    to R's intersect()).

    The ordering must match R's behavior so that it is easy to compare the outputs of
    the two packages.
    """
    return pd.Index(x).intersection(y).tolist()
