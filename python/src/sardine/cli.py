# =============================================================================
# python/src/sardine/cli.py — the `sardine` command-line interface
# =============================================================================

import argparse
import os
import sys
import warnings

from ._version import __version__
from .sardine import example, sardine
from .utils import SardineError

DESCRIPTION = """\
Estimate the association between an epigenetic clock and a phenotype from EWAS
summary statistics, clock weights and a reference DNAm covariance matrix. No 
individual summary statistics needed.

Writes <out>_results.tsv (one row per estimated parameter) and
<out>_metadata.tsv (CpG overlap counts and Var(Y) diagnostics), or results.tsv
and metadata.tsv inside <out> if it is a folder. Progress and warnings go to
stderr.
"""

EPILOG = """\
Option values (ewas_direction, ewas_model, ewas_DNAm_scale, ewas_pheno_type)
are matched ignoring case and whitespace. See
https://github.com/qlu-lab/SARDINE/blob/main/inputs.md for every input's
accepted values and file formats.

example:
  sardine --ewas_sumstats ewas.tsv --clock clock.tsv \\
          --DNAm_cov_matrix sigma.tsv --out results/my_trait
"""


class _ExampleAction(argparse.Action):
    """--example: run sardine() on the bundled example data, print the result and
    exit, without needing any other arguments (like --version)."""

    def __init__(
        self,
        option_strings,
        dest=argparse.SUPPRESS,
        default=argparse.SUPPRESS,
        help=None,
    ):
        super().__init__(
            option_strings=option_strings,
            dest=dest,
            default=default,
            nargs=0,
            help=help,
        )

    def __call__(self, parser, namespace, values, option_string=None):
        print(example())
        parser.exit()


def build_parser():
    parser = argparse.ArgumentParser(
        prog="sardine",
        description=DESCRIPTION,
        epilog=EPILOG,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    required = parser.add_argument_group("required arguments")
    required.add_argument(
        "--ewas_sumstats",
        required=True,
        metavar="FILE",
        help="EWAS summary statistics, with columns for CpG, BETA, SE and N",
    )
    required.add_argument(
        "--clock",
        required=True,
        metavar="FILE",
        help="Clock weights, with columns for CpG and CoefficientTraining",
    )
    required.add_argument(
        "--DNAm_cov_matrix",
        required=True,
        metavar="FILE",
        help="Reference DNAm covariance matrix: CpG IDs in the first column, and a header "
        "listing the same CpG IDs in the same order",
    )
    required.add_argument(
        "--out",
        required=True,
        metavar="PREFIX",
        help="Output prefix: writes PREFIX_results.tsv and PREFIX_metadata.tsv. If PREFIX ends "
        "with / or is an existing folder, writes results.tsv and metadata.tsv inside it",
    )

    model = parser.add_argument_group("EWAS model")
    model.add_argument(
        "--ewas_direction",
        default="M~Y",
        metavar="{M~Y,Y~M}",
        help="Direction of the EWAS regression, outcome ~ exposure (default: %(default)s)",
    )
    model.add_argument(
        "--ewas_model",
        default="linear",
        metavar="{linear,logistic,cox}",
        help="Regression model the EWAS used (default: %(default)s)",
    )
    model.add_argument(
        "--ewas_DNAm_scale",
        default="beta",
        metavar="{beta,m}",
        help="Methylation scale of the EWAS: beta values or M-values (default: %(default)s)",
    )
    model.add_argument(
        "--ewas_pheno_type",
        default="continuous",
        metavar="{continuous,binary}",
        help="Phenotype type; binary also estimates gamma_logistic (default: %(default)s)",
    )

    optional = parser.add_argument_group("optional inputs")
    optional.add_argument(
        "--mean_beta",
        metavar="FILE",
        help="Mean beta value per CpG; required when --ewas_DNAm_scale m",
    )
    optional.add_argument(
        "--varY",
        type=float,
        metavar="NUMBER",
        help="Phenotype variance, with the same covariate adjustment as the EWAS "
        "(default: estimated from the summary statistics)",
    )
    optional.add_argument(
        "--N",
        type=float,
        metavar="NUMBER",
        help="EWAS sample size (default: the median of the summary statistics' N column)",
    )
    optional.add_argument(
        "--cox_pi",
        type=float,
        metavar="NUMBER",
        help="Cox only: proportion of cases at the end of the study",
    )
    optional.add_argument(
        "--cox_E",
        type=float,
        metavar="NUMBER",
        help="Cox only: number of cases at the end of the study",
    )
    optional.add_argument(
        "--scale",
        action="store_true",
        help="Rescale estimates (necessary for Horvath and PedBE clocks)",
    )
    optional.add_argument(
        "--age_adult",
        type=float,
        default=20,
        metavar="NUMBER",
        help="Adult age used in the --scale rescaling (default: %(default)s)",
    )

    general = parser.add_argument_group("general")
    general.add_argument(
        "--quiet",
        action="store_true",
        help="Don't print progress messages (warnings still print)",
    )
    general.add_argument(
        "--example",
        action=_ExampleAction,
        help="Run SARDINE on the small example dataset bundled with the package, print the "
        "result and exit (a quick check that the installation works)",
    )
    general.add_argument(
        "--version", action="version", version=f"sardine {__version__}"
    )
    return parser


def _format_warning(message, category, filename, lineno, line=None):
    # Show SARDINE's warnings as "Warning: [WARN] ...", without Python's file/line details
    return f"Warning: {message}\n"


def output_paths(out):
    """
    The results and metadata file paths for --out.

    If out ends with a path separator or is an existing folder, the files go inside it
    as results.tsv and metadata.tsv; otherwise out is a prefix for <out>_results.tsv and
    <out>_metadata.tsv.
    """
    separators = (os.sep, os.altsep) if os.altsep else (os.sep,)
    if out.endswith(separators) or os.path.isdir(out):
        return os.path.join(out, "results.tsv"), os.path.join(out, "metadata.tsv")
    return f"{out}_results.tsv", f"{out}_metadata.tsv"


def main(argv=None):
    """
    Run the CLI.

    Returns the exit code: 0 on success, 1 if the inputs are invalid.
    """
    args = build_parser().parse_args(argv)
    warnings.formatwarning = _format_warning

    try:
        result = sardine(
            args.ewas_sumstats,
            ewas_direction=args.ewas_direction,
            ewas_model=args.ewas_model,
            ewas_DNAm_scale=args.ewas_DNAm_scale,
            ewas_pheno_type=args.ewas_pheno_type,
            clock=args.clock,
            DNAm_cov_matrix=args.DNAm_cov_matrix,
            mean_beta=args.mean_beta,
            varY=args.varY,
            scale=args.scale,
            age_adult=args.age_adult,
            N=args.N,
            cox_pi=args.cox_pi,
            cox_E=args.cox_E,
            verbose=not args.quiet,
        )
    except SardineError as e:
        print(str(e).strip(), file=sys.stderr)
        return 1

    results_path, metadata_path = output_paths(args.out)
    out_dir = os.path.dirname(results_path)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    result.results.to_csv(results_path, sep="\t", index=False)
    result.metadata.to_csv(metadata_path, sep="\t", index=False)

    if not args.quiet:
        print(f"Wrote {results_path} and {metadata_path}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
