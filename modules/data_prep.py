### modules/data_prep.py
### Gelman (2 SD) normalization utilities shared by every analysis in this repository.
###
### Scope
###   This repository distributes an *already-prepared* analysis dataframe
###   (data/final/df_public_*.csv).  All raw-REDCap ingestion, behavioural-task
###   merging and derived-variable construction were performed upstream and are
###   not reproduced here, so this module contains only the normalization layer.
###
### What lives here
###   gelman_standardize(...)    Low-level (x - mean) / (2 * sd) transform.
###   is_gelman_normalized(...)  Predicate: is this column already on the 2 SD scale?
###   normalize_analysis_df(...) The single entry point that applies the project's
###                              normalization conventions to an analysis df before
###                              it is written out for brms.  Which convention a
###                              variable falls under is declared in
###                              VARIABLE_REGISTRY (modules/master_config.py).
###   verify_normalization(...)  Post-hoc assertion that every registry-declared
###                              column ended up on the scale it was supposed to.
### Used by
###   03_hpc/generate_hpc_jobs.py — builds df_foranalysis_master.csv for the cluster.
###
### Gotcha
###   "Gelman normalization" divides by 2 SD, not 1 SD.  A correctly normalized
###   column therefore has mean 0 and sd 0.5, not sd 1.

import os, sys
from pathlib import Path
sys.path.insert(0, os.path.dirname(__file__))

import pandas as pd
import numpy as np

from master_config import PUBLIC_DF_PREFIX


def most_recent_public_df(data_dir):
    """Return the path of the most recently modified shipped analysis dataframe.

    Looks for ``{PUBLIC_DF_PREFIX}*.csv`` directly inside *data_dir* (not
    recursively, so ``data/final/backups/`` is ignored) and returns the newest by
    modification time.  This mirrors how ``load_public_wide_df()`` in
    ``02_hgf_modeling/hgf_pipeline.py`` selects its input, so both halves of the
    repository pick up a new data export without any code edit.

    Parameters
    ----------
    data_dir : str or pathlib.Path
        The ``data/final`` directory.

    Returns
    -------
    pathlib.Path
        Path to the selected CSV.  The choice is printed so it always appears in
        the run log — worth checking when results look unexpected.

    Raises
    ------
    FileNotFoundError
        If *data_dir* does not exist or contains no matching CSV.
    """
    data_dir = Path(data_dir)
    if not data_dir.is_dir():
        raise FileNotFoundError(f'Data directory not found: {data_dir}')

    matches = [p for p in data_dir.glob(f'{PUBLIC_DF_PREFIX}*.csv')
               if 'LABELS' not in p.name]
    if not matches:
        raise FileNotFoundError(
            f'No {PUBLIC_DF_PREFIX}*.csv found in {data_dir}. '
            f'See the repository README for where data/final/ comes from.'
        )

    chosen = max(matches, key=lambda p: p.stat().st_mtime)
    if len(matches) > 1:
        print(f'  [data] {len(matches)} candidates in {data_dir}; '
              f'using most recent: {chosen.name}')
    else:
        print(f'  [data] using {chosen.name}')
    return chosen

# Binary covariates that should be modeled as factors in brms, not Gelman-normalized.
BINARY_FACTOR_VARS = [
    'amph_lifetime',
    'race_bipoc',
    'race_asian',
    'sex_v2',
    'mental_illness2_v2',
    'inhalants_lifetime',
    'coke_lifetime',
]



# ── Hardware (display-class) control covariate ───────────────────────────────


def gelman_standardize(df, cols, skip_binary=True, binary_vars=None):
    """
    Gelman 2SD normalization: (x - mean) / (2 * sd).

    Parameters
    ----------
    df : pd.DataFrame
        Input dataframe.
    cols : iterable
        Columns to normalize.
    skip_binary : bool
        If True, do not normalize columns listed in binary_vars.
    binary_vars : iterable | None
        Binary columns to leave untouched. Defaults to BINARY_FACTOR_VARS.
    """
    out = df.copy()
    binary_set = set(binary_vars if binary_vars is not None else BINARY_FACTOR_VARS)

    for col in cols:
        if col not in out.columns:
            print(f'  [gelman_standardize] Skipping {col} (missing column)')
            continue

        if skip_binary and col in binary_set:
            # Keep BRMS factor-style binaries unscaled.
            print(f'  [gelman_standardize] Leaving binary column unscaled: {col}')
            continue

        mean_val = out[col].mean()
        std_val = out[col].std()
        if pd.isna(std_val) or pd.isna(mean_val) or std_val == 0:
            print(f'  [gelman_standardize] Skipping {col} (NaN mean/std or zero variance)')
            continue

        out[col] = (out[col] - mean_val) / (2 * std_val)

    return out


# A Gelman-normalized column has mean 0 and sd 0.5.  That pair is exactly the fixed
# point of gelman_standardize(), so "normalize only if not already at it" produces the
# same values as normalizing unconditionally, without re-running the transform and
# accumulating floating-point drift.
GELMAN_ATOL = 1e-9


def is_gelman_normalized(series, atol=GELMAN_ATOL):
    """
    True if `series` is already Gelman 2SD normalized (mean 0, sd 0.5).

    Returns False for all-NaN or zero-variance columns: those cannot be normalized,
    so reporting them as normalized would hide the reason they were skipped.
    """
    s = series.dropna()
    if len(s) < 2:
        return False
    sd = s.std()
    if pd.isna(sd) or sd == 0:
        return False
    return abs(s.mean()) < atol and abs(sd - 0.5) < atol


def normalize_analysis_df(
    df,
    normalize_in_place=(),
    ensure_normalized_copy=(),
    categorical=(),
    monotonic=(),
    verbose=True,
):
    """
    Apply the project's normalization conventions to an analysis dataframe.

    This is the single place where an analysis df is Gelman-normalized before being
    written out for brms.  Which convention a variable falls under is declared in
    `VARIABLE_REGISTRY` (modules/master_config.py) and read here via the derived
    `NEED_NON_NORMALIZED` / `INPLACE_NORMALIZED` lists.

    Four rules, applied in order:

    1. `NEED_NON_NORMALIZED` — the brms family constrains the raw scale (gamma,
       beta, zero_inflated_beta, hurdle/negbinomial).  The raw column is left
       untouched and a separate `{col}_normalized` column is created.
    2. `INPLACE_NORMALIZED` — the family is an unbounded real (student_t/gaussian).
       The raw column is normalized and `{col}_normalized` is an alias of it, so
       the same name is valid in both the mediator formula and the DV formula.
    3. `normalize_in_place` — remaining continuous covariates and main predictors
       are normalized in place.  `categorical` variables are excluded (R applies
       as.factor()) as are `monotonic` ones (must stay raw integers for mo()).
    4. `ensure_normalized_copy` — every job predictor is guaranteed a
       `{col}_normalized` column.  Predictors already normalized by rules 1-3 are
       aliased; any others keep their raw column and receive a normalized copy.

    Parameters
    ----------
    df : pd.DataFrame
        Analysis dataframe. Not modified; a copy is returned.
    normalize_in_place : iterable
        Columns to Gelman-normalize in place (typically all covariates across every
        model type, plus the main predictors).
    ensure_normalized_copy : iterable
        Columns that must end up with a `{col}_normalized` counterpart (typically
        every predictor referenced by a job array).
    categorical : iterable
        Columns handled as R factors — never normalized.
    monotonic : iterable
        Ordinal columns wrapped in mo() — never normalized.
    verbose : bool
        Print a per-rule summary.

    Returns
    -------
    pd.DataFrame
        A new dataframe with normalization applied.

    Notes
    -----
    Gelman 2SD normalization is idempotent: a column already at mean 0 / sd 0.5
    maps to (x - 0) / (2 * 0.5) = x.  Calling this function on an already-
    normalized dataframe is therefore safe.
    """
    from master_config import NEED_NON_NORMALIZED, INPLACE_NORMALIZED

    def _std(d, cols):
        return gelman_standardize(d, cols, skip_binary=True, binary_vars=BINARY_FACTOR_VARS)

    def _norm_if_needed(d, cols):
        """Gelman-normalize only the columns not already at mean 0 / sd 0.5."""
        todo = [c for c in cols if c in d.columns and not is_gelman_normalized(d[c])]
        return (_std(d, todo) if todo else d), todo

    out = df.copy()
    stats = {}

    # Rule 1 — raw column preserved, normalized values in a separate column.
    need_present = [c for c in NEED_NON_NORMALIZED if c in out.columns]
    for col in need_present:
        if f'{col}_normalized' not in out.columns:
            out[f'{col}_normalized'] = out[col].copy()
    out, stats['rule1'] = _norm_if_needed(out, [f'{c}_normalized' for c in need_present])

    # Rule 2 — raw column normalized; `{col}_normalized` aliases it.
    inplace_present = [c for c in INPLACE_NORMALIZED if c in out.columns]
    out, stats['rule2'] = _norm_if_needed(out, inplace_present)
    for col in inplace_present:
        out[f'{col}_normalized'] = out[col].copy()

    # Rule 3 — everything else brms will see as a continuous term.
    in_place = sorted(
        set(normalize_in_place)
        - set(categorical) - set(monotonic)
        - set(NEED_NON_NORMALIZED) - set(INPLACE_NORMALIZED)
    )
    out, stats['rule3'] = _norm_if_needed(out, in_place)

    # Rule 4 — guarantee a `{col}_normalized` column for every job predictor.
    in_place_set = set(in_place)
    copies_made = []
    for col in ensure_normalized_copy:
        if col not in out.columns:
            raise KeyError(f'Predictor {col!r} not found in dataframe')
        if not pd.api.types.is_numeric_dtype(out[col]):
            raise TypeError(f'Non-numeric predictor {col!r}')
        norm_col = f'{col}_normalized'
        if norm_col in out.columns:
            continue
        out[norm_col] = out[col].copy()
        if col not in in_place_set:
            out, _ = _norm_if_needed(out, [norm_col])
        copies_made.append(col)

    verify_normalization(out, verbose=verbose)

    if verbose:
        print(f'  [normalize] rule 1 need_non_normalized : {len(need_present)} cols '
              f'({len(stats["rule1"])} normalized, rest already were)')
        print(f'  [normalize] rule 2 inplace_normalized  : {len(inplace_present)} cols '
              f'({len(stats["rule2"])} normalized, rest already were)')
        print(f'  [normalize] rule 3 in-place            : {len(in_place)} cols '
              f'({len(stats["rule3"])} normalized, rest already were)')
        print(f'  [normalize] rule 4 predictor copies    : {len(copies_made)} cols')

    return out


def verify_normalization(df, atol=GELMAN_ATOL, verbose=True):
    """
    Assert the normalization invariants declared in VARIABLE_REGISTRY.

    For every NEED_NON_NORMALIZED column present in `df`:
      - `{col}_normalized` exists and is Gelman-normalized (mean 0, sd 0.5)
      - the raw column is NOT normalized

    For every INPLACE_NORMALIZED column present in `df`:
      - the raw column is Gelman-normalized
      - `{col}_normalized` exists and is Gelman-normalized

    Raises ValueError listing every violation.  A raw NEED_NON_NORMALIZED column
    that looks normalized is the serious case: those variables carry a gamma, beta,
    or hurdle/negbinomial family whose support excludes the negative values that
    normalization introduces, so the column is unusable as a brms response and the
    original values cannot be recovered from the dataframe.
    """
    from master_config import NEED_NON_NORMALIZED, INPLACE_NORMALIZED

    problems = []
    checked = 0

    for col in NEED_NON_NORMALIZED:
        if col not in df.columns:
            continue
        checked += 1
        norm_col = f'{col}_normalized'
        if norm_col not in df.columns:
            problems.append(f'{norm_col}: missing (need_non_normalized requires it)')
        elif not is_gelman_normalized(df[norm_col], atol):
            s = df[norm_col].dropna()
            problems.append(
                f'{norm_col}: not normalized '
                f'(mean={s.mean():.6g}, sd={s.std():.6g}; expected mean 0, sd 0.5)'
            )
        if is_gelman_normalized(df[col], atol):
            problems.append(
                f'{col}: RAW COLUMN IS NORMALIZED — its brms family requires the raw '
                f'scale, and the raw values are no longer present in the dataframe'
            )

    for col in INPLACE_NORMALIZED:
        if col not in df.columns:
            continue
        checked += 1
        for c, role in ((col, 'raw'), (f'{col}_normalized', 'alias')):
            if c not in df.columns:
                problems.append(f'{c}: missing (inplace_normalized requires it)')
                continue
            if not is_gelman_normalized(df[c], atol):
                s = df[c].dropna()
                problems.append(
                    f'{c}: {role} column not normalized '
                    f'(mean={s.mean():.6g}, sd={s.std():.6g}; expected mean 0, sd 0.5)'
                )

    if problems:
        raise ValueError(
            'Normalization verification failed:\n  - ' + '\n  - '.join(problems)
        )
    if verbose:
        print(f'  [normalize] verified {checked} registry-declared columns '
              f'(atol={atol})')


