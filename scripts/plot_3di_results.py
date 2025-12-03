import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from typing import Optional
from pathlib import Path


def _add_variant_order(df, variant_col="variant", descending=False):
    """정렬용 숫자 컬럼(예: '8f','9f','10f' → 8,9,10)을 추가.

    기본 정렬은 descending=False 로 하여 8f,9f,10f 순서가 되게 함.
    """
    df = df.copy()
    # If the variant column is already a categorical with an ordered categories list,
    # respect that order and compute _variant_order from the category codes.
    try:
        import pandas as _pd
        if hasattr(df[variant_col].dtype, 'categories') and getattr(df[variant_col].dtype, 'ordered', False):
            # map categories to ordinal positions
            cat_to_idx = {c: i for i, c in enumerate(df[variant_col].dtype.categories)}
            df["_variant_order"] = df[variant_col].astype(str).map(cat_to_idx).fillna(-1).astype(int)
            return df.sort_values("_variant_order", ascending=True)
    except Exception:
        pass

    try:
        df["_variant_order"] = df[variant_col].astype(str).str.extract(r"(\\d+)").astype(float)
    except Exception:
        df["_variant_order"] = df[variant_col].astype(str)

    # build explicit ordered categories so seaborn/matplotlib respects the x-axis order
    try:
        ordered = (
            df[[variant_col, "_variant_order"]].drop_duplicates().sort_values("_variant_order", ascending=not bool(descending))[variant_col].tolist()
        )
        df[variant_col] = _pd.Categorical(df[variant_col].astype(str), categories=ordered, ordered=True)
    except Exception:
        pass

    return df.sort_values("_variant_order", ascending=not bool(descending))


def plot_pair_metrics_over_pairs(
    df: pd.DataFrame,
    metrics: list,
    pair_col: str = "pair",
    agg: str = "mean",
    figsize=(10, 5),
    save_dir: Optional[str] = None,
    show: bool = True,
):
    """
    For each metric in `metrics`, aggregate by `pair` (default mean) and draw a bar plot comparing pairs.

    Produces one PNG per metric if `save_dir` is set (files named metric.png).
    """
    try:
        import seaborn as sns
        import matplotlib.pyplot as plt
        import pandas as pd
    except Exception:
        raise RuntimeError('pandas/seaborn/matplotlib are required for plotting')

    dfp = df.copy()
    for metric in metrics:
        if metric not in dfp.columns:
            print(f"Warning: metric '{metric}' not found in data -> skipping")
            continue
        if agg == 'mean':
            agg_df = dfp.groupby(pair_col, as_index=False)[metric].mean()
        elif agg == 'median':
            agg_df = dfp.groupby(pair_col, as_index=False)[metric].median()
        else:
            agg_df = dfp.groupby(pair_col, as_index=False)[metric].agg(agg)

        plt.figure(figsize=figsize)
        sns.barplot(data=agg_df, x=pair_col, y=metric, palette='Set3')
        plt.title(f"{metric} by pair ({agg})")
        plt.ylabel(metric)
        plt.xlabel('pair')
        plt.xticks(rotation=45, ha='right')
        plt.tight_layout()
        if save_dir:
            out = Path(save_dir) / f"{metric}_by_pair.png"
            plt.savefig(out, dpi=300)
        if show:
            plt.show()
        else:
            plt.close()
        if show:
            plt.show()
        else:
            plt.close()


def plot_metric_by_variant_pairs(
    df: pd.DataFrame,
    metric: str,
    variant_col: str = "variant",
    pair_col: str = "pair",
    w_col: str = "w",
    w_value: Optional[int] = 1,
    figsize=(8, 4),
    save_path: Optional[str] = None,
    show: bool = True,
):
    """
    Grouped bar chart: x=variant (ordered), one bar per pair for the given metric.

    By default uses w == w_value rows; pass w_value=None to use all rows but then
    the function will aggregate by pair/variant (mean) which the user said they
    didn't want — so default is w=1.
    """
    try:
        import seaborn as sns
        import matplotlib.pyplot as plt
        import pandas as _pd
    except Exception:
        raise RuntimeError('pandas/seaborn/matplotlib are required for plotting')

    dfp = df.copy()
    if metric not in dfp.columns:
        raise ValueError(f"Metric '{metric}' not in dataframe")

    if w_value is not None:
        dfp = dfp[dfp[w_col] == w_value]
    if dfp.empty:
        raise ValueError(f"No rows available for w={w_value}")

    # --- if metric is an alignment length, add TM baseline as 'tm' variant
    # detect tm column candidates
    if 'alignment_len' in metric or 'alignment_len' in metric.lower():
        # candidate TM length columns
        tm_candidates = ['len_Ttm', 'len_Qtm', 'len_target', 'len_T3']
        found_tm = next((c for c in tm_candidates if c in dfp.columns), None)
        if found_tm is not None:
            # build one-row per pair with variant='tm' and metric value equal to found_tm
            tm_rows = dfp.groupby(pair_col, as_index=False)[found_tm].max()
            tm_rows = tm_rows.rename(columns={found_tm: metric})
            tm_rows[variant_col] = 'tm'
            # append and continue
            dfp = pd.concat([dfp, tm_rows], ignore_index=True, sort=False)

    # build explicit variant ordering and ensure 'tm' is placed just left of the smallest numeric variant (typically 8f)
    # collect present variants
    variants_present = dfp[variant_col].astype(str).unique().tolist()
    # extract numeric variants and sort by numeric part descending
    import re
    numeric_variants = [v for v in variants_present if re.search(r"(\d+)", v)]
    def _num(v):
        m = re.search(r"(\d+)", v)
        return int(m.group(1)) if m else -1
    ordered_numeric = sorted(numeric_variants, key=_num)
    ordered = ordered_numeric.copy()
    # insert 'tm' at the leftmost position so it appears first on the x-axis
    if 'tm' in variants_present:
        ordered.insert(0, 'tm')

    try:
        dfp[variant_col] = pd.Categorical(dfp[variant_col].astype(str), categories=ordered, ordered=True)
    except Exception:
        pass

    # now plot without mentioning w in title (user requested no 'w=1' text)
    plt.figure(figsize=figsize)
    # draw all variants (tm is included and will appear leftmost because we inserted it at index 0)
    ax = sns.barplot(data=dfp, x=variant_col, y=metric, hue=pair_col)
    plt.title(f"{metric} per pair by variant")
    plt.ylabel(metric)
    plt.xlabel('variant')
    plt.ylim(bottom=0)
    plt.legend(title='pair', bbox_to_anchor=(1.05,1), loc='upper left')
    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=300)
    if show:
        plt.show()
    else:
        plt.close()


def plot_sequence_lengths_by_pair(
    df: pd.DataFrame,
    pair_col: str = "pair",
    target_len_cols: list = None,
    query_len_cols: list = None,
    agg: str = "max",
    figsize=(10, 5),
    save_path: Optional[str] = None,
    show: bool = True,
):
    """
    Draw a side-by-side barplot of sequence lengths per pair (target vs query).

    The function will try to find sensible length columns automatically. By default
    it looks for ['len_target', 'len_T3', 'len_Ttm'] for the target and
    ['len_Q3', 'len_Qtm', 'len_query'] for the query. If a column is missing it
    will be skipped.

    Aggregation default is 'max' (sequence lengths should be identical across
    variants/w, use max to be robust). Other supported agg: 'mean', 'median'.
    """
    try:
        import seaborn as sns
        import matplotlib.pyplot as plt
        import pandas as _pd
    except Exception:
        raise RuntimeError('pandas/seaborn/matplotlib are required for plotting')

    dfp = df.copy()

    # candidate columns
    if target_len_cols is None:
        target_len_cols = ['len_target', 'len_T3', 'len_Ttm']
    if query_len_cols is None:
        query_len_cols = ['len_Q3', 'len_Qtm', 'len_query', 'len_Q']

    found_target = next((c for c in target_len_cols if c in dfp.columns), None)
    found_query = next((c for c in query_len_cols if c in dfp.columns), None)

    if found_target is None and found_query is None:
        raise ValueError('No suitable length columns found in dataframe')

    agg_func = 'max' if agg == 'max' else ('mean' if agg == 'mean' else 'median')

    # build agg mapping
    agg_map = {}
    if found_target:
        agg_map[found_target] = agg_func
    if found_query:
        agg_map[found_query] = agg_func

    grp = dfp.groupby([pair_col], as_index=False).agg(agg_map)

    # normalize column names for plotting
    plot_df = grp.rename(columns={found_target: 'target_len'} if found_target else {})
    if found_query:
        plot_df = plot_df.rename(columns={found_query: 'query_len'})

    # melt for seaborn barplot
    melt_cols = [c for c in ['target_len', 'query_len'] if c in plot_df.columns]
    if not melt_cols:
        raise ValueError('No length columns available to plot after detection')

    df_melt = plot_df.melt(id_vars=[pair_col], value_vars=melt_cols, var_name='side', value_name='length')

    plt.figure(figsize=figsize)
    sns.barplot(data=df_melt, x=pair_col, y='length', hue='side', palette='Set2')
    plt.title('Sequence lengths by pair (target vs query)')
    plt.ylabel('length (aa/tokens)')
    plt.xlabel('pair')
    plt.xticks(rotation=45, ha='right')
    plt.legend(title='side')
    plt.tight_layout()

    if save_path:
        out = Path(save_path)
        plt.savefig(out, dpi=300)
    if show:
        plt.show()
    else:
        plt.close()

    return df_melt

def plot_metric_variant_lines(
    df: pd.DataFrame,
    metric: str,
    variant_col: str = "variant",
    pair_col: str = "pair",
    w_col: str = "w",
    figsize=(8, 4),
    save_path: Optional[str] = None,
    show: bool = True,
):
    """
    Line plot across variants: one line per pair (hue=pair). If multiple w exist
    the marker/style will reflect w values (style=w) so you can see per-w
    variation. This avoids crowded bar groups.
    """
    try:
        import seaborn as sns
        import matplotlib.pyplot as plt
        import pandas as _pd
    except Exception:
        raise RuntimeError('pandas/seaborn/matplotlib are required for plotting')

    dfp = df.copy()
    if metric not in dfp.columns:
        raise ValueError(f"Metric '{metric}' not in dataframe")

    dfp = _add_variant_order(dfp, variant_col=variant_col, descending=False)

    plt.figure(figsize=figsize)
    ax = sns.lineplot(
        data=dfp,
        x=variant_col,
        y=metric,
        hue=pair_col,
        style=w_col if w_col in dfp.columns else None,
        markers=True,
        dashes=False,
    )
    plt.title(f"{metric} by variant (lines per pair)")
    plt.ylabel(metric)
    plt.xlabel('variant')
    plt.ylim(bottom=0)
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles=handles, labels=labels, title='pair / w', bbox_to_anchor=(1.05,1), loc='upper left')
    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=300)
    if show:
        plt.show()
    else:
        plt.close()


def plot_pair_f1_all(
    df: pd.DataFrame,
    f1_col: str = "f1",
    variant_col: str = "variant",
    pair_col: str = "pair",
    w_col: str = "w",
    aggregate: str = "mean",   # "mean", "median", or None
    figsize=(7, 5),
    title: str = "Pair-level F1 by variant (color=pair, marker=w)",
    save_path: Optional[str] = None,
    show: bool = True,
    per_w: bool = False,
):
    """
    한 플롯에 모든 pair와 w를 동시에 그리는 pair-level F1 plot.

    - x축: variant (예: 8f, 9f, 10f)
    - y축: pair-level F1 (f1_col)
    - 색(hue): pair
    - marker(style): w (0,1,2,...)
    - CSV가 여러 개여도 df에 다 concat만 되어 있으면 OK.

    df에 필요한 최소 컬럼:
        [variant_col, pair_col, w_col, f1_col]
    """

    df_plot = _add_variant_order(df, variant_col=variant_col, descending=False)

    if aggregate is not None:
        agg_func = "mean" if aggregate == "mean" else "median"
        # use observed=True to avoid including unused categorical combinations
        df_plot = (
            df_plot
            .groupby([variant_col, pair_col, w_col, "_variant_order"], as_index=False, observed=True)
            .agg({f1_col: agg_func})
        )

    def _draw(df_to_plot, out_path=None, title_suffix=None):
        plt.figure(figsize=figsize)
        ax = sns.lineplot(
            data=df_to_plot,
            x=variant_col,
            y=f1_col,
            hue=pair_col,    # 색: pair
            style=w_col,     # marker/linestyle: w
            markers=True,
            dashes=False,
        )
        ax.set_ylim(0.0, 1.0)
        ax.set_ylabel("Pair-level F1")
        ax.set_xlabel("3Di variant")
        title_full = title if not title_suffix else f"{title} ({title_suffix})"
        ax.set_title(title_full)

        handles, labels = ax.get_legend_handles_labels()
        ax.legend(handles=handles, labels=labels, title="pair / w", bbox_to_anchor=(1.05, 1), loc="upper left")

        plt.tight_layout()
        if out_path is not None:
            plt.savefig(out_path, dpi=300)
        if show:
            plt.show()
        else:
            plt.close()
        return ax

    # if per_w requested, plot one figure per w value
    if per_w:
        w_vals = sorted(df_plot[w_col].dropna().unique())
        last_ax = None
        for w in w_vals:
            df_w = df_plot[df_plot[w_col] == w]
            out = None
            if save_path is not None:
                # insert _w{w} before extension
                if save_path.lower().endswith('.png'):
                    out = save_path[:-4] + f"_w{w}.png"
                else:
                    out = save_path + f"_w{w}.png"
            last_ax = _draw(df_w, out_path=out, title_suffix=f"w={w}")
        return last_ax
    else:
        return _draw(df_plot, out_path=save_path)


def plot_block_f1_w1_all(
    df: pd.DataFrame,
    range_target_col: str = "range_target_f1",
    range_query_col: str = "range_query_f1",
    variant_col: str = "variant",
    pair_col: str = "pair",
    w_col: str = "w",
    w_value: int = 1,
    aggregate: str = "mean",
    per_w: bool = False,
    figsize=(7, 5),
    title: str = "Block-level F1 by variant",
    save_path: Optional[str] = None,
    show: bool = True,
):
    """
    Block-level F1 plot.

    - per_w=False:
        * df에서 w_col == w_value 인 행만 사용
        * x: variant, y: block_f1, hue: pair (한 장)

    - per_w=True:
        * df[w_col].unique() 각각에 대해
          df[df[w_col] == w_val] 로 필터링해서
          w별로 다른 그림을 그림 (파일명 뒤에 _w{w_val}.png 붙임)
    """

    # block_f1를 먼저 row 단위로 계산
    df = df.copy()
    df["block_f1"] = 0.5 * (df[range_target_col] + df[range_query_col])

    # -----------------------------
    # 1) per_w = False : 특정 w만 사용
    # -----------------------------
    if not per_w:
        df_w = df[df[w_col] == w_value].copy()
        if df_w.empty:
            raise ValueError(f"No rows with {w_col} == {w_value}")

        df_w = _add_variant_order(df_w, variant_col=variant_col)

        if aggregate is not None:
            agg_func = "mean" if aggregate == "mean" else "median"
            df_w = (
                df_w
                .groupby([variant_col, pair_col, "_variant_order"], as_index=False, observed=True)
                .agg(block_f1=("block_f1", agg_func))
            )

        plt.figure(figsize=figsize)
        ax = sns.lineplot(
            data=df_w,
            x=variant_col,
            y="block_f1",
            hue=pair_col,
            marker="o",
        )
        ax.set_ylim(0.0, 1.0)
        ax.set_ylabel(f"Block-level F1 (w={w_value})")
        ax.set_xlabel("3Di variant")
        ax.set_title(title if title else f"Block-level F1 (w={w_value}) by variant")

        handles, labels = ax.get_legend_handles_labels()
        ax.legend(handles=handles, labels=labels, title="pair", bbox_to_anchor=(1.05, 1), loc="upper left")

        plt.tight_layout()

        if save_path is not None:
            plt.savefig(save_path, dpi=300)
        if show:
            plt.show()
        else:
            plt.close()

        return ax

    # -----------------------------
    # 2) per_w = True : 모든 w에 대해 각각 그림
    # -----------------------------
    axes = {}
    unique_ws = sorted(df[w_col].unique())
    for w_val in unique_ws:
        df_w = df[df[w_col] == w_val].copy()
        if df_w.empty:
            continue

        df_w = _add_variant_order(df_w, variant_col=variant_col)

        if aggregate is not None:
            agg_func = "mean" if aggregate == "mean" else "median"
            df_w = (
                df_w
                .groupby([variant_col, pair_col, "_variant_order"], as_index=False, observed=True)
                .agg(block_f1=("block_f1", agg_func))
            )

        plt.figure(figsize=figsize)
        ax = sns.lineplot(
            data=df_w,
            x=variant_col,
            y="block_f1",
            hue=pair_col,
            marker="o",
        )
        ax.set_ylim(0.0, 1.0)
        ax.set_ylabel(f"Block-level F1 (w={w_val})")
        ax.set_xlabel("3Di variant")
        this_title = f"Block-level F1 (w={w_val}) by variant"
        ax.set_title(this_title)

        handles, labels = ax.get_legend_handles_labels()
        ax.legend(handles=handles, labels=labels, title="pair", bbox_to_anchor=(1.05, 1), loc="upper left")

        plt.tight_layout()

        if save_path is not None:
            # per_w=True 이므로 파일명을 w별로 분리
            out_path = save_path
            if out_path.endswith(".png"):
                out_path = out_path.replace(".png", f"_w{w_val}.png")
            else:
                out_path = out_path + f"_w{w_val}.png"
            plt.savefig(out_path, dpi=300)
        if show:
            plt.show()
        else:
            plt.close()

        axes[w_val] = ax

    return axes
"""
plot_3di_results.py

Plotting helpers for tmalign_3di_match_pipeline results.
This module provides three plotting functions (using pandas/seaborn/matplotlib):
 - plot_pair_f1
 - plot_block_f1_vs_w
 - plot_pair_f1_mean_via_w

These functions accept either a pandas.DataFrame or a path to a CSV produced by
`tmalign_3di_match_pipeline.py` (the variant summary CSV).

They are intentionally self-contained so consumers can call them from a shell
script or notebook.
"""
from pathlib import Path
from typing import Union


def _load_dataframe(df_or_path):
    """Accept either a pandas DataFrame or a path to CSV and return a DataFrame."""
    try:
        import pandas as pd
    except Exception:
        raise RuntimeError('pandas is required for plotting functions')
    if isinstance(df_or_path, pd.DataFrame):
        return df_or_path.copy()
    p = df_or_path
    if isinstance(p, (str, Path)):
        return pd.read_csv(str(p))
    raise ValueError('df_or_path must be a pandas.DataFrame or path to CSV')


def plot_pair_f1(df_or_path, pair=None, save_path=None, show=True):
    """Barplot of F1 per variant for a single pair.

    Parameters:
    - df_or_path: pandas DataFrame or path to CSV with at least columns ['pair','variant','f1']
    - pair: pair id string like 'A_vs_B'. If None and CSV contains a single pair, that pair is used.
    - save_path: if provided, save figure to this path.
    - show: whether to plt.show()
    """
    df = _load_dataframe(df_or_path)
    try:
        import seaborn as sns
        import matplotlib.pyplot as plt
    except Exception:
        raise RuntimeError('seaborn and matplotlib are required for plotting')
    if 'pair' in df.columns:
        if pair is None:
            pairs = df['pair'].unique()
            if len(pairs) == 1:
                pair = pairs[0]
            else:
                raise ValueError('Multiple pairs found in data; pass pair parameter')
        dfp = df[df['pair'] == pair]
    else:
        dfp = df.copy()
    if dfp.empty:
        raise ValueError('No data for requested pair')
    plt.figure(figsize=(6,4))
    sns.barplot(data=dfp, x='variant', y='f1', palette='Set2')
    plt.title(f'F1 by variant ({pair})' if pair else 'F1 by variant')
    plt.ylabel('F1')
    plt.xlabel('Variant')
    plt.ylim(0,1)
    plt.tight_layout()
    if save_path:
        plt.savefig(save_path)
    if show:
        plt.show()
    plt.close()


def plot_block_f1_vs_w(df_or_path, w_col='w', variant=None, save_path=None, show=True):
    """Line plot of block-level metric (block_recall_tau) vs w.

    Expects a CSV/DataFrame containing columns ['w','block_recall_tau'] and optionally 'variant' or 'pair'.
    """
    df = _load_dataframe(df_or_path)
    try:
        import seaborn as sns
        import matplotlib.pyplot as plt
    except Exception:
        raise RuntimeError('seaborn and matplotlib are required for plotting')
    if w_col not in df.columns or 'block_recall_tau' not in df.columns:
        raise ValueError(f"Data must contain columns '{w_col}' and 'block_recall_tau'.")
    if variant and 'variant' in df.columns:
        df = df[df['variant'] == variant]
    plt.figure(figsize=(6,4))
    hue = 'variant' if 'variant' in df.columns else ('pair' if 'pair' in df.columns else None)
    import seaborn as sns
    sns.lineplot(data=df, x=w_col, y='block_recall_tau', hue=hue, marker='o')
    plt.title('Block recall (tau) vs window w')
    plt.ylabel('block_recall_tau')
    plt.xlabel(w_col)
    plt.ylim(0,1)
    plt.tight_layout()
    if save_path:
        plt.savefig(save_path)
    if show:
        plt.show()
    plt.close()


def plot_pair_f1_mean_via_w(df_or_path, w_col='w', save_path=None, show=True):
    """Plot mean F1 (across variants) for each pair as function of w.

    Expects columns ['w','pair','f1'] in the CSV/DataFrame. Aggregates mean f1 per (pair,w) and plots.
    """
    df = _load_dataframe(df_or_path)
    try:
        import seaborn as sns
        import matplotlib.pyplot as plt
        import pandas as pd
    except Exception:
        raise RuntimeError('pandas, seaborn and matplotlib are required for plotting')
    if w_col not in df.columns or 'pair' not in df.columns or 'f1' not in df.columns:
        raise ValueError("Data must contain columns 'w', 'pair', and 'f1'")
    grp = df.groupby([w_col, 'pair'], as_index=False)['f1'].mean()
    plt.figure(figsize=(8,5))
    sns.lineplot(data=grp, x=w_col, y='f1', hue='pair', marker='o')
    plt.title('Mean F1 per pair vs window w')
    plt.ylabel('mean F1')
    plt.xlabel(w_col)
    plt.ylim(0,1)
    plt.tight_layout()
    if save_path:
        plt.savefig(save_path)
    if show:
        plt.show()
    plt.close()
