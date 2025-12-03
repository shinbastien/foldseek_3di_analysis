"""
plot_pair_block.py

Two plotting functions requested by the user:
- plot_pair_f1_by_variant_and_w
- plot_block_f1_w1_by_variant

Usage:
    import pandas as pd
    from scripts.plot_pair_block import plot_pair_f1_by_variant_and_w, plot_block_f1_w1_by_variant

    df = pd.read_csv('combined_variant_summary.csv')
    plot_pair_f1_by_variant_and_w(df)
    plot_block_f1_w1_by_variant(df, w_value=1)

"""
from typing import Optional
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns


def plot_pair_f1_by_variant_and_w(
    df: pd.DataFrame,
    figsize=(6, 4),
    aggregate: Optional[str] = "mean",
    title: str = "Pair-level F1 by variant and tolerance w",
):
    """
    Line plot of pair-level F1.
    - x축: variant (8f / 9f / 10f ...)
    - y축: pair-level F1
    - 색(hue): tolerance w (0 / 1 / 2 / ...)
    - 여러 pair가 있으면 variant×w 별로 평균(or median)을 그려줌.

    Required columns in df:
        ['pair', 'variant', 'w', 'pair_f1']
    """

    df_plot = df.copy()

    # 숫자 순서 정렬을 위해 variant 정렬 (원하는 순서: 8f, 9f, 10f -> 오름차순)
    # If variant is already an ordered categorical (set globally), keep it
    try:
        import pandas as _pd
        if hasattr(df_plot['variant'].dtype, 'categories') and getattr(df_plot['variant'].dtype, 'ordered', False):
            df_plot['variant_order'] = df_plot['variant'].astype(str).map({c:i for i,c in enumerate(df_plot['variant'].dtype.categories)}).fillna(-1).astype(float)
        else:
            df_plot["variant_order"] = df_plot["variant"].astype(str).str.extract(r"(\d+)")[0].astype(float)
            # build ordered categories ascending (8,9,10...)
            ordered_variants = (
                df_plot[['variant','variant_order']].drop_duplicates().sort_values('variant_order', ascending=True)['variant'].tolist()
            )
            df_plot['variant'] = _pd.Categorical(df_plot['variant'], categories=ordered_variants, ordered=True)
    except Exception:
        df_plot["variant_order"] = df_plot["variant"].astype(str)

    # 여러 pair가 들어왔을 때는 variant×w 단위로 aggregate
    if aggregate is not None:
        if aggregate == "mean":
            df_agg = (
                df_plot.groupby(["variant", "variant_order", "w"], as_index=False)
                .agg(pair_f1=("pair_f1", "mean"))
            )
        elif aggregate == "median":
            df_agg = (
                df_plot.groupby(["variant", "variant_order", "w"], as_index=False)
                .agg(pair_f1=("pair_f1", "median"))
            )
        else:
            raise ValueError("aggregate must be 'mean', 'median', or None")
    else:
        df_agg = df_plot

    plt.figure(figsize=figsize)
    ax = sns.lineplot(
        data=df_agg,
        x="variant",
        y="pair_f1",
        hue="w",
        marker="o",
    )
    ax.set_ylim(0.0, 1.0)
    ax.set_ylabel("Pair-level F1")
    ax.set_xlabel("3Di variant")
    ax.set_title(title)
    ax.legend(title="w (tolerance)")
    plt.tight_layout()
    return ax


def plot_block_f1_w1_by_variant(
    df: pd.DataFrame,
    w_value: int = 1,
    figsize=(5, 4),
    aggregate: Optional[str] = "mean",
    title: str = "Block-level F1 (w=1) by variant",
):
    """
    Bar plot of block-level F1 at a fixed w (default w=1).
    - x축: variant
    - y축: block-level F1
    - 여러 pair가 있으면 variant별로 평균(or median)을 그림.

    Required columns in df:
        ['pair', 'variant', 'w', 'block_f1']
    """

    df_w = df[df["w"] == w_value].copy()
    if df_w.empty:
        raise ValueError(f"No rows with w == {w_value}")

    # 정렬용 (8f,9f,10f 순서)
    try:
        import pandas as _pd
        if hasattr(df_w['variant'].dtype, 'categories') and getattr(df_w['variant'].dtype, 'ordered', False):
            df_w['variant_order'] = df_w['variant'].astype(str).map({c:i for i,c in enumerate(df_w['variant'].dtype.categories)}).fillna(-1).astype(float)
        else:
            df_w["variant_order"] = df_w["variant"].astype(str).str.extract(r"(\d+)")[0].astype(float)
            ordered_variants = (
                df_w[['variant','variant_order']].drop_duplicates().sort_values('variant_order', ascending=True)['variant'].tolist()
            )
            df_w['variant'] = _pd.Categorical(df_w['variant'], categories=ordered_variants, ordered=True)
    except Exception:
        df_w["variant_order"] = df_w["variant"].astype(str)

    plt.figure(figsize=figsize)

    if aggregate is not None:
        # seaborn이 estimator로 집계하게 하는 방식
        estimator = np.mean if aggregate == "mean" else np.median
        ax = sns.barplot(
            data=df_w,
            x="variant",
            y="block_f1",
            estimator=estimator,
            ci="sd",            # error bar: 표준편차 (원하면 "sd" → None)
        )
    else:
        # aggregate 안 하고 raw 값들을 그대로 찍고 싶으면 stripplot을 섞을 수도 있음
        ax = sns.barplot(
            data=df_w,
            x="variant",
            y="block_f1",
            estimator=np.mean,
            ci="sd",
            color="lightgray",
        )
        sns.stripplot(
            data=df_w,
            x="variant",
            y="block_f1",
            dodge=True,
            alpha=0.7,
            ax=ax,
        )

    ax.set_ylim(0.0, 1.0)
    ax.set_ylabel(f"Block-level F1 (w={w_value})")
    ax.set_xlabel("3Di variant")
    ax.set_title(title)
    plt.tight_layout()
    return ax
