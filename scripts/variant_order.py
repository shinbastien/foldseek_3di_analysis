from typing import List
import pandas as pd

# Canonical variant order (left-to-right on x-axis). 'tm' first, then numeric variants ascending.
CANONICAL_ORDER: List[str] = ["tm", "8f", "9f", "10f"]


def apply_variant_order(df: pd.DataFrame, variant_col: str = "variant") -> pd.DataFrame:
    """Return a copy of df where `variant_col` is converted to an ordered Categorical
    using the canonical order. Any variants not in CANONICAL_ORDER are appended
    after the canonical list in their appearance order.

    This is safe to call multiple times; if the column is already an ordered
    categorical with the same categories, it will be left unchanged.
    """
    df = df.copy()
    vals = list(dict.fromkeys(df[variant_col].astype(str).tolist()))
    # build ordered list: include canonical if present, then any extras
    ordered = [v for v in CANONICAL_ORDER if v in vals]
    extras = [v for v in vals if v not in ordered]
    ordered.extend(extras)

    try:
        df[variant_col] = pd.Categorical(df[variant_col].astype(str), categories=ordered, ordered=True)
    except Exception:
        pass
    return df
