import re
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Publication-quality TeX-compatible styling
plt.rcParams.update(
    {
        "font.family": "serif",
        "font.size": 11,
        "axes.labelsize": 12,
        "axes.titlesize": 12,
        "xtick.labelsize": 10,
        "ytick.labelsize": 10,
        "legend.fontsize": 9,
        "lines.linewidth": 1.5,
        "text.usetex": False,  # Switch to True if running in a TeX environment
    }
)

# 1. Read CSV, treating Excel errors like #NUM!, #N/A, blanks as NaN
df = pd.read_csv(
    "l10_m0_lp10_mp0_subse-mpmath-Figure.csv",
    na_values=["#NUM!", "#N/A", "NaN", "nan", ""],
)

# 2. Identify the parameter column 'a' and the 7 series columns
a_col = df.columns[0]
val_cols = df.columns[1:8]

# Sort sequentially by parameter 'a'
df[a_col] = pd.to_numeric(df[a_col], errors="coerce")
df = df.dropna(subset=[a_col]).sort_values(a_col)


# Helper function to format k_max labels cleanly for LaTeX legends
def format_kmax_label(col_name):
    # Extracts numerical value if column is named like 'kmax10', 'k_max_10', '10', etc.
    digits = re.findall(r"\d+", str(col_name))
    if digits:
        return rf"$k_{{\mathrm{{max}}}} = {digits[0]}$"
    return str(col_name)


# 3. Plotting setup
fig, ax = plt.subplots(figsize=(6.8, 4.4), dpi=300)
colors = plt.cm.plasma(np.linspace(0.05, 0.85, len(val_cols)))

for col, color in zip(val_cols, colors):
    series = pd.to_numeric(df[col], errors="coerce")
    valid_mask = series.notna()

    ax.plot(
        df.loc[valid_mask, a_col],
        series[valid_mask],
        label=format_kmax_label(col),
        color=color,
    )

# 4. Axes & formatting
ax.set_xlabel(r"$a$")
ax.set_ylabel(r"$\log_{10} |A_1 - M_0|$")
ax.grid(True, linestyle="--", linewidth=0.5, alpha=0.7)

ax.legend(
    loc="best",
    frameon=True,
    fancybox=False,
    edgecolor="black",
    framealpha=0.9,
    ncol=2,
)

plt.tight_layout()

# Save vector PDF directly for pdflatex / lualatex inclusion
plt.savefig("a1_m0_kmax_comparison.pdf", bbox_inches="tight")
plt.savefig("a1_m0_kmax_comparison.png", bbox_inches="tight", dpi=300)
plt.show()
