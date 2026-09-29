import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Load CSV, stripping BOM and whitespace from headers
csv_path = "delta_l3_m1_lp2_mp0-Figure.csv"
df = pd.read_csv(csv_path, encoding="utf-8-sig")
df.columns = df.columns.astype(str).str.strip().str.replace("\ufeff", "")

# Fallback: identify 'a' by name, or take column 0 if unmatched
a_matches = [c for c in df.columns if c.lower() in ["a", "param_a", "parama"]]
col_a = a_matches[0] if a_matches else df.columns[0]

# Identify A1 and A2 comparison columns
col_a1_matches = [
    c
    for c in df.columns
    if "a1" in c.lower() or ("error" in c.lower() and "1" in c)
]
col_a2_matches = [
    c
    for c in df.columns
    if "a2" in c.lower() or ("error" in c.lower() and "2" in c)
]

# Fall back to positional columns (col 1 and col 2) if header patterns differ
col_a1 = col_a1_matches[0] if col_a1_matches else df.columns[1]
col_a2 = col_a2_matches[0] if col_a2_matches else df.columns[2]

# Cast to numeric and sort
df[col_a] = pd.to_numeric(df[col_a], errors="coerce")
df[col_a1] = pd.to_numeric(df[col_a1], errors="coerce")
df[col_a2] = pd.to_numeric(df[col_a2], errors="coerce")
df = df.dropna(subset=[col_a]).sort_values(col_a)

# Convert from log10 if values are negative exponents, else keep raw errors
y1 = 10.0 ** df[col_a1] if (df[col_a1].dropna() < 0).all() else df[col_a1]
y2 = 10.0 ** df[col_a2] if (df[col_a2].dropna() < 0).all() else df[col_a2]

# Plotting setup
fig, ax = plt.subplots(figsize=(10, 5), dpi=300)

ax.plot(
    df[col_a],
    y1,
    label=r"$|\mathrm{A1} - \mathrm{Num}|\ (k_{\mathrm{max}} = 50)$",
    color="#1f77b4",
    linewidth=1.8,
    linestyle="-",
)

ax.plot(
    df[col_a],
    y2,
    label=r"$|\mathrm{A1} - \mathrm{Num}|\ (k_{\mathrm{max}} = 100)$",
    color="#ff7f0e",
    linewidth=1.8,
    linestyle="--",
)

ax.set_yscale("log")
ax.set_xlim(-1.5, 51.5)
ax.set_ylim(5e-11, 1.2e-4)

ax.grid(True, which="major", linestyle="--", linewidth=0.7, color="#d3d3d3")
ax.grid(True, which="minor", linestyle=":", linewidth=0.5, color="#e5e5e5")

ax.set_xlabel("Parameter $a$", fontsize=10, labelpad=6)
ax.set_ylabel("Absolute Error (Log Scale)", fontsize=10, labelpad=6)
ax.set_title(
    r"Absolute Error of Summation Methods relative to Ground Truth (l=3, m=1, lp=2, mp=0)",
    fontsize=12,
    pad=8,
)

ax.legend(
    loc="upper left",
    frameon=True,
    fancybox=False,
    edgecolor="#e0e0e0",
    facecolor="white",
    framealpha=0.95,
    fontsize=9,
)

plt.tight_layout()
plt.savefig("error_overview.png", dpi=300, bbox_inches="tight")
plt.show()
