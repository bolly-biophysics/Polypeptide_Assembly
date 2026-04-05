import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap, LogNorm
from matplotlib.ticker import FuncFormatter

# Using the new absolute transition count matrix
transitions = np.array([
    [7228460,      176,        0,        0,        0,        0,        0,        0,        0],
    [    140,   382187,      402,        1,        0,        0,        0,        0,        0],
    [      0,      374,   587305,      257,        0,        0,        0,        0,        0],
    [      0,        1,      237,   588042,      142,        0,        0,        0,        0],
    [      0,        0,        0,      137,   139225,       55,        0,        0,        0],
    [      0,        0,        0,        0,       53,    72626,        0,        0,        0],
    [      0,        0,        0,        0,        0,        0,        0,        0,        0],
    [      0,        0,        0,        0,        0,        0,        0,        0,        0],
    [      0,        0,        0,        0,        0,        0,        0,        0,        0]
])

# State definitions (extended to S8)
STATE_NAMES = {
    0: "S0 (8H₂O)",
    1: "S1 (7H₂O+1COO⁻)",
    2: "S2 (6H₂O+2COO⁻)",
    3: "S3 (5H₂O+3COO⁻)",
    4: "S4 (4H₂O+4COO⁻)",
    5: "S5 (3H₂O+5COO⁻)",
    6: "S6 (2H₂O+6COO⁻)",
    7: "S7 (1H₂O+7COO⁻)",
    8: "S8 (0H₂O+8COO⁻)"
}

# Create DataFrame
df_transitions = pd.DataFrame(
    transitions,
    index=[f"From {STATE_NAMES[i]}" for i in range(9)],
    columns=[f"To {STATE_NAMES[i]}" for i in range(9)]
)

# Custom color gradient (white to dark green)
colors = ["#FFFFFF", "#CC4E3E"]  # white → dark green
cmap = LinearSegmentedColormap.from_list("custom_green", colors)

# Figure settings (adjusted for 9x9 matrix)
plt.figure(figsize=(14, 12))
sns.set_style("white")
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['axes.labelpad'] = 15  # increased axis label spacing

# Adjust logarithmic color scale based on data range
norm = LogNorm(vmin=1, vmax=10**7)

# Smart annotation color: auto-select black/white for contrast
def auto_text_color(val):
    return "black" if val < transitions.max()/2 else "white"

# Create heatmap
heatmap = sns.heatmap(
    df_transitions,
    annot=True,
    fmt="d",
    cmap=cmap,
    norm=norm,
    linewidths=0.5,
    linecolor="lightgray",
    cbar_kws={
        'label': 'Transition Count (log scale)',
        'ticks': [1, 10, 100, 1000, 10000, 100000, 1000000, 10000000],
        'format': FuncFormatter(lambda x, _: f"{x/1e6:.1f}M" if x >= 1e6 else f"{x/1e3:.0f}k" if x >= 1e3 else f"{x:.0f}")
    },
    annot_kws={"fontsize": 10},
    square=True
)

# Dynamically adjust annotation colors
for text in heatmap.texts:
    val = int(text.get_text())
    text.set_color(auto_text_color(val))

# Graph styling
plt.title("Ca²⁺ Hydration State Transition Matrix\n(9×9)", pad=25, fontsize=18, fontweight='bold')
plt.xlabel("To State", fontsize=14)
plt.ylabel("From State", fontsize=14)

# Adjust label rotation and font size
heatmap.set_xticklabels(heatmap.get_xticklabels(), rotation=45, ha='right', fontsize=10)
heatmap.set_yticklabels(heatmap.get_yticklabels(), rotation=0, fontsize=10)

# Save as EPS
plt.tight_layout()
plt.savefig(
    "Ca_9x9_state_transition_heatmap.eps",
    format='eps',
    dpi=300,
    bbox_inches='tight',
    facecolor='white'
)
print("Heatmap saved as Ca_9x9_state_transition_heatmap.eps")
plt.show()
