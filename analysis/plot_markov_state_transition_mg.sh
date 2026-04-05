import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap, LogNorm
from matplotlib.ticker import FuncFormatter

# Updated 7x7 transfer matrix data
transitions = np.array([
    [41106525, 23, 0, 0, 0, 0, 0],
    [2, 3517551, 6, 0, 0, 0, 0],
    [0, 2, 375711, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0]
])

# State definitions
STATE_NAMES = {
    0: "S0 (6H₂O)",
    1: "S1 (5H₂O+1COO⁻)",
    2: "S2 (4H₂O+2COO⁻)",
    3: "S3 (3H₂O+3COO⁻)",
    4: "S4 (2H₂O+4COO⁻)",
    5: "S5 (1H₂O+5COO⁻)",
    6: "S6 (0H₂O+6COO⁻)"
}

# Create DataFrame
df_transitions = pd.DataFrame(
    transitions,
    index=[f"From {STATE_NAMES[i]}" for i in range(7)],
    columns=[f"To {STATE_NAMES[i]}" for i in range(7)]
)

# Custom dark green → white gradient
colors = ["#FFFFFF", "#007239"]  # pure white → dark green
cmap = LinearSegmentedColormap.from_list("custom_green", colors)

# Figure settings
plt.figure(figsize=(12, 10))  # larger canvas for 7x7 matrix
sns.set_style("white")
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['axes.labelpad'] = 15  # increase axis label spacing

# Log color scale (1 to 10^8 to accommodate the larger values)
norm = LogNorm(vmin=1, vmax=10**8)

# Smart label coloring: auto-select black/white for contrast
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
        'ticks': [1, 10, 100, 1000, 10000, 100000, 1000000, 10000000, 100000000],
        'format': FuncFormatter(lambda x, _: f"{x/1e6:.1f}M" if x >= 1e6 else f"{x/1e3:.0f}k" if x >= 1e3 else f"{x:.0f}")
    },
    annot_kws={"fontsize": 10},
    square=True
)

# Dynamic label color adjustment
for text in heatmap.texts:
    val = int(text.get_text())
    text.set_color(auto_text_color(val))

# Style improvements
plt.title("Mg²⁺ Hydration State Transition Matrix\n(7×7)", pad=25, fontsize=18, fontweight='bold')
plt.xlabel("To State", fontsize=14)
plt.ylabel("From State", fontsize=14)
heatmap.set_xticklabels(heatmap.get_xticklabels(), rotation=45, ha='right', fontsize=11)
heatmap.set_yticklabels(heatmap.get_yticklabels(), rotation=0, fontsize=11)

# Save as EPS
plt.tight_layout()
plt.savefig(
    "7x7_state_transition_heatmap.eps",
    format='eps',
    dpi=300,
    bbox_inches='tight',
    facecolor='white'
)
print("Heatmap saved as 7x7_state_transition_heatmap.eps")
plt.show()
