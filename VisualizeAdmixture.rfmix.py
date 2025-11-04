

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# -------------------------
# User-defined file paths
# -------------------------
clumpp_file = "MergeDEMI_OKA.filtered.ref_sample_2.rfmix.allKhwe.chrAll.Q.txt"  # Averaged CLUMPP file (one row per individual)
pop_file = "MergeDEMI_OKA.filtered.ref_sample_2.rfmix.allKhwe.chrAll.pop.txt"          # Population file: one population per line, in order
output_img = "Figure_1_rfmix_ancestry_plot.grey_line_solid.pdf"     # Output image file

# -------------------------
# Step 1: Read the CLUMPP file
# -------------------------
# The CLUMPP file is assumed to have no header and K+5 columns per row.
# Columns: [FID, IndID, IID, dummy1, dummy2, Cluster1, Cluster2, ..., ClusterK]
df = pd.read_csv(clumpp_file, delim_whitespace=True, header=None)

# Determine the number of clusters (K)
num_cols = df.shape[1]


# Rename columns for clarity
df.columns = ["sample_name", "Bantu", "Khoisan", "EA"]
cluster_cols = ["Bantu", "Khoisan", "EA"]
# -------------------------
# Step 2: Read the Population File
# -------------------------
# The population file should have one population name per line.
pop_df = pd.read_csv(pop_file, header=None, names=["Population"])

print(df.shape[0])
print(pop_df.shape[0])
# Verify that the number of individuals matches
if df.shape[0] != pop_df.shape[0]:
    raise ValueError("The number of rows in the CLUMPP file does not match the number of lines in the population file.")

# Assign population labels by row order
df["Population"] = pop_df["Population"]

# -------------------------
# Step 3: Sort by Population for Grouping
# -------------------------
#df.sort_values(by="Population", inplace=True)
df.reset_index(drop=True, inplace=True)

population_averages = df.groupby("Population")[["Bantu"]].mean()


print(population_averages)

population_averages = df.groupby("Population")[["Khoisan"]].mean()


print(population_averages)

population_averages = df.groupby("Population")[["EA"]].mean()


print(population_averages)


# -------------------------
# Step 4: Create the Stacked Bar Plot
# -------------------------
plt.figure(figsize=(12, 6))
# Remove the legend by setting legend=False
ax = df[cluster_cols].plot(kind="bar", stacked=True, width=1.1, edgecolor="none", legend=False, figsize=(12, 6), color = ("red","orange","darkblue"))

for p in ax.patches:
    p.set_linewidth(0)
    p.set_edgecolor((0, 0, 0, 0))  # fully transparent

#plt.xlabel("Individuals grouped by Population label")
plt.ylabel("Individual ancestry proportion")
#plt.title("Rfmix K=3 Juhoan/Kwangali/Somali")

# Remove individual tick labels
ax.set_xticks([])

# -------------------------
# Step 5: Add Population Labels & Vertical Separators
# -------------------------
# Group by population to determine where to place labels and separators.
grouped = df.groupby("Population", sort=False)

tick_positions = []
tick_labels = []
boundaries = []

# For each population, compute the average index (for label) and the maximum index (for separator)
for pop_name, group in grouped:
    indices = group.index.to_numpy()
    tick_positions.append(np.mean(indices))
    tick_labels.append(pop_name)
    boundaries.append(indices.max())

# Set one tick label per population group at the average index position.
ax.set_xticks(tick_positions)
ax.set_xticklabels(tick_labels, rotation=0, ha='center', fontsize=12)

# Draw vertical dashed lines between groups (i.e. after each population except the last).
for boundary in boundaries[:-1]:
    ax.axvline(x=boundary + 0.5, color='grey', linestyle='-', linewidth=1)

# Increase bottom margin to ensure labels have enough room.
plt.subplots_adjust(bottom=0.25)
plt.tight_layout()

# -------------------------
# Step 6: Save the Plot to a File and Display It
# -------------------------
plt.savefig(output_img, dpi=300, bbox_inches="tight")
print(f"Plot saved to {output_img}")
plt.show()