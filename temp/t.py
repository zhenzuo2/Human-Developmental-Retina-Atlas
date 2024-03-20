sc.pp.normalize_total(adata, target_sum=1e6)
sc.pp.log1p(adata)


markers = ['ALDH1A1', 'ALDH1A2','ALDH1A3']
df = sc.get.obs_df(adata, markers + ["Region","Days"])
grouped_df = df.groupby(["Region","Days"])
result = grouped_df.mean().reset_index()

# Set the 'Date' column as the index
result.set_index('Days', inplace=True)

# Plot multiple lines using plt
plt.clf()
fig = plt.gcf()
fig.set_size_inches(8, 8)
plt.plot(result.index, result['ALDH1A1'], marker='o', linestyle='-', label='ALDH1A1')
plt.plot(result.index, result['ALDH1A2'], marker='s', linestyle='-', label='ALDH1A2')
plt.plot(result.index, result['ALDH1A3'], marker='^', linestyle='-', label='ALDH1A3')

# Customize the plot (optional)
plt.title('')
plt.xlabel('')
plt.ylabel('log CPM')
plt.grid(True)
plt.legend()  # Add legend to distinguish lines
plt.ylim(0, 1.5)
plt.tight_layout()
plt.savefig(
    "path to save",
    dpi=600,
    bbox_inches="tight",
    transparent=True,
)