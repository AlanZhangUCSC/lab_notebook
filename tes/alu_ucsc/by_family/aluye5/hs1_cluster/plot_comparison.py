import numpy as np
import pandas as pd
from scipy import stats
import matplotlib.pyplot as plt

samples = {
  'hg38': pd.read_csv('hg38.stats.tsv', sep='\t'),
  'hs1_outcluster': pd.read_csv('hs1_outcluster.stats.tsv', sep='\t'),
  'hs1_incluster': pd.read_csv('hs1_incluster.stats.tsv', sep='\t'),
}

names = list(samples.keys())

for policy, col in [('A (clips ignored)', 'id_A'),
                    ('B (clips penalized)', 'id_B')]:
  print(f"\n{'=' * 60}\nPolicy {policy}\n{'=' * 60}")

  for name, df in samples.items():
    x = df[col].values
    q = np.quantile(x, [0.05, 0.25, 0.5, 0.75, 0.95])
    print(f"{name}: n={len(x)}, mean={x.mean():.4f}, "
          f"q05={q[0]:.4f}, median={q[2]:.4f}, q95={q[4]:.4f}")

  arrays = [samples[n][col].values for n in names]

  ad = stats.anderson_ksamp(arrays)
  print(f"\nAnderson-Darling k-sample: stat={ad.statistic:.4f}, p={ad.pvalue:.2e}")

  print("Pairwise:")
  for i in range(len(names)):
    for j in range(i + 1, len(names)):
      ks = stats.ks_2samp(arrays[i], arrays[j])
      w = stats.wasserstein_distance(arrays[i], arrays[j])
      print(f"  {names[i]} vs {names[j]}: D={ks.statistic:.4f}, "
            f"p={ks.pvalue:.2e}, W={w:.5f}")

print(f"\n{'=' * 60}\nClip rate diagnostics\n{'=' * 60}")
for name, df in samples.items():
  clip_rate = (df['soft'] + df['hard']) / df['query_len']
  q = np.quantile(clip_rate, [0.5, 0.75, 0.95, 0.99])
  print(f"{name}: median={q[0]:.4f}, q75={q[1]:.4f}, "
        f"q95={q[2]:.4f}, q99={q[3]:.4f}")

clip_arrays = [(samples[n]['soft'] + samples[n]['hard']) / samples[n]['query_len']
               for n in names]
ad_clip = stats.anderson_ksamp(clip_arrays)
print(f"\nAnderson-Darling on clip rates: stat={ad_clip.statistic:.4f}, "
      f"p={ad_clip.pvalue:.2e}")

fig, axes = plt.subplots(2, 2, figsize=(12, 10))

for n in names:
  axes[0, 0].hist(samples[n]['id_A'], bins=100, density=True,
                   histtype='step', label=n)
  axes[0, 1].hist(samples[n]['id_B'], bins=100, density=True,
                   histtype='step', label=n)
  xs_a = np.sort(samples[n]['id_A'].values)
  xs_b = np.sort(samples[n]['id_B'].values)
  axes[1, 0].plot(xs_a, np.arange(1, len(xs_a) + 1) / len(xs_a), label=n)
  axes[1, 1].plot(xs_b, np.arange(1, len(xs_b) + 1) / len(xs_b), label=n)

axes[0, 0].set_title('Clips ignored')
axes[0, 1].set_title('Clips penalized')
axes[0, 0].set_xlabel('Identity'); axes[0, 0].set_ylabel('Density')
axes[0, 1].set_xlabel('Identity'); axes[0, 1].set_ylabel('Density')
axes[1, 0].set_xlabel('Identity'); axes[1, 0].set_ylabel('ECDF')
axes[1, 1].set_xlabel('Identity'); axes[1, 1].set_ylabel('ECDF')
for ax in axes.flat:
  ax.legend()

plt.tight_layout()
plt.savefig('identity_comparison.png', dpi=150)