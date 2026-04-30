import matplotlib.pyplot as plt
import numpy as np

cpus = np.array([1, 2, 4, 8, 16, 32, 64])
sars_runtimes = np.array([112.342, 54.034, 26.366, 14.222, 7.696, 4.201, 3.68])
rsv_runtimes = np.array([71.958, 35.169, 16.766, 8.281, 4.767, 2.751, 2.325])
hiv_runtimes = np.array([645.849, 298.386, 140.827, 64.137, 32.982, 17.682, 11.287])

sars_relative = sars_runtimes[0] / sars_runtimes
rsv_relative = rsv_runtimes[0] / rsv_runtimes
hiv_relative = hiv_runtimes[0] / hiv_runtimes

fig, ax = plt.subplots(figsize=(12, 8))

ax.plot(cpus, sars_relative, marker='o', linewidth=2, markersize=8, label='SARS_20K', color='#e74c3c')
ax.plot(cpus, rsv_relative, marker='s', linewidth=2, markersize=8, label='RSV_4K', color='#3498db')
ax.plot(cpus, hiv_relative, marker='^', linewidth=2, markersize=8, label='HIV_20K', color='#2ecc71')

ax.plot(cpus, cpus, '--', color='gray', alpha=0.7, linewidth=1.5, label='Linear speedup')

for i, cpu in enumerate(cpus):
  ax.annotate(f'{sars_relative[i]:.1f}X', 
              (cpu, sars_relative[i]), 
              textcoords="offset points", 
              xytext=(0, 10), 
              ha='center', 
              fontsize=9,
              color='#e74c3c')
  
  ax.annotate(f'{rsv_relative[i]:.1f}X', 
              (cpu, rsv_relative[i]), 
              textcoords="offset points", 
              xytext=(0, -15), 
              ha='center', 
              fontsize=9,
              color='#3498db')
  
  ax.annotate(f'{hiv_relative[i]:.1f}X', 
              (cpu, hiv_relative[i]), 
              textcoords="offset points", 
              xytext=(0, 10), 
              ha='center', 
              fontsize=9,
              color='#2ecc71')

ax.set_xlabel('Number of CPUs', fontsize=14, fontweight='bold')
ax.set_ylabel('Speedup (relative to 1 CPU)', fontsize=14, fontweight='bold')
ax.set_title('Thread Efficiency: Speedup vs Number of CPUs', fontsize=16, fontweight='bold')
ax.set_xscale('log', base=2)
ax.set_yscale('log', base=2)
ax.set_xticks(cpus)
ax.set_xticklabels(cpus)
ax.grid(True, alpha=0.3, linestyle='--')
ax.legend(fontsize=12, loc='upper left')

plt.tight_layout()
plt.savefig('thread_efficiency.png', dpi=300, bbox_inches='tight')
plt.show()

print("Plot saved to thread_efficiency.png")
print("\nSpeedup summary:")
print(f"SARS:  1 CPU -> 64 CPUs: {sars_relative[-1]:.2f}X speedup (efficiency: {sars_relative[-1]/64*100:.1f}%)")
print(f"RSV:   1 CPU -> 64 CPUs: {rsv_relative[-1]:.2f}X speedup (efficiency: {rsv_relative[-1]/64*100:.1f}%)")
print(f"HIV:   1 CPU -> 64 CPUs: {hiv_relative[-1]:.2f}X speedup (efficiency: {hiv_relative[-1]/64*100:.1f}%)")