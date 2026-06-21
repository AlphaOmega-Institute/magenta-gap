import numpy as np
import matplotlib.pyplot as plt

theta    = np.load("theta_uniform.npy")
ratio    = np.load("ratio_lhs.npy")
beta     = np.load("beta_uniform.npy")
alpha    = np.load("alpha_uniform.npy")
margin   = beta - ratio
deg      = np.degrees(theta)

# Clip to non-degenerate region only
BETA_MAX = 10.0
mask = beta < BETA_MAX

print(f"Non-degenerate angles: {mask.sum()} / {len(mask)}")
print(f"|α'/α| in region: max={ratio[mask].max():.4f}")
print(f"β       in region: min={beta[mask].min():.4f}, "
      f"max={beta[mask].max():.4f}")
print(f"Bound satisfied:   {(margin[mask]>0).sum()} / {mask.sum()}")
print(f"Tightest margin:   {margin[mask].min():.4f}")

fig, axes = plt.subplots(2, 1, figsize=(10, 7), sharex=True)

ax = axes[0]
ax.plot(deg[mask], ratio[mask], 'b-', lw=1.8, label="|α'/α|")
ax.plot(deg[mask], beta[mask],  'r-', lw=1.8, label='β (Fisher)')
ax.axhline(1.0, color='gray', lw=0.8, ls=':', label='β=1 flat')
ax.set_ylabel('ratio')
ax.set_title('Non-degenerate region (β < 10): |α\'(θ)/α(θ)| vs β(θ)')
ax.legend()
ax.grid(True, alpha=0.3)

ax = axes[1]
m = margin[mask]
d = deg[mask]
ax.plot(d, m, 'g-', lw=1.8)
ax.axhline(0, color='r', lw=1.2, ls='--', label='bound boundary')
ax.fill_between(d, m, 0, where=(m < 0),
                color='red',   alpha=0.35, label='violation')
ax.fill_between(d, m, 0, where=(m >= 0),
                color='green', alpha=0.15, label='satisfied')
ax.set_xlabel('hue angle θ (degrees)')
ax.set_ylabel('β − |α\'/α|')
ax.set_title('Margin in non-degenerate region (positive = bound satisfied)')
ax.legend()
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig("step5c_clipped.png", dpi=150)
plt.close()
print("Saved step5c_clipped.png")
