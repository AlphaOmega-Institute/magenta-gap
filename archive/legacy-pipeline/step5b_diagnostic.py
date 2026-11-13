import numpy as np

theta    = np.load("theta_uniform.npy")
alpha    = np.load("alpha_uniform.npy")
ratio    = np.load("ratio_lhs.npy")
beta     = np.load("beta_uniform.npy")
margin   = beta - ratio

deg = np.degrees(theta)

# Where is β smallest? That's where the bound is tightest.
idx_bmin = np.argmin(beta)
print(f"β minimum: {beta[idx_bmin]:.4f} at θ={deg[idx_bmin]:.1f}°")
print(f"  |α'/α| there: {ratio[idx_bmin]:.4f}")
print(f"  margin there: {margin[idx_bmin]:.4f}")

# Where is |α'/α| largest vs local β?
tightest = np.argmin(margin)
print(f"\nTightest margin: {margin[tightest]:.4f} at θ={deg[tightest]:.1f}°")
print(f"  |α'/α|={ratio[tightest]:.4f}, β={beta[tightest]:.4f}")

# Region where β < 10 (non-degenerate Fisher)
nondegenerate = beta < 10
print(f"\nAngles with β < 10 (non-degenerate): {nondegenerate.sum()}")
if nondegenerate.sum() > 0:
    print(f"  θ range: {deg[nondegenerate].min():.1f}° to "
          f"{deg[nondegenerate].max():.1f}°")
    print(f"  |α'/α| in that region: max={ratio[nondegenerate].max():.4f}, "
          f"mean={ratio[nondegenerate].mean():.4f}")
    print(f"  β in that region:      min={beta[nondegenerate].min():.4f}, "
          f"mean={beta[nondegenerate].mean():.4f}")
    print(f"  Bound satisfied there: "
          f"{(margin[nondegenerate]>0).mean()*100:.1f}%")

# Region where |α'/α| > 1 (flat-space violations)
violations_flat = ratio > 1.0
print(f"\nFlat-space violation angles ({violations_flat.sum()} total):")
print(f"  β at those angles: min={beta[violations_flat].min():.3f}, "
      f"max={beta[violations_flat].max():.1f}, "
      f"mean={beta[violations_flat].mean():.1f}")
