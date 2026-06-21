import numpy as np
import matplotlib.pyplot as plt

# --- Step 4: Differentiate α(θ), compute |α'/α| ---

theta = np.load("theta_uniform.npy")   # 360 points, uniform in [-π, π]
alpha = np.load("alpha_uniform.npy")

dtheta = theta[1] - theta[0]           # uniform spacing ~0.01745 rad

# Central differences, periodic wrap
alpha_prime = np.gradient(alpha, dtheta)
# np.gradient uses central differences internally;
# handle periodicity manually at endpoints
alpha_prime[0]  = (alpha[1] - alpha[-1]) / (2 * dtheta)
alpha_prime[-1] = (alpha[0] - alpha[-2]) / (2 * dtheta)

ratio = np.abs(alpha_prime) / alpha    # |α'/α| — left side of bound

print("--- α(θ) ---")
print(f"  min:  {alpha.min():.4f}  max: {alpha.max():.4f}")
print("\n--- |α'(θ)| ---")
print(f"  min:  {np.abs(alpha_prime).min():.4f}  max: {np.abs(alpha_prime).max():.4f}")
print("\n--- |α'/α| ---")
print(f"  min:  {ratio.min():.4f}")
print(f"  max:  {ratio.max():.4f}")
print(f"  mean: {ratio.mean():.4f}")
print(f"  fraction > 1.0: {(ratio > 1.0).mean()*100:.1f}%")

# Where is the ratio largest?
peak_idx = np.argmax(ratio)
print(f"\n  Peak |α'/α| = {ratio[peak_idx]:.4f} "
      f"at θ = {np.degrees(theta[peak_idx]):.1f}°")

np.save("alpha_prime.npy", alpha_prime)
np.save("ratio_lhs.npy", ratio)

# --- Plot ---
fig, axes = plt.subplots(2, 1, figsize=(10, 7), sharex=True)

ax = axes[0]
ax.plot(np.degrees(theta), alpha, 'k-', lw=1.5)
ax.set_ylabel('α(θ)')
ax.set_title('Gamut boundary and its log-derivative')
ax.grid(True, alpha=0.3)

ax = axes[1]
ax.plot(np.degrees(theta), ratio, 'b-', lw=1.5, label='|α\'(θ)/α(θ)|')
ax.axhline(1.0, color='r', lw=1.2, ls='--', label='β = 1  (flat-space bound)')
ax.set_xlabel('hue angle θ (degrees)')
ax.set_ylabel("|α'/α|")
ax.legend()
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig("step4_ratio.png", dpi=150)
plt.close()
print("\nSaved step4_ratio.png")
