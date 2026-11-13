import numpy as np
import matplotlib.pyplot as plt

# --- Step 2: Spectral locus → gamut boundary α(θ) ---

wavelengths = np.load("wavelengths.npy")
LMS = np.load("LMS.npy")

# Restrict to visible range where cone responses are non-negligible
vis_mask = (wavelengths >= 390) & (wavelengths <= 700)
wl_vis = wavelengths[vis_mask]
LMS_vis = LMS[vis_mask]          # shape (N_vis, 3)

# --- Chromaticity: project to (l, m) = (L, M) / (L+M+S) ---
total = LMS_vis.sum(axis=1, keepdims=True)
total = np.where(total < 1e-12, 1e-12, total)   # avoid /0 at wings
lms_chroma = LMS_vis / total                      # shape (N_vis, 3)
l_c = lms_chroma[:, 0]
m_c = lms_chroma[:, 1]

# White point: equal-energy illuminant E in LMS chromaticity
# (sum of each cone row / total, integrated across spectrum → equal weights)
# For illuminant E: L=M=S equal, so chromaticity = (1/3, 1/3, 1/3)
l_w, m_w = 1/3, 1/3

# Centre on white point → opponent-like chromatic coordinates
x_c = l_c - l_w
y_c = m_c - m_w

# --- Polar coordinates about white ---
r_locus = np.sqrt(x_c**2 + y_c**2)
theta_locus = np.arctan2(y_c, x_c)   # radians, -π to π

# Sort by angle for a clean curve
sort_idx = np.argsort(theta_locus)
theta_sorted = theta_locus[sort_idx]
r_sorted = r_locus[sort_idx]
wl_sorted = wl_vis[sort_idx]

# --- α(θ): max chroma at each angle ---
# The spectral locus IS the boundary for self-luminous colours,
# so r_sorted is already α(θ) along the locus.
alpha_theta = r_sorted
theta_alpha = theta_sorted

print(f"Spectral locus: {len(wl_vis)} points, "
      f"θ range {np.degrees(theta_sorted[0]):.1f}° to "
      f"{np.degrees(theta_sorted[-1]):.1f}°")
print(f"α range: {alpha_theta.min():.4f} – {alpha_theta.max():.4f}")
print(f"α mean:  {alpha_theta.mean():.4f}")

# Save for step 3
np.save("theta_alpha.npy", theta_alpha)
np.save("alpha_theta.npy", alpha_theta)
np.save("wl_sorted.npy", wl_sorted)

# Save wavelength-ordered Cartesian chromaticity coords for fix1
# (x_c[0], y_c[0]) is the 390 nm violet endpoint and
# (x_c[-1], y_c[-1]) is the 700 nm red endpoint of the purple chord.
np.save("x_c_vis.npy", x_c)
np.save("y_c_vis.npy", y_c)

# --- Plot ---
fig, axes = plt.subplots(1, 2, figsize=(12, 5))

# Left: spectral locus in chromatic plane
ax = axes[0]
sc = ax.scatter(x_c, y_c, c=wl_vis, cmap='nipy_spectral', s=8,
                vmin=390, vmax=700)
ax.plot(0, 0, 'k+', markersize=10, label='white (E)')
ax.set_xlabel('l − l_w')
ax.set_ylabel('m − m_w')
ax.set_title('Spectral locus (LMS chromaticity, centred on E)')
ax.legend()
ax.set_aspect('equal')
plt.colorbar(sc, ax=ax, label='wavelength (nm)')

# Right: α(θ) around the hue circle
ax = axes[1]
ax.plot(np.degrees(theta_alpha), alpha_theta, 'k-', lw=1.5)
ax.set_xlabel('hue angle θ (degrees)')
ax.set_ylabel('α(θ)  [max chromaticity radius]')
ax.set_title('Gamut boundary α(θ) — spectral locus')
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig("step2_spectral_locus.png", dpi=150)
plt.close()
print("Saved step2_spectral_locus.png")
