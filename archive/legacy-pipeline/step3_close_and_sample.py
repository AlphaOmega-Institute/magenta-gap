import numpy as np
import matplotlib.pyplot as plt
import colour

# --- Recompute chromatic coords (same as step 2) ---
cmfs = colour.colorimetry.MSDS_CMFS_LMS["Stockman & Sharpe 2 Degree Cone Fundamentals"]
wavelengths = cmfs.wavelengths
LMS = cmfs.values
vis_mask = (wavelengths >= 390) & (wavelengths <= 700)
wl_vis  = wavelengths[vis_mask]
LMS_vis = LMS[vis_mask]
total   = LMS_vis.sum(axis=1, keepdims=True)
total   = np.where(total < 1e-12, 1e-12, total)
lms_c   = LMS_vis / total
l_w, m_w = 1/3, 1/3
x_c = lms_c[:, 0] - l_w
y_c = lms_c[:, 1] - m_w

# Endpoints: 390 nm (violet end) and 700 nm (red end)
x_violet, y_violet = x_c[0],  y_c[0]
x_red,    y_red    = x_c[-1], y_c[-1]

print(f"Violet endpoint (390nm): ({x_violet:.4f}, {y_violet:.4f}), "
      f"θ={np.degrees(np.arctan2(y_violet,x_violet)):.1f}°")
print(f"Red endpoint   (700nm): ({x_red:.4f}, {y_red:.4f}), "
      f"θ={np.degrees(np.arctan2(y_red,x_red)):.1f}°")

# --- Purple line: linear interpolation between red and violet endpoints ---
N_purple = 60
t = np.linspace(0, 1, N_purple, endpoint=False)
x_purple = x_red + t * (x_violet - x_red)
y_purple = y_red + t * (y_violet - y_red)
r_purple = np.sqrt(x_purple**2 + y_purple**2)
theta_purple = np.arctan2(y_purple, x_purple)

# --- Combine spectral locus + purple line ---
theta_raw = np.load("theta_alpha.npy")
alpha_raw = np.load("alpha_theta.npy")

theta_all = np.concatenate([theta_raw, theta_purple])
alpha_all = np.concatenate([alpha_raw, r_purple])

# Sort combined set by angle
idx = np.argsort(theta_all)
theta_all = theta_all[idx]
alpha_all = alpha_all[idx]

# --- Resample uniformly: 360 points, one per degree ---
N = 360
theta_uniform = np.linspace(-np.pi, np.pi, N, endpoint=False)
alpha_uniform = np.interp(theta_uniform,
                          theta_all,
                          alpha_all,
                          period=2*np.pi)

print(f"\nResampled α(θ): {N} points")
print(f"α range: {alpha_uniform.min():.4f} – {alpha_uniform.max():.4f}")

np.save("theta_uniform.npy", theta_uniform)
np.save("alpha_uniform.npy", alpha_uniform)

# --- Plot ---
fig, ax = plt.subplots(figsize=(8, 4))
ax.plot(np.degrees(theta_uniform), alpha_uniform, 'k-', lw=1.5,
        label='spectral locus + purple line')
ax.axvspan(np.degrees(theta_purple.min()),
           np.degrees(theta_purple.max()),
           alpha=0.12, color='purple', label='purple line region')
ax.set_xlabel('hue angle θ (degrees)')
ax.set_ylabel('α(θ)')
ax.set_title('Closed gamut boundary α(θ) — 360 uniform samples')
ax.legend()
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig("step3_alpha_closed.png", dpi=150)
plt.close()
print("Saved step3_alpha_closed.png")