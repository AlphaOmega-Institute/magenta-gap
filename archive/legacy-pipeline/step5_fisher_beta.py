import numpy as np
import matplotlib.pyplot as plt
import colour

# --- Step 5: Fisher information metric → β(θ) ---
#
# Model: quantal catches are Poisson, so Fisher metric on
# chromaticity = J_ij = Σ_k (∂a_k/∂c_i)(∂a_k/∂c_j) / a_k
# where a_k = ∫ l_k(λ) Φ(λ) dλ  (cone k catch for stimulus Φ)
# and c = (l, m) are chromaticity coords.
#
# For monochromatic stimuli Φ = δ(λ - λ0), a_k = l_k(λ0).
# We compute J at each point on the spectral locus, then
# extract the angular component B = g_θθ and form β = √B / R.

cmfs = colour.colorimetry.MSDS_CMFS_LMS[
    "Stockman & Sharpe 2 Degree Cone Fundamentals"]
wavelengths = cmfs.wavelengths
LMS_raw = cmfs.values                  # shape (441, 3)

vis_mask = (wavelengths >= 390) & (wavelengths <= 700)
wl_vis   = wavelengths[vis_mask]
LMS_vis  = LMS_raw[vis_mask]           # shape (N, 3), columns L M S

# Chromaticity l = L/(L+M+S), m = M/(L+M+S)
total = LMS_vis.sum(axis=1)
total = np.where(total < 1e-12, 1e-12, total)
l_c   = LMS_vis[:, 0] / total
m_c   = LMS_vis[:, 1] / total

l_w, m_w = 1/3, 1/3
x_c = l_c - l_w
y_c = m_c - m_w

# --- Fisher metric at each spectral locus point ---
# For Poisson model with 3 cones and chromaticity parameterisation:
# J_ij = Σ_k [ (∂L_k/∂c_i) / L_k - (∂L_k/∂c_j) / L_k ] * L_k
#       = Σ_k (∂ log L_k/∂c_i)(∂ log L_k/∂c_j) * L_k
#
# where c_1 = l, c_2 = m and L_k are the raw catches.
# ∂L/∂l = ∂(l·T)/∂l etc. — use numerical differentiation on the
# chromaticity map, which is smoother than symbolic.

dl = 1e-5   # finite difference step in chromaticity

# We need J in (x_c, y_c) coords = shifted (l, m).
# Shift doesn't change the metric (translation invariant).

N = len(wl_vis)
J11 = np.zeros(N)
J12 = np.zeros(N)
J22 = np.zeros(N)

for i in range(N):
    L = LMS_vis[i, 0]
    M = LMS_vis[i, 1]
    S = LMS_vis[i, 2]
    T = L + M + S
    if T < 1e-12:
        continue

    # l = L/T, m = M/T, s = S/T, l+m+s = 1
    # LMS are fixed (monochromatic); l,m are coords on unit simplex
    # Derivatives of log-catches w.r.t. l and m:
    # log L = log(l·T) ... but T itself doesn't depend on coord choice
    # For Fisher on the chromaticity simplex:
    # ∂log(L_k)/∂l_i = δ_{ki}/l_k - 1  (standard simplex Fisher)

    # Fisher on the 2-simplex (l, m, s=1-l-m):
    # g_ij = Σ_k δ_{ki}δ_{kj}/l_k - sum issues...
    # Clean form: g_ij = 1/l_i δ_ij + 1/s  (i,j ∈ {l,m})

    s_c = 1 - l_c[i] - m_c[i]
    s_c = max(s_c, 1e-12)

    # Fisher metric on 2D chromaticity simplex
    # (Rao 1945, standard result):
    # g_ll = 1/l + 1/s
    # g_mm = 1/m + 1/s
    # g_lm = 1/s
    J11[i] = 1.0/l_c[i] + 1.0/s_c
    J22[i] = 1.0/m_c[i] + 1.0/s_c
    J12[i] = 1.0/s_c

# --- Convert J to polar coords centred on white ---
# In (x,y) = (l-lw, m-mw) coords, the metric tensor is J.
# Angular component: g_θθ = r² (J transformed to polar)
# g_θθ = J11 sin²θ - 2 J12 sinθ cosθ + J22 cos²θ  ... wait:
# x = r cosθ, y = r sinθ
# ∂x/∂θ = -r sinθ,  ∂y/∂θ = r cosθ
# g_θθ = J11 (r sinθ)² - 2 J12 (r sinθ)(r cosθ) + J22 (r cosθ)²
#       = r² [J11 sin²θ + J22 cos²θ - 2 J12 sinθ cosθ]

r_locus   = np.sqrt(x_c**2 + y_c**2)
theta_loc = np.arctan2(y_c, x_c)
sin_t = np.sin(theta_loc)
cos_t = np.cos(theta_loc)

B_locus = r_locus**2 * (J11 * sin_t**2 + J22 * cos_t**2
                        - 2 * J12 * sin_t * cos_t)

# β = √B / R  (R = r at boundary = r_locus here)
beta_locus = np.sqrt(np.maximum(B_locus, 0)) / np.maximum(r_locus, 1e-12)

print(f"β on spectral locus:")
print(f"  min:  {beta_locus.min():.3f}")
print(f"  max:  {beta_locus.max():.3f}")
print(f"  mean: {beta_locus.mean():.3f}")

# --- Resample β onto same uniform θ grid as α ---
theta_uniform = np.load("theta_uniform.npy")
alpha_uniform = np.load("alpha_uniform.npy")
ratio_lhs     = np.load("ratio_lhs.npy")

# Sort locus by angle for interpolation
sort_idx    = np.argsort(theta_loc)
theta_s     = theta_loc[sort_idx]
beta_s      = beta_locus[sort_idx]

beta_uniform = np.interp(theta_uniform, theta_s, beta_s)

np.save("beta_uniform.npy", beta_uniform)

# --- THE TEST ---
margin   = beta_uniform - ratio_lhs   # positive = bound satisfied
fraction_satisfied = (margin > 0).mean()
print(f"\n|α'/α| < β  satisfied at {fraction_satisfied*100:.1f}% of angles")
print(f"Max violation (|α'/α| - β): {(-margin).max():.4f}")
print(f"Mean margin:                 {margin.mean():.4f}")

# --- Plot ---
fig, axes = plt.subplots(2, 1, figsize=(10, 7), sharex=True)

ax = axes[0]
ax.plot(np.degrees(theta_uniform), ratio_lhs,    'b-',  lw=1.5,
        label="|α'/α|")
ax.plot(np.degrees(theta_uniform), beta_uniform,  'r-',  lw=1.5,
        label='β (Fisher)')
ax.axhline(1.0, color='gray', lw=0.8, ls=':', label='β=1 flat')
ax.set_ylabel('ratio')
ax.set_title('Test: |α\'(θ)/α(θ)| vs β(θ) — Fisher metric bound')
ax.legend()
ax.grid(True, alpha=0.3)

ax = axes[1]
ax.plot(np.degrees(theta_uniform), margin, 'g-', lw=1.5)
ax.axhline(0, color='r', lw=1.2, ls='--', label='bound boundary')
ax.fill_between(np.degrees(theta_uniform), margin, 0,
                where=(margin < 0), color='red', alpha=0.3,
                label='violation')
ax.fill_between(np.degrees(theta_uniform), margin, 0,
                where=(margin >= 0), color='green', alpha=0.15,
                label='satisfied')
ax.set_xlabel('hue angle θ (degrees)')
ax.set_ylabel('β − |α\'/α|')
ax.set_title('Margin (positive = bound satisfied)')
ax.legend()
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig("step5_fisher_test.png", dpi=150)
plt.close()
print("Saved step5_fisher_test.png")
