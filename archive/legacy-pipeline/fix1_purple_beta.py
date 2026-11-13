import numpy as np
import colour

def fisher_metric_simplex(l, m):
    s = 1.0 - l - m
    s = np.maximum(s, 1e-12)
    l = np.maximum(l, 1e-12)
    m = np.maximum(m, 1e-12)
    J11 = 1.0/l + 1.0/s
    J22 = 1.0/m + 1.0/s
    J12 = 1.0/s
    return J11, J12, J22

def main():
    # Load everything from previous steps
    theta_uniform = np.load("theta_uniform.npy")
    alpha_uniform = np.load("alpha_uniform.npy")
    ratio_lhs     = np.load("ratio_lhs.npy")
    beta_old      = np.load("beta_uniform.npy")
    x_c           = np.load("x_c_vis.npy")
    y_c           = np.load("y_c_vis.npy")

    # Purple chord endpoints (390nm violet, 700nm red)
    x_violet, y_violet = x_c[0],  y_c[0]
    x_red,    y_red    = x_c[-1], y_c[-1]

    # White point
    l_w, m_w = 1/3, 1/3

    # For each uniform grid point, decide: spectral or purple?
    # Purple region: θ between red endpoint and violet endpoint
    # going through the bottom (i.e., -134.6° to -24.3°)
    theta_red    = np.arctan2(y_red,    x_red)
    theta_violet = np.arctan2(y_violet, x_violet)

    deg = np.degrees(theta_uniform)
    deg_red    = np.degrees(theta_red)
    deg_violet = np.degrees(theta_violet)

    print(f"Red endpoint:    θ = {deg_red:.2f}°")
    print(f"Violet endpoint: θ = {deg_violet:.2f}°")

    # Purple mask: angles strictly between violet and red endpoints
    # (the arc that has no spectral representation)
    purple_mask = (deg > deg_violet) & (deg < deg_red)
    print(f"Purple region: {purple_mask.sum()} grid angles")

    # For each purple-region grid point, find the Cartesian point on the
    # chord at that angle, then compute its (l, m, s) chromaticity,
    # then compute β from the Fisher metric at that point.
    beta_fixed = beta_old.copy()

    for i in np.where(purple_mask)[0]:
        th = theta_uniform[i]
        r  = alpha_uniform[i]   # radius of purple chord at this angle

        # Cartesian coords of this point on the purple chord
        x = r * np.cos(th)
        y = r * np.sin(th)

        # Chromaticity (un-center)
        l = x + l_w
        m = y + m_w
        s = 1.0 - l - m

        # Fisher metric at this chromaticity
        J11, J12, J22 = fisher_metric_simplex(l, m)

        # β = sqrt(J11 sin²θ + J22 cos²θ - 2 J12 sinθ cosθ)
        sin_t = np.sin(th)
        cos_t = np.cos(th)
        expr  = J11*sin_t**2 + J22*cos_t**2 - 2*J12*sin_t*cos_t
        beta_fixed[i] = np.sqrt(np.maximum(expr, 0))

    # Report changes
    changed = np.where(purple_mask)[0]
    print(f"\nβ changes in purple region:")
    print(f"  Old β (interpolated): min={beta_old[changed].min():.3f}, max={beta_old[changed].max():.1f}")
    print(f"  New β (from mixture chromaticity): min={beta_fixed[changed].min():.3f}, max={beta_fixed[changed].max():.3f}")

    # Re-run bound check with fixed β
    margin_fixed = beta_fixed - ratio_lhs
    nd_fixed     = beta_fixed < 10.0

    print(f"\n=== RESULTS WITH FIXED β ===")
    print(f"Non-degenerate (β<10): {nd_fixed.sum()} angles")
    print(f"Bound satisfied there: {(margin_fixed[nd_fixed]>0).sum()}/{nd_fixed.sum()}")
    print(f"Tightest margin: {margin_fixed[nd_fixed].min():.4f}")
    tightest = np.where(nd_fixed)[0][np.argmin(margin_fixed[nd_fixed])]
    deg_u = np.degrees(theta_uniform)
    print(f"  at θ={deg_u[tightest]:.1f}°: |α'/α|={ratio_lhs[tightest]:.4f}, β={beta_fixed[tightest]:.4f}")

    # Save fixed β
    np.save("beta_uniform_fixed.npy", beta_fixed)
    print("\nSaved: beta_uniform_fixed.npy")

if __name__ == "__main__":
    main()
