import numpy as np
import colour

def main():
    data    = colour.DATA_MACADAM_1942_ELLIPSES
    x0      = data[:, 0]
    y0      = data[:, 1]
    a_calc  = data[:, 5] * 1e-3
    b_calc  = data[:, 6] * 1e-3
    phi     = np.radians(data[:, 7])

    print(f"MacAdam ellipses: {len(data)} centers in CIE xy")

    # ------------------------------------------------------------------
    # Invert ellipses → metric tensor in CIE xy
    # ------------------------------------------------------------------
    c = np.cos(phi);  s = np.sin(phi)
    G_xx = c**2/a_calc**2 + s**2/b_calc**2
    G_yy = s**2/a_calc**2 + c**2/b_calc**2
    G_xy = s*c*(1/a_calc**2 - 1/b_calc**2)

    # ------------------------------------------------------------------
    # Stockman-Sharpe XYZ → LMS matrix (Stockman & Sharpe 2000)
    # ------------------------------------------------------------------
    M_xyz2lms = np.array([
        [ 0.210576,  0.855098, -0.0396983],
        [-0.417076,  1.17726,   0.0786957],
        [ 0.0,       0.0102995, 0.723625 ]
    ])

    def xy_to_lms_chroma(x, y):
        """CIE xy (Y=1) → LMS chromaticity (normalized to unit sum)."""
        XYZ = np.array([x/y, 1.0, (1-x-y)/y])
        LMS = M_xyz2lms @ XYZ
        total = LMS.sum()
        return LMS / max(total, 1e-12)

    # ------------------------------------------------------------------
    # Ellipse centers in LMS chromaticity
    # ------------------------------------------------------------------
    lms0 = np.array([xy_to_lms_chroma(x0[i], y0[i]) for i in range(25)])
    l0, m0 = lms0[:, 0], lms0[:, 1]

    l_w, m_w = 1/3, 1/3
    x_mac     = l0 - l_w
    y_mac     = m0 - m_w
    r_mac     = np.sqrt(x_mac**2 + y_mac**2)
    theta_mac = np.arctan2(y_mac, x_mac)
    deg_mac   = np.degrees(theta_mac)

    print(f"Ellipse centers: θ range {deg_mac.min():.1f}° – {deg_mac.max():.1f}°")
    print(f"                 r range {r_mac.min():.4f} – {r_mac.max():.4f}")

    # ------------------------------------------------------------------
    # Jacobian d(l,m)/d(x,y) via central differences → pull back metric
    # G_lms = J^{-T} G_xy J^{-1}
    # ------------------------------------------------------------------
    eps = 1e-6
    G_ll = np.zeros(25);  G_mm = np.zeros(25);  G_lm = np.zeros(25)

    for i in range(25):
        xi, yi = x0[i], y0[i]
        lms_xp = xy_to_lms_chroma(xi+eps, yi)
        lms_xm = xy_to_lms_chroma(xi-eps, yi)
        lms_yp = xy_to_lms_chroma(xi, yi+eps)
        lms_ym = xy_to_lms_chroma(xi, yi-eps)

        J = np.array([
            [(lms_xp[0]-lms_xm[0])/(2*eps), (lms_yp[0]-lms_ym[0])/(2*eps)],
            [(lms_xp[1]-lms_xm[1])/(2*eps), (lms_yp[1]-lms_ym[1])/(2*eps)]
        ])
        G_xy_mat = np.array([[G_xx[i], G_xy[i]],
                             [G_xy[i], G_yy[i]]])
        J_inv        = np.linalg.inv(J)
        G_lms_mat    = J_inv.T @ G_xy_mat @ J_inv
        G_ll[i]      = G_lms_mat[0, 0]
        G_mm[i]      = G_lms_mat[1, 1]
        G_lm[i]      = G_lms_mat[0, 1]

    # ------------------------------------------------------------------
    # β_MacAdam at each center
    # β = sqrt(G_ll sin²θ + G_mm cos²θ - 2 G_lm sinθ cosθ)
    # ------------------------------------------------------------------
    sin_t    = np.sin(theta_mac)
    cos_t    = np.cos(theta_mac)
    expr_mac = G_ll*sin_t**2 + G_mm*cos_t**2 - 2*G_lm*sin_t*cos_t
    beta_mac = np.sqrt(np.maximum(expr_mac, 0))

    print(f"β_MacAdam: min={beta_mac.min():.3f}, max={beta_mac.max():.3f}, mean={beta_mac.mean():.3f}")

    # ------------------------------------------------------------------
    # Interpolate |α'/α| to ellipse center angles
    # ------------------------------------------------------------------
    theta_uniform = np.load("theta_uniform.npy")
    ratio_lhs     = np.load("ratio_lhs.npy")
    ratio_at_mac  = np.interp(theta_mac, theta_uniform, ratio_lhs, period=2*np.pi)

    # ------------------------------------------------------------------
    # Bound test
    # ------------------------------------------------------------------
    margin_mac = beta_mac - ratio_at_mac

    print(f"\n{'i':>3} {'θ':>8} {'|α\'/α|':>8} {'β_Mac':>10} {'margin':>8} {'pass':>5}")
    print("-" * 55)
    for i in range(25):
        print(f"{i:>3} {deg_mac[i]:>8.1f}° {ratio_at_mac[i]:>8.4f} "
              f"{beta_mac[i]:>10.4f} {margin_mac[i]:>8.4f} "
              f"{'✓' if margin_mac[i]>0 else '✗':>5}")

    print(f"\nBound satisfied: {(margin_mac>0).sum()} / 25")
    print(f"Tightest margin: {margin_mac.min():.4f}")
    worst = np.argmin(margin_mac)
    print(f"Worst case: θ={deg_mac[worst]:.1f}°, |α'/α|={ratio_at_mac[worst]:.4f}, β={beta_mac[worst]:.4f}")

    np.save("beta_macadam.npy",  beta_mac)
    np.save("theta_macadam.npy", theta_mac)
    np.save("ratio_macadam.npy", ratio_at_mac)
    print("\nSaved: beta_macadam.npy, theta_macadam.npy, ratio_macadam.npy")

if __name__ == "__main__":
    main()
