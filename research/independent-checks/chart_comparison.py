"""The bound's scale in the appendix chart versus the opponent chart of the July manuscript.

Inputs   Stockman-Sharpe table; hull of the 390-700 nm chromaticities, in both normalizations.
Checks   max (u/beta)^2 over the retained open-edge samples (beta < 10 in the appendix chart), first in the chart
         (l - l_w, m - m_w), then in the opponent chart x = (l - m)/2, y = (l + m)/4 - s/2 with the Fisher tensor transported.
Status   Numerical.  Both reviewers obtained 0.376081907 / 0.430243986 (white-balanced) and 0.402831775 / 0.485244688 (raw).
Scope    (agreed in review) The calculation re-applies the manuscript's calibration and angular identification in the new
         chart, so the change exposes a missing identification between Sec. 5.1 and the appendix; it does not contradict
         tensor covariance or the conditional theorem.  Applying the opponent formulas to normalized catches at fixed
         total is itself a step the manuscript does not justify.
"""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from common import load_table, REPO
import sys, numpy as np
from scipy.spatial import ConvexHull
wl_all, LMS_all = load_table()
Eint = np.trapezoid(LMS_all, wl_all, axis=0); m = (wl_all >= 390) & (wl_all <= 700)
M2 = np.array([[0.5, -0.5], [0.75, 0.75]]); M2i = np.linalg.inv(M2)
for name, X, w3 in (("white-balanced", LMS_all[m]/Eint, np.ones(3)), ("raw, E-centred", LMS_all[m], Eint)):
    c = X / X.sum(1, keepdims=True); w = (w3/w3.sum())[:2]; pts = c[:, :2] - w
    V = pts[ConvexHull(pts).vertices]; t = np.linspace(1e-5, 1 - 1e-5, 101); a_, b_ = [], []
    for a, b in zip(V, np.roll(V, -1, axis=0)):
        d = b - a
        for tt in t:
            p = a + tt*d; l, mm = p + w; s = 1 - l - mm
            if min(l, mm, s) <= 0: continue
            G = np.array([[1/l + 1/s, 1/s], [1/s, 1/mm + 1/s]])
            e = np.array([-p[1], p[0]])/np.hypot(*p); beta = np.sqrt(e @ G @ e)
            if beta >= 10: continue
            u = abs(p @ d)/abs(p[0]*d[1] - p[1]*d[0])
            p2, d2 = M2 @ p, M2 @ d; G2 = M2i.T @ G @ M2i; e2 = np.array([-p2[1], p2[0]])/np.hypot(*p2)
            a_.append((u/beta)**2); b_.append((abs(p2 @ d2)/abs(p2[0]*d2[1] - p2[1]*d2[0]))**2/(e2 @ G2 @ e2))
    print(f"{name}: retained {len(a_)}; appendix chart {max(a_):.9f}; opponent chart {max(b_):.9f}")
