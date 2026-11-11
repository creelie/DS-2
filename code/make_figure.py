"""Draws fig_supports.png: supports of x, x^2, x^3, x^2x^2, x^3x for x = L_1 + L_{1+tau}."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import io, contextlib
with contextlib.redirect_stdout(io.StringIO()):
    from verify_jordan import prod
from fractions import Fraction as Fr
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
tau = (1 + 5 ** 0.5) / 2
x = {(1, 0): Fr(1), (1, 1): Fr(1)}
x2 = prod(x, x); x3 = prod(x2, x)
rows = [(r"$x$", x), (r"$x^2$", x2), (r"$x^3$", x3), (r"$x^2x^2$", prod(x2, x2)), (r"$x^3x$", prod(x3, x))]
fig = plt.figure(figsize=(8.6, 5.2))
ax = fig.add_subplot(111, projection="3d")
cols = ["#4C72B0", "#55A868", "#8172B2", "#C44E52", "#DD8452"]
for i, (lab, u) in enumerate(rows):
    xs = [p + q * tau for (p, q) in u]; hs = [float(c) for c in u.values()]
    ax.bar3d(xs, [i] * len(xs), [0] * len(xs), 0.55, 0.45, hs, color=cols[i], alpha=0.85, shade=True)
top = 9 + 13 * tau
ax.bar3d([top], [4], [0], 0.55, 0.45, [0.25], color="black", shade=True)
ax.set_yticks(range(len(rows))); ax.set_yticklabels([r[0] for r in rows])
ax.set_xlabel("position on the real line", labelpad=8); ax.set_zlabel("coefficient", labelpad=6)
ax.view_init(elev=22, azim=-78)
plt.tight_layout()
out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "paper", "fig_supports.png")
plt.savefig(out, dpi=200, bbox_inches="tight", pad_inches=0.15)
print("wrote", os.path.abspath(out))
