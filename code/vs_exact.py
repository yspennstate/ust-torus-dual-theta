"""Compare sampled class sets (ust_torus_shear.py / ust_metal output) with the exact law from exact_small.py.
Usage: python vs_exact.py exact.json samples.tsv"""
import sys, json, math
from collections import Counter

ex = json.load(open(sys.argv[1], encoding="utf-8"))
p = {k: v / ex["ntrees"] for k, v in ex["counts"].items()}
cnt = Counter()
tot = 0
for line in open(sys.argv[2], encoding="utf-8"):
    f = line.rstrip("\n").split("\t")
    cls = sorted(tuple(map(int, x.split(","))) for x in f[4].split(";")) if f[4] else []
    cnt[str(tuple(cls))] += 1
    tot += 1
chi2 = 0.0
worst = 0.0
print(f"{sys.argv[2]}: {tot} samples against {ex['ntrees']} trees")
for k in sorted(set(p) | set(cnt), key=lambda k: -p.get(k, 0)):
    e = p.get(k, 0.0)
    o = cnt.get(k, 0) / tot
    se = math.sqrt(max(e * (1 - e), 1e-300) / tot)
    z = (o - e) / se
    worst = max(worst, abs(z))
    if e > 0:
        chi2 += (cnt.get(k, 0) - e * tot) ** 2 / (e * tot)
    print(f"  {e:.6f} {o:.6f} z={z:+.2f} {k}")
print(f"chi2 {chi2:.1f} on {len(p) - 1} dof, worst |z| {worst:.2f}, unseen-in-exact classes "
      f"{sum(v for k, v in cnt.items() if k not in p)}")
