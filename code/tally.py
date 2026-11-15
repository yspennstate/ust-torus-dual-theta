"""Independent tally of ust_metal / ust_torus_shear output (written for the paper, separate from classify.py).
For each file: samples, figure-eight frequency, and the superbases in decreasing probability with
Q(S) = (|v1|^2 + |v2|^2 + |v3|^2) / Im tau, v = a + b tau, tau = (s + i n2) / n1.
Usage: python tally.py [--top K] file.tsv.gz [...]"""
import sys, gzip, math
from collections import Counter

args = sys.argv[1:]
top = 8
if args and args[0] == "--top":
    top = int(args[1]); args = args[2:]
for fn in args:
    op = gzip.open if fn.endswith(".gz") else open
    cnt = Counter(); tot = 0; eight = 0; other = 0; tau = None
    with op(fn, "rt", encoding="utf-8") as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            n1, n2, s, m = int(p[0]), int(p[1]), int(p[2]), int(p[3])
            tau = complex(s, n2) / n1
            tot += 1
            cls = frozenset(tuple(map(int, x.split(","))) for x in p[4].split(";")) if p[4] else frozenset()
            if m == 2:
                eight += 1
            elif m == 3:
                cnt[cls] += 1
            else:
                other += 1
    print(f"{fn.split('/')[-1]}: tau {tau.real:.4f}+{tau.imag:.4f}i  samples {tot}  figure-eight {eight / tot:.5f}"
          f"  other {other}")
    for cls, v in cnt.most_common(top):
        vs = sorted(cls, key=lambda c: abs(c[0] + c[1] * tau))
        Q = sum(abs(a + b * tau) ** 2 for a, b in vs) / tau.imag
        pr = v / tot
        print(f"   {pr:.5f} +- {math.sqrt(pr * (1 - pr) / tot):.5f}  Q {Q:7.3f}  {vs}")
