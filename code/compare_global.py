"""Two-sample test of the superbase laws of two conformally equivalent tori, conditioned on a theta core.

Usage: python compare_global.py A.tsv.gz B.tsv.gz m11 m12 m21 m22 [--pool 200]
The matrix sends a class (a, b) of B to the class (m11*a + m12*b, m21*a + m22*b) of A, as in compare.py. Superbases
with fewer than --pool trees in the two files together are pooled into one cell. Prints the Pearson chi-square of the
2 x K table, its degrees of freedom and the figure-eight counts of both files."""
import gzip
import math
import sys
from collections import Counter


def load(fn, M):
    cnt, tot, eight = Counter(), 0, 0
    op = gzip.open if fn.endswith(".gz") else open
    with op(fn, "rt", encoding="utf-8") as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            tot += 1
            if int(p[3]) != 3:
                eight += 1
                continue
            cls = []
            for x in p[4].split(";"):
                a, b = map(int, x.split(","))
                a, b = M[0] * a + M[1] * b, M[2] * a + M[3] * b
                if a < 0 or (a == 0 and b < 0):
                    a, b = -a, -b
                cls.append((a, b))
            cnt[tuple(sorted(cls))] += 1
    return cnt, tot, eight


def main():
    args = sys.argv[1:]
    pool = 200
    if "--pool" in args:
        i = args.index("--pool")
        pool = int(args[i + 1])
        args = args[:i] + args[i + 2:]
    A, tA, eA = load(args[0], (1, 0, 0, 1))
    B, tB, eB = load(args[1], tuple(int(x) for x in args[2:6]))
    keys = set(A) | set(B)
    big = [k for k in keys if A[k] + B[k] >= pool]
    rows = [(A[k], B[k]) for k in big]
    rest = (sum(A[k] for k in keys if k not in big), sum(B[k] for k in keys if k not in big))
    if rest[0] + rest[1] > 0:
        rows.append(rest)
    nA, nB = sum(r[0] for r in rows), sum(r[1] for r in rows)
    chi = 0.0
    for a, b in rows:
        t = a + b
        for o, n in ((a, nA), (b, nB)):
            e = t * n / (nA + nB)
            chi += (o - e) ** 2 / e
    dof = len(rows) - 1
    z8 = (eA / tA - eB / tB) / math.sqrt(eA / tA ** 2 + eB / tB ** 2) if eA + eB else 0.0
    print(f"theta law: chi2 {chi:.3f} on {dof} dof ({len(big)} superbases, pooled rest {rest})")
    print(f"figure-eight: A {eA} of {tA}, B {eB} of {tB}, difference {z8:+.2f} se")


if __name__ == "__main__":
    main()
