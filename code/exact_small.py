"""Exact law of the set of homology classes of fundamental cycles of a uniform spanning tree on a small torus
Z^2 / <(n1, 0), (s, n2)>, by listing every spanning tree. Independent of the samplers (no random walks).
Usage: python exact_small.py n1 n2 s [out.json]"""
import sys, json, itertools
from collections import Counter

n1, n2, s = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
N = n1 * n2


def vid(x, y):
    # reduce (x, y) modulo the lattice: first y modulo n2 (moving by (s, n2)), then x modulo n1
    k = y // n2
    x -= k * s
    y -= k * n2
    return (x % n1) + n1 * y


edges = []   # (u, v, dx, dy): v is the neighbour of u in direction (dx, dy)
for y in range(n2):
    for x in range(n1):
        u = x + n1 * y
        edges.append((u, vid(x + 1, y), 1, 0))
        edges.append((u, vid(x, y + 1), 0, 1))
E = len(edges)


def find(p, a):
    while p[a] != a:
        p[a] = p[p[a]]
        a = p[a]
    return a


def classes(tree):
    adj = [[] for _ in range(N)]
    for i in tree:
        u, v, dx, dy = edges[i]
        adj[u].append((v, dx, dy))
        adj[v].append((u, -dx, -dy))
    pos = [None] * N
    pos[0] = (0, 0)
    stack = [0]
    while stack:
        u = stack.pop()
        for v, dx, dy in adj[u]:
            if pos[v] is None:
                pos[v] = (pos[u][0] + dx, pos[u][1] + dy)
                stack.append(v)
    ts = set(tree)
    out = set()
    for i in range(E):
        if i in ts:
            continue
        u, v, dx, dy = edges[i]
        DX = pos[u][0] + dx - pos[v][0]
        DY = pos[u][1] + dy - pos[v][1]
        if DX == 0 and DY == 0:
            continue
        assert DY % n2 == 0
        b = DY // n2
        assert (DX - b * s) % n1 == 0
        a = (DX - b * s) // n1
        if a < 0 or (a == 0 and b < 0):
            a, b = -a, -b
        out.add((a, b))
    return tuple(sorted(out))


cnt = Counter()
ntree = 0
for sub in itertools.combinations(range(E), N - 1):
    p = list(range(N))
    ok = True
    for i in sub:
        u, v = edges[i][0], edges[i][1]
        ru, rv = find(p, u), find(p, v)
        if ru == rv:
            ok = False
            break
        p[ru] = rv
    if not ok:
        continue
    ntree += 1
    cnt[classes(sub)] += 1
print(f"{n1}x{n2} s{s}: {ntree} spanning trees")
for k, v in cnt.most_common():
    print(f"{v:8d} {v / ntree:.6f} {k}")
if len(sys.argv) > 4:
    json.dump({"n1": n1, "n2": n2, "s": s, "ntrees": ntree, "counts": {str(k): v for k, v in cnt.items()}},
              open(sys.argv[4], "w"), indent=1)
