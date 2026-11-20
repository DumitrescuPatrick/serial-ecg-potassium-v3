import numpy as np

def effect(c, n, seed=20261119):
    rng = np.random.default_rng(seed)
    cells = []
    for k in sorted(set(c.stratum.dropna()) & set(n.stratum.dropna())):
        a = c.loc[c.stratum == k, 'abs_change'].to_numpy()
        b = n.loc[n.stratum == k, 'abs_change'].to_numpy()
        if min(len(a), len(b)) >= 3:
            cells.append((k, a, b))
    if not cells:
        return dict(status='insufficient_common_support', cases=len(c), controls=len(n))
    total = sum((len(b) for _, a, b in cells))
    v = 0.0
    boot = np.zeros(2000)
    for k, a, b in cells:
        w = len(b) / total
        v += w * (a.mean() - b.mean())
        boot += w * (rng.choice(a, (2000, len(a)), replace=True).mean(1) - rng.choice(b, (2000, len(b)), replace=True).mean(1))
    return dict(status='estimated', effect=v, CI95=np.quantile(boot, [0.025, 0.975]).tolist(), cases=sum((len(a) for _, a, b in cells)), controls=total, strata=len(cells), cells=[dict(stratum=str(k), cases=len(a), controls=len(b), case_mean=float(a.mean()), control_mean=float(b.mean())) for k, a, b in cells])
