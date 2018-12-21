import numpy as np
xs = np.linspace(0.0, 2.0, 10)
ys = np.linspace(0.0, 1.0, 10)
ts = np.linspace(0.0, 10.0, 10)
x, y, t = np.meshgrid(xs, ys, ts, indexing='ij')
for k in range(10):
    for j in range(10):
        for i in range(10):
            px = x[i, j, k]
            py = y[i, j, k]
            pt = t[i, j, k]
            print(px + i)
            print(py + j)
            print(pt + k)