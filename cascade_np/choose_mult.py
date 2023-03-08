import numpy as np

v1 = np.array([2, 5])
v2 = np.array([3, 4])
v3 = np.array([4, 3])
v4 = np.array([5, 1])

aa = np.array([v1, v2, v3, v4])
bb = np.amin(aa, axis = 0)

print(aa.min(0))
print(aa.argmin(0))
