import numpy as np
import matplotlib.pyplot as plt
import sys

data = np.loadtxt(sys.argv[1])
t    = data[:, 0]
x1   = data[:, 1]
x2   = data[:, 2]

plt.plot(t, x1, label='x1')
plt.plot(t, x2, label='x2')

plt.legend()
plt.show()

