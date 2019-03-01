import matplot.pyplot as plt
from scipy import special
import numpy as np


ai, aip, bi, bip = special.airy(x)
x = linspace(-15, 5, 201)
plt.plot(x, ai, 'r', label='Ai(x)')
plt.plot(x, bi, 'b--', label='Bi(x)')
plt.ylim(-0.5,1.0)
plt.grid
plt.legend(loc='upper left')
plt.show()
