import matplotlib
from matplotlib.patches import Wedge
import matplotlib.pyplot as plt
import numpy as np
import math


"""
fig=plt.figure()
ax=fig.add_subplot(111) 

s1 = Wedge((0,0), 0.6, -30, 30, color="b",hatch='-', alpha=0.5)
s2 = Wedge((0,0), 0.6, 90, 120, color="b",hatch='-', alpha=0.5)
s3 = Wedge((0,0), 0.6, -120, -90, color="b",hatch='-', alpha=0.5)

ax.add_artist(s1)
ax.add_artist(s2)
ax.add_artist(s3)

plt.autoscale(enable=True)

plt.show()
"""
x = np.linspace(-20.0, 20, 201)
x1 = np.linspace(1,20,501)
x2 = np.linspace(-20.0,0.0,101)
x3 = np.linspace(0,20,101)
def y1(x):
    return list(map(lambda x1:  math.tan(math.pi/3.0)*x1,x))

def y2(x):
    return list(map(lambda x1:  -math.tan(math.pi/3.0)*x1,x))

def y3(x):
    return list(map(lambda x1:  math.tan(math.pi)*x1,x))

def c1a(x):
    return list(map(lambda x1: math.tan(math.pi/3)*math.sqrt(x1**2-1),x))

def c1b(x):
    return list(map(lambda x1: -math.tan(math.pi/3)*math.sqrt(x1**2-1),x))

def c2a_x(x,theta):
    return list(map(lambda x1: x1*math.cos(theta) - math.tan(math.pi/3)*math.sqrt(x1**2-1)*math.sin(theta),x))

def c2a_y(x,theta):
    return list(map(lambda x1: x1*math.sin(theta) + math.tan(math.pi/3)*math.sqrt(x1**2-1)*math.cos(theta),x))

def c2b_x(x,theta):
    return list(map(lambda x1: -x1*math.cos(theta) + math.tan(math.pi/3)*math.sqrt(x1**2-1)*math.sin(theta),x))

def c2b_y(x,theta):
    return list(map(lambda x1: -x1*math.sin(theta) - math.tan(math.pi/3)*math.sqrt(x1**2-1)*math.cos(theta),x))

plt.plot(x3,y1(x3),color='k')
plt.plot(x3,y2(x3),color='k')
plt.plot(x2,y3(x2),color='k')

plt.text(3, 0.0, "$S_0$", fontdict=None, withdash=False)
plt.text(-2.25, 2.25, "$S_1$", fontdict=None, withdash=False)
plt.text(-2.25, -2.25, "$S_{-1}$", fontdict=None, withdash=False)


plt.plot(x1,c1a(x1),color='k', linestyle='dashed')
plt.plot(x1,c1b(x1),color='k', linestyle='dashed')

plt.plot(c2a_x(x1,2*math.pi/3),c2a_y(x1,2*math.pi/3),color='k', linestyle='dashed')
plt.plot(c2a_x(x1,-2*math.pi/3),c2a_y(x1,-2*math.pi/3),color='k', linestyle='dashed')

plt.plot(c2a_x(x1,2*math.pi/3),c2b_y(x1,2*math.pi/3),color='k', linestyle='dashed')
plt.plot(c2a_x(x1,-2*math.pi/3),c2b_y(x1,-2*math.pi/3),color='k', linestyle='dashed')

plt.arrow(1.0, -0.3, 0.0, 0.2, head_width=0.2, head_length=0.2, fc='k', ec='k')
plt.arrow(-0.2, 1.1, -0.1, -0.1, head_width=0.2, head_length=0.2, fc='k', ec='k')
plt.arrow(-0.2, -1.03, +0.1, -0.1, head_width=0.2, head_length=0.2, fc='k', ec='k')


plt.ylim(-5,5)
plt.xlim(-5,5)
plt.grid(linestyle='-',linewidth=1.0)
plt.legend(loc='upper left')
plt.show()
