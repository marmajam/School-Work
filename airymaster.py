import matplotlib.pyplot as plt
from scipy import special
from scipy.special import gamma
import numpy as np
import math

x = np.linspace(-15, 15, 201)
c1 = 0.355028053887817
c2 = 0.258819403792807
ai, aip, bi, bip = special.airy(x)

def factorial(n):
    if n <= 0:
        return 1
    else:
        return n*factorial(n-1)

def pochhammer(a,k):
    if k == 0:
        return 1
    else:
        return gamma(a+k)/gamma(a)
    
def f(x):
    return 1 + (1/factorial(3))*(x**3) + ((1*4)/factorial(6))*x**(6) + ((1*4*7)/factorial(9))*(x**9)
     
def g(x):
    return x + (2/factorial(4))*(x**4) + ((2*5)/factorial(7))*(x**7) + ((2*5*8)/factorial(10))*(x**10)

def f_k(x,k):
    return ((3**k)*pochhammer(1/3,k)*(x**(3*k)))/(factorial(3*k))

def g_k(x,k):
    return ((3**k)*pochhammer(2/3,k)*(x**(3*k+1)))/(factorial(3*k +1))


def f_n(x,n):
    m = 0
    s = 0
    while s <= n:
        m = m + f_k(x,s)
        s = s + 1
    return m

def g_n(x,n):
    m=0
    s=0
    while s <= n:
        m = m + g_k(x,s)
        s = s +1
    return m

def ai_taylor(x):
    return c1*f(x) - c2*g(x)
        
def bi_taylor(x):
    return math.sqrt(3) * (c1*f(x) + c2*g(x))

def ai_taylor_n(x,n):
    return c1*f_n(x,n) - c2*g_n(x,n)
        
def bi_taylor_n(x,n):
    return math.sqrt(3) * (c1*f_n(x,n) + c2*g_n(x,n))


    
"""
plt.plot(x, ai, 'r', label='Ai(x)')
plt.plot(x, bi, 'b--', label='Bi(x)')
"""

plt.plot(x, bi_taylor_n(x,50), 'r', label='Bi_taylor(x)')
plt.plot(x, ai_taylor_n(x,50), 'b--', label='Ai_taylor(x)')

plt.ylim(-0.5,1.0)
plt.xlim(-20.0,20.0)
plt.grid
plt.legend(loc='upper left')
plt.show()
