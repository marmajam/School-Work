import matplotlib.pyplot as plt
from scipy import special
from scipy.special import gamma
import numpy as np
import math
import cmath
from fractions import Fraction
from decimal import *

getcontext().prec = 50

x = np.linspace(-6, 4, 501)
x1 = np.linspace(-20, -6, 501)
x2 = np.linspace(4, 20, 501)

c1 = (3**(-2/3))/gamma(2/3)
c2 = (3**(-1/3))/gamma(1/3)
ai, aip, bi, bip = special.airy(x)
exp1 = Decimal(-0.25)
exp2 = Decimal(-0.5)

def cos_n(x):
    return Decimal(math.cos(zeta_n(x) +Decimal(math.pi/4.0)))

def sin_n(x):
    return Decimal(math.sin(zeta_n(x) +Decimal(math.pi/4.0)))

def factorial(n):
    if n <= 0:
        return 1
    else:
        return n*factorial(n-1)

def pochhammer(a,k):
    if k == 0:
        return 1
    else:
        return gamma(a +k)/gamma(a)

def c_k(k):
    if k == 0:
        return Decimal(1)
    else:
        return Decimal(gamma(3*k+ 1/2)/((54**k)*factorial(k)*gamma(k+ 1/2)))

def d_k(k):
    if k ==0:
        return 1
    else:
        return -((6*k +1)/(6*k -1))*c_k(k)

def zeta_p(x):
    return (2.0/3.0)*(x**1.5)

def zeta_n(x):
    return -Decimal(2.0/3.0)*(-Decimal(x)**Decimal(1.5))

def f_k(x,k):
    return ((3**k) * pochhammer(1/3, k)*(x**(3*k)))/(factorial(3*k))

def g_k(x,k):
    return ((3**k) * pochhammer(2/3, k)*(x**(3*k +1)))/(factorial(3*k +1))

def f_n(x,n):
    m = 0
    s = 0
    while s <= n:
        m = m + f_k(x,s)
        s += 1
    return m

def g_n(x,n):
    m=0
    s=0
    while s <= n:
        m = m + g_k(x,s)
        s += 1
    return m

def ai_asymp_pos_sum(x,k):
    return pow(Decimal(-1),Decimal(k))*c_k(k)*pow(zeta_p(x),-Decimal(k))

def ai_asymp_neg_sum1(x,k):
    return pow(Decimal(-1),Decimal(k))*c_k(2*k)*pow(zeta_n(x),Decimal(-2*k))

def ai_asymp_neg_sum2(x,k):
    return pow(Decimal(-1),Decimal(k))*c_k(2*k +1)*pow(zeta_n(x),Decimal(-(2*k)-1))

def ai_asymp_pos_x(x,n):
        m=0
        s=0
        while s <= n:
            pow_x = Decimal((x**exp1).real)
            pow_pi = Decimal((math.pi**exp2).real)
            m = m + Decimal(0.5)*pow_pi*pow_x*math.exp(-zeta_p(x))*ai_asymp_pos_sum(x,s)
            s += 1
        return m

def ai_asymp_neg_x(x,n):
    m=0
    s=0
    while s <= n:
        pow_x = Decimal((x**exp1).real)
        pow_pi = Decimal((math.pi**exp2).real)
        m = m + pow_pi*pow_x*(sin_n(x)*ai_asymp_neg_sum1(x,s) -cos_n(x)*ai_asymp_neg_sum2(x,s))
        s += 1
    return m

def ai_asymp_pos(x,n):
    return list(map(lambda x1:  ai_asymp_pos_x(x1,n),x))

def ai_asymp_neg(x,n):
    return list(map(lambda x1:  ai_asymp_neg_x(x1,n),x))

def ai_taylor_n(x,n):
    return c1*f_n(x,n) - c2*g_n(x,n)
        
def bi_taylor_n(x,n):
    return math.sqrt(3) * (c1*f_n(x,n) + c2*g_n(x,n))
    
"""
plt.plot(x, ai, 'r', label='Ai(x)')
plt.plot(x, bi, 'b--', label='Bi(x)')
"""
"""
plt.plot(x, bi_taylor_n(x,50), 'r', label='Bi_taylor(x)')
plt.plot(x, ai_taylor_n(x,50), 'b--', label='Ai_taylor(x)')
"""
plt.plot(x2, ai_asymp_pos(x2,30), 'b--', label='Ai_taylor(x)')
plt.plot(x1, ai_asymp_neg(x1,30), 'b--', label='Ai_taylor(x)')
plt.plot(x, ai_taylor_n(x,50), 'r--', label='Ai_taylor(x)')

plt.ylim(-5,5)
plt.xlim(-20,20)
plt.grid
plt.legend(loc='upper left')
plt.show()
