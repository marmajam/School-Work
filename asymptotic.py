import matplotlib.pyplot as plt
from scipy import special
from scipy.special import gamma
import numpy as np
import math
import cmath
from fractions import Fraction
import operator

x = np.linspace(-12.0, 4.0, 801)
x1 = np.linspace(4.0, -12.0, 801)


c1 = (3**(-2/3))/gamma(2/3)
c2 = (3**(-1/3))/gamma(1/3)
ai, aip, bi, bip = special.airy(x)
en14 = Fraction('-1/4')
en12 = Fraction('-1/2')
fr23 = Fraction('2/3')
fr32 = Fraction('3/2')
fr13 = Fraction('1/3')

def factorial(n):
    if n <= 0:
        return 1
    else:
        return n*factorial(n-1)

def gamma_k(k):
    if k <= 1:
        return 1
    else:
        return factorial(k-1)


def pochhammer(x,k):
    if k == 0:
        return 1
    else:
        return gamma_k(x +k)/gamma_k(x)


def c_k(k):
    if k == 0:
        return 1
    else:
        return gamma_k(3*k+ 1/2)/((54**k)*factorial(k)*gamma_k(k+ 1/2))


def d_k(k):
    if k ==0:
        return 1
    else:
        return -(6*k +1)*c_k(k)/(6*k -1)
  

def zeta(x):
    if x<=0:
        return -fr23*((-x)**fr32)
    else:
        return fr23*(x**fr32)





def cos_ai(x):
    return math.cos(zeta(x) +(3.0*math.pi/4.0))


def sin_ai(x):
    return math.sin(zeta(x) +(3.0*math.pi/4.0))

def cos_bi(x):
    return math.cos(zeta(x) -(1.0*math.pi/4.0))

def sin_bi(x):
    return math.sin(zeta(x) -(1.0*math.pi/4.0))


def ai_asymp_pos_sum(x,k):
    return pow(-1,k)*c_k(k)*pow(zeta(x),-k)

def bi_asymp_pos_sum(x,k):
    return c_k(k)*pow(zeta(x),-k)

def ai_asymp_neg_sum1(x,k):
    return pow(-1,k)*c_k(2*k)*pow(zeta(x),(-2*k))

def ai_asymp_neg_sum2(x,k):
    return pow(-1,k)*c_k(2*k +1)*pow(zeta(x),-(2*k+1))


def ai_asymp_pos_x(x,n):
        m=0
        s=0
        while s <= n:
            pow_x = abs(x**en14)
            pow_pi = abs(math.pi**en12)
            m = m + 0.5*pow_pi*pow_x*math.exp(-zeta(x))*ai_asymp_pos_sum(x,s)
            s += 1
        return m


def bi_asymp_pos_x(x,n):
        m=0
        s=0
        while s <= n:
            pow_x = abs(x**en14)
            pow_pi = abs(math.pi**en12)
            m = m +pow_pi*pow_x*math.exp(zeta(x))*bi_asymp_pos_sum(x,s)
            s += 1
        return m


def ai_asymp_neg_x(x,n):
    m=0
    s=0
    while s <= n:
        pow_x = abs(x**en14)
        pow_pi = abs(math.pi**en12)
        m = m + pow_pi*pow_x* (sin_ai(x)*ai_asymp_neg_sum1(x,s) -cos_ai(x)*ai_asymp_neg_sum2(x,s))
        s += 1
    return m

def bi_asymp_neg_x(x,n):
    m=0
    s=0
    while s <= n:
        pow_x = abs(x**en14)
        pow_pi = abs(math.pi**en12)
        m = m + pow_pi*pow_x*(cos_bi(x)*ai_asymp_neg_sum1(x,s) +sin_bi(x)*ai_asymp_neg_sum2(x,s))
        s += 1
    return m


def ai_asymp_pos(x,n):
    return list(map(lambda x1:  ai_asymp_pos_x(x1,n),x))


def ai_asymp_neg(x,n):
    return list(map(lambda x1:  ai_asymp_neg_x(x1,n),x))

def bi_asymp_pos(x,n):
    return list(map(lambda x1:  bi_asymp_pos_x(x1,n),x))


def bi_asymp_neg(x,n):
    return list(map(lambda x1:  bi_asymp_neg_x(x1,n),x))

plt.plot(x,ai_asymp_neg(x,10), 'b',label='Ai_Asymp')
plt.plot(x,bi_asymp_neg(x,10), 'r',label='Bi_Asymp')

plt.ylim(-5,5)
plt.xlim(-12,-2)
plt.grid(linestyle='-',linewidth=1.0)
plt.legend(loc='upper left')
plt.show()

plt.plot(x1,ai_asymp_pos(x1,10), 'b', label='Ai_Asymp')
plt.plot(x1,bi_asymp_pos(x1,10), 'r', label='Bi_Asymp')

plt.ylim(-5,5)
plt.xlim(2,12)
plt.grid(linestyle='-',linewidth=1.0)
plt.legend(loc='upper left')
plt.show()
