import matplotlib.pyplot as plt
from scipy import special
from scipy.special import gamma
import numpy as np
import math


def factorial(n):
    if n <= 0:
        return 1
    else:
        return n*factorial(n-1)

def pochhammer(a1,k1):
    if k1 == 0:
        return 1
    else:
        return (gamma(a1+k1))*(1/(gamma(a1)))

def f_k(z,k):
    return ((3**k)*pochhammer(1/3,k)*(z**(3*k)))/(factorial(3*k))

def g_k(z,k):
    return ((3**k)*pochhammer(2/3,k)*(z**(3*k+1)))/(factorial(3*k +1))