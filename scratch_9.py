import numpy as np
from numpy import array, identity, diagonal
from math import sqrt
import time

a = np.array([[4., 2., 2., 1.], [2., -3., 1., 1.], [2., 1., 3., 1.], [1., 1., 1., 2.]])
print ("Матриця\n", a)

def jacobi(ain, tol=1.0e-9):

    def maxElem(a):  # Знаходження максимального недіагонального елемента a[k,1]
        n = len(a)
        aMax = 0.0
        for i in range(n - 1):
            for j in range(i + 1, n):
                if abs(a[i, j]) >= aMax:
                    aMax = abs(a[i, j])
                    k = i;
                    l = j
        return aMax, k, l

    def rotate(a, p, k, l):  # Ротація, щоб a[k,l] = 0
        n = len(a)
        aDiff = a[l, l] - a[k, k]
        if abs(a[k, l]) < abs(aDiff) * 1.0e-36:
            t = a[k, l] / aDiff
        else:
            phi = aDiff / (2.0 * a[k, l])
            t = 1.0 / (abs(phi) + sqrt(phi ** 2 + 1.0))
            if phi < 0.0: t = -t
        c = 1.0 / sqrt(t ** 2 + 1.0);
        s = t * c
        tau = s / (1.0 + c)
        temp = a[k, l]
        a[k, l] = 0.0
        a[k, k] = a[k, k] - t * temp
        a[l, l] = a[l, l] + t * temp
        for i in range(k):  # Випадок i < k
            temp = a[i, k]
            a[i, k] = temp - s * (a[i, l] + tau * temp)
            a[i, l] = a[i, l] + s * (temp - tau * a[i, l])
        for i in range(k + 1, l):  # Випадок k < i < l
            temp = a[k, i]
            a[k, i] = temp - s * (a[i, l] + tau * a[k, i])
            a[i, l] = a[i, l] + s * (temp - tau * a[i, l])
        for i in range(l + 1, n):  # Випадок i > l
            temp = a[k, i]
            a[k, i] = temp - s * (a[l, i] + tau * temp)
            a[l, i] = a[l, i] + s * (temp - tau * a[l, i])
        for i in range(n):  # Оновлення матриці
            temp = p[i, k]
            p[i, k] = temp - s * (p[i, l] + tau * p[i, k])
            p[i, l] = p[i, l] + s * (temp - tau * p[i, l])

    a = np.copy(ain)
    n = len(a)
    maxRot = 5 * (n ** 2)  # Ліміт кількості ротацій
    p = identity(n) * 1.0  # Ініціалізація транспонованої матриці
    for i in range(maxRot):  # Ротація Якобі
        aMax, k, l = maxElem(a)
        if aMax < tol: return diagonal(a), p
        rotate(a, p, k, l)
    print ('Метод Якобі не сходиться')


print ("\n---Метод Якобі:---\n")

tic = time.perf_counter()
wj, vj = jacobi(a)
tic2 = time.perf_counter()
print ("Час", tic2 - tic)
print ("Власні числа:\n", wj)
print ("Власні вектори:\n", vj)






