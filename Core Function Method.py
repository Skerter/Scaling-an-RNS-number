import numpy as np
from timeit import default_timer


def mul_inv(a, b):
    b0 = b
    x0, x1 = 0, 1
    if b == 1: return 1
    while a > 1:
        q = a // b
        a, b = b, a % b
        x0, x1 = x1 - q * x0, x0
    if x1 < 0: x1 += b0
    return x1


def calculate_P(p, n):
    P = 1
    for i in range(n):
        P *= p[0, i]
    return P


def calculate_P_i(p, n, P):
    P_i = np.zeros((1, n))
    for i in range(n):
        P_i[0, i] = P / p[0, i]
    return P_i


def orthogonal_bases(p, n, P_i):
    m = np.zeros((1, n))
    B = np.zeros((n, n))
    for i in range(n):
        m[0, i] = mul_inv(P_i[0, i], p[0, i])
        B[i, i] = m[0, i] * P_i[0, i]
    return B


def calculate_x_rns(A, n, p):
    x_rns = np.zeros((1, n))
    for i in range(n):
        x_rns[0, i] = A % p[0, i]
    return x_rns


def define_w_i(n):
    w_i = np.zeros((1, n))
    for i in range(n):
        w_i[0, i] = 1
    return w_i


def calculate_C_Bi(n, B, w_i, P_i, P, p, S_p):
    C_Bi = np.zeros((1, n))
    for i in range(n):
        for j in range(n):
            C_Bi[0, i] += w_i[0, j]*(B[i,i]//p[0,j])
    return C_Bi


def calculate_S_p(w_i, n, P_i):
    S_p = 0
    for i in range(n):
        S_p += w_i[0, i]*P_i[0, i]
    return S_p


def calculate_S_x(C_Bi, X, n, S_p):
    S_x = 0
    for i in range(n):
        S_x += X[0, i]*C_Bi[0, i]
    S_x = S_x % S_p
    return S_x


def sign_X(S_x, P_half):   
    if S_x < P_half:
        return 1
    elif S_x >= P_half:
        return 0


p = np.array([[15, 8, 17]], np.float64)
n = len(p[0])
P = calculate_P(p, n)
P_half = P / 2
P_i = calculate_P_i(p, n, P)
w_i = np.array([[0, 0, 1]], np.float64)

A = 47
X = calculate_x_rns(A, n, p)


t = default_timer()
B = orthogonal_bases(p, n, P_i)
S_p = calculate_S_p(w_i, n, P_i)
C_Bi = calculate_C_Bi(n, B, w_i, P_i, P, p, S_p)

S_x = calculate_S_x(C_Bi, X, n, S_p)
signX = sign_X(S_x, P_half)
print('{:.9f}'.format(default_timer() - t))
