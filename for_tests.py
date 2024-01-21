import numpy as np


def calculate_x_rns(A, n, p):
    """
    Перевод числа в СОК
    :param A: Само число
    :param n: Количество модулей
    :param p: Модули СОКа
    :return: Число в СОК (массив)
    """
    x_rns = np.zeros((1, n))
    for i in range(n):
        x_rns[0, i] = A % p[0, i]
    return x_rns


def calculate_P(p, n):
    """
    Вычисление динамического диапазона
    :param p: Модули СОКа
    :param n: Количество модулей
    :return: P (динамический диапазон)
    """
    P = 1
    for i in range(n):
        P *= p[0, i]
    return P


def mul_inv(a, b):
    """
    Вычисление мультипликативной инверсии P_i^-1
    :param a: <НАЗВАНИЕ P_i> P/p_i
    :param b: p_i (i-й модуль СОКа)
    :return: P_i^-1 (i-ая мультипликативная инверсия, НЕ МАССИВ, А ОДИН ЭЛЕМЕНТ)
    TODO: переделать вычисление мультипликативной инверсии
    """
    b0 = b
    x0, x1 = 0, 1
    if b == 1: return 1
    while a > 1:
        q = a // b
        a, b = b, a % b
        x0, x1 = x1 - q * x0, x0
    if x1 < 0: x1 += b0
    return x1


def orthogonal_bases(p, n, P_i, P_i_1):
    """
    Вычисление ортогональных базисов
    :param p: Модули СОКа
    :param n: Количество модулей
    :param P_i: Название не знаю, это массив длиной n, где каждый элемент равен P/p_i
    :param P_i_1: Мультипликативные инверсии
    :return: массив B (ортогонольные базисы), матрица
    TODO: переделать принцип хранения ортогональных базисов
    """
    m = np.zeros((1, n))
    B = np.zeros((n, n))
    for i in range(n):
        m[0, i] = mul_inv(P_i[0, i], p[0, i])
        B[i, i] = P_i_1[0, i] * P_i[0, i]
    return B


def calculate_P_i(p, n, P):
    """
    Вычисление <НАЗВАНИЕ P_i>
    :param p: Модули СОКа
    :param n: Количество модулей
    :param P: Динамический диапазон
    :return: Массив P_i (<НАЗВАНИЕ P_i>)
    TODO: узнать название P_i
    """
    P_i = np.zeros((1, n))
    for i in range(n):
        P_i[0, i] = P / p[0, i]
    return P_i


def calculate_C_P(w_i, P_i, n):
    """
    Вычисление значения функции ядра Акушского
    :param w_i: Веса
    :param P_i: <НАЗВАНИЕ P_i>
    :param n: Количество модулей
    :return: C(P)
    """
    C_P = 0
    for i in range(n):
        C_P += w_i[0, i] * P_i[0, i]
    return C_P


def calculate_P_i_1_mod_p_i(P_i, p_i):
    """
    Вычисление значения мультипликативной инверсии по модулю p_i
    :param P_i: <НАЗВАНИЕ P_i> P/p_i
    :param p_i: i-й модуль СОКа
    :return: y (мультипликативная инверсия)
    """
    for y in range(9999):
        if (y * P_i - 1) % p_i == 0:
            return y


def calculate_C_Bi(P, p, w_i, n, C_P):
    """
    Вычисление формулы (12) из статьи
    :param P: Динамический диапазон
    :param p: Модули СОКа
    :param w_i: Веса
    :param n: Количество модулей
    :param C_P: Значение функции ядра Акушского
    :return: C(Bi)
    """
    C_Bi = np.zeros((1, n))
    for i in range(n):
        C_Bi[0, i] = (C_P * calculate_P_i_1_mod_p_i(P // p[0, i], p[0, i]) - w_i[0, i]) / p[0, i]
    return C_Bi


p = np.array([[7, 11, 13, 17, 19, 23]], np.float64)  # Модули СОКа
n = len(p[0])  # Количество модулей
w_i = np.array([[1, 1, 1, 1, 1, 1]], np.float64)  # Веса, которые я взял рандомно

P = calculate_P(p, n)         # Динамический диапазон
print("P:", P)
P_i = calculate_P_i(p, n, P)  # <НАЗВАНИЕ P_i>
print("P_i:", P_i)
P_i_1 = np.zeros((1, n))
for i in range(n):  # Вычисление мультипликативных инверсий
    P_i_1[0, i] = mul_inv(P_i[0, i], p[0, i])
print("P_i^-1:", P_i_1)
C_P = calculate_C_P(w_i, P_i, n)
B = orthogonal_bases(p, n, P_i, P_i_1)  # Вычисление ортогональных базисов
C_Bi = calculate_C_Bi(P, p, w_i, n, C_P)

P_J = np.array([[7, 17, 23], [0, 3, 5]])  # TODO: Понять принцип разбиения динамического диапазона на 2 набора
P_K = np.array([[11, 13, 19], [1, 2, 4]])
C_J_P = 2737
w_i1 = np.array([[0, -2, 1, 0, 2, 0]], np.float64)   # Веса для C_J(P)
C_K_P = 2717
w_i2 = np.array([[-1, 0, 0, -2, 0, 6]], np.float64)  # Веса для C_K(P) TODO: Понять, как берутся эти веса

C_J_Bi = calculate_C_Bi(P, p, w_i1, n, C_J_P)  # C_J(Bi)
C_K_Bi = calculate_C_Bi(P, p, w_i2, n, C_K_P)  # C_K(Bi)
print("C_J(Bi):", *C_J_Bi)
print("C_K(Bi):", *C_K_Bi)

dC_Bi = C_J_Bi - C_K_Bi  # dC(Bi)
print("dC(Bi):", *dC_Bi)

# Пример 1 (2.3)
print("-----Пример 1 (2.3)-----")
X = 1859107  # Это число в примере масштабируется на 2717
print("X:", X)
X_rns = calculate_x_rns(X, n, p)  # Перевод числа Х в СОК
print("X in RNS:", *X_rns)

C_J_X_rns = np.zeros((1, 6))
for i in P_J[1]:
    elem = X_rns * C_J_Bi
    print(f"C_J(X)mod{int(p[0, i])}: {elem.sum() % p[0, i]}")
    C_J_X_rns[0, i] = elem.sum() % p[0, i]

dC_X_mod20 = X_rns * dC_Bi
dC_X_mod20 = dC_X_mod20.sum() % 20
print("dC(X)mod20:", dC_X_mod20)

for i in P_K[1]:
    elem = X_rns * C_K_Bi
    print(f"C_K(X)mod{int(p[0, i])}: {elem.sum() % p[0, i]}")
    elem2 = (elem.sum() % p[0, i]) + (dC_X_mod20 % p[0, i])
    print(f"C_J(X)mod{int(p[0, i])}: {elem2.sum() % p[0, i]}")
    C_J_X_rns[0, i] = elem2.sum() % p[0, i]

print('C_J(X) in RNS:', *C_J_X_rns)

C_J_X = C_J_X_rns * P_i * P_i_1
C_J_X = C_J_X.sum()
print("C_J(X):", C_J_X % P)
