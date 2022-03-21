from cmath import exp
import numpy as np
from matplotlib import pyplot as plt
import math


def create_data(num):

    m = [0 for x in range(0, 1 << (num-1))]
    Bel = [0 for x in range(0, 1 << (num-1))]
    Pl = [0 for x in range(0, 1 << (num-1))]
    Q = [0 for x in range(0, 1 << (num-1))]

    for iter in range(1, num):
        m[1 << (iter-1)] = 1/num

    for i in range(0, 1 << (num-1)):
        print(num, bin(i), 'Bel', end=' ')
        for j in range(0, i+1):
            if i | j == i:
                Bel[i] += m[j]
        print(Bel[i])

    for i in range(0, 1 << (num-1)):
        print(num, bin(i), 'Pl', end=' ')
        Pl[i] = 1-Bel[~i]
        print(Pl[i])

    for i in range(0, 1 << (num-1)):
        print(num, bin(i), 'Q', end=' ')
        for j in range(i, 1 << (num-1)):
            if i | j == j:
                Q[i] += m[j]
        print(Q[i])

    return m, Bel, Pl, Q


def calculate_capacity(num):
    ans = 0
    for iter in range(0, 15):
        if num & (1 << iter):
            ans += 1
    return ans


def definition1(m, Bel, Pl, Q):
    C_H = 0
    for iter in range(0, 1 << 15):
        if m[iter] != 0:
            if Bel[iter] != 0:
                C_H -= m[iter]*np.log2(Bel[iter])
    return C_H


def definition2(m, Bel, Pl, Q):
    H_S = 0
    for iter in range(0, 1 << 15):
        if m[iter] != 0:
            if Q[iter] != 0:
                H_S -= m[iter]*np.log2(Q[iter])
    return H_S


def definition3(m, Bel, Pl, Q):
    E_Y = 0
    for iter in range(0, 1 << 15):
        if m[iter] != 0:
            if Pl[iter] != 0:
                E_Y -= m[iter]*np.log2(Pl[iter])
    return E_Y


def definition4(m, Bel, Pl, Q):
    l_DP = 0
    for iter in range(0, 1 << 15):
        if m[iter] != 0:
            l_DP += m[iter]*np.log2(calculate_capacity(iter))
    return l_DP


def definition6(m, Bel, Pl, Q):
    D_KR = 0
    for iter in range(0, 1 << 15):
        D_KR_TEMP = 0
        for j in range(1, 1 << 15):
            if m[j] != 0:
                D_KR_TEMP += m[j]*calculate_capacity(iter & j)/calculate_capacity(j)
        print('6', D_KR_TEMP)
        if D_KR_TEMP != 0:
            D_KR -= m[iter]*np.log2(D_KR_TEMP)
    return D_KR


def definition7(m, Bel, Pl, Q):
    D_KP = 0
    for iter in range(1, 1 << 15):
        D_KP_TEMP = 0
        for j in range(0, 1 << 15):
            if m[j] != 0:
                D_KP_TEMP += m[j]*calculate_capacity(iter & j)/calculate_capacity(iter)
        print('7', D_KP_TEMP)
        if D_KP_TEMP != 0:
            D_KP -= m[iter]*np.log2(D_KP_TEMP)
    return D_KP


def definition8(m, Bel, Pl, Q):
    H_P = 0
    for iter in range(0, 1 << 15):
        if m[iter] != 0:
            H_P -= m[iter]*np.log2(m[iter]/calculate_capacity(iter))
    return H_P


def definition9(m, Bel, Pl, Q):
    D_GP = 0
    for iter in range(1, 1 << 15):
        D_GP_TEMP = 0
        for j in range(0, 1 << 15):
            if (m[j] != 0) and iter != 0:
                D_GP_TEMP += m[j]*(1-calculate_capacity(iter & j)/calculate_capacity(iter | j))
        print('9', D_GP_TEMP)
        if D_GP_TEMP != 0:
            D_GP -= m[iter]*np.log2(D_GP_TEMP)
    return D_GP


def definition10(m, Bel, Pl, Q):
    H_J = 0
    for jter in range(1, 16):
        iter = 1 << (jter-1)
        H_J_TEMP = 0
        for j in range(0, 1 << 15):
            if m[j] != 0:
                H_J_TEMP += m[j]/calculate_capacity(j)
        print('10', H_J_TEMP)
        if H_J_TEMP != 0:
            H_J -= H_J_TEMP*np.log2(H_J_TEMP)
    return H_J


def definition11(m, Bel, Pl, Q):
    H_JS = 0
    K = 0
    for iter in range(1, 16):
        K += Pl[1 << (iter-1)]
    K = 1/K
    for iter in range(1, 16):
        if Pl[1 << (iter-1)] != 0:
            H_JS -= K*Pl[1 << (iter-1)]*np.log2(K*Pl[1 << (iter-1)])
    # for iter in range(1, 1 << 15):
    #     if m[iter] != 0:
    #         H_JS += m[iter]*np.log2(calculate_capacity(iter))
    return H_JS


def definition12(m, Bel, Pl, Q):
    H_PQ = 0
    K = 0
    for iter in range(1, 16):
        K += Pl[1 << (iter-1)]
    K = 1/K
    for iter in range(1, 1 << 15):
        Pm = 0
        for j in range(1, 16):
            temp = 1 << (j-1)
            if temp | iter == iter:
                Pm += K*Pl[iter]
        if m[iter] != 0 and Pm != 0:
            H_PQ -= m[iter]*Pm
    return H_PQ


def definition13(m, Bel, Pl, Q):
    e = math.e
    U_exp = e
    Temp = e-(1/15**2)*exp(1/15**2)
    C = 0
    for iter in range(1, 1 << 15):
        if m[iter] > 0:
            C = C | iter
    C = calculate_capacity(C)
    for iter in range(1, 1 << 15):
        if calculate_capacity(iter) == 1:
            U_exp -= m[iter]*exp(m[iter])
        else:
            if calculate_capacity(iter)*C != 0:
                temp = m[iter]/(calculate_capacity(iter)*C)
                U_exp -= temp*exp(temp)
    U_exp = U_exp/Temp
    return U_exp


def definition14(m, Bel, Pl, Q):
    SU_SW = 0
    for jter in range(1, 16):
        iter = 1 << (jter-1)
        SU_SW += -(Bel[iter]+Pl[iter])/2*np.log2((Bel[iter]+Pl[iter])/2)+(-Bel[iter]+Pl[iter])/2
    return SU_SW


def d_I(a1, b1, a2, b2):
    ans = np.sqrt(((a1+b1-a2-b2)/2)**2+1/3*((b1-a1-b2+a2)/2)**2)
    return ans


def definition15(m, Bel, Pl, Q):
    TU = 1
    X = np.sqrt(3)/15
    for jter in range(1, 16):
        iter = 1 << (jter-1)
        TU -= X*d_I(Bel[iter], Pl[iter], 0, 1)
    return TU


def d_E(a1, b1, a2, b2):
    ans = np.sqrt((a1-a2)**2+(b1-b2)**2)
    return ans


def definition16(m, Bel, Pl, Q):
    TU_E = 0
    for jter in range(1, 16):
        iter = 1 << (jter-1)
        TU_E += 1-d_E(Bel[iter], Pl[iter], 0, 1)
    return TU_E


def Deng_f(m, Bel, Pl, Q, num):
    E_D = 0
    for iter in range(1, 1 << (num-1)):
        if m[iter] != 0:
            E_D -= m[iter]*np.log2(m[iter]/(2**calculate_capacity(iter)-1))
    return E_D


if __name__ == '__main__':
    # print('熵')
    A = [0 for x in range(1, 21)]
    Uncertainly1 = [0 for x in range(1,  21)]
    Uncertainly2 = [0 for x in range(1,  21)]
    Uncertainly3 = [0 for x in range(1,  21)]
    Uncertainly4 = [0 for x in range(1,  21)]
    Uncertainly5 = [0 for x in range(1,  21)]
    Uncertainly6 = [0 for x in range(1,  21)]
    Uncertainly7 = [0 for x in range(1,  21)]
    Uncertainly8 = [0 for x in range(1,  21)]
    Uncertainly9 = [0 for x in range(1,  21)]
    Uncertainly10 = [0 for x in range(1, 21)]
    Uncertainly11 = [0 for x in range(1, 21)]
    Uncertainly12 = [0 for x in range(1, 21)]
    Uncertainly13 = [0 for x in range(1, 21)]
    Uncertainly14 = [0 for x in range(1, 21)]
    Uncertainly15 = [0 for x in range(1, 21)]
    Uncertainly16 = [0 for x in range(1, 21)]
    Deng = [0 for x in range(1, 21)]
    for iter in range(1, 21):
        A[iter-1] = iter
        m, Bel, Pl, Q = create_data(iter)
        # Uncertainly1[iter-1] = definition1(m, Bel, Pl, Q)
        # Uncertainly2[iter-1] = definition2(m, Bel, Pl, Q)
        # Uncertainly3[iter-1] = definition3(m, Bel, Pl, Q)
        # Uncertainly4[iter-1] = definition4(m, Bel, Pl, Q)
        # Uncertainly5[iter-1] = Uncertainly3[iter-1]+Uncertainly4[iter-1]
        # Uncertainly6[iter-1] = definition6(m, Bel, Pl, Q)
        # Uncertainly7[iter-1] = definition7(m, Bel, Pl, Q)
        # Uncertainly8[iter-1] = definition8(m, Bel, Pl, Q)
        # Uncertainly9[iter-1] = definition9(m, Bel, Pl, Q)
        # Uncertainly10[iter-1] = definition10(m, Bel, Pl, Q)
        # Uncertainly11[iter-1] = definition11(m, Bel, Pl, Q) + Uncertainly4[iter-1]
        # Uncertainly12[iter-1] = definition12(m, Bel, Pl, Q) + Uncertainly4[iter-1]
        # Uncertainly13[iter-1] = definition13(m, Bel, Pl, Q)
        # Uncertainly14[iter-1] = definition14(m, Bel, Pl, Q)
        # Uncertainly15[iter-1] = definition15(m, Bel, Pl, Q)
        # Uncertainly16[iter-1] = definition16(m, Bel, Pl, Q)
        Deng[iter-1] = Deng_f(m, Bel, Pl, Q, iter)

    plt.figure(figsize=(12, 8), dpi=72)
    plt.ylim(0, 20)
    plt.xlim(0, 5)
    plt.xlabel("A")
    plt.ylabel("Uncertainly")
    plt.tick_params(axis='both')
    # plt.plot(A, Uncertainly1, color='r', label='definition1')
    # plt.plot(A, Uncertainly2, color='g', label='definition2')
    # plt.plot(A, Uncertainly3, color='b', label='definition3')
    # plt.plot(A, Uncertainly4, color='y', label='definition4')
    # plt.plot(A, Uncertainly5, color='k', label='definition5')
    # plt.plot(A, Uncertainly6, '-.', color='r', label='definition6')
    # plt.plot(A, Uncertainly7, '-.', color='g', label='definition7')
    # plt.plot(A, Uncertainly8, '-.', color='b', label='definition8')
    # plt.plot(A, Uncertainly9, '-.', color='y', label='definition9')
    # plt.plot(A, Uncertainly10, '-.', color='k', label='definition10')
    # plt.plot(A, Uncertainly11, '--', color='r', label='definition11')
    # plt.plot(A, Uncertainly12, '--', color='g', label='definition12')
    # plt.plot(A, Uncertainly13, '--', color='b', label='definition13')
    # plt.plot(A, Uncertainly14, '--', color='y', label='definition14')
    # plt.plot(A, Uncertainly15, '--', color='k', label='definition15')
    # plt.plot(A, Uncertainly16, ':', color='g', label='definition16')
    plt.plot(A, Deng, ':', color='r', label='Deng')

    plt.legend()
    plt.show()
