import math
from math import exp
from math import e
import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


def create_data(a, b):
    m = [0 for x in range(0, 1 << 2)]
    Bel = [0 for x in range(0, 1 << 2)]
    Pl = [0 for x in range(0, 1 << 2)]
    Q = [0 for x in range(0, 1 << 2)]

    m[2] = a
    m[1] = b
    m[3] = 1-a-b

    Bel[1] = m[1]
    Bel[2] = m[2]
    Bel[3] = 1

    Pl[0] = 1-Bel[3]
    Pl[1] = 1-Bel[2]
    Pl[2] = 1-Bel[1]
    Pl[3] = 1-Bel[0]

    Q[0] = 1
    Q[1] = Bel[3]
    Q[2] = Bel[3]
    Q[3] = Bel[3]

    return m, Bel, Pl, Q


def calculate_capacity(num):
    ans = 0
    for iter in range(0, 20):
        if num & (1 << iter):
            ans += 1
    return ans


def definition1(m, Bel, Pl, Q, num):
    C_H = 0
    for iter in range(0, 1 << num):
        if m[iter] != 0 and Bel[iter] != 0:
            C_H -= m[iter]*np.log2(Bel[iter])
    return C_H


def definition2(m, Bel, Pl, Q, num):
    H_S = 0
    for iter in range(0, 1 << num):
        if m[iter] != 0:
            if Q[iter] != 0:
                H_S -= m[iter]*np.log2(Q[iter])
    return H_S


def definition3(m, Bel, Pl, Q, num):
    E_Y = 0
    for iter in range(0, 1 << num):
        if m[iter] != 0 and Pl[iter] != 0:
            E_Y -= m[iter]*np.log2(Pl[iter])
    return E_Y


def definition4(m, Bel, Pl, Q, num):
    l_DP = 0
    for iter in range(0, 1 << num):
        if m[iter] != 0:
            l_DP += m[iter]*np.log2(calculate_capacity(iter))
    return l_DP


def definition6(m, Bel, Pl, Q, num):
    D_KR = 0
    for iter in range(0, 1 << num):
        D_KR_TEMP = 0
        if m[iter] != 0:
            for j in range(1, 1 << 15):
                if m[j] != 0:
                    D_KR_TEMP += m[j]*calculate_capacity(iter & j)/calculate_capacity(j)
            if D_KR_TEMP != 0:
                D_KR -= m[iter]*np.log2(D_KR_TEMP)
    return D_KR


def definition7(m, Bel, Pl, Q, num):
    D_KP = 0
    for iter in range(1, 1 << num):
        D_KP_TEMP = 0
        if m[iter] != 0:
            for j in range(0, 1 << num):
                if m[j] != 0:
                    D_KP_TEMP += m[j]*calculate_capacity(iter & j)/calculate_capacity(iter)
            if D_KP_TEMP != 0:
                D_KP -= m[iter]*np.log2(D_KP_TEMP)
    return D_KP


def definition8(m, Bel, Pl, Q, num):
    H_P = 0
    for iter in range(0, 1 << num):
        if m[iter] != 0:
            H_P -= m[iter]*np.log2(m[iter]/calculate_capacity(iter))
    return H_P


def definition9(m, Bel, Pl, Q, num):
    D_GP = 0
    for iter in range(1, 1 << num):
        D_GP_TEMP = 0
        if m[iter] != 0:
            for j in range(0, 1 << num):
                if (m[j] != 0) and iter != 0:
                    D_GP_TEMP += m[j]*(1-calculate_capacity(iter & j)/calculate_capacity(iter | j))
            if D_GP_TEMP != 0:
                D_GP -= m[iter]*np.log2(D_GP_TEMP)
    return D_GP


def definition10(m, Bel, Pl, Q, num):
    H_J = 0
    for jter in range(1, num + 1):
        iter = 1 << (jter-1)
        Bet_m = 0
        for j in range(0, 1 << num):
            if m[j] != 0 and iter | j == j:
                Bet_m += m[j]/calculate_capacity(j)
        if Bet_m != 0:
            H_J -= Bet_m*np.log2(Bet_m)
    return H_J


def definition11(m, Bel, Pl, Q, num):
    H_JS = 0
    K = 0
    for iter in range(1, num + 1):
        K += Pl[1 << (iter-1)]
    K = 1/K
    for iter in range(1, num+1):
        if Pl[1 << (iter-1)] != 0:
            H_JS -= K*Pl[1 << (iter-1)]*np.log2(K*Pl[1 << (iter-1)])
    return H_JS


def definition12(m, Bel, Pl, Q, num):
    H_PQ = 0
    K = 0
    for iter in range(1, num+1):
        K += Pl[1 << (iter-1)]
    K = 1/K
    for iter in range(1, 1 << num):
        Pm = 0
        for j in range(1, num+1):
            temp = 1 << (j-1)
            if temp | iter == iter:
                Pm += K*Pl[iter]
        if m[iter] != 0 and Pm != 0:
            H_PQ -= m[iter]*Pm
    return H_PQ


def definition13(m, Bel, Pl, Q, num):
    e = math.e
    U_exp = e
    Temp = e-(1/num**2)*exp(1/num**2)
    C = 0
    for iter in range(1, 1 << num):
        if m[iter] > 0:
            C = C | iter
    C = calculate_capacity(C)
    for iter in range(1, 1 << num):
        if m[iter] != 0:
            if calculate_capacity(iter) == 1:
                U_exp -= m[iter]*exp(m[iter])
            else:
                if calculate_capacity(iter)*C != 0:
                    temp = m[iter]/(calculate_capacity(iter)*C)
                    U_exp -= temp*exp(temp)
    U_exp = U_exp/Temp
    return U_exp


def definition14(m, Bel, Pl, Q, num):
    SU_SW = 0
    for jter in range(1, num+1):
        iter = 1 << (jter-1)
        SU_SW += -((Bel[iter]+Pl[iter])/2)*(np.log2((Bel[iter]+Pl[iter])/2))+(Pl[iter]-Bel[iter])/2
    return SU_SW


def d_I(a1, b1, a2, b2):
    ans = np.sqrt(((a1+b1-a2-b2)/2)**2+1/3*((b1-a1-b2+a2)/2)**2)
    return ans


def definition15(m, Bel, Pl, Q, num):
    TU = 1
    X = np.sqrt(3)/num
    for jter in range(1, num+1):
        iter = 1 << (jter-1)
        TU -= X*d_I(Bel[iter], Pl[iter], 0, 1)
    return TU


def d_E(a1, b1, a2, b2):
    ans = np.sqrt((a1-a2)**2+(b1-b2)**2)
    return ans


def definition16(m, Bel, Pl, Q, num):
    TU_E = 0
    for jter in range(1, num+1):
        iter = 1 << (jter-1)
        TU_E += (1-d_E(Bel[iter], Pl[iter], 0, 1))
    return TU_E


def Deng_f(m, Bel, Pl, Q, num):
    E_D = 0
    for iter in range(1, 1 << num+1):
        if m[iter] != 0:
            E_D -= m[iter]*np.log2(m[iter]/(2**calculate_capacity(iter)-1))
    return E_D


# def get(num):

#     m = [0 for x in range(0, 1 << num)]
#     Bel = [0 for x in range(0, 1 << num)]
#     Pl = [0 for x in range(0, 1 << num)]
#     Q = [0 for x in range(0, 1 << num)]

#     for iter in range(1, num + 1):
#         m[1 << (iter-1)] = 1/num

#     for iter in range(0, 1 << num):
#         number, fre, Type, value = input().split(' ')
#         Bel[eval(fre)] = eval(value)
#     for iter in range(0, 1 << num):
#         number, fre, Type, value = input().split(' ')
#         Pl[eval(fre)] = eval(value)
#     for iter in range(0, 1 << num):
#         number, fre, Type, value = input().split(' ')
#         Q[eval(fre)] = eval(value)

#     return m, Bel, Pl, Q


if __name__ == '__main__':
    Z = np.zeros((100, 100))
    # Z = np.array(Z)
    for i in range(0, 100):
        for j in range(0, 100):
            m, Bel, Pl, Q = create_data(i*0.05, j*0.05)
            Z[i][j] = definition1(m, Bel, Pl, Q, 2)
    fig = plt.figure()
    ax = Axes3D(fig)
    X = np.arange(0, 0.5, 0.005)
    Y = np.arange(0, 0.5, 0.005)
    X, Y = np.meshgrid(X, Y)
    ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap='rainbow')

    plt.draw()
    plt.show()
    # ax = plt.subplot(projection='3d')  # 创建一个三维的绘图工程

    # ax.scatter(X, Y, Z, c='r')  # 绘制数据点,颜色是红色

    # ax.set_zlabel('Z')  # 坐标轴
    # ax.set_ylabel('Y')
    # ax.set_xlabel('X')

    # plt.draw()
    # plt.show()
