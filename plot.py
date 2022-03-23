import os
import sys
import getopt
import numpy as np
from math import e
from math import exp
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


def create_data1(num):
    N1 = 0b000000000011100
    N2 = 0b000000000100000
    N3 = 0b111111111111111
    N4 = 0
    for i in range(0, num):
        N4 = N4 + (1 << i)

    m = [0 for x in range(0, 1 << 15)]
    Bel = [0 for x in range(0, 1 << 15)]
    Pl = [0 for x in range(0, 1 << 15)]
    Q = [0 for x in range(0, 1 << 15)]

    m[N1] = 0.05
    m[N2] = 0.05
    m[N3] = 0.1
    m[N4] = 0.8

    for i in range(0, 1 << 15):
        print(bin(N4), bin(i), 'Bel', end=' ')
        for j in range(0, i+1):
            if i | j == i:
                Bel[i] += m[j]
        print(Bel[i])

    for i in range(0, 1 << 15):
        print(bin(N4), bin(i), 'Pl', end=' ')
        Pl[i] = 1-Bel[~i]
        print(Pl[i])

    for i in range(0, 1 << 15):
        print(bin(N4), bin(i), 'Q', end=' ')
        for j in range(i, 1 << 15):
            if i | j == j:
                Q[i] += m[j]
        print(Q[i])

    return m, Bel, Pl, Q


def create_data2(num):

    m = [0 for x in range(0, 1 << num)]
    Bel = [0 for x in range(0, 1 << num)]
    Pl = [0 for x in range(0, 1 << num)]
    Q = [0 for x in range(0, 1 << num)]

    for iter in range(1, num + 1):
        m[1 << (iter-1)] = 1/num

    for jter in range(1, num + 1):
        i = 1 << (jter-1)
        Bel[i] = 1/num
        i = ~i
        Bel[i] = (num-1)/num

    for jter in range(1, num+1):
        iter = 1 << (jter-1)
        Pl[iter] = 1-Bel[~iter]

    for jter in range(1, num+1):
        iter = 1 << (jter-1)
        Q[iter] = m[iter]
    return m, Bel, Pl, Q


def create_data3(a, b):
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

    Q[0] = m[1]+m[2]+m[3]+m[0]
    Q[1] = m[1]+m[3]
    Q[2] = m[2]+m[3]
    Q[3] = m[3]

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


def definition5(m, Bel, Pl, Q, num):
    E_Y = 0
    for iter in range(0, 1 << num):
        if m[iter] != 0 and Pl[iter] != 0:
            E_Y -= m[iter]*np.log2(Pl[iter])
    l_DP = 0
    for iter in range(0, 1 << num):
        if m[iter] != 0:
            l_DP += m[iter]*np.log2(calculate_capacity(iter))
    return E_Y+l_DP


def definition6(m, Bel, Pl, Q, num):
    D_KR = 0
    for iter in range(0, 1 << num):
        D_KR_TEMP = 0
        if m[iter] != 0:
            for j in range(1, 1 << num):
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
    l_DP = 0
    for iter in range(0, 1 << num):
        if m[iter] != 0:
            l_DP += m[iter]*np.log2(calculate_capacity(iter))
    return H_JS+l_DP


def definition12(m, Bel, Pl, Q, num):
    H_PQ = 0
    K = 0
    for iter in range(1, num + 1):
        K += Pl[1 << (iter-1)]
    K = 1/K
    for iter in range(1, 1 << num):
        Pm = 0
        for j in range(1, num + 1):
            temp = 1 << (j-1)
            if temp | iter == iter:
                Pm += K*Pl[temp]
        if m[iter] != 0 and Pm != 0:
            H_PQ -= m[iter]*Pm
    l_DP = 0
    for iter in range(0, 1 << num):
        if m[iter] != 0:
            l_DP += m[iter]*np.log2(calculate_capacity(iter))
    return H_PQ+l_DP


def definition13(m, Bel, Pl, Q, num):
    if num == 1:
        return 0
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
    for iter in range(1, 1 << num):
        if m[iter] != 0:
            E_D -= m[iter]*np.log2(m[iter]/(2**calculate_capacity(iter)-1))
    return E_D


def plot3(fun, num):
    print('case3', num)
    Z = np.zeros((100, 100))
    for i in range(0, 100):
        for j in range(0, 100):
            m, Bel, Pl, Q = create_data3(i*0.005, j*0.005)
            Z[i][j] = fun(m, Bel, Pl, Q, 2)
    fig = plt.figure()
    ax = Axes3D(fig)
    X = np.arange(0, 0.5, 0.005)
    Y = np.arange(0, 0.5, 0.005)
    X, Y = np.meshgrid(X, Y)
    ax.invert_yaxis()
    ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap='rainbow')
    ax.set_xlabel('a')
    ax.set_ylabel('b')
    ax.set_zlabel('Uncertainty')

    plt.draw()
    try:
        plt.savefig('./result/case3_pic/' + num + '.png')
    except:
        os.system('mkdir -p ./result/case3_pic')
        plt.savefig('./result/case3_pic/' + num + '.png')
    # plt.show()
    plt.close()


def get(num, f):

    N1 = 0b000000000011100
    N2 = 0b000000000100000
    N3 = 0b111111111111111
    N4 = 0
    for i in range(0, num):
        N4 = N4 + (1 << i)

    m = [0 for x in range(0, 1 << 15)]
    Bel = [0 for x in range(0, 1 << 15)]
    Pl = [0 for x in range(0, 1 << 15)]
    Q = [0 for x in range(0, 1 << 15)]

    m[N1] = 0.05
    m[N2] = 0.05
    m[N3] = 0.1
    m[N4] = 0.8

    for iter in range(0, 1 << 15):
        number, fre, Type, value = f.readline().strip('\n').split(' ')
        Bel[eval(fre)] = eval(value)
    for iter in range(0, 1 << 15):
        number, fre, Type, value = f.readline().strip('\n').split(' ')
        Pl[eval(fre)] = eval(value)
    for iter in range(0, 1 << 15):
        number, fre, Type, value = f.readline().strip('\n').split(' ')
        Q[eval(fre)] = eval(value)

    return m, Bel, Pl, Q


def plot_case1():
    A = [0 for x in range(1, 15)]
    Uncertainty1 = [0 for x in range(1, 15)]
    Uncertainty2 = [0 for x in range(1, 15)]
    Uncertainty3 = [0 for x in range(1, 15)]
    Uncertainty4 = [0 for x in range(1, 15)]
    Uncertainty5 = [0 for x in range(1, 15)]
    Uncertainty6 = [0 for x in range(1, 15)]
    Uncertainty7 = [0 for x in range(1, 15)]
    Uncertainty8 = [0 for x in range(1, 15)]
    Uncertainty9 = [0 for x in range(1, 15)]
    Uncertainty10 = [0 for x in range(1, 15)]
    Uncertainty11 = [0 for x in range(1, 15)]
    Uncertainty12 = [0 for x in range(1, 15)]
    Uncertainty13 = [0 for x in range(1, 15)]
    Uncertainty14 = [0 for x in range(1, 15)]
    Uncertainty15 = [0 for x in range(1, 15)]
    Uncertainty16 = [0 for x in range(1, 15)]
    Deng = [0 for x in range(1, 15)]
    f = open("./src/case1.txt", "r", encoding='utf-8')
    for iter in range(1, 15):
        print('case1', 'Size of A =', iter)
        A[iter-1] = iter
        m, Bel, Pl, Q = get(iter, f)
        Uncertainty1[iter-1] = definition1(m, Bel, Pl, Q, 15)
        Uncertainty2[iter-1] = definition2(m, Bel, Pl, Q, 15)
        Uncertainty3[iter-1] = definition3(m, Bel, Pl, Q, 15)
        Uncertainty4[iter-1] = definition4(m, Bel, Pl, Q, 15)
        Uncertainty5[iter-1] = Uncertainty3[iter-1] + Uncertainty4[iter-1]
        Uncertainty6[iter-1] = definition6(m, Bel, Pl, Q, 15)
        Uncertainty7[iter-1] = definition7(m, Bel, Pl, Q, 15)
        Uncertainty8[iter-1] = definition8(m, Bel, Pl, Q, 15)
        Uncertainty9[iter-1] = definition9(m, Bel, Pl, Q, 15)
        Uncertainty10[iter-1] = definition10(m, Bel, Pl, Q, 15)
        Uncertainty11[iter-1] = definition11(m, Bel, Pl, Q, 15) + Uncertainty4[iter-1]
        Uncertainty12[iter-1] = definition12(m, Bel, Pl, Q, 15) + Uncertainty4[iter-1]
        Uncertainty13[iter-1] = definition13(m, Bel, Pl, Q, 15)
        Uncertainty14[iter-1] = definition14(m, Bel, Pl, Q, 15)
        Uncertainty15[iter-1] = definition15(m, Bel, Pl, Q, 15)
        Uncertainty16[iter-1] = definition16(m, Bel, Pl, Q, 15)
        Deng[iter-1] = Deng_f(m, Bel, Pl, Q, 15)
    f.close()
    plt.figure(figsize=(12, 8), dpi=72)
    plt.ylim(0, 15)
    plt.xlim(0, 15)
    plt.xlabel("A")
    plt.ylabel("Uncertainty")
    plt.tick_params(axis='both')
    plt.plot(A, Uncertainty1, color='r', label='definition1')
    plt.plot(A, Uncertainty2, color='g', label='definition2')
    plt.plot(A, Uncertainty3, color='b', label='definition3')
    plt.plot(A, Uncertainty4, color='y', label='definition4')
    plt.plot(A, Uncertainty5, color='k', label='definition5')
    plt.plot(A, Uncertainty6, '-.', color='r', label='definition6')
    plt.plot(A, Uncertainty7, '-.', color='g', label='definition7')
    plt.plot(A, Uncertainty8, '-.', color='b', label='definition8')
    plt.plot(A, Uncertainty9, '-.', color='y', label='definition9')
    plt.plot(A, Uncertainty10, '-.', color='k', label='definition10')
    plt.plot(A, Uncertainty11, '--', color='r', label='definition11')
    plt.plot(A, Uncertainty12, '--', color='g', label='definition12')
    plt.plot(A, Uncertainty13, '--', color='b', label='definition13')
    plt.plot(A, Uncertainty14, '--', color='y', label='definition14')
    plt.plot(A, Uncertainty15, '--', color='k', label='definition15')
    plt.plot(A, Uncertainty16, ':', color='g', label='definition16')
    plt.plot(A, Deng, ':', color='r', label='Deng')
    plt.legend()
    try:
        plt.savefig('./result/case1.png')
    except:
        os.system('mkdir result')
        plt.savefig('./result/case1.png')
    # plt.show()
    plt.close()


def case1():
    plot_case1()
    # create_data1()


def case2():
    A = [0 for x in range(1, 21)]
    Uncertainty1 = [0 for x in range(1, 21)]
    Uncertainty2 = [0 for x in range(1, 21)]
    Uncertainty3 = [0 for x in range(1, 21)]
    Uncertainty4 = [0 for x in range(1, 21)]
    Uncertainty5 = [0 for x in range(1, 21)]
    Uncertainty6 = [0 for x in range(1, 21)]
    Uncertainty7 = [0 for x in range(1, 21)]
    Uncertainty8 = [0 for x in range(1, 21)]
    Uncertainty9 = [0 for x in range(1, 21)]
    Uncertainty10 = [0 for x in range(1, 21)]
    Uncertainty11 = [0 for x in range(1, 21)]
    Uncertainty12 = [0 for x in range(1, 21)]
    Uncertainty13 = [0 for x in range(1, 21)]
    Uncertainty14 = [0 for x in range(1, 21)]
    Uncertainty15 = [0 for x in range(1, 21)]
    Uncertainty16 = [0 for x in range(1, 21)]
    Deng = [0 for x in range(1, 21)]
    for iter in range(1, 21):
        print('case2', 'N =', iter)
        A[iter-1] = iter
        m, Bel, Pl, Q = create_data2(iter)
        # Uncertainty1[iter-1] = definition1(m, Bel, Pl, Q, iter)
        # Uncertainty2[iter-1] = definition2(m, Bel, Pl, Q, iter)
        # Uncertainty3[iter-1] = definition3(m, Bel, Pl, Q, iter)
        Uncertainty4[iter-1] = definition4(m, Bel, Pl, Q, iter)
        # Uncertainty5[iter-1] = Uncertainty3[iter-1] + Uncertainty4[iter-1]
        # Uncertainty6[iter-1] = definition6(m, Bel, Pl, Q, iter)
        # Uncertainty7[iter-1] = definition7(m, Bel, Pl, Q, iter)
        # Uncertainty8[iter-1] = definition8(m, Bel, Pl, Q, iter)
        Uncertainty9[iter-1] = definition9(m, Bel, Pl, Q, iter)
        # Uncertainty10[iter-1] = definition10(m, Bel, Pl, Q, iter)
        # Uncertainty11[iter-1] = definition11(m, Bel, Pl, Q, iter) + Uncertainty4[iter-1]
        # Uncertainty12[iter-1] = definition12(m, Bel, Pl, Q, iter) + Uncertainty4[iter-1]
        Uncertainty13[iter-1] = definition13(m, Bel, Pl, Q, iter)
        # Uncertainty14[iter-1] = definition14(m, Bel, Pl, Q, iter)
        Uncertainty15[iter-1] = definition15(m, Bel, Pl, Q, iter)
        Uncertainty16[iter-1] = definition16(m, Bel, Pl, Q, iter)
        Deng[iter-1] = Deng_f(m, Bel, Pl, Q, iter)
    plt.figure(figsize=(12, 8), dpi=72)
    plt.ylim(0, 5)
    plt.xlim(0, 20)
    plt.xlabel("A")
    plt.ylabel("Uncertainty")
    plt.tick_params(axis='both')
    plt.plot(A, Uncertainty1, color='r', label='definition1')
    plt.plot(A, Uncertainty2, color='g', label='definition2')
    plt.plot(A, Uncertainty3, color='b', label='definition3')
    plt.plot(A, Uncertainty4, color='y', label='definition4')
    plt.plot(A, Uncertainty5, color='k', label='definition5')
    plt.plot(A, Uncertainty6, '-.', color='r', label='definition6')
    plt.plot(A, Uncertainty7, '-.', color='g', label='definition7')
    plt.plot(A, Uncertainty8, '-.', color='b', label='definition8')
    plt.plot(A, Uncertainty9, '-.', color='y', label='definition9')
    plt.plot(A, Uncertainty10, '-.', color='k', label='definition10')
    plt.plot(A, Uncertainty11, '--', color='r', label='definition11')
    plt.plot(A, Uncertainty12, '--', color='g', label='definition12')
    plt.plot(A, Uncertainty13, '--', color='b', label='definition13')
    plt.plot(A, Uncertainty14, '--', color='y', label='definition14')
    plt.plot(A, Uncertainty15, '--', color='k', label='definition15')
    plt.plot(A, Uncertainty16, ':', color='g', label='definition16')
    plt.plot(A, Deng, ':', color='r', label='Deng')
    plt.legend()
    try:
        plt.savefig('./result/case2.png')
    except:
        os.system('mkdir result')
        plt.savefig('./result/case2.png')
    # plt.show()
    plt.close()


def case3():
    plot3(definition1, 'definition1')
    plot3(definition2, 'definition2')
    plot3(definition3, 'definition3')
    plot3(definition4, 'definition4')
    plot3(definition5, 'definition5')
    plot3(definition6, 'definition6')
    plot3(definition7, 'definition7')
    plot3(definition8, 'definition8')
    plot3(definition9, 'definition9')
    plot3(definition10, 'definition10')
    plot3(definition11, 'definition11')
    plot3(definition12, 'definition12')
    plot3(definition13, 'definition13')
    plot3(definition14, 'definition14')
    plot3(definition15, 'definition15')
    plot3(definition16, 'definition16')
    plot3(Deng_f, 'Deng')


if __name__ == '__main__':
    for i in range(1, len(sys.argv)):
        op = sys.argv[i]
        if op == "-case1":
            case1()
        elif op == "-case2":
            case2()
        elif op == "-case3":
            case3()
        elif op == "-h":
            print('\n[Options]')
            print('\n[-case1] --------------------------------------------------------- case1')
            print('[-case2] --------------------------------------------------------- case2')
            print('[-case3] --------------------------------------------------------- case3')
            print('[-h]     --------------------------------------------------------- help\n')
            sys.exit(0)
    print('FINISHED !')
