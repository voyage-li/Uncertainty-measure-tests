def create_data(num):
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
        Bel[i] = ((Bel[i]*100)//1)/100
        print(Bel[i])

    for i in range(0, 1 << 15):
        print(bin(N4), bin(i), 'Pl', end=' ')
        for j in range(0, 1 << 15):
            if m[j] != 0 and i & j != 0:
                Pl[i] += m[j]
        Pl[i] = ((Pl[i]*100)//1)/100
        print(Pl[i])

    for i in range(0, 1 << 15):
        print(bin(N4), bin(i), 'Q', end=' ')
        for j in range(i, 1 << 15):
            if i | j == j:
                Q[i] += m[j]
        Q[i] = ((Q[i]*100)//1)/100
        print(Q[i])

    return m, Bel, Pl, Q


if __name__ == '__main__':
    for iter in range(1, 15):
        m, Bel, Pl, Q = create_data(iter)
