import numpy as np
import matplotlib.pyplot as plt

# заполняем входной файл
def write(filename, x, p):
    f = open(filename, "w")
    f.write(str(len(x)) + "\n")
    for i in x:
        f.write(str(i) + "\n")
    f.write(str(p[0]))
    for i in p[1:]:
        f.write("\n" + str(i))
    f.close()


# создание матрицы информационной матрицы М
def makeM(x, p, delta):
    M = np.zeros((len(func(x[0], delta)), len(func(x[0], delta))))
    for i in range(len(x)):
        M += p[i] * make_partM(func(x[i], delta))
    return M


# создание части информационной матрицы матрицы, для каждого элемента спектра плана
def make_partM(fx):
    M = np.zeros((len(fx), len(fx)))
    for i in range(len(fx)):
        for j in range(len(fx)):
            M[i][j] = fx[i] * fx[j]
    return M


# возвращает вектор со значениями нашей функции
def func(x, delta):
    return np.array([1, x, mu1(x, delta), x * mu1(x, delta)])


def mu1(dot, delta):
    if -delta >= dot >= -1:
        return 1.
    elif delta >= dot >= -delta:
        return (delta - dot) / 2 * delta
    elif dot > delta:
        return 0


def mu2(dot, delta):
    if -delta >= dot >= -1:
        return 0.
    elif delta >= dot >= -delta:
        return (delta + dot) / 2 * delta
    elif dot > delta:
        return 1.


# создание начального плана
def makeUPlan(N):
    return list(np.linspace(-1, 1, N)), [1.0 / N for _ in range(N)]


# создание матрицы D из матрицы М
def makeD(M):
    return np.linalg.inv(M)


# Добавление новой точки
def addNewPoint(x, p, grid, delta, flag, usereps=0.01):
    M = makeM(x, p, delta)
    D = makeD(M)
    curA = A_optim(D)
    max, maxdot = findMaxFi(grid, D, delta)
    if flag:
        x.append(maxdot)
        p = addToP(curA, p, x, delta)
    eps = max * usereps
    delta = abs(max - np.trace(D @ D @ M))
    print(np.trace(D))
    return delta, eps, x, p


# Проверка на D-оптимальность ( вычисление критерия для заданной матрицы)
def A_optim(D):
    return -np.trace(D)


# функция нахождения максимума на сетке, для добавления новой точки
def findMaxFi(grid, D, delta):
    max = func(grid[0], delta) @ D @ D @ func(grid[0], delta).T
    maxdot = grid[0]
    for i in grid:
        f = func(i, delta) @ D @ D @ func(i, delta).T
        if f > max:
            max = f
            maxdot = i
    return max, maxdot


# Подбор веса для новой точки
def addToP(curA, p, x, delta):
    newD = curA - 1
    ksy = 1
    while curA > newD:
        newP = p.copy()
        for i in range(len(newP)):
            newP[i] = (1.0 - ksy / len(newP)) * newP[i]
        newP.append(ksy / len(newP))
        newM = makeM(x, newP, delta)
        newD = makeD(newM)
        newD = A_optim(newD)
        ksy /= 2
    return newP


# объединение близких точек
def unionCloseDots(x, p):
    newX = [x[0]]
    newP = [p[0]]
    for i in range(1, len(x)):
        index = findClose(x[i], newX)
        if index == -1:
            newX.append(x[i])
            newP.append(p[i])
        else:
            newP[index] += p[i]
    return newX, newP


# нахождение в массиве X индекса элемента, который близок к элементу x
def findClose(x, X):
    for i in range(len(X)):
        vec = np.sqrt((x - X[i]) ** 2)
        if vec < 0.025:
            return i
    return -1


# удаление точек с малыми весами
def removeDotsWithSmallWeitgh(x, p):
    sum = 0
    index = 0
    for i in range(len(p)):
        if p[i] < 0.01:
            sum += p[i]
            p[i] = 0
            x[i] = 0
            index += 1
    for i in range(index):
        p.remove(0)
        x.remove(0)
    sum /= len(p)
    for i in range(len(p)):
        p[i] += sum
    return x, p


def Plot(x, grid, D, label, delta):
    Dfunc = list(map(lambda x: func(x, delta) @ D @ D @ func(x, delta), grid))
    maxD = max(Dfunc)
    for i in x:
        plt.scatter(i, 0)
    for i in x:
        plt.plot([i, i], [0, maxD], color='grey')
    plt.plot([x[0], x[0]], [0, maxD], color='grey', label='План')
    plt.plot(grid, Dfunc, color='red', label='d(x,e)')
    F = [mu1, mu2]
    Color = ['black', 'green']
    Type = ['--', '.-']
    for f, t, c in zip(F, Type, Color):
        pltmu = list(map(lambda x: f(x, delta), grid))
        plt.plot(grid, pltmu, t, color=c, label="{i}".format(i=f.__name__))
    plt.legend(loc='best')
    plt.title(label)
    plt.xlabel('x')
    plt.ylabel('f(x)')
    plt.show()


def main():
    Delta = 0.5
    N = 25
    delta = 1
    eps = 0.01
    iteration = 1
    sumiteration = 0
    ito = 0

    x, p = makeUPlan(N)
    write("input.txt", x, p)
    D = makeD(makeM(x, p, Delta))
    grid = np.linspace(-1, 1, 101)
    Plot(x, grid, D, "График начального плана", Delta)

    while abs(delta) > eps:
        smalliteration = 0
        while delta > eps:
            delta, eps, x, p = addNewPoint(x, p, grid, Delta, True)
            smalliteration += 1
            if smalliteration > 10000:
                break
        x, p = unionCloseDots(x, p)

        x, p = removeDotsWithSmallWeitgh(x, p)
        ito += 1
        sumiteration += smalliteration
        delta, eps, x, p = addNewPoint(x, p, grid, Delta, False)

        iteration += 1
        print(f'Номер итерации {iteration}')
        if iteration > 1000:
            break

    write("output.txt ", x, p)
    D = makeD(makeM(x, p, Delta))
    Plot(x, grid, D, 'График оптимального плана', Delta)


if __name__ == '__main__':
    main()
