import numpy as np
import matplotlib.pyplot as plt
import matplotlib

matplotlib.use("TKAgg")

# Линейная однофакторная модель





def model(x, delta):
    return [[1], [x], [mu1(x, delta)], [mu1(x, delta) * x]]

# Вычисление информационных матриц для каждого плана


def inf_matrix(x, p, delta):
    M = np.zeros((6, 6))
    for i in range(len(x)):
        M += p[i] * np.dot(model(x[i], delta),
                           np.transpose(model(x[i], delta)))
    return M


# Вычисление дисперсионной матрицы

def disp_matrix(M):
    return np.linalg.inv(M)

# Создание плана


def make_plan():
    x, p = [], []
    grid = np.arange(-1, 1.01, 0.25)
    
    for i in range(len(grid)):
        x.append([grid[i]])
        p.append(1/n)

    return x, p

# Проверка оптимальности


def A_optimal(x_plan, p_plan, delta):
    n = 81
    max = 0
    grid = np.arange(-1, 1.01, 0.25)
    x = []

    for i in range(len(grid)):
        for j in range(len(grid)):
            x.append([grid[j], grid[i]])

    for i in range(n):
        condition = np.trace(np.dot(np.linalg.matrix_power(disp_matrix(inf_matrix(
            x_plan, p_plan)), 2), np.dot(model(x[i], delta), np.transpose(model(x[i], delta)))))
        if condition > max:
            max = condition

    trd = np.trace(disp_matrix(inf_matrix(x_plan, p_plan)))

    return max, trd


def plot(x_plan, p_plan):
    fig, ax = plt.subplots()
    cmap = plt.cm.coolwarm
    for i in range(len(x_plan)):
        ax.scatter(x_plan[i][0], x_plan[i][1], s=p_plan[i]*1000, cmap=cmap)
    plt.grid()
    plt.show()


def main():
    flag = True                     # Флаг выхода из градиентного спуска
    n = 81                          # Кол-во точек на сетке
    grad = np.zeros(n)              # Вектор градиента
    count = 0                       # Cчётчик итераций
    delta = 0.5
    proj = np.zeros(n)

    x_plan, p_plan = make_plan()    # Начальный план

    while flag:
        flag = False
        M = inf_matrix(x_plan, p_plan)
        D = disp_matrix(M)

        lambdas = 0.01  # Шаг по весам
        q = n - np.count_nonzero(p_plan)

        for i in range(81):     # Вычисление градиента
            grad[i] = np.trace(np.dot(np.linalg.matrix_power(disp_matrix(inf_matrix(
                x_plan, p_plan, delta)), 2), np.dot(model(x_plan[i], delta), np.transpose(model(x_plan[i], delta)))))

        grad /= np.linalg.norm(grad)

        avg = 0.0
        for i in range(n):
            if p_plan[i] != 0:
                avg += grad[i]
        avg /= n - q

        for i in range(n):
            if p_plan[i] != 0:
                if abs(grad[i] - avg) > 1e-10:
                    flag = True

        for j in range(0, n):
            proj[j] = grad[j] - avg
            if p_plan[j] == 0:
                if proj[j] > 0:
                    flag = True
                else:
                    proj[j] = 0

        if count % 200 == 0:
            print(np.trace(D))
            # print(np.linalg.norm(proj))

        if flag:
            for i in range(0, n):
                if proj[i] < 0 and lambdas > - p_plan[i] / proj[i]:
                    lambdas = - p_plan[i] / proj[i]

            for i in range(0, n):
                p_plan[i] += lambdas * proj[i]

        count += 1

    n_plan, np_plan = [], []
    for i in range(len(p_plan)):
        if p_plan[i] != 0:
            n_plan.append(x_plan[i])
            np_plan.append((p_plan[i]))

    print(count)
    a, b = A_optimal(n_plan, np_plan, delta)
    print(np.trace(disp_matrix(inf_matrix(n_plan, np_plan, delta))))


if __name__ == '__main__':
    main()
