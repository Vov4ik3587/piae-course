from tkinter import *
from tkinter.ttk import *
import numpy as np
import matplotlib.pyplot as plt

plt.rcParams['font.family'] = 'fantasy'
plt.rcParams['font.fantasy'] = 'Arial'


def inserterErorr(value):
    output.delete("1.0", "end")
    output.insert("1.0", value)


def inserter(value):
    output.insert("1.0", value)


def clear(event):
    """ Clears entry form """
    caller = event.widget
    caller.delete("0", "end")


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
def makeM(x, p, linear):
    M = np.zeros((len(func(x[0], linear)), len(func(x[0], linear))))
    for i in range(len(x)):
        M += p[i] * make_partM(func(x[i], linear))
    return M


# создание части информационной матрицы матрицы, для каждого элемента спектра плана
def make_partM(fx):
    M = np.zeros((len(fx), len(fx)))
    for i in range(len(fx)):
        for j in range(len(fx)):
            M[i][j] = fx[i] * fx[j]
    return M


# возвращает вектор со значениями нашей функции
def func(x, linear):
    if linear:
        return np.array([mu1(x), x * mu1(x), mu2(x), x * mu2(x)])
    else:
        return np.array([mu1(x), x * mu1(x), x ** 2 * mu1(x), mu2(x), x * mu2(x), x ** 2 * mu2(x)])


def mu1(dot):
    if -0.5 > dot >= -1.0:
        return 1.0
    elif 0.5 >= dot >= -0.5:
        return 1 / 2 - dot
    elif 1.0 >= dot > 0.5:
        return 0


def mu2(dot):
    if -0.5 > dot >= -1.0:
        return 0
    elif 0.5 >= dot >= -0.5:
        return 1 / 2 + dot
    elif 1.0 >= dot > 0.5:
        return 1


# создание начального плана
def makeUPlan(N):
    return list(np.linspace(-1, 1, N)), [1.0 / N for _ in range(N)]


# создание матрицы D из матрицы М
def makeD(M):
    return np.linalg.inv(M)


# Добавление новой точки
def addNewPoint(x, p, grid, linear, flag, usereps):
    M = makeM(x, p, linear)
    D = makeD(M)
    curientD = D_optim(M)
    max, maxdot = findMaxFi(grid, D, linear)
    if flag:
        x.append(maxdot)
        p = addToP(curientD, p, x, linear)
    eps = max * usereps
    delta = abs(max - np.trace(np.dot(M, D)))
    # np.around((np.linalg.det(M)), 9))
    inserter("%04.14e \n" % np.linalg.det(M))
    print(np.linalg.det(M))
    return delta, eps, x, p


# Проверка на D-оптимальность ( вычисление критерия для заданной матрицы)
def D_optim(M):
    return np.linalg.det(M)


# функция нахождения максимума на сетке, для добавления новой точки
def findMaxFi(grid, D, linear):
    max = np.dot(np.dot(func(grid[0], linear), D), func(grid[0], linear).T)
    maxdot = grid[0]
    for i in grid:
        f = np.dot(np.dot(func(i, linear), D), func(i, linear).T)
        if f > max:
            max = f
            maxdot = i
    return max, maxdot


# Подбор веса для новой точки
def addToP(curientD, p, x, linear):
    newD = curientD - 1
    ksy = 1
    while curientD > newD:
        newP = p.copy()
        for i in range(len(newP)):
            newP[i] = (1.0 - ksy / len(newP)) * newP[i]
        newP.append(ksy / len(newP))
        newM = makeM(x, newP, linear)
        newD = D_optim(newM)
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


def Plot(x, grid, D, label, linear):
    Dfunc = list(map(lambda x: np.dot(
        np.dot(func(x, linear), D), func(x, linear)), grid))
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
        pltmu = list(map(lambda x: f(x), grid))
        plt.plot(grid, pltmu, t, color=c, label="{i}".format(i=f.__name__))
    plt.legend(loc='best')
    plt.title(label, family="verdana")
    plt.xlabel('x')
    plt.ylabel('f(x)')
    plt.show()


def main():
    output.delete("1.0", "end")
    try:
        e_val = int(e.get())
    except ValueError:
        inserterErorr("Количество элементов должно быть целым!")
    try:
        r_val = float(r.get())
    except ValueError:
        inserterErorr("Формат ввода точности: 0.1 \n Повторите ввод!")
    if combobox.current() == 0:
        linear = True
    else:
        linear = False

    x, p = makeUPlan(e_val)

    write("input.txt", x, p)
    D = makeD(makeM(x, p, linear))
    grid = np.linspace(-1, 1, 101)
    Plot(x, grid, D, "График начального плана", linear)
    delta = 1
    eps = 0.001
    iteration = 1
    sumiteration = 0
    ito = 0
    while abs(delta) > eps:
        smalliteration = 0
        while delta > eps:
            delta, eps, x, p = addNewPoint(x, p, grid, linear, True, r_val)
            smalliteration += 1
            if smalliteration > 10000:
                break
        x, p = unionCloseDots(x, p)

        x, p = removeDotsWithSmallWeitgh(x, p)
        ito += 1
        inserter("\nКоличество добавленных точек: %s \n\n" % smalliteration)
        sumiteration += smalliteration
        inserter("\nОчистка № %s \n" % ito)
        delta, eps, x, p = addNewPoint(x, p, grid, linear, False, r_val)

        iteration += 1
        print(iteration)
        if iteration > 1000:
            break

    write("output.txt ", x, p)
    D = makeD(makeM(x, p, linear))
    Plot(x, grid, D, "График оптимального плана", linear)
    inserter("Оптимальный план построен за %s итераций \n" % sumiteration)
    allw = output.get('1.0', END)
    output.delete("1.0", "end")
    inserter("%s " % revers(allw))


def revers(s):
    s = s.split("\n")
    s.reverse()
    a = "Результаты: \n\nОпределитель матрицы M:"
    for i in s:
        a += i + '\n'
    return a


root = Tk()
# установить расположения окна по центру
x = (root.winfo_screenwidth() - root.winfo_reqwidth()) / 8
y = (root.winfo_screenheight() - root.winfo_reqheight()) / 8
root.wm_geometry("+%d+%d" % (x, y))
root.title("Синтез непрерывных D-планов эксперимента для нечетких однофакторных моделей с двумя подобластями определения")
root.minsize(900, 600)
root.resizable(width=False, height=False)

frame = Frame(root)
frame.grid()

# выбор модели
lab_mod = Label(frame, text="Выберите модель:").grid(
    column=1, row=0, pady=(10, 0), padx=(30, 0))
combobox = Combobox(frame, values=[u"Линейная", u"Квадратичная"], height=2)
combobox.set(u"Линейная")
combobox.grid(column=2, row=0, pady=(10, 0), padx=(5, 0))

# точность
lab_t = Label(frame, text="Введите точность:").grid(
    column=3, row=0, pady=(10, 0), padx=(30, 0))
vt = StringVar()
vt.set('0.01')
r = Entry(frame, width=23, textvariable=vt)
r.grid(row=0, column=4, padx=(5, 0), pady=(10, 0))
r.bind("<FocusIn>", clear)

# количество
lab_n = Label(frame, text="Введите количество элементов:").grid(
    column=5, row=0, padx=(30, 0), pady=(10, 0))
v = StringVar()
v.set('25')
e = Entry(frame, width=23, textvariable=v)
e.grid(row=0, column=6, padx=(5, 0), pady=(10, 0))
e.bind("<FocusIn>", clear)

but = Button(frame, text="Решить", command=main)
but.grid(row=0, column=7, padx=(30, 0), pady=(10, 0))

output = Text(frame, bg="white", font="Arial 12", width=100, height=30)
output.grid(row=1, column=1, columnspan=7, rowspan=4, pady=(10, 0))

scrollbar = Scrollbar(command=output.yview, orient=VERTICAL)
scrollbar.grid(row=0, column=1, sticky='ns', padx=(30, 0))
output.configure(yscrollcommand=scrollbar.set)

root.mainloop()
