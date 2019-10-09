import numpy as np
import matplotlib.pyplot as plt


def func_polynomial(v):
    return 0.567307 - 1.73076 * v ** 2 - 1e-6 * v ** 3 + 1.201921 * v ** 4 + 1e-6 * v ** 5


def func_original(v):
    return 1.0 / (1 + 25 * v ** 2)


def func_spline(v):
    a = [1.720641, -3.314769, -0.000000, 3.314769, -1.720641, 0.000000]
    b = [0.000000, 2.064770, -1.912954, -1.912954, 2.064770, 0.000000]
    c = [-0.121453, 0.704455, 0.765181, -0.765182, -0.704455, 0.121453]
    left, right = 0, 6
    while left < right:
        mid = (left + right) // 2
        if x[mid] < v:
            left = mid + 1
        elif x[mid] > v:
            right = mid - 1
        else:
            return y[mid]
    left = min(left, 5)
    left = left if x[left] < v else left - 1
    left = max(0, left)
    h = v - x[left]
    if left < 6:
        return ((a[left] * h + b[left]) * h + c[left]) * h + y[left]
    else:
        return b[-1] * h + c[-1] * h + y[-1]


def plot(f, x, label):
    funcx = np.arange(x[0], x[-1], 0.01)
    funcy = []
    for i in funcx:
        funcy.append(f(i))
    plt.plot(funcx, funcy, label=label)


def do_plot():
    plt.xlabel = 'x'
    plt.ylabel = 'y'
    plt.legend(loc='best')
    plt.grid(linestyle='-.')
    plt.savefig('img/cubic_spline.png')
    plt.show()


if __name__ == '__main__':
    x = np.array([-1.0, -0.6, -0.2, 0.2, 0.6, 1.0])
    y = func_original(x)
    plt.scatter(x, y)
    plot(func_polynomial, x, 'polynomial')
    plot(func_original, x, 'origin')
    plot(func_spline, x, 'cubic_spline')
    do_plot()
