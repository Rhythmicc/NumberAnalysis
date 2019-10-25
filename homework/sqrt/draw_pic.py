#!/usr/bin/env python3
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import axes3d


def main():
    x, y = np.meshgrid(np.arange(-5, 5, 1e-2), np.arange(-5, 5, 1e-2))
    z = np.sqrt(x ** 2 + y ** 2)
    ax = plt.gca(projection='3d')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    ax.plot_surface(x, y, np.sin(z), cmap='rainbow')
    ax.grid(linestyle='-.')
    plt.title('$sin(\sqrt{x^2+y^2})$')
    plt.savefig('sqrt_homework.png')
    plt.show()


if __name__ == '__main__':
    main()
