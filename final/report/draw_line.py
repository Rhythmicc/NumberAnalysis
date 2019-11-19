import matplotlib.pyplot as plt
import numpy as np
import sys


def draw_line(xx, yy, label):
    plt.scatter(xx, yy)
    plt.plot(xx, yy, label=label)


if __name__ == '__main__':
    data = []
    with open(sys.argv[1], 'r') as f:
        for line in f.readlines():
            if line[0] == 'U' or line[0] == 'V' or line.startswith('Total'):
                line = set(line.split())
                for x in line:
                    try:
                        x = float(x)
                        data.append(x)
                    except:
                        continue
    x = np.arange(1, 21)
    y_uX = []
    y_uY = []
    y_v = []
    y_T = []
    for i in range(0, len(data), 4):
        y_uX.append(data[i])
        y_uY.append(data[i + 1])
        y_v.append(data[i + 2])
        y_T.append(data[i + 3])
    print(y_uX)
    print(y_uY)
    print(y_v)
    print(y_T)
    title = sys.argv[1].split('.')[0]
    plt.title(title)
    draw_line(x, y_uX, 'update X')
    draw_line(x, y_uY, 'update Y')
    draw_line(x, y_v, 'Validate')
    draw_line(x, y_T, 'Total')
    plt.xlabel('times')
    plt.ylabel('ms')
    plt.legend(loc='best')
    plt.grid(linestyle='-.')
    plt.savefig(title+'.png')
    plt.show()
