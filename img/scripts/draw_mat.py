#!/usr/bin/env python3
import matplotlib.pyplot as plt
import numpy as np
from scipy import sparse
import sys


def check_symmetry(mtx):
    a, b = mtx.shape
    if a != b:
        print("NO!")
    else:
        for i in range(a):
            for j in range(i):
                if mtx[i][j] != mtx[j][i]:
                    print("NO")
                    return
        print("YES")


if __name__ == "__main__":
    with open("input/"+sys.argv[1]+".txt", 'r') as f:
        n = int(f.readline())
        mat = np.zeros((n, n))
        for i in range(n):
            for j in range(n):
                mat[i][j] = float(f.readline())
    mat = sparse.coo_matrix(mat)
    x = mat.col.tolist()
    y = mat.row.tolist()
    T = mat.data.tolist()
    plt.scatter(x, y, c=T)
#plt.savefig('img/Ac_mat.png')
    plt.show()
