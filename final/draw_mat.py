#!/usr/bin/env python3
import matplotlib.pyplot as plt
import numpy as np
from scipy import sparse, io


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


def load_mtx(filepath: str, toType: str):
    try:
        mtx = io.mmread(filepath)
        if toType == 'bsr':
            return sparse.bsr_matrix(mtx)
        elif toType == 'coo':
            return sparse.coo_matrix(mtx)
        elif toType == 'csr':
            return sparse.csr_matrix(mtx)
        elif toType == 'csc':
            return sparse.csc_matrix(mtx)
        elif toType == 'array':
            return np.array(mtx)
        elif toType == 'dok':
            return sparse.dok_matrix(mtx)
        elif toType == 'dia':
            return sparse.dia_matrix(mtx)
        else:
            print("Unknown type:%s" % toType)
    except IOError:
        print('No such file: %s' % filepath)
    except ValueError:
        print('error:%s' % filepath)
    return None


if __name__ == "__main__":
    mat = load_mtx('sample.mtx', 'coo')
    x = mat.col.tolist()
    y = mat.row.tolist()
    T = mat.data.tolist()
    plt.figure(figsize=(4, 6))
    plt.scatter(x, y, c=T)
    plt.show()
