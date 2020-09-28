import matplotlib.pyplot as plt
import numpy as np
from scipy import sparse


def load_matrix(path):
    with open(path, 'r') as f:
        l = f.read().split()
        n = int(l[0])
        l = l[1:]
        mat = np.array([float(i) for i in l[:n*n]]).reshape(n, n)
        l = l[n*n:]
        vec = np.array([float(i) for i in l]).reshape(n, 1)
    return n, mat, vec


def load_sparse_matrix(path):
    with open(path, 'r') as f:
        l = f.read().split()
        n = int(l[0])
        l = l[1:]
        mat = sparse.coo_matrix(np.array([float(i) for i in l[:n*n]]).reshape(n, n), shape=(n, n))
        l = l[n*n:]
        vec = sparse.coo_matrix(np.array([float(i) for i in l]).reshape(n, 1), shape=(n, 1))
    return n, mat, vec


def draw_mat(n, mat: sparse.coo_matrix):
    row = mat.row.tolist()
    col = mat.col.tolist()
    fig = plt.figure()
    plt.scatter(col, row)
    plt.show()


if __name__ == '__main__':
    n, mat, vec = load_matrix('input/Ac.txt')
    plt.imshow(mat)
    plt.show()
