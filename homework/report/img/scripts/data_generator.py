#!/usr/bin/env python3
from scipy import sparse
import numpy as np
import tqdm


def random_vector(n):
    p = np.random.random() * 0.1 + 0.9
    rand_coo = sparse.rand(1, n, density=p, format='coo')
    rand_coo *= np.random.randint(100, 1000)
    rand_coo = rand_coo.toarray()[0]
    rand_coo = [round(i, 3) for i in rand_coo]
    return rand_coo


def random_matrix(n):
    p = np.random.random() * 0.1 + 0.9
    _mtx = sparse.rand(n, n, density=p, format='coo') * np.random.randint(100, 1000)
    _mtx = _mtx.toarray().tolist()
    for i in range(n):
        for j in range(n):
            _mtx[i][j] = round(_mtx[i][j], 3)
    return _mtx


def make_symmetric(_mtx):
    n = len(_mtx)
    for i in range(n):
        for j in range(i):
            _mtx[i][j] = _mtx[j][i]
    return _mtx


if __name__ == '__main__':
    with open('./cmake-build-debug/input.txt', 'w') as f:
        f.write('4\n')
        bar = tqdm.tqdm(total=15)
        for i in range(4):
            sz = (1 << i) * 1000
            f.write(str(sz) + '\n')
            mtx = make_symmetric(random_matrix(sz))
            for line in mtx:
                f.write(' '.join([str(i) for i in line])+'\n')
            vec = random_vector(sz)
            f.write(' '.join([str(i) for i in vec])+'\n')
            bar.update(1 << i)
