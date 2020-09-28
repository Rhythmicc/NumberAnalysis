import matplotlib.pyplot as plt
import numpy as np

config = {
    'Matrix': 'Af',
    'row_size': 1220,
    'T': 20
}

if __name__ == "__main__":
    with open('%s_time' % config['Matrix'], 'r') as f:
        data = [float(i) for i in f.read().split()]
    doo = data[:config['T']]
    for i in range(config['T']):
        doo[i] = 2 * config['row_size'] ** 2 / doo[i] / 1e6
    x = np.arange(1, config['T'] + 1, 1)
    aver_doo = sum(doo) / config['T']

    max_doo = max(doo)
    max_doo_index = doo.index(max_doo) + 1
    min_doo = min(doo)
    min_doo_index = doo.index(min_doo) + 1

    plt.title('%s.txt aver_conjugate: %f' % (config['Matrix'], aver_doo))
    plt.plot(x, doo, label='conjugate')
    plt.scatter(x, doo, alpha=0.5, edgecolors='white')
    plt.scatter([max_doo_index], [max_doo], label='max_cg')
    plt.scatter([min_doo_index], [min_doo], label='min_cg')
    plt.xlabel('Times')
    plt.ylabel('MFlop/s')
    plt.legend(loc='best')
    plt.grid(linestyle='-.')
    plt.savefig('./img/%s.png' % config['Matrix'])
    plt.show()
