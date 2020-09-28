import matplotlib.pyplot as plt
import numpy as np

config = {
    'Matrix': 'Ac',
    'row_size': 1182,
    'T': 20
}

if __name__ == "__main__":
    with open('%s_time' % config['Matrix'], 'r') as f:
        data = [float(i) for i in f.read().split()]
    doo = data[:config['T']]
    cro = data[config['T']:]
    for i in range(config['T']):
        doo[i] = 2 * config['row_size'] ** 2 / doo[i] / 1e6
        cro[i] = 2 * config['row_size'] ** 2 / cro[i] / 1e6
    x = np.arange(1, config['T'] + 1, 1)
    aver_doo = sum(doo) / config['T']
    aver_cro = sum(cro) / config['T']

    max_doo = max(doo)
    max_doo_index = doo.index(max_doo) + 1
    min_doo = min(doo)
    min_doo_index = doo.index(min_doo) + 1

    max_cro = max(cro)
    max_cro_index = cro.index(max_cro) + 1
    min_cro = min(cro)
    min_cro_index = cro.index(min_cro) + 1
    print(max_doo, min_doo, max_cro, min_cro)

    plt.title('%s.txt aver_cro: %f' % (config['Matrix'], aver_cro))
    '''
    plt.plot(x, doo, label='Lu_doolittle')
    plt.scatter(x, doo, alpha=0.5, edgecolors='white')
    plt.scatter([max_doo_index], [max_doo], label='max_doo')
    plt.scatter([min_doo_index], [min_doo], label='min_doo')
    '''
    plt.plot(x, cro, label='Lu_crout')
    plt.scatter(x, cro, alpha=0.5, edgecolors='white')
    plt.scatter([max_cro_index], [max_cro], label='max_cro')
    plt.scatter([min_cro_index], [min_cro], label='min_cro')

    plt.xlabel('Times')
    plt.ylabel('MFlop/s')
    plt.legend(loc='best')
    plt.grid(linestyle='-.')
    plt.savefig('./img/%s_cro.png' % config['Matrix'])
    plt.show()
