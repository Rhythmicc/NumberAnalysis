import matplotlib.pyplot as plt
import numpy as np

config = {
    'Matrix': 'Am',
    'row_size': 3534,
    'T': 100
}

if __name__ == "__main__":
    with open('%s_time' % config['Matrix'], 'r') as f:
        data = [float(i) for i in f.read().split()]
    GS = data[:config['T']]
    sor = data[config['T']:config['T']*2]
    con = data[-config['T']:]
    for i in range(config['T']):
        GS[i] = 2 * config['row_size'] ** 2 / GS[i] / 1e6
        sor[i] = 2 * config['row_size'] ** 2 / sor[i] / 1e6
        con[i] = 2 * config['row_size'] ** 2 / con[i] / 1e6
    x = np.arange(1, config['T'] + 1, 1)
    aver_gs = sum(GS) / config['T']
    aver_sor = sum(sor) / config['T']
    aver_con = sum(con) / config['T']

    max_gs = max(GS)
    max_gs_index = GS.index(max_gs) + 1
    min_gs = min(GS)
    min_gs_index = GS.index(min_gs) + 1

    max_sor = max(sor)
    max_sor_index = sor.index(max_sor) + 1
    min_sor = min(sor)
    min_sor_index = sor.index(min_sor) + 1

    max_con = max(con)
    max_con_index = con.index(max_con)
    min_con = min(con)
    min_con_index = con.index(min_con)

    print(max_gs, min_gs, max_sor, min_sor, max_con, min_con)

    plt.title('%s.txt aver_GS: %f' % (config['Matrix'], aver_gs))

    plt.plot(x, GS, label='GS')
    plt.scatter(x, GS, alpha=0.5, edgecolors='white')
    plt.scatter([max_gs_index], [max_gs], label='max_GS')
    plt.scatter([min_gs_index], [min_gs], label='min_GS')
    '''
    plt.plot(x, sor, label='sor')
    plt.scatter(x, sor, alpha=0.5, edgecolors='white')
    plt.scatter([max_sor_index], [max_sor], label='max_sor')
    plt.scatter([min_sor_index], [min_sor], label='min_sor')

    plt.plot(x, con, label='cg')
    plt.scatter(x, con, alpha=0.5, edgecolors='white')
    plt.scatter([max_con_index], [max_con], label='max_cg')
    plt.scatter([min_con_index], [min_con], label='min_cg')
    '''
    plt.xlabel('Times')
    plt.ylabel('MFlop/s')
    plt.legend(loc='best')
    plt.grid(linestyle='-.')
    plt.savefig('./img/%s_gs.png' % config['Matrix'])
    plt.show()
