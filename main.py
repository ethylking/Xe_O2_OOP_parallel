from classes import *

if __name__ == '__main__':
    n = 5_000_000
    pool = multiprocessing.Pool()
    exp = pool.map(Calculation,[n for i in range(5)])
    FG = []
    for exp_i in exp:
        for i in range(n):
            exp_i.points[i].Fg()
            FG += [exp_i.points[i].fg]
        exp_i.sum_points()
        print(exp_i.S, exp_i.sigma)
