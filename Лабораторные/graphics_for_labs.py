import matplotlib.pyplot as plt
import numpy as np

class Graphic:


    def _calc_Us(k, u):    
        ulist = u[k]
        return ulist
  
    def _calc_U(x, t, k, U, dx, l, r):
        xarr = np.arange(l, r + dx, dx) 
        Ulist = [U(x_, t[k]) for x_ in xarr]
        return xarr, Ulist

    def _draw_level(xarr, x, ulist, Ulist, axes, m, n):
        axes[m, n].plot(xarr, Ulist)
        axes[m, n].scatter(x, ulist)
        axes[m, n].grid()

    def draw_levels(x, t, u, U, l, r, dx=0.1):
        fig, axes = plt.subplots(2, 2)
        fig.set_figheight(10)
        fig.set_figwidth(15)
        dt = len(t) // 4
        k = 1
        for i in range(2):
            for j in range(2):
                ulist = Graphic._calc_Us(k * dt, u) 
                xarr, Ulist = Graphic._calc_U(x, t, k * dt, U, dx, l, r)
                Graphic._draw_level(xarr, x, ulist, Ulist, axes, i, j)
                k += 1
            

    def _calc_vars_t(x, t, u, U):
        epss = []
        for i, t_ in enumerate(t):
            s = 0
            for j, x_ in enumerate(x):
                s += (U(x_, t_) - u[i][j]) ** 2
            epss.append(s ** 0.5)
        return epss

    def _calc_vars_x(x, t, u, U):
        epss = []
        for i, x_ in enumerate(x):
            s = 0
            for j, t_ in enumerate(t):
                s += (U(x_, t_) - u[j][i]) ** 2
            epss.append(s ** 0.5)
        return epss

    def draw_variance(x, t, u, U):
        plt.figure(figsize=(10,7))
        plt.plot(t, Graphic._calc_vars_t(x, t, u, U))
        plt.plot(x, Graphic._calc_vars_x(x, t, u, U))
        plt.legend(["Ошибка по t", "Ошибка по x"])
        plt.grid()
        plt.show()
    
    
    def draw_for_compare(x, t, u_1, u_2, U, l, r, level, dx=0.1):
        xarr, Ulist = Graphic._calc_U(x, t, level, U, dx, l, r)
        ulist_1 = Graphic._calc_Us(level, u_1)
        ulist_2 = Graphic._calc_Us(level, u_2)
    
        plt.figure(figsize=(10,7))
        plt.plot(xarr, Ulist)
        plt.scatter(x, ulist_1, color="blue")
        plt.scatter(x, ulist_2, color="red")
        plt.legend(["Аналитическое решение", "u1", "u2"])
        plt.grid()
        plt.show()
