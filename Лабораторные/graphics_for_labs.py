import matplotlib.pyplot as plt
import numpy as np
import plotly
import plotly.graph_objs as go
from plotly.subplots import make_subplots
from mpl_toolkits.mplot3d import Axes3D

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
    
    def _calc_vars_xy(x, y, u, U):
        epss = []
        for i, x_ in enumerate(x):
            s = 0
            for j, y_ in enumerate(y):
                s += (U(x_, y_) - u[j][i]) ** 2
            epss.append(s ** 0.5)
        return epss
    
    def _calc_vars_yx(x, y, u, U):
        epss = []
        for i, y_ in enumerate(y):
            s = 0
            for j, x_ in enumerate(x):
                s += (U(x_, y_) - u[i][j]) ** 2
            epss.append(s ** 0.5)
        return epss

    def draw_variance(x, t, u, U):
        plt.figure(figsize=(10,7))
        plt.plot(t, Graphic._calc_vars_t(x, t, u, U))
        plt.plot(x, Graphic._calc_vars_x(x, t, u, U))
        plt.legend(["Ошибка по t", "Ошибка по x"])
        plt.grid()
        plt.show()
        
    def draw_variance_xy(x, y, u, U):
        plt.figure(figsize=(10,7))
        plt.plot(x, Graphic._calc_vars_xy(x, y, u, U))
        plt.plot(y, Graphic._calc_vars_yx(y, x, u, U))
        plt.legend(["Ошибка по x", "Ошибка по y"])
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
        
        
    def _calc_U_2d(U, t, m1, m2, lx, rx, ly, ry, dx, dy):
        xs = np.arange(lx, rx * m1 + dx, dx)
        ys = np.arange(ly, ry * m2 + dy, dy)
        U_ = []

        for x_ in xs:
            row = []
            for y_ in ys:
                row.append(U(x_, y_, t, m1, m2))
            U_.append(row)
        return xs, ys, U_
    
    def _calc_U_xy(U, lx, rx, ly, ry, dx, dy):
        xs = np.arange(lx, rx + dx, dx)
        ys = np.arange(ly, ry + dy, dy)
        U_ = []

        for x_ in xs:
            row = []
            for y_ in ys:
                row.append(U(x_, y_))
            U_.append(row)
        return xs, ys, U_
    
    def _calc_Us_2d(k, u):
        return u[k]
    
    def draw_surface(x, y, u, U, lx, rx, ly, ry):
        fig = plt.figure(figsize=(10, 10))
        ax = Axes3D(fig)
    
        dx = (rx - lx) / 1000
        dy = (ry - ly) / 1000
        xarr, yarr, Ulist = Graphic._calc_U_xy(U, lx, rx, ly, ry, dx, dy)
        ax.plot_surface(np.array(xarr), np.array(yarr), np.array(Ulist))
        ax.plot_wireframe(x, y, u, color="black")
        ax.set(xlabel='$x$', ylabel='$y$', zlabel='$U$')
        fig.tight_layout()
        
    def draw_3d(x, y, u, U, lx, rx, ly, ry, label):
        plotly.offline.init_notebook_mode()
    
        dx = (rx - lx) / 1000
        dy = (ry - ly) / 1000
        xs, ys, U_ = Graphic._calc_U_xy(U, lx, rx, ly, ry, dx, dy)

        fig = make_subplots(rows=1, cols=2,
                    specs=[[{'is_3d': True}, {'is_3d': True}]],
                    subplot_titles=[label, 'Solution'],
                    )

        fig.add_trace(go.Surface(x=x, y=y, z=u), 1, 1)
        fig.add_trace(go.Surface(x=xs, y=ys, z=U_), 1, 2)
        fig.show()

    def draw_surfaces_levels(x, y, t, m1, m2, u, U, lx, rx, ly, ry):
        fig, axes = plt.subplots(2, 2)
        fig.set_figheight(10)
        fig.set_figwidth(15)
    
        dx = (rx * m1 - lx) / 1000
        dy = (ry * m2 - ly) / 1000
        dt = len(t) // 4
        for k in range(4):
            ulist = Graphic._calc_Us_2d(k * dt, u) 
            xarr, yarr, Ulist = Graphic._calc_U_2d(U, k * dt, m1, m2, lx, rx, ly, ry, dx, dy)
            ax = fig.add_subplot(2, 2, k + 1, projection='3d')
            ax.plot_surface(np.array(xarr), np.array(yarr), np.array(Ulist))
            ax.plot_wireframe(x, y, ulist, color="black")
            ax.set(xlabel='$x$', ylabel='$y$', zlabel='$U$')
            fig.tight_layout()
    
    def draw_level_3d(x, y, k, m1, m2, u, U, lx, rx, ly, ry, label):
        plotly.offline.init_notebook_mode()
    
        dx = (rx * m1 - lx) / 1000
        dy = (ry * m2 - ly) / 1000
        xs, ys, U_ = Graphic._calc_U_2d(U, k, m1, m2, lx, rx, ly, ry, dx, dy)

        fig = make_subplots(rows=1, cols=2,
                    specs=[[{'is_3d': True}, {'is_3d': True}]],
                    subplot_titles=[label, 'Solution'],
                    )

        fig.add_trace(go.Surface(x=x, y=y, z=u[k]), 1, 1)
        fig.add_trace(go.Surface(x=xs, y=ys, z=U_), 1, 2)
        fig.update_layout(title_text="Sulutions at t = {}".format(k))
        fig.show()