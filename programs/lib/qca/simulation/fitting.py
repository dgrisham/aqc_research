#!/usr/bin/python3

from math import sqrt, e, pi
import numpy as np
import scipy.odr.odrpack as odrpack
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy import interpolate
from os import environ
import simulation.plotting as pt

font = {'family':'normal', 'weight':'bold', 'size':16}
mpl.rc('font',**font)



# Orthogonal distance regression
# ------------------------------
def flin(B, x):
    return B[0]*x + B[1]

def fpow(B, x):
    return B[1]*x**B[0] + B[2]

def gaussian(B, x):
    return 1.0/(B[0]*sqrt(2*pi)) * e ** ( - ( x - B[1] )**2 / (2 * B[0]**2) )

def ftwo_gaussian(B, x):
    return 1.0/(2*B[0]*sqrt(2*pi)) * e ** ( - ( x - B[1] )**2 / (2 * B[0]**2) ) +\
           1.0/(2*B[2]*sqrt(2*pi)) * e ** ( - ( x - B[3] )**2 / (2 * B[2]**2) )

def fpoly(B, x):
    return sum([b*x**m for m, b in enumerate(B)])

def fexp(B, x):
    return B[1]*B[0]**(x)

def fnexp(B, x, base=e):
    return B[0] + sum([B[i]*base**(x*B[i+1]) for i in range(1, len(B) - 1)])

def fsexp(B, x):
    return B[0]*B[1]**(x**B[2]) + B[3]

def do_odr(f, x, xe, y, ye, estimates):
    model = odrpack.Model(f)
    data = odrpack.RealData(x, y, sx=xe, sy=ye)
    odr = odrpack.ODR(data, model, estimates)
    output = odr.run()
    return output

def chi2_calc(f, betas, x, y):
    chi2 = 0
    for xval, yval in zip(x,y):
        chi2 += (yval - f(betas, xval))**2
    return chi2

def f_fits(func, beta_est, x_list, y_list, xerr = None, yerr = None):
    xerr = [0.0000001] * len(x_list) if xerr is None else xerr
    yerr = [0.0000001] * len(y_list) if yerr is None else yerr
    fit = do_odr(func, x_list, xerr, y_list, yerr, beta_est)
    chi2 = chi2_calc(func, fit.beta, x_list, y_list)
    return fit.beta, chi2, fit.sd_beta

def plot_f_fits(func, beta_est, x_list, y_list, ax, label, color, xerr =
    None, yerr = None, kwargs={}):
    xerr = [0.0000001] * len(x_list) if xerr is None else xerr
    yerr = [0.0000001] * len(y_list) if yerr is None else yerr

    fit_beta, chi2,  = \
            f_fits(func, beta_est, x_list, y_list, xerr = xerr, yerr = yerr)

    xs = np.linspace(min(x_list), max(x_list), 100)
    ax.plot(x_list, y_list, 'o', label = 'tmax = ' + str(label), color = color)
    ax.legend(loc = 'lower right')
    ax.set_yscale('log', basey=2)
    return










if __name__ == '__main__':
    import simulation.fio as io
    import simulation.measures as ms
    import simulation.time_evolve as time_evolve
    import matplotlib.pyplot as plt
    import simulation.plotting as ptt
    import scipy.stats as stats
    import h5py
    

    params =  {
                    'output_dir' : 'Hphase',

                    'L'    : 19,
                    'T'    : 60,
                    'mode' : 'alt',
                    'S'    :  6,
                    'V'    : 'HP_90',
                    'IC'   : 'f0',
                    'BC'   : '1'
                                    }


    res = h5py.File(io.default_file_name(params, 'data', '.hdf5'))
    z_grid = ms.get_diag_vecs(res['zz'][::])

    fig = plt.figure(1)
    ax = fig.add_subplot(111)
    img = (1.0-z_grid)/2
    plt.imshow(img, origin='lower', aspect=1)
    for row in img:
        row = row/sum(row)
        plt.plot(row)
        jbar = sum(j*c for j,c in enumerate(row))
        m, s = stats.norm.fit(row)
        plt.scatter(jbar, row[int(jbar)])
        p = stats.norm.pdf(row, m, s)
        plt.plot(row, p)
        plt.show()

    '''

    def fit_speed(grid):
        fig = plt.figure(1)
        ax = fig.add_subplot(111)

        img = (1.0-grid)/2

        # block S6 12 H
        #ks = [0, 11, 21, 31, 41]

        # sweep S6 12 H
        ks = [0, 7, 21, 27, 40]
        for p in range(len(ks)-1):
            pimg = img[ks[p]:ks[p+1], 0:L]
            j = np.argmax(pimg, axis=1)
            pt = range(ks[p], ks[p+1])
            z = img[(pt, j)]
            dz = 1/z
            ax.errorbar(j, pt, xerr=dz, yerr=dz)
            js = np.linspace(0, L-1, 150)
            Bs, chi2 = f_fits(flin, [1.0, 0.0], j, pt, xerr=dz, yerr=dz)
            print(1/Bs[0])
            ax.plot(js, flin(Bs, js),  color='k')
        ax.set_ylim([-0.5, 75.5])
        ax.set_xlim([-0.5, (L-1)+0.5])

        plt.savefig(environ['HOME'] +
                '/documents/qca/notebook_figs/fit_sweep_L12_S6_VH.pdf',
            format='pdf', dpi=300, bbox_inches='tight')

    def fit_measure(x, y, func, Bsest):
        xs = np.linspace(min(x), max(x), 100)
        Bs, chi2 = f_fits(func, Bsest, x, y)
        print(Bs, chi2)
        plt.plot(xd, y, 'o', label = 'data')
        plt.plot(xs, func(Bs, xs))
        plt.show()

    data = np.array([x**2-.3*x for x in np.linspace(-2, 3, 30)])
    data = data+np.random.rand(len(data))
    fit_measure(data, fpoly, [1.0, 1.0, 1.0,1.0,1.0,1.0])

    xd = np.linspace(-10, 10, 20)
    data = np.array([e**(.2*x) + 0.1*e**(-.12*x) for x in xd])
    data = data+np.random.rand(len(data))
    fit_measure(xd, data, fnexp, [0.0, 1.0, 0.2, 0.1, -0.12])

    #fit_speed(z_grid)
    '''
