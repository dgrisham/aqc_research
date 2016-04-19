#!/usr/bin/python3

import matplotlib.pyplot as plt
import matplotlib as mpl
import simulation.plotting as pt
import matplotlib.gridspec as gridspec
import simulation.fio as io
import simulation.measures as ms
import numpy as np
import h5py

from matplotlib import animation

font = {'size':12, 'weight' : 'normal'}
mpl.rcParams['mathtext.fontset'] = 'stix'
mpl.rcParams['font.family'] = 'STIXGeneral'
mpl.rc('font',**font)



# First set up the figure, the axis, and the plot element we want to animate

# initialization function: plot the background of each frame
def init():
        grd.set_data([[],[]])
        return grd,

# animation function.  This is called sequentially
def animate(i):
        page = grid_list[i]
        grd.set_data(page)
        return grd,


def make_V_name(V):
    if V in ['X', 'Y', 'Z']:
        name = '\sigma^{}'.format(V.lower())
    else:
        name = V
    if len(name.split('_')) == 2:
        name = '('.join(name.split('_'))+'^\circ)'
    else:
        name = '('.join(name.split('_'))
    return name

def make_mode_name(mode):
    if mode is 'sweep':
        name = '\mathrm{SWP}'
    elif mode is 'alt':
        name = '\mathrm{ALT}'
    elif mode is 'block':
        name = '\mathrm{BLK}'
    else:
        name =''
    return name

def make_U_name(mode, S, V):
    S = str(S)
    name = r'$U^{%s}_{%s}(%s)$' % (make_mode_name(mode), S, make_V_name(V))

    return name




output_dir = 'fock_IC'
data_repo = '/mnt/ext0/qca_output/'+output_dir+'/data/'
#data_repo = None

deg_list = range(0, 105, 15)
deg_list = [0]
fixed_params_dict = {
            'output_dir' : [output_dir],
            'L' : [17],
            'T' : [1000],
            'IC': ['c3_f1'],
            'BC': ['1_00'],
            'mode': ['alt'],
            'S' : [14]
             }

var_params_dict = {
            'V' : ['HP_' + str(deg) for deg in deg_list]
             }

params_list_list = io.make_params_list_list(fixed_params_dict, var_params_dict)


def plot_grid(grid, ax, span=[0,100], n_xticks = 4, n_yticks = 6):
    im = ax.imshow( grid,
                    origin = 'lower',
                    vmin = 0.0,
                    vmax = 1.0,
                    interpolation = 'none',
                    aspect = '1',
                    rasterized = True)

    x_tick_labels = range(len(grid[0]))
    y_tick_labels = range(*span)

    ylabel = 'Iteration'
    xlabel = 'Site'
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    delta = max(1, int(len(x_tick_labels)/n_xticks))
    ax.set_xticks(range(0, len(x_tick_labels), delta ))
    ax.set_xticklabels(x_tick_labels[::delta])

    delta = max(1, int(len(y_tick_labels)/n_yticks))
    ax.set_yticks(range(0, len(y_tick_labels), delta ))
    ax.set_yticklabels(y_tick_labels[::delta])

    box = ax.get_position()
    cax = plt.axes([box.x1-0.1, box.y0+0.1, 0.04, box.height - 0.19])
    cb = plt.colorbar(im, cax = cax, ticks = [-1.0, -0.5, 0.0, 0.5, 1.0])
    cb.ax.tick_params(labelsize=12)
    cb.set_label(r'$\langle \sigma_j^z \rangle$', rotation=0, labelpad = -22,
            y=1.12)

    return im


grid_list = []
title_list = []

n_xticks = 4
n_yticks = 5
span = [100, 200]
for params_list in params_list_list:
    for params in params_list:
        output_dir = params['output_dir']
        mode = params['mode']
        S = params['S']
        V = params['V']
        T = params['T']
        L = params['L']
        IC = params['IC']

        title = make_U_name(mode, S, V)

        if data_repo is not None:
            sname = io.sim_name(params)
            res_path = data_repo + sname + '_v0.hdf5'
        else:
            res_path = io.default_file_name(params, 'data', '.hdf5')
        
        res = h5py.File(res_path)

        exp = ms.get_diag_vecs(res['zz'][::])
        s = res['s'][::]
        grid = exp[span[0]:span[1], 0:L]

        grid_list.append(grid)

        title_list.append(title)






fig = plt.figure(figsize=(2,3))
ax = fig.add_subplot(111)
im=plot_grid(grid_list[0], ax, span=span)
title = ax.set_title(title_list[0], fontsize=10)
fig.subplots_adjust(bottom=0.2)

def init():
    im.set_data(grid_list[0])
    title.set_text(title_list[0])
    return (im, title)

# animation function.  This is called sequentially
def animate(i):
    a=im.get_array()
    a=grid_list[i]
    im.set_array(a)
    title.set_text(title_list[i])
    return (im, title)


anim = animation.FuncAnimation(fig, animate, init_func=init,
                                       frames=len(grid_list), interval=150)

bn = io.base_name(output_dir, 'plots')
#anim.save(bn + 'L21_alt_zavg.mp4', fps=8, extra_args=['-vcodec', 'libx264'])

plt.show()
