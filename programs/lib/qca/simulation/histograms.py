#!/usr/bin/python3
# =============================================================================
# Modification of load_data.py to make histograms
# By Logan Hillberry
# =============================================================================


import h5py
import numpy as np
from itertools import product
import matplotlib.pyplot as plt
import simulation.fio as io
import simulation.plotting as pt

def main():



    # locate the data you're interested in
    output_dir = 'fock_IC'
    mount_point = '/mnt/ext0/'
    data_repo = mount_point + 'qca_output/' + output_dir + '/data/'
    #data_repo = None


    # describe the simulations you'd like to load

    # params looped through at outer layer
    fixed_params_dict = {
                'output_dir' : [output_dir],
                'T'          : [1000],
                'BC'         : ['1_00'],
                'mode'       : ['alt'],
                'S'          : [1,2,3,4,5,6,7,8,9,10,11,12,13,14]
                 }

    # params looped through at inner layer
    var_params_dict = {
                'L'   : [20],
                'V'   : ['HP_'+str(deg) for deg in [0, 45,90]],
                'IC'  : ['c3_f0-2'],
                 }

    # histogram parameters
    partition_size = 10
    n_partitions = 99
    include_remaining = True
    n_bins = 20

    barmode = False
    aggregate_mode = True

    measure_name = "$ND$"

    colormap = plt.cm.get_cmap('jet')

    # function to extract data to put into histogram
    def obtain_hist(data_set,i,j):
        # helper
        def delta(values):
            return np.abs(values[i:j-1] - values[i+1:j])

        #return data_set['Y'][i:j] # Y
        #return delta(data_set['Y']) # \Delta Y

        return data_set['ND'][i:j] # ND
        #return delta(data_set['ND']) # \Delta ND

        #return data_set['CC'][i:j] # ND
        #return delta(data_set['CC']) # \Delta ND

        #return data_set['m'][i:j] # Mutual information
        #return delta(data_set['m']) # Mutual information change

        #return data_set['s'][i:j] # s
        #return delta(data_set['s']) # \Delta s

        #return get_diag_vecs(data_set['zz'][i:j]) # Z
        #return np.abs(get_diag_vecs(data_set['zz'][i:j-1]) -
        #        get_diag_vecs(data_set['zz'][i+1:j]))  # \Delta Z


    #aggregate data

    params_list_list = io.make_params_list_list(fixed_params_dict, var_params_dict)
    fignum = 0
    for m, params_list in enumerate(params_list_list):

        all_means = []
        all_stds = []
        names = []
        labels = []
        print(fignum)
        fig = plt.figure(fignum,figsize=(4,3))
        ax = fig.add_subplot(1,1,1)
        for n, params in enumerate(params_list):
            if data_repo is not None:
                sname = io.sim_name(params)
                data_path = data_repo + sname + '_v0.hdf5'
            else:
                sname = io.sim_name(params)
                data_path = io.default_file_name(params, 'data', '.hdf5')

            data_set = h5py.File(data_path)
            names.append(sname)

            simulation_hists = []

            #aggregate
            means = []
            stds = []

            #bin range
            maximum = 0
            minimum = 0

            for i in range(0, partition_size*n_partitions, partition_size):

                # get data from helper defined above
                partition_hist = obtain_hist(data_set,i,i+partition_size)

                # dataset may be multi-dimensional. reduce to one dimension
                while (len(partition_hist.shape) > 1):
                    partition_hist = np.concatenate(partition_hist)

                simulation_hists.append((partition_hist,i,i+partition_size))

                means.append(np.mean(partition_hist))
                stds.append(np.std(partition_hist))

                # update minima and maxima
                if (np.max(partition_hist) > maximum): maximum = np.max(partition_hist)
                if (np.min(partition_hist) > minimum): minimum = np.min(partition_hist)

            # same thing for the remaining iterations
            if include_remaining:
                low = partition_size*n_partitions
                partition_hist = obtain_hist(data_set,low,1000)
                while (len(partition_hist.shape) > 1):
                    partition_hist = np.concatenate(partition_hist)
                simulation_hists.append((partition_hist,low,1000))
                means.append(np.mean(partition_hist))
                stds.append(np.std(partition_hist))
                if (np.max(partition_hist) > maximum): maximum = np.max(partition_hist)
                if (np.min(partition_hist) > minimum): minimum = np.min(partition_hist)

            means = np.array(means)
            stds = np.array(stds)

            all_means.append(means)
            all_stds.append(stds)
            labels.append(r'$\theta = '+str(pt.get_th(params['V']))+'^\circ$')


            #compute bin range
            if (minimum == maximum):
                raise ValueError("Values are constant throughout simulation.")
            bins = np.arange(minimum,maximum,(maximum-minimum)/n_bins)
            maxbin = 0


            if  aggregate_mode: continue

            plt.clf()

            idx = 0
            for (values, i, j)  in simulation_hists:
                #make histogram
                hist, _ = np.histogram(values,bins=bins)
                hist = hist/(j-i) #normalize

                # help decide y axis maximum
                if (np.max(hist) > maxbin): maxbin = np.max(hist)

                fraction = idx/len(simulation_hists)
                color  = colormap(fraction)

                if barmode:
                    #plot as bar graph
                    width = (0.5 + 0.5*fraction) * (bins[1] - bins[0])
                    center = (bins[:-1] + bins[1:]) / 2
                    plt.bar(center, hist, align='center', width=width, color=color, alpha=0.5, label="%d to %d" % (i,j))
                else:
                    #plot as filled region
                    centers = (bins[:-1] + bins[1:]) / 2
                    plt.fill_between(centers,0,hist,color=color,alpha=0.5,label="%d to %d" % (i,j))
                idx += 1


            #configure plot
            plt.xlim(bins[0],bins[-1])
            plt.ylim(0,maxbin+1)

            plt.title("Histograms of %s in %s in groups of %d iterations" % (measure_name, sname, partition_size))
            plt.ylabel("Frequency (normalized by partition size)")
            plt.xlabel(measure_name)

            plt.grid(True)
            plt.legend(loc='best')
            plt.show()

        if aggregate_mode:
            if include_remaining:
                agg_range = np.arange(n_partitions+1)
            else:
                agg_range = np.arange(n_partitions)

            for i in range(len(all_means)):

                fraction = i/len(all_means)
                means = all_means[i]
                stds = all_stds[i]
                label = labels[i]
                ax.plot(agg_range,means, color=colormap(fraction), label=label)
                ax.fill_between(agg_range, means+stds, means-stds, alpha=0.5,
                        color=colormap(fraction))
                #plt.plot(agg_range, [1.0/20]*len(means))

            ax.set_title(r'$S='+str(params['S'])+'$,'+r'  $\tau=' +
                    str(partition_size) + '$')
            ax.set_ylabel(measure_name)
            ax.set_xlabel("Partition index")
            ax.grid(True)
            ax.legend(loc='best')
            ax.set_xscale("log", nonposx='clip')
            ax.set_yscale("log", nonposy='clip')
        fignum += 1
    io.multipage('./../output/fock_IC/plots/L20_ND_avgs.pdf')
        #plt.show()

# set default behavior of the file
if __name__ == '__main__':
    main()
