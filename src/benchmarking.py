# -*- coding: utf-8 -*-
# # 4.4. Benchmarking

# # 4.4.1. Kulmala and Su
#
# Calculating evaporation rate over RH

if __name__ == '__main__':
    Kulmala_data = np.load('src/kulmala_data.npy', allow_pickle = True)
    Kulmala_RH_list = np.load('src/kulmala_rh_list.npy', allow_pickle = True)

    Su_data = np.load('src/Su_data.npy', allow_pickle = True)
    Su_RH_list = np.load('src/Su_rh_list.npy', allow_pickle = True)

if __name__ == '__main__':
    solution = aqueous_NaCl
    R0 = 25e-6 # metres
    T = 298.15 # Kelvin
    mfs = 0
    time = 25 # seconds

    history_list = []

    #RH_range = np.linspace(0, 1, 100)**0.5
    RH_range = np.linspace(0, 1, 21)

    for RH in RH_range:
        print('\r'+str(round(RH,3)), end = '')
        droplet, trajectory = simulate(time, solution, T, RH, R0, T, mfs)

        # Obtain a table giving a history of *all* droplet parameters.
        history = droplet.complete_trajectory(trajectory)
        history_list.append(history)

    sadkat_data = history_list
    sadkat_RH_list = RH_range

if __name__ == '__main__':
    def trim_data(time_data, property_data):
        '''
        Clips off first 10% and last 10% of data for calcualtion of evaporation rate.
        '''

        time_start = time_data.max() * 0.1
        time_end = time_data.max() * 0.9

        selection = (time_data >= time_start) & (time_data <= time_end)

        time_data = time_data[selection]

        property_data = property_data[selection]

        return time_data, property_data

    def get_evaporation_rate(t, r2):
        '''
        Takes time data and r^2 data and calculated evaporation data.

        '''
        return -np.gradient(r2, t)

    def get_temperature_suppression(T_droplet, T_environment = 298.15):
        '''
        Takes time data and r^2 data and calculated evaporation data.

        '''
        return T_environment - T_droplet

if __name__ == '__main__':
    #preparing Kulmala data
    for dataset, RH in zip(Kulmala_data, Kulmala_RH_list):
        dataset['k_um2pers'] = get_evaporation_rate(dataset.t_s, dataset.R2_um2)
        dataset['radius2_um2'] = dataset.R2_um2
        dataset['radius_um'] = np.sqrt(dataset.R2_um2)
        dataset['time_s'] = dataset.t_s
        dataset['temperature_K'] = dataset.T_K
        dataset['delta_T_K'] = get_temperature_suppression(dataset.temperature_K)

    Kulmala_k_means = []
    Kulmala_T_suppression = []
    for dataset, RH in zip(Kulmala_data, Kulmala_RH_list):
        Kulmala_k_means.append(trim_data(dataset.t_s, dataset.k_um2pers)[1].mean())
        Kulmala_T_suppression.append(trim_data(dataset.time_s, dataset.delta_T_K)[1].mean())

    #4th order polynomial fit of Kulmala k vs RH
    Kulmala_poly_model = np.poly1d(np.polyfit(Kulmala_RH_list,
                                              Kulmala_k_means,
                                              4))
    #4th order polynomial fit of Su ∆T vs RH
    Kulmala_T_poly_model = np.poly1d(np.polyfit(Kulmala_RH_list,
                                                Kulmala_T_suppression,
                                                4))
    #preparing Su data
    for dataset, RH in zip(Su_data, Su_RH_list):
        dataset['radius_um'] = dataset.radius_m * 1e6
        dataset['radius2_um2'] = dataset.radius_um ** 2
        dataset['k_um2pers'] = get_evaporation_rate(dataset.time_s, dataset.radius2_um2)
        dataset['temperature_K'] = dataset.delta_T + 298.15
        dataset['delta_T_K'] = - dataset.delta_T

    Su_k_means = []
    Su_T_suppression = []
    for dataset, RH in zip(Su_data, Su_RH_list):
        Su_k_means.append(trim_data(dataset.time_s, dataset.k_um2pers)[1].mean())
        Su_T_suppression.append(trim_data(dataset.time_s, dataset.delta_T_K)[1].mean())

    #4th order polynomial fit of Su k vs RH
    Su_poly_model = np.poly1d(np.polyfit(Su_RH_list,
                                               Su_k_means,
                                               4))

    #4th order polynomial fit of Su k vs RH
    Su_T_poly_model = np.poly1d(np.polyfit(Su_RH_list,
                                           Su_T_suppression,
                                           4))
    #preparing SADKAT data
    for dataset, RH in zip(sadkat_data, sadkat_RH_list):
        dataset['time_s'] = dataset.time
        dataset['radius_um'] = dataset.radius * 1e6
        dataset['radius2_um2'] = dataset.radius_um ** 2
        dataset['k_um2pers'] = get_evaporation_rate(dataset.time, dataset.radius2_um2)
        dataset['temperature_K'] = dataset.temperature
        dataset['delta_T_K'] = get_temperature_suppression(dataset.temperature_K)

    sadkat_k_means = []
    sadkat_T_suppression = []
    for dataset, RH in zip(sadkat_data, sadkat_RH_list):
        sadkat_k_means.append(trim_data(dataset.time, dataset.k_um2pers)[1].mean())
        sadkat_T_suppression.append(trim_data(dataset.time_s, dataset.delta_T_K)[1].mean())

    #4th order polynomial fit of sadkat k vs RH
    sadkat_poly_model = np.poly1d(np.polyfit(sadkat_RH_list,
                                               sadkat_k_means,
                                               4))

    #4th order polynomial fit of Su ∆T vs RH
    sadkat_T_poly_model = np.poly1d(np.polyfit(sadkat_RH_list,
                                           sadkat_T_suppression,
                                           4))

if __name__ == '__main__':
    def plot_rh_dependent_data(data, RH_list,
                               k_means, k_poly_model,
                               delta_T_means, delta_t_poly_model,
                               cmap = mpl.cm.viridis):

        norm = mpl.colors.Normalize(vmin= RH_list.min(), vmax=RH_list.max())
        colours = cmap(RH_list)

        #plot radius
        for dataset, RH, colour in zip(data, RH_list, colours):
            plt.plot(dataset.time_s.values,
                     dataset.radius_um.values,
                     color = colour)
        plt.xlabel('Time / s')
        plt.ylabel('Radius / µm')
        plt.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap), label = '% RH' )
        plt.show()

        #plot r^2
        for dataset, RH, colour in zip(data, RH_list, colours):
            plt.plot(dataset.time_s.values,
                     dataset.radius2_um2.values,
                     color = colour)
        plt.xlabel('Time / s')
        plt.ylabel('Radius$^2$ / µm$^2$')
        plt.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap), label = '% RH' )
        plt.show()

        # plot evaporation rate
        for dataset, RH, colour in zip(data, RH_list, colours):
            plt.plot(dataset.time_s.values,
                     dataset.k_um2pers.values,
                     color = colour, ls = ':')
            plt.plot(trim_data(dataset.time_s, dataset.k_um2pers)[0],
                     trim_data(dataset.time_s, dataset.k_um2pers)[1],
                     color = colour)
        plt.xlabel('RH')
        plt.ylabel(r'Evaporation Rate / µm$^2$s$^{-1}$')
        plt.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap), label = '% RH' )
        plt.show()

        #plot k vs RH
        for k, RH, colour in zip(k_means, RH_list, colours):
            plt.scatter(RH,
                        k,
                        color = colour, s = 100)
        plt.plot(np.linspace(0,1,100),
                 k_poly_model(np.linspace(0,1,100)),
                 lw = 0.5,
                 color = 'k')
        plt.xlabel('RH')
        plt.ylabel(r'Evaporation Rate / µm$^2$s$^{-1}$')
        plt.show()

        # plot temperature
        for dataset, RH, colour in zip(data, RH_list, colours):
            plt.plot(dataset.time_s.values,
                     dataset.temperature_K.values,
                     color = colour, ls = ':')
            plt.plot(trim_data(dataset.time_s, dataset.temperature_K)[0],
                     trim_data(dataset.time_s, dataset.temperature_K)[1],
                     color = colour)
        plt.xlabel('RH')
        plt.ylabel(r'Droplet Temperature / K')
        plt.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap), label = '% RH' )
        plt.show()

        #plot T_suppression vs RH
        for T, RH, colour in zip(delta_T_means, RH_list, colours):
            plt.scatter(RH,
                        T,
                        color = colour, s = 100)
        plt.plot(np.linspace(0,1,100),
                 delta_t_poly_model(np.linspace(0,1,100)),
                 lw = 0.5,
                 color = 'k')
        plt.xlabel('RH')
        plt.ylabel(r'Temperature Suppression / K')
        plt.show()

        return

if __name__ == '__main__':
    plot_rh_dependent_data(Kulmala_data, Kulmala_RH_list,
                           Kulmala_k_means, Kulmala_poly_model,
                           Kulmala_T_suppression, Kulmala_T_poly_model,
                           mpl.cm.winter_r)

    plot_rh_dependent_data(Su_data, Su_RH_list,
                           Su_k_means, Su_poly_model,
                           Su_T_suppression, Su_T_poly_model,
                           mpl.cm.plasma_r)

    plot_rh_dependent_data(sadkat_data, sadkat_RH_list,
                           sadkat_k_means, sadkat_poly_model,
                           sadkat_T_suppression, sadkat_T_poly_model,
                           mpl.cm.cool_r)

if __name__ == '__main__':

    from matplotlib.patches import Patch
    from matplotlib.lines import Line2D

    #comparing evaporation rate data

    saturation_to_percent_convert = 100
    model_rh_range = np.linspace(0,100,101)
    model_saturation_range = np.linspace(0,1,101)

    poly_line_width = 1.5

    fig, ax = plt.subplots(dpi = resolution)

    for k, RH in zip(Kulmala_k_means, Kulmala_RH_list):
        ax.scatter(RH * saturation_to_percent_convert,
                    k,
                    color = '#004D40', s = 50, marker = 'P')
    ax.plot(model_rh_range,
             Kulmala_poly_model(model_saturation_range),
             poly_line_width,
             color = '#004D40')

    for k, RH in zip(Su_k_means, Su_RH_list):
        ax.scatter(RH * saturation_to_percent_convert,
                    k,
                    color = '#1E88E5', s = 50, marker = 'x')
    ax.plot(model_rh_range,
             Su_poly_model(model_saturation_range),
             poly_line_width,
             color = '#1E88E5')

    for k, RH in zip(sadkat_k_means, sadkat_RH_list):
        ax.scatter(RH * saturation_to_percent_convert,
                    k,
                    color = '#D81B60', s = 50, marker = 'o')
    ax.plot(model_rh_range,
             sadkat_poly_model(model_saturation_range),
             poly_line_width,
             color = '#D81B60')

    ax.set_xlabel('RH / %')
    ax.set_ylabel(r'Evaporation Rate / µm$^2$s$^{-1}$')

    legend_elements = [Line2D([0], [0], marker='P', color='#004D40', label='Kulmala',
                              markerfacecolor='#004D40', markersize=11, lw = 1),
                       Line2D([0], [0], marker='x', color='#1E88E5', label='Su',
                              markerfacecolor='#1E88E5', markersize=11, lw = 1),
                       Line2D([0], [0], marker='o', color='#D81B60', label='SADKAT',
                              markerfacecolor='#D81B60', markersize=11, lw = 1),]

    ax.legend(handles=legend_elements)

    plt.show()

    fig, ax1 = plt.subplots(dpi = resolution)

    # comparing temperature suppresssion data
    for k, RH in zip(Kulmala_T_suppression, Kulmala_RH_list):
        ax1.scatter(RH * saturation_to_percent_convert,
                    k,
                    color = '#004D40', s = 50, marker = 'P')
    ax1.plot(model_rh_range,
             Kulmala_T_poly_model(model_saturation_range),
             poly_line_width,
             color = '#004D40')

    for T, RH in zip(Su_T_suppression, Su_RH_list):
        ax1.scatter(RH * saturation_to_percent_convert,
                    T,
                    color = '#1E88E5', s = 50, marker = 'x')
    ax1.plot(model_rh_range,
             Su_T_poly_model(model_saturation_range),
             poly_line_width,
             color = '#1E88E5')
    ax1.set_xlabel('RH / %')

    for T, RH in zip(sadkat_T_suppression, sadkat_RH_list):
        ax1.scatter(RH * saturation_to_percent_convert,
                    T,
                    color = '#D81B60', s = 50)
    ax1.plot(model_rh_range,
             sadkat_T_poly_model(model_saturation_range),
             poly_line_width,
             color = '#D81B60')
    ax1.set_xlabel('RH / %')
    ax1.set_ylabel(r'Temperature Suppression / K')
    legend_elements = [Line2D([0], [0], marker='P', color='#004D40', label='Kulmala',
                              markerfacecolor='#004D40', markersize=11, lw = 1),
                       Line2D([0], [0], marker='x', color='#1E88E5', label='Su',
                              markerfacecolor='#1E88E5', markersize=11, lw = 1),
                       Line2D([0], [0], marker='o', color='#D81B60', label='SADKAT',
                              markerfacecolor='#D81B60', markersize=11, lw = 1),
                       Line2D([0], [0], ls = ':', lw = 2, color = 'k', label = '∆T = 3 K'),]

    ax1.axhline(3, ls = ':', lw = 2, color = 'k', label = '∆T = 3 K' )

    ax2 = inset_axes(ax1,
                     width="100%", height="100%",
                     bbox_to_anchor=(.5, .5, .45, .45),bbox_transform=ax1.transAxes, loc=3)#fig.add_axes([0.18, 0.2, 0.3, 0.3],)
    ax2.patch.set_alpha(0)

    for k, RH in zip(Kulmala_T_suppression, Kulmala_RH_list):
        ax2.scatter(RH * saturation_to_percent_convert,
                    k,
                    color = '#004D40', s = 50, marker = 'P')
    ax2.plot(model_rh_range,
             Kulmala_T_poly_model(model_saturation_range),
             poly_line_width,
             color = '#004D40')

    for T, RH in zip(Su_T_suppression, Su_RH_list):
        ax2.scatter(RH * saturation_to_percent_convert,
                    T,
                    color = '#1E88E5', s = 50, marker = 'x')
    ax2.plot(model_rh_range,
             Su_T_poly_model(model_saturation_range),
             poly_line_width,
             color = '#1E88E5')

    for T, RH in zip(sadkat_T_suppression, sadkat_RH_list):
        ax2.scatter(RH * saturation_to_percent_convert,
                    T,
                    color = '#D81B60', s = 50)
    ax2.plot(model_rh_range,
             sadkat_T_poly_model(model_saturation_range),
             poly_line_width,
             color = '#D81B60')
    ax2.axhline(3, xmax = 0.65, ls = ':', lw = 2, color = 'k', label = '∆T = 3 K' )

    ax2.set_xlim(0.75 * saturation_to_percent_convert, 1 * saturation_to_percent_convert)
    ax2.set_ylim(0,3.5)


    ax1.legend(handles=legend_elements)
    plt.show()

# # 4.4.2. Loading experimental data for benchmarking

if __name__ == '__main__':
    import os

    def load_experimental_data():

        # function to walk though files to find edb data and conditions
        def list_files(dir):                                                                                                  
            r = []                                                                                                            
            subdirs = [x[0] for x in os.walk(dir)]                                                                            
            for subdir in subdirs:  
                print(subdir)
                if 'p' in subdir:
                    print('probe')
                if 's' in subdir:
                    print('sample')
                files = os.walk(subdir).__next__()[2]                                                                             
                if (len(files) > 0):                                                                                          
                    for file in files:                                                                                        
                        r.append(os.path.join(subdir, file))
                        print(file)
            return r 


        #get list of experiments and make a dict
        experiments_list = []


        for i in os.walk('src/EDB data for benchmarking'):

            experiment_folder = os.path.basename(os.path.normpath(i[0]))
            #clip out the id part of the folder
            if experiment_folder == 'EDB data for benchmarking':
                pass
            else:
                experiment_id = experiment_folder[0:13]
                experiments_list.append(experiment_id)

        #list
        experiments_list = list(set(experiments_list))
        # dict
        experiment_ids = {exp_id: {} for exp_id in experiments_list}

        #go through again and get smaples and probes
        for root, dirs, files in os.walk('src/EDB data for benchmarking'):

            experiment_folder = os.path.basename(os.path.normpath(root))
            #clip out the id part of the folder
            if experiment_folder == 'EDB data for benchmarking':
                pass
            else:
                experiment_id = experiment_folder[0:13]
                experiments_list.append(experiment_id)

            if experiment_folder.find('p') < 14 and experiment_folder.find('p') > 0:
                experiment_ids[experiment_id]['p'] = root

                #extracting temperature while considering probe file
                experiment_ids[experiment_id]['T'] = experiment_folder[experiment_folder.find('p')+1:experiment_folder.find('deg')]

            if experiment_folder.find('s') < 14 and experiment_folder.find('s') > 0:
                experiment_ids[experiment_id]['s'] = root

        for key in experiment_ids:
            #getting RH values, and renumber experiements to match droplet files
            rh_filepath = experiment_ids[key]['p'] + '/RH results file.txt'
            experiment_ids[key]['exp_props'] = pd.read_csv(rh_filepath, sep = '\t')
            experiment_ids[key]['exp_props'].rename(columns = {'Droplet':'droplet', 'Best RH':'RH'},inplace=1)
            experiment_ids[key]['exp_props'].droplet = experiment_ids[key]['exp_props'].droplet -1

            #adding temp to df
            experiment_ids[key]['exp_props']['T'] = float(experiment_ids[key]['T'])

            filenames = next(os.walk(experiment_ids[key]['s']), (None, None, []))[2]  # [] if no file. get list of files in sample folder
            filenames.remove('All data')
            filenames.sort()

            filenames = [experiment_ids[key]['s'] +'/' + element for element in filenames]
            experiment_ids[key]['exp_props']['droplet_path'] = filenames

            experiment_ids[key]['exp_props']['k'] = np.nan

            for i in experiment_ids[key]['exp_props']['droplet']:
                droplet_data = np.loadtxt(experiment_ids[key]['exp_props'].droplet_path[0]).T
                droplet_data[1] = droplet_data[1]/1e6

                '''plt.plot(droplet_data[0], droplet_data[1]**2)
                plt.show()
                plt.plot(droplet_data[0], np.gradient(droplet_data[1]**2, droplet_data[0]))
                plt.axhline(np.gradient(droplet_data[1]**2, droplet_data[0]).mean(), c= 'k')
                plt.show()'''

                experiment_ids[key]['exp_props']['k'][i] = -np.gradient(droplet_data[1]**2, droplet_data[0]).mean()

        # sort by temperature
        experiment_ids = dict(sorted(experiment_ids.items(), key=lambda item: (item[1]['T'], item[1]['exp_props'].RH.mean())))

        for key in experiment_ids:
            plt.scatter(experiment_ids[key]['exp_props']['RH'], experiment_ids[key]['exp_props']['k'], label = (experiment_ids[key]['exp_props']['T'][0],experiment_ids[key]['exp_props']['RH'][0]))
        plt.legend()
        plt.ylabel('Evaporation Rate / m$^2$s$^{-1}$')
        plt.xlabel('Relative Humidity')
        plt.show()

        T_list = np.array(list(set(list(float(experiment_ids[key]['T']) for key in experiment_ids ))))
        T_list.sort()
        print("Temperatures: ", T_list)

        return experiment_ids, T_list + T_freezing

if __name__ == '__main__':
    experiments, experiment_T_range = load_experimental_data()

# # 4.4.3 Iterating over a 2 dimensional parameter space
#
# A sample range of T and RH

if __name__ == '__main__':
    solution = aqueous_NaCl
    R0 = 25e-6 # metres
    mfs = 0
    time = 25 # seconds

    history_list = []
    T_range = np.linspace(14.85 + T_freezing, 29.85 + T_freezing,16) # Kelvin
    #T_range = np.append(T_range, np.array(T_list) + T_freezing) 
    #T_range.sort()

    RH_range = np.linspace(0,1, 11) # % RH

    print(T_range)
    print(RH_range)

    for T in T_range:
        T_history_list = []
        for RH in RH_range:
            print('\r'+str(round(RH,3)), ' ',str(round(T,3)),  end = '')
            droplet, trajectory = simulate(time, solution, T, RH, R0, T, mfs)

            # Obtain a table giving a history of *all* droplet parameters.
            history = droplet.complete_trajectory(trajectory)
            #history_list.append(history)
            T_history_list.append(history)
        history_list.append(T_history_list) 

    sadkat_2d_data = history_list
    sadkat_RH_range = RH_range
    sadkat_T_range = T_range

# The range of T used in experiemts

if __name__ == '__main__':
    solution = aqueous_NaCl
    R0 = 25e-6 # metres
    mfs = 0
    time = 25 # seconds

    history_list = []

    #RH_range = np.linspace(0,1, 11) # % RH

    for T in experiment_T_range:
        T_history_list = []
        for RH in RH_range:

            print('\r'+str(round(RH,3)), ' ',str(round(T,3)),  end = '')
            droplet, trajectory = simulate(time, solution, T, RH, R0, T, mfs)

            # Obtain a table giving a history of *all* droplet parameters.
            history = droplet.complete_trajectory(trajectory)
            #history_list.append(history)
            T_history_list.append(history)
        history_list.append(T_history_list) 

    sadkat_experiemtal_comparison_data = history_list

if __name__ == '__main__':
    def get_k_surface(data_2d, temp_range, rel_humidity_range):
        '''
        Takes list of histories generated by sadkat in format [T0[RH0...],T1[RH0...]] and respective ranges of T and RH.
        Returns evaporation rates in transposed form (easier for surface plots).
        '''

        k = []

        for i, T in enumerate(temp_range):

            T_k = []
            for dataset, RH in zip(data_2d[i], rel_humidity_range):

                evaporation_rate = get_evaporation_rate(dataset.time, (dataset.radius * 1e6) ** 2)
                _, evaporation_rate = trim_data(dataset.time, evaporation_rate)

                T_k.append(evaporation_rate.mean())

            k.append(T_k)

        k = np.array(k).T
        return k

    def get_delta_T_surface(data_2d, temp_range, rel_humidity_range):
        '''
        Takes list of histories generated by sadkat in format [T0[RH0...],T1[RH0...]] and respective ranges of T and RH.
        Returns evaporation rates in transposed form (easier for surface plots).
        '''

        T_supp = []

        for i, T in enumerate(temp_range):

            T_delta_T = []
            for dataset, RH in zip(data_2d[i], rel_humidity_range):

                delta_T = get_temperature_suppression(dataset.temperature, T)
                _, delta_T = trim_data(dataset.time, delta_T)

                T_delta_T.append(delta_T.mean())
            T_supp.append(T_delta_T)

        T_supp = np.array(T_supp).T
        return T_supp

if __name__ == '__main__':
    sadkat_2d_evaporation_rates = get_k_surface(sadkat_2d_data,T_range, RH_range)
    sadkat_2d_delta_T = get_delta_T_surface(sadkat_2d_data,T_range, RH_range)
    sadkat_experimental_evaporation_rates = get_k_surface(sadkat_experiemtal_comparison_data, experiment_T_range, RH_range)
    sadkat_experimental_delta_T = get_delta_T_surface(sadkat_experiemtal_comparison_data, experiment_T_range, RH_range)

# # 4.4.4. EDB Comparison
#

if __name__ == '__main__':
    import matplotlib.colors as colors

    '''
    Source: https://www.wolframalpha.com/input/?i=a+%2B+b*x+%2B+c*y+%2B++d*x*y+%2B++e*x%5E2+%2B+f*y%5E2+%3D+0+solve+for+y

    '''

    def RH_sol1(T, k, a, b, c, d, e, f):
        return - (np.sqrt((c+d*T)**2 - 4*f*((a-k) + T*(b + e*T))) + c + d*T) / (2*f)

    def RH_sol2(T, k, a, b, c, d, e, f):
        return (np.sqrt((c+d*T)**2 - 4*f*((a-k) + T*(b + e*T))) + c + d*T) / (2*f)

    def RH_sol3(T, k, a, b, c, d, e, f):
        return - ((a-k) + T*(b + e*T)) / (c + d*T)


    '''
    Source: https://www.wolframalpha.com/input/?i=a+%2B+b*x+%2B+c*y+%2B++d*x*y+%2B++e*x%5E2+%2B+f*y%5E2+%3D+0+solve+for+x
    '''

    def T_sol1(T, k, a, b, c, d, e, f):
        return - (np.sqrt((b+d*T)**2 - 4*e*((a-k) + T*(c + f*T))) + b + d*T) / (2*e)

    def T_sol2(T, k, a, b, c, d, e, f):
        return (np.sqrt((b+d*T)**2 - 4*e*((a-k) + T*(c + f*T))) + b + d*T) / (2*e)

    def T_sol3(T, k, a, b, c, d, e, f):
        return - ((a-k) + T*(c + f*T)) / (b + d*T)

    def get_constant_k_line(Temps, k, params):
        '''Uses inverse function of surface fit, RH_sol1, to find the values of RH corresponding to given T and evaporation rate'''
        return Temps[((RH_sol1(Temps, k, *params)) >= 0) & ((RH_sol1(Temps, k, *params)) <= 1)], RH_sol1(Temps, k, *params)[((RH_sol1(Temps, k, *params)) >= 0) & ((RH_sol1(Temps, k, *params)) <= 1)]

    def fit_evaporation_surface(T_range, RH_range, k_surface):
        '''
        using curve fit to fit to the values of evaporation rate from sadkat over temp and rh

        code taken from:
        https://gist.github.com/silgon/24b56f8ae857ff4ab397
        https://scipython.com/blog/non-linear-least-squares-fitting-of-a-two-dimensional-data/
        '''


        def _poly_func(x, a, b, c, d, e, f):
            '''
            A 2nd order binary polynomial surface function for fitting

            ''' 
            return a + b*x[0] + c*x[1] + d*x[0]*x[1]+ e*x[0]**2 + f*x[1]**2

        # using curve fit to get polynomial surface fit
        X,Y = np.meshgrid(T_range,RH_range)
        size = X.shape
        x1_1d = X.reshape((1, np.prod(size)))
        x2_1d = Y.reshape((1, np.prod(size)))
        xdata = np.vstack((x1_1d, x2_1d)) #creating 1d x data for curvefit

        ydata = k_surface.flatten() # using evaporation rates as ydata

        popt, pcov = curve_fit(_poly_func, xdata, ydata,)
        z_fit = _poly_func(xdata, *popt)
        Z_fit = z_fit.reshape(size)


        return X, Y, Z_fit, popt, np.sqrt(np.diag(pcov))
    
    def myround(x, base=5):
        return base * round(float(x)/base)

    def myround_down(x, base=5):
        return base * (round(float(x)/base) - 1)
    def myround_up(x, base=5):
        return base * (round(float(x)/base) + 1)

    def k_surface_plot(T, RH, k_data,
                       fit_k_surface = False,
                       k_ref_value = None,
                       experiment_T = None, experiment_data = None, sadkat_experiment_data = None,
                       delta_T_data = None):
    
        Rh_convert_to_percent = 100

        T_mesh, RH_mesh, k_fit_surface_data, k_surface_parameters, k_surface_parameters_errors = fit_evaporation_surface(T, Rh_convert_to_percent * RH, k_data)

        #angles for viewing 3d plot
        alpha_ang = 30
        beta_ang = 135

        import matplotlib as mpl
        from mpl_toolkits.mplot3d import Axes3D

        # some markers and styles for use in plotting
        marker_list = ['x', 'o', 'v', '^', 's', 'P', 'p', '*', 'D']
        ls_list = [':','-.','--','-']

        colour_ticks = 25
        cmap = mpl.cm.jet
        boundaries = np.arange(myround_down(k_data.T.min(), colour_ticks),myround_up(k_data.T.max(), colour_ticks),colour_ticks)
        norm = mpl.colors.BoundaryNorm(boundaries=boundaries, ncolors=256)

        #3d plot of evaporation rate over temperature and rh
        fig = plt.figure(figsize=(21,13), dpi = resolution)
        ax = fig.gca(projection='3d')

        #surface of data
        ax.plot_surface(T_mesh, RH_mesh, k_data, alpha = 0.2, cmap = cmap, zorder = 0)
        #white contours to show levels
        ax.contour(T_mesh, RH_mesh, k_data, linewidths = 2, colors = 'w', levels = boundaries)
        #ax.contour(X, Y, Z, zdir='z', offset=-100, cmap=cm.coolwarm)
        #ax.contour(T_mesh, RH_mesh, k_data, zdir = 'x', offset = T.max(), linewidths = 2, cmap = cmap,)# levels = boundaries)
        #ax.contour(T_mesh, RH_mesh, k_data, zdir = 'y', offset = 0, linewidths = 2, cmap = cmap,)# levels = boundaries)
        #ax.contour(T_mesh, RH_mesh, k_data, zdir = 'z', offset = -10, linewidths = 2, cmap = cmap, levels = boundaries)
        #scatter data
        #for i, Temp in enumerate(T):
        #    ax.scatter(Temp, RH, k_data.T[i], alpha = 0.9, c = k_data.T[i], cmap = cmap, norm=norm)


        #ax.contour(T_mesh, RH_mesh, k_data, linewidths = 2, levels = np.arange(myround_down(k_data.T.min(), colour_ticks),myround_up(k_data.T.max(), colour_ticks),colour_ticks), cmap = cmap, norm = norm, zorder = 0)

        #fitted surface
        if fit_k_surface != False:
            ax.plot_wireframe(T_mesh, RH_mesh, k_fit_surface_data, alpha = 1, color = 'k', lw = 0.5)

        color_Ts = ['#FE6100', '#DC267F', '#004D40', '#648FFF']
        color_Ts = dict((el,col) for el, col in zip(experiment_T,color_Ts))

        if isinstance(experiment_T, np.ndarray):
                if isinstance(sadkat_experiment_data, np.ndarray):

                    #experiemental temperatures, model data
                    for i, (temp, style) in enumerate(zip(experiment_T, ls_list)):
                        ax.plot(np.full_like(RH,temp), Rh_convert_to_percent * RH, sadkat_experiment_data.T[i],
                                color = color_Ts[temp],
                                ls = '-', alpha = 1, lw = 2, label = str(round(float(temp))) + ' K', zorder = 1)

        if isinstance(experiment_T, np.ndarray):
                if isinstance(experiment_data, dict):    

                    #plotting experimental data

                    for i, (key) in enumerate(experiment_data):
                        #plotting errors
                        #temperature uncertainty
                        #ax.plot(xs = [experiment_data[key]['exp_props']['T'][0] + T_freezing - 0.5,
                        #              experiment_data[key]['exp_props']['T'][0] + T_freezing + 0.5],
                        #            ys = [experiment_data[key]['exp_props']['RH'].mean(),
                        #                  experiment_data[key]['exp_props']['RH'].mean()],
                        #            zs = [(experiment_data[key]['exp_props']['k'].mean()) * 1e12,
                        #                  (experiment_data[key]['exp_props']['k'].mean()) * 1e12],
                        #        color = 'k')
                        # RH std dev
                        #ax.plot(xs = [experiment_data[key]['exp_props']['T'][0] + T_freezing,
                        #              experiment_data[key]['exp_props']['T'][0] + T_freezing],
                        #            ys = [experiment_data[key]['exp_props']['RH'].mean() - experiment_data[key]['exp_props']['RH'].std(),
                        #                  experiment_data[key]['exp_props']['RH'].mean() + experiment_data[key]['exp_props']['RH'].std()],
                        #            zs = [(experiment_data[key]['exp_props']['k'].mean()) * 1e12,
                        #                  (experiment_data[key]['exp_props']['k'].mean()) * 1e12],
                        #        color = 'k')
                        # k std dev
                        #ax.plot(xs = [experiment_data[key]['exp_props']['T'][0] + T_freezing,
                        #              experiment_data[key]['exp_props']['T'][0] + T_freezing],
                        #            ys = [experiment_data[key]['exp_props']['RH'].mean(),
                        #                  experiment_data[key]['exp_props']['RH'].mean()],
                        #            zs = [(experiment_data[key]['exp_props']['k'].mean() - experiment_data[key]['exp_props']['k'].std()) * 1e12,
                        #                  (experiment_data[key]['exp_props']['k'].mean() + experiment_data[key]['exp_props']['k'].std()) * 1e12],
                        #        color = 'k')
                        #plotting data
                        ax.scatter(xs = experiment_data[key]['exp_props']['T'][0] + T_freezing,
                                    ys = experiment_data[key]['exp_props']['RH'].mean() * Rh_convert_to_percent,
                                    zs = experiment_data[key]['exp_props']['k'].mean() * 1e12,
                                    color = color_Ts[experiment_data[key]['exp_props']['T'][0] + T_freezing] ,marker = 's', s = 100,
                                    label = str(round(experiment_data[key]['exp_props']['T'][0] + T_freezing,3)) + ' K, ' + str(round(100* experiment_data[key]['exp_props']['RH'].mean(),1)) + ' % RH')

        if k_ref_value != None:

            #k line for given evaporation rate:
            ax.plot3D(*get_constant_k_line(T, k_ref_value, k_surface_parameters), zs = k_ref_value, lw = 5, color = 'k', label = str(round(k_ref_value,3)) + ' µm$^2$s$^{-1}$')

        ax.view_init(alpha_ang, beta_ang)

        ax.set_ylim(Rh_convert_to_percent * RH.min(),Rh_convert_to_percent * RH.max())
        ax.set_zlim(k_data.min())

        ax.set_xlabel('T / K', labelpad = 25)
        ax.set_ylabel('RH / %', labelpad = 25)
        ax.set_zlabel(r'Evaporation Rate / µm$^2$s$^{-1}$', labelpad = 25)

        #plt.legend(bbox_to_anchor = (1.5,1))
        #plt.colorbar(mpl.cm.ScalarMappable(norm, cmap), shrink = 0.5, pad = 0.1, alpha = 0.3)

        plt.show()

        ### Contour plot

        fig, ax = plt.subplots(figsize = figure_size, dpi = resolution)
        colour_ticks = 25
        cmap = mpl.cm.jet
        boundaries = np.arange(myround_down(k_data.T.min(), colour_ticks),myround_up(k_data.T.max() + colour_ticks, colour_ticks),colour_ticks)
        norm = mpl.colors.BoundaryNorm(boundaries=boundaries, ncolors=256)
        ax.contourf(T, RH * 100, k_data, levels = boundaries, cmap = cmap, norm = norm)
        CS = ax.contour(T, RH * 100, k_data, linewidths = 2, linestyles = '-', levels = boundaries, colors = 'k')
        ax.clabel(CS, levels = boundaries, inline = True, inline_spacing = 50, fmt = '%1.0f')
        if isinstance(experiment_T, np.ndarray):
                if isinstance(experiment_data, dict):    

                    #plotting experimental data
                    for i, key, in enumerate(experiment_data):
                        ax.scatter(x = experiment_data[key]['exp_props']['T'][0] + T_freezing,
                                   y = experiment_data[key]['exp_props']['RH'].mean() * 100,
                                   c = 'w',
                                   cmap = cmap, norm = norm,
                                   marker = 's', s = 130)
                        ax.scatter(x = experiment_data[key]['exp_props']['T'][0] + T_freezing,
                                   y = experiment_data[key]['exp_props']['RH'].mean() * 100,
                                   c = experiment_data[key]['exp_props']['k'].mean() * 1e12,
                                   cmap = cmap, norm = norm,
                                   marker = 's', s = 100)

        plt.colorbar(mpl.cm.ScalarMappable(norm,cmap), label = r'Evaporation Rate / µm$^2$s$^{-1}$')
        ax.set_xlabel('T / K')
        ax.set_ylabel('RH / %')
        plt.show()

        ### Temperature supression

        if isinstance(delta_T_data, np.ndarray):

            cmap_RedBlue = colors.LinearSegmentedColormap.from_list("RedBlue",['r', 'b'])
            temp_colour_ticks = 5
            temp_boundaries = np.arange(delta_T_data.T.min(),delta_T_data.T.max(),temp_colour_ticks)
            temp_norm = mpl.colors.BoundaryNorm(boundaries=temp_boundaries, ncolors=256)

            #3d plot of evaporation rate over temperature and rh
            fig = plt.figure(figsize=(21,13), dpi = resolution)
            ax = fig.gca(projection='3d')

            #white contours to show levels
            ax.contour(T_mesh, RH_mesh, delta_T_data, linewidths = 2, colors = 'w', levels = temp_boundaries)
            #surface of data
            ax.plot_surface(T_mesh, RH_mesh, delta_T_data, alpha = 0.3, cmap = cmap_RedBlue)

            #scatter data
            #for i, Temp in enumerate(T):
            #    ax.scatter(Temp, RH, delta_T_data.T[i], alpha = 0.9, c = delta_T_data.T[i], cmap = cmap_RedBlue, norm = temp_norm)

            ax.view_init(alpha_ang, beta_ang)

            ax.set_ylim(Rh_convert_to_percent * RH.min(),Rh_convert_to_percent * RH.max())

            ax.set_xlabel('T / K', labelpad = 25)
            ax.set_ylabel('RH / %', labelpad = 25)
            ax.set_zlabel(r'Temperature Suppression / K', labelpad = 25)

            #plt.legend(bbox_to_anchor = (1.5,1))

            plt.show()
        #fitted surface
        if fit_k_surface != False:

            residuals = k_data - k_fit_surface_data

            residuals_max = myround_up(np.abs(residuals).max())

            fig, ax = plt.subplots(figsize = figure_size, dpi = resolution)

            colour_ticks = 50
            cmap = mpl.cm.jet
            boundaries = np.arange(myround_down(k_data.T.min(), colour_ticks),myround_up(k_data.T.max(), colour_ticks),colour_ticks)
            norm = mpl.colors.BoundaryNorm(boundaries=boundaries, ncolors=256)

            CS = ax.contour(T, RH, k_data, linewidths = 0.5, linestyles = '-', levels = boundaries, colors = 'k')
            ax.clabel(CS, levels = boundaries, inline = True, inline_spacing = 50, fmt = '%1.0f')
            ax.contour(T, RH, k_fit_surface_data, linewidths = 0.5, linestyles = ':', levels = boundaries, colors = 'k')


            colour_ticks = 9
            cmap = mpl.cm.seismic
            boundaries_resid = np.linspace(-residuals_max,residuals_max,colour_ticks,)
            norm = mpl.colors.BoundaryNorm(boundaries=boundaries_resid, ncolors=256)

            ax.contourf(T, RH, k_data - k_fit_surface_data, levels = boundaries_resid, cmap = cmap)
            ax.set_ylabel('RH')
            ax.set_xlabel('Temperature / K')
            plt.colorbar(mpl.cm.ScalarMappable(norm,cmap), label = 'Evaporation Rate, Simulation - Fit / µm$^2$s$^{-1}$')
            plt.show()


    return

if __name__ == '__main__':
    k_surface_plot(sadkat_T_range, sadkat_RH_range,
                   sadkat_2d_evaporation_rates,
                   True,
                   k_ref_value = 200,
                   experiment_T = experiment_T_range,
                   experiment_data = experiments,
                   sadkat_experiment_data = sadkat_experimental_evaporation_rates,
                   delta_T_data = sadkat_2d_delta_T)

if __name__ == '__main__':
    #sample data available from https://doi.org/10.5281/zenodo.6396580

    #sample_2d_data = np.load('src/sample_2d_data.npy',allow_pickle = True)
    #sample_RH_range = np.load('src/sample_RH_range.npy',allow_pickle = True)
    #sample_T_range = np.load('src/sample_T_range.npy',allow_pickle = True)
    #sample_2d_evaporation_rates = get_k_surface(sample_2d_data,sample_T_range, sample_RH_range)
    #sample_2d_delta_T = get_delta_T_surface(sample_2d_data,sample_T_range, sample_RH_range)

    #k_surface_plot(sample_T_range, sample_RH_range,
    #               sample_2d_evaporation_rates,
    #               False,
    #               k_ref_value = None,
    #               experiment_T = None,
    #               experiment_data = None,
    #               sadkat_experiment_data = None,
    #               delta_T_data = sample_2d_delta_T)

if __name__ == '__main__':
    from mpl_toolkits.axes_grid1.inset_locator import inset_axes

    def EDB_comparison(experiment_number, T_err = 1, RH_err = None, sample_number = 0, solution = aqueous_NaCl, mfs = 0, time = 25):


        RH = experiments[experiment_number]['exp_props']['RH'].mean()
        if isinstance(RH_err, float):
            pass
        else:
            RH_err = experiments[experiment_number]['exp_props']['RH'].std()
            if RH_err < 1e-6:
                RH_err = 0.025

        RH_upper = RH + RH_err
        RH_lower = RH - RH_err

        T = experiments[experiment_number]['exp_props']['T'].mean() + T_freezing
        T_upper = T + T_err
        T_lower = T - T_err

        R0s= []
        t0s = []

        for exp in experiments[experiment_number]['exp_props'].droplet_path:
            R0s.append(pd.read_csv(exp, sep = '\t', header = None, names=['time_s', 'radius_um', 'RI', None]).radius_um[0])
            t0s.append(pd.read_csv(exp, sep = '\t', header = None, names=['time_s', 'radius_um', 'RI', None]).time_s[0])

        R0s = np.array(R0s)
        t0s = np.array(t0s)

        R0 = R0s.mean()/1e6 # metres
        t0 = t0s.mean() # seconds

        R0_err = R0s.std()/1e6 # metres
        t0_err = t0.std() # seconds

        print('RH: ', RH, ' ± ', RH_err)
        print('T: ', T, ' ± ', T_err)
        print('R0: ', R0, ' ± ', R0_err)
        print('t0: ', t0, ' ± ', t0_err)

        print('simulating experiement', end = '\r')
        specific_droplet, specific_trajectory = simulate(time, solution, T, RH, R0, T, mfs)
        result = specific_droplet.complete_trajectory(specific_trajectory)

        print('simulating RH range', end = '\r')
        specific_droplet, specific_trajectory = simulate(time, solution, T, RH_upper, R0, T, mfs)
        result_posRHerr = specific_droplet.complete_trajectory(specific_trajectory)
        specific_droplet, specific_trajectory = simulate(time, solution, T, RH_lower, R0, T, mfs)
        result_negRHerr = specific_droplet.complete_trajectory(specific_trajectory)

        print('simulating T range', end = '\r')
        specific_droplet, specific_trajectory = simulate(time, solution, T_upper, RH, R0, T_upper, mfs)
        result_posTerr = specific_droplet.complete_trajectory(specific_trajectory)
        specific_droplet, specific_trajectory = simulate(time, solution, T_lower, RH, R0, T_lower, mfs)
        result_negTerr = specific_droplet.complete_trajectory(specific_trajectory)


        fig, ax1 = plt.subplots(dpi = resolution)

        #plotting RH error range
        ax1.fill(np.append(result_posRHerr.time + t0,
                           result_negRHerr.time[::-1] + t0),
                 np.append(result_posRHerr.radius/1e-6,
                           result_negRHerr.radius[::-1]/1e-6),
                 alpha = 1, color = '#FFC107', zorder = 0, label = str(round(100 * RH,3)) + '± ' + str(round(100 * RH_err,3)) + ' % RH')

        #plotting T error range
        ax1.fill(np.append(result_posTerr.time + t0,
                           result_negTerr.time[::-1] + t0),
                 np.append(result_posTerr.radius/1e-6,
                           result_negTerr.radius[::-1]/1e-6),
                 alpha = 0.5, color = '#1E88E5', zorder = 0, label = str(round(T,3)) + '± ' + str(round(T_err,3)) + ' K')

        #Plotting results
        ax1.plot(result['time'] + t0,
                 result['radius'] / 1e-6,
                 ':', lw = 5, color = '#D81B60', label = str(round(100 * RH,3)) + ' % RH, 298 K')

        #scattering EDB data
        df = pd.DataFrame(columns=['time_s', 'radius_um', 'RI', None])
        for exp in experiments[experiment_number]['exp_props'].droplet_path:
            #df = pd.read_csv(exp, sep = '\t', header = None, names=['time_s', 'radius_um', 'RI', None])
            df = pd.concat([df , pd.read_csv(exp, sep = '\t', header = None, names=['time_s', 'radius_um', 'RI', None])])

        ax1.scatter(df.time_s,
                    df.radius_um,
                    s = 1, color = 'k', zorder = 1, label = 'All EDB data')

        ax1.set_xlabel('Time / s')
        ax1.set_ylabel('Radius / µm')
        ax1.set_xlim(0)

        plt.legend(loc='best')


        ax2 = inset_axes(ax1,
                         width="100%", height="100%",
                         bbox_to_anchor=(.1, .1, .6, .5),bbox_transform=ax1.transAxes, loc=3)#fig.add_axes([0.18, 0.2, 0.3, 0.3],)
        ax2.patch.set_alpha(0)

        #plotting RH error range
        ax2.fill(np.append(result_posRHerr.time + t0,
                           result_negRHerr.time[::-1] + t0),
                 np.append(result_posRHerr.radius/1e-6,
                           result_negRHerr.radius[::-1]/1e-6),
                 alpha = 1, color = '#FFC107', zorder = 0, label = str(round(100 * RH,3)) + '± ' + str(round(100 * RH_err,3)) + ' % RH')

        #plotting T error range
        ax2.fill(np.append(result_posTerr.time + t0,
                           result_negTerr.time[::-1] + t0),
                 np.append(result_posTerr.radius/1e-6,
                           result_negTerr.radius[::-1]/1e-6),
                 alpha = 0.5, color = '#1E88E5', zorder = 0, label = str(round(T,3)) + '± ' + str(round(T_err,3)) + ' K')

        #Plotting results
        ax2.plot(result['time'] + t0,
                 result['radius'] / 1e-6,
                 ':', lw = 5, color = '#D81B60', label = str(round(100 * RH,3)) + ' % RH, 298 K')

        #scattering EDB data
        ax2.scatter(df.time_s,
                    df.radius_um,
                    s = 1, color = 'k', zorder = 1, label = 'All EDB data')

        ax2.set_ylim(df.loc[df.time_s <0.5].radius_um.min(),
                     df.loc[df.time_s <0.5].radius_um.max())
        ax2.set_xlim(0,0.5)

        plt.show()
        return

if __name__ == '__main__':
    EDB_comparison('21-12-07_1618', RH_err= 0.025, )
    EDB_comparison('21-12-15_1604', )

# # 4.4.5. FDC comparison

if __name__ == '__main__':
    def aerodyamic_diameter(his, gasflow = 0):
        return his.time[his.reynolds_number < 0.1], np.sqrt(( 18 * 1.81e-5 * ((his.vz[his.reynolds_number < 0.1] - gasflow)) ) / (1000 * 9.81))


if __name__ == '__main__':

    # load data
    NaCl_data = np.load('src/FDC data for benchmarking/NaCl_FDC_data.npy', allow_pickle=True)
    NaCl_data_parameters = np.load('src/FDC data for benchmarking/NaCl_FDC_data_parameters.npy', allow_pickle=True)
    NaCl_data_RHs = np.load('src/FDC data for benchmarking/NaCl_FDC_data_RHs.npy', allow_pickle=True)
    NaCl_data_RHs = NaCl_data_RHs/100

    initial_points = 5

    NaCl_v_x0s = []
    NaCl_v_y0s = []
    NaCl_R0s = []
    NaCl_t_finals = []
    NaCl_gas_flows = []

    for result, RH, param in zip(NaCl_data, NaCl_data_RHs, NaCl_data_parameters):
        NaCl_v_x0s.append((result['all'].head(initial_points).x_velocity_mm_per_s.mean() / 1e3,
                          result['all'].head(initial_points).x_velocity_mm_per_s.std() / 1e3))
        NaCl_v_y0s.append((result['all'].head(initial_points).y_velocity_mm_per_s.mean() / 1e3,
                          result['all'].head(initial_points).y_velocity_mm_per_s.std() / 1e3))
        NaCl_R0s.append((result['all'].head(initial_points).d_e_um_mean.mean() / 2e6,
                        result['all'].head(initial_points).d_e_um_mean.std() / 2e6))
        NaCl_t_finals.append(result['all'].Time_s.max())
        NaCl_gas_flows.append(param['CalculatedGasFlow']*(1/60e6)/(0.02 **2))

# +
if __name__ == '__main__':

    solution = aqueous_NaCl
    T = 294.15 # Kelvin
    mfs = 0.05

    #time = 5 # seconds

    NaCl_simulations = []
    NaCl_simulations_w_eff = []

    for vx, vy, R0, RH, t,gas_flow in zip(NaCl_v_x0s, NaCl_v_y0s, NaCl_R0s, NaCl_data_RHs, NaCl_t_finals, NaCl_gas_flows):

        droplet, trajectory = simulate(t, solution, T, RH, R0[0], T, mfs, initial_position = np.array([0.1/1e3, 0, 1/1e3]), initial_velocity = np.array([vx[0], 0, vy[0]]), gravity = np.array([0, 0, 9.80665]), gas_velocity = np.array([0,0,gas_flow]), terminate_on_equilibration=False)
        # Obtain a table giving a history of *all* droplet parameters.
        history = droplet.complete_trajectory(trajectory)
        NaCl_simulations.append(history)
        droplet, trajectory = simulate(t, solution, T, RH, R0[0], T, mfs, initial_position = np.array([0.1/1e3, 0, 1/1e3]), initial_velocity = np.array([vx[0], 0, vy[0]]), gravity = np.array([0, 0, 9.80665]), gas_velocity = np.array([0,0,gas_flow]), terminate_on_efflorescence=True, eff_threshold=0.238)
        history = droplet.complete_trajectory(trajectory)
        NaCl_simulations_w_eff.append(history)



# -

if __name__ == '__main__':
    import matplotlib as mpl

    def plot_trajectories_and_evaporation(exp_data, exp_RHs, exp_v_x0s, exp_v_y0s, exp_R0s, exp_t_finals, exp_gas_flows, sim_data, sim_data_w_eff):

        #get colourmap and norm 
        bounds = np.power(10.0, np.arange(-4, 2))
        ncolors = len(bounds) -1
        cmap = mpl.cm.get_cmap('plasma_r', ncolors) # Colour map (there are many others)
        norm = mpl.colors.BoundaryNorm(boundaries=bounds, ncolors=ncolors)
        alpha = 0.5
        mapper = mpl.cm.ScalarMappable(norm=norm, cmap= cmap)


        fig, axes = plt.subplots(len(exp_data), 2, sharex='col', sharey=False, figsize = (10,4 * len(exp_data)), dpi = resolution)
        plt.subplots_adjust(hspace=0, wspace=0.25)
        for i, (axis_pair, result, simulation, simulation_eff, R0, RH) in enumerate(zip(axes, exp_data, sim_data, sim_data_w_eff, exp_R0s, exp_RHs)):

            #Trajectory
            axis_pair[0].set_xscale('log')
            axis_pair[0].set_yscale('log')
            axis_pair[0].set_ylim(150,0.3)
            axis_pair[0].set_title(str(100 * RH) + ' % RH', x = 0.5, y = 0.85,fontsize = 20)
            axis_pair[0].set_ylabel('Vertical Position / mm')

            #experiment
            axis_pair[0].scatter(result['all'].x_position_mm_mean, result['all'].y_position_mm_mean,
                                 marker = 's', s = 7,
                                 c = result['all'].Time_s, cmap = cmap, norm = norm, label = 'FDC Measurement')
            axis_pair[0].plot(result['all'].x_position_mm_mean, result['all'].y_position_mm_mean,
                              lw = 0.5, color = 'k', zorder = 0)
            #simulation
            axis_pair[0].scatter(simulation.x / 1e-3, simulation.z / 1e-3,
                                 marker = 'o', s = 7, lw = 1,
                                 c = simulation.time, cmap = cmap, norm = norm, label = 'SADKAT Simulation')
            axis_pair[0].plot(simulation.x / 1e-3, simulation.z / 1e-3,
                              lw = 0.5, color = 'r', zorder = 0)


            axis_pair[1].set_ylim(0,1)

            # evaporation (normalised)
            axis_pair[1].set_ylabel('R$^2$ / R$_0^2$')

            axis_pair[1].errorbar(result['all'].Time_s / (R0[0]/1e-6) ** 2,
                                  ((result['all'].d_e_um_mean/2) / (R0[0]/1e-6)) ** 2,
                                  ((result['all'].d_e_um_mean/2 + result['all'].d_e_um_std/2) / (R0[0]/1e-6)) ** 2 - ((result['all'].d_e_um_mean/2) / (R0[0]/1e-6)) ** 2,
                                  fmt = 'o', ms = 1.5,
                                  elinewidth = 1, capsize = 1, color = 'k',
                                  label = 'FDC Measurement')

            axis_pair[1].plot(simulation[simulation.time <= simulation_eff.tail(1).time.values[0]].time / (R0[0]/1e-6) ** 2,
                              (simulation[simulation.time <= simulation_eff.tail(1).time.values[0]].radius/1e-6) ** 2 / (R0[0]/1e-6) ** 2,
                              lw = 1, c = 'r',
                              label = 'SADKAT Simulation')

            axis_pair[1].plot(simulation[simulation.time > simulation_eff.tail(1).time.values[0]].time / (R0[0]/1e-6) ** 2,
                              (simulation[simulation.time > simulation_eff.tail(1).time.values[0]].radius/1e-6) ** 2 / (R0[0]/1e-6) ** 2,
                              lw = 1, c = 'r', ls = '--')

            if simulation_eff.time.max() < simulation.time.max():

                axis_pair[0].scatter(simulation_eff.tail(1).x / 1e-3,
                                     simulation_eff.tail(1).z / 1e-3,
                                     marker = '*', color = 'r', s = 50, zorder = 3,
                                     label = 'SADKAT Efflorescence')

                axis_pair[1].scatter(simulation_eff.tail(1).time / (R0[0]/1e-6) ** 2,
                                     (simulation_eff.tail(1).radius/1e-6) ** 2 / (R0[0]/1e-6) ** 2,
                                     marker = '*', color = 'r', s = 50, zorder = 3,
                                     label = 'SADKAT Efflorescence')
            if i == 0:
                axis_pair[0].legend(loc = 'lower left', fontsize = 'small')
                axis_pair[1].legend(loc = 'upper right', fontsize = 'small')
                leg = axis_pair[0].get_legend()
                leg.legendHandles[0].set_color('orange')
                leg.legendHandles[1].set_color('orange')
            else:
                yticks = axis_pair[1].yaxis.get_major_ticks()
                yticks[5].label1.set_visible(False)

            if i == len(exp_data) - 1:
                axis_pair[0].set_xlabel('Horizontal Position / mm')

                axis_pair[1].set_xlim(0,0.0075)
                axis_pair[1].set_xlabel('(t / R$_0^2$) / (s / µm$^2$s$^{-1}$)')
        return

if __name__ == '__main__':
    plot_trajectories_and_evaporation(NaCl_data, NaCl_data_RHs, NaCl_v_x0s, NaCl_v_y0s, NaCl_R0s, NaCl_t_finals, NaCl_gas_flows, NaCl_simulations, NaCl_simulations_w_eff)


if __name__ == '__main__':
    def plot_trajectories_and_evaporation_single_graph(exp_data, exp_RHs, exp_v_x0s, exp_v_y0s, exp_R0s, exp_t_finals, exp_gas_flows, sim_data, sim_data_w_eff, colours):

        #get colourmap and norm 
        bounds = np.power(10.0, np.arange(-4, 2))
        ncolors = len(bounds) -1
        cmap = mpl.cm.get_cmap('plasma_r', ncolors) # Colour map (there are many others)
        norm = mpl.colors.BoundaryNorm(boundaries=bounds, ncolors=ncolors)
        alpha = 0.5
        mapper = mpl.cm.ScalarMappable(norm=norm, cmap= cmap)


        fig, axis_pair = plt.subplots(1, 2, sharex='col', sharey=False, figsize = (13,8), dpi = resolution)
        plt.subplots_adjust(hspace=0, wspace=0.25)
        for i, (result, simulation, simulation_eff, R0, RH, colour) in enumerate(zip(exp_data, sim_data, sim_data_w_eff, exp_R0s, exp_RHs, colours)):

            #Trajectory
            axis_pair[0].set_xscale('log')
            axis_pair[0].set_yscale('log')
            axis_pair[0].set_ylim(150,0.3)
            axis_pair[0].set_ylabel('Vertical Position / mm')

            #experiment
            axis_pair[0].scatter(result['all'].x_position_mm_mean, result['all'].y_position_mm_mean,
                                 marker = 's', s = 7,
                                 c = result['all'].Time_s, cmap = cmap, norm = norm, label = 'FDC Data')
            axis_pair[0].plot(result['all'].x_position_mm_mean, result['all'].y_position_mm_mean,
                              lw = 0.5, color = 'k', zorder = 0)
            #simulation
            axis_pair[0].scatter(simulation.x / 1e-3, simulation.z / 1e-3,
                                 marker = 'o', s = 7, lw = 1,
                                 c = simulation.time, cmap = cmap, norm = norm, label = 'SADKAT Simulation')
            axis_pair[0].plot(simulation.x / 1e-3, simulation.z / 1e-3,
                              lw = 0.5, color = 'r', zorder = 0)


            axis_pair[1].set_ylim(0,1)
            axis_pair[1].set_xlim(0,0.0075)

            # evaporation (normalised)
            axis_pair[1].set_ylabel('R$^2$ / R$_0^2$')

            axis_pair[1].errorbar(result['all'].Time_s / (R0[0]/1e-6) ** 2,
                                  ((result['all'].d_e_um_mean/2) / (R0[0]/1e-6)) ** 2,
                                  ((result['all'].d_e_um_mean/2 + result['all'].d_e_um_std/2) / (R0[0]/1e-6)) ** 2 - ((result['all'].d_e_um_mean/2) / (R0[0]/1e-6)) ** 2,
                                  fmt = 'o', ms = 1.5,
                                  elinewidth = 1, capsize = 1, color = colour)

            axis_pair[1].plot(simulation[simulation.time <= simulation_eff.tail(1).time.values[0]].time / (R0[0]/1e-6) ** 2,
                              (simulation[simulation.time <= simulation_eff.tail(1).time.values[0]].radius/1e-6) ** 2 / (R0[0]/1e-6) ** 2,
                              lw = 1, c = colour,
                              label = str(round(RH,2)*100) + ' % RH')

            axis_pair[1].plot(simulation[simulation.time > simulation_eff.tail(1).time.values[0]].time / (R0[0]/1e-6) ** 2,
                              (simulation[simulation.time > simulation_eff.tail(1).time.values[0]].radius/1e-6) ** 2 / (R0[0]/1e-6) ** 2,
                              lw = 1, c = colour, ls = '--')

            if simulation_eff.time.max() < simulation.time.max():

                axis_pair[0].scatter(simulation_eff.tail(1).x / 1e-3,
                                     simulation_eff.tail(1).z / 1e-3,
                                     marker = '*', color = 'k', s = 100, zorder = 3)

                axis_pair[1].scatter(simulation_eff.tail(1).time / (R0[0]/1e-6) ** 2,
                                     (simulation_eff.tail(1).radius/1e-6) ** 2 / (R0[0]/1e-6) ** 2,
                                     marker = '*', color = 'k', s = 100, zorder = 3)
            if i == 0:
                axis_pair[0].legend(loc = 'lower left', fontsize = 'small')
                leg = axis_pair[0].get_legend()
                leg.legendHandles[0].set_color('orange')
                leg.legendHandles[1].set_color('orange')
            else:
                yticks = axis_pair[1].yaxis.get_major_ticks()
                yticks[5].label1.set_visible(False)

            if i == len(exp_data) - 1:
                axis_pair[0].set_xlabel('Horizontal Position / mm')

                axis_pair[1].set_xlabel('(t / R$_0^2$) / (s / µm$^2$s$^{-1}$)')
        axis_pair[1].legend(loc = 'upper right', fontsize = 'small')
        return

if __name__ == '__main__':
    indices = [0,2,3,4,6]
    
    plot_trajectories_and_evaporation_single_graph([NaCl_data[index] for index in indices],
                                                   [NaCl_data_RHs[index] for index in indices],
                                                   [NaCl_v_x0s[index] for index in indices],
                                                   [NaCl_v_y0s[index] for index in indices],
                                                   [NaCl_R0s[index] for index in indices],
                                                   [NaCl_t_finals[index] for index in indices],
                                                   [NaCl_gas_flows[index] for index in indices],
                                                   [NaCl_simulations[index] for index in indices],
                                                   [NaCl_simulations_w_eff[index] for index in indices],
                                                   ['#FFB000',  '#FE6100', '#DC267F', '#785EF0', '#648FFF'])

# ## Comparing Evaporation Rate

if __name__ == '__main__':
    rolling_window = 7
    colours = ['#E69F00',  '#56B4E9', '#009E73', '#F0E442', '#0072B2', '#D55E00', '#CC79A7']

    fig,ax = plt.subplots(dpi = resolution)

    for result, RH, simulation, R0, colour in zip(NaCl_data, NaCl_data_RHs, NaCl_simulations, NaCl_R0s, colours):
        data = simulation[simulation.time >= 0.1]
        line, = ax.plot(data.time,
                        - np.gradient((data.radius / 1e-6) ** 2, data.time), ls = '--', lw = 3,
                        label = str(round(100 * RH)) + ' % RH', color = colour)

        exp_data = result['all'].rolling(rolling_window, center = True,).mean()
        exp_err = result['all'].rolling(rolling_window, center = True).std()
        exp_data = exp_data[exp_data.Time_s > 0.1]
        exp_err = exp_err[exp_err.Time_s > 0.1]


        ax.scatter(exp_data.Time_s,
                   - np.gradient((exp_data.d_e_um_mean/2) ** 2, exp_data.Time_s),
                   color = line.get_color(), s = 2, marker = 'o')

        ax.fill_between(exp_data.Time_s,
                        - np.gradient(((exp_data.d_e_um_mean + exp_data.d_e_um_std)/2) ** 2, exp_data.Time_s),
                        - np.gradient(((exp_data.d_e_um_mean - exp_data.d_e_um_std)/2) ** 2, exp_data.Time_s),
                        color = line.get_color(), alpha = 0.1,)



    ax.axhline(0, color = 'k' , lw = 0.5)
    ax.legend()
    ax.set_ylim(-50, 400)
    ax.set_xlabel('Time / s')
    ax.set_ylabel(r'Evaporation Rate / µm$^2$s$^{-1}$')
    plt.show()

    fig,ax = plt.subplots(dpi = resolution)

    for result, RH, simulation, R0, colour in zip(NaCl_data, NaCl_data_RHs, NaCl_simulations, NaCl_R0s, colours):
        data = simulation[simulation.time >= 0.1]
        line, = ax.plot(data.time / R0[0] ** 2 / 1e12,
                        - np.gradient((data.radius / 1e-6) ** 2, data.time), ls = '--', lw = 3,
                        label = str(round(100 * RH)) + ' % RH', color = colour)

        exp_data = result['all'].rolling(rolling_window, center = True,).mean()
        exp_err = result['all'].rolling(rolling_window, center = True).std()
        exp_data = exp_data[exp_data.Time_s > 0.1]
        exp_err = exp_err[exp_err.Time_s > 0.1]


        ax.scatter(exp_data.Time_s / R0[0] ** 2 / 1e12,
                   - np.gradient((exp_data.d_e_um_mean/2) ** 2, exp_data.Time_s),
                   color = line.get_color(), s = 2, marker = 'o')

        ax.fill_between(exp_data.Time_s / R0[0] ** 2 / 1e12,
                        - np.gradient(((exp_data.d_e_um_mean + exp_data.d_e_um_std)/2) ** 2, exp_data.Time_s),
                        - np.gradient(((exp_data.d_e_um_mean - exp_data.d_e_um_std)/2) ** 2, exp_data.Time_s),
                        color = line.get_color(), alpha = 0.1,)



    ax.axhline(0, color = 'k' , lw = 0.5)
    ax.legend()
    ax.set_ylim(-50, 400)
    ax.set_xlim(0, 0.008)
    ax.set_xlabel('(Time / R$_0^2$) / (s / µm$^2$)')
    ax.set_ylabel(r'Evaporation Rate / µm$^2$s$^{-1}$')
    plt.show()

# # 4.5.1. Sensitivity Analysis

if __name__ == '__main__':
    NaCl_10pctRH_data = NaCl_data[2]
    NaCl_10pctRH_v_x0 = NaCl_v_x0s[2]
    NaCl_10pctRH_v_y0 = NaCl_v_y0s[2]
    NaCl_10pctRH_R_0 = NaCl_R0s[2]
    NaCl_10pctRH_t = NaCl_t_finals[2]
    NaCl_10pctRH_gas_flow = NaCl_gas_flows[2]

if __name__ == '__main__':

    def plot_NaCl_sensitivity(his, his_pos, his_neg, df, chart_label):
        fig, (ax0, ax1) = plt.subplots(1,2, figsize = figure_size, dpi = 400)
        ax0.plot(his.x/1e-3, his.z/1e-3, 'red', lw = 2, label = 'Simulation')
        ax0.plot(his_pos.x/1e-3, his_pos.z/1e-3, '--', color = 'red', lw = 2, label = '+ve error')
        ax0.plot(his_neg.x/1e-3, his_neg.z/1e-3, ':', color = 'red', lw = 2, label = '-ve error')
        ax0.errorbar(df.x_position_mm_mean,
                     df.y_position_mm_mean,
                     df.x_position_mm_std,
                     df.y_position_mm_std,
                     fmt = 'o', color = 'k', elinewidth = 1, capsize = 2, ms = 3, label = 'FDC data')

        ax0.set_xscale('log')
        ax0.set_yscale('log')
        ax0.set_ylim(100,0.9)
        ax0.set_ylabel('Vertical Position / mm')
        ax0.set_xlabel('Horizontal Position / mm')

        ax1.plot(his.time, his.radius / 1e-6, 'red', lw = 2, label = 'Simulation')
        ax1.plot(his_pos.time, his_pos.radius / 1e-6, '--', color = 'red', lw = 2, label = '+ve error')
        ax1.plot(his_neg.time, his_neg.radius / 1e-6, ':', color = 'red', lw = 2, label = '-ve error')
        ax1.errorbar(df.Time_s,
                     df.d_e_um_mean / 2,
                     df.d_e_um_std / 2,
                     fmt = 'o', color = 'k', elinewidth = 1, capsize = 2, ms = 3, label = 'FDC data')

        ax1.set_ylabel('Radius / µm')
        ax1.set_xlabel('Time / s')
        ax1.legend()

        fig.suptitle(chart_label, fontsize = 20)
        plt.show()

if __name__ == '__main__':
    T = T_freezing + 21
    mfs = 0.05
    RH = 0.1

    droplet, trajectory = simulate(NaCl_10pctRH_t,
                                   solution,
                                   T,
                                   RH,
                                   NaCl_10pctRH_R_0[0],
                                   T,
                                   mfs,
                                   initial_position = np.array([0.1/1e3, 0, 1/1e3]),
                                   initial_velocity = np.array([NaCl_10pctRH_v_x0[0],
                                                                0,
                                                                NaCl_10pctRH_v_y0[0]]),
                                   gravity = np.array([0, 0, 9.80665]),
                                   gas_velocity = np.array([0,0,NaCl_10pctRH_gas_flow]),
                                   terminate_on_equilibration=False)
    # Obtain a table giving a history of *all* droplet parameters.
    history = droplet.complete_trajectory(trajectory)

    droplet, trajectory = simulate(NaCl_10pctRH_t,
                                   solution,
                                   T,
                                   RH,
                                   NaCl_10pctRH_R_0[0] + NaCl_10pctRH_R_0[1],
                                   T,
                                   mfs,
                                   initial_position = np.array([0.1/1e3, 0, 1/1e3]),
                                   initial_velocity = np.array([NaCl_10pctRH_v_x0[0],
                                                                0,
                                                                NaCl_10pctRH_v_y0[0]]),
                                   gravity = np.array([0, 0, 9.80665]),
                                   gas_velocity = np.array([0,0,NaCl_10pctRH_gas_flow]),
                                   terminate_on_equilibration=False)
    # Obtain a table giving a history of *all* droplet parameters.
    history_pos_R0 = droplet.complete_trajectory(trajectory)

    droplet, trajectory = simulate(NaCl_10pctRH_t,
                                   solution,
                                   T,
                                   RH,
                                   NaCl_10pctRH_R_0[0] - NaCl_10pctRH_R_0[1],
                                   T,
                                   mfs,
                                   initial_position = np.array([0.1/1e3, 0, 1/1e3]),
                                   initial_velocity = np.array([NaCl_10pctRH_v_x0[0],
                                                                0,
                                                                NaCl_10pctRH_v_y0[0]]),
                                   gravity = np.array([0, 0, 9.80665]),
                                   gas_velocity = np.array([0,0,NaCl_10pctRH_gas_flow]),
                                   terminate_on_equilibration=False)
    # Obtain a table giving a history of *all* droplet parameters.
    history_neg_R0 = droplet.complete_trajectory(trajectory)

    #plot_NaCl_sensitivity(history, history_pos_R0, history_neg_R0, NaCl_10pctRH_data['all'], 'Radius: 18.6 µm ± 0.4 µm')

    ################################

    droplet, trajectory = simulate(NaCl_10pctRH_t,
                                   solution,
                                   T,
                                   RH,
                                   NaCl_10pctRH_R_0[0],
                                   T,
                                   mfs,
                                   initial_position = np.array([0.1/1e3, 0, 1/1e3]),
                                   initial_velocity = np.array([NaCl_10pctRH_v_x0[0] + NaCl_10pctRH_v_x0[1],
                                                                0,
                                                                NaCl_10pctRH_v_y0[0]]),
                                   gravity = np.array([0, 0, 9.80665]),
                                   gas_velocity = np.array([0,0,NaCl_10pctRH_gas_flow]),
                                   terminate_on_equilibration=False)
    # Obtain a table giving a history of *all* droplet parameters.
    history_pos_Vx = droplet.complete_trajectory(trajectory)

    droplet, trajectory = simulate(NaCl_10pctRH_t,
                                   solution,
                                   T,
                                   RH,
                                   NaCl_10pctRH_R_0[0],
                                   T,
                                   mfs,
                                   initial_position = np.array([0.1/1e3, 0, 1/1e3]),
                                   initial_velocity = np.array([NaCl_10pctRH_v_x0[0] - NaCl_10pctRH_v_x0[1],
                                                                0,
                                                                NaCl_10pctRH_v_y0[0]]),
                                   gravity = np.array([0, 0, 9.80665]),
                                   gas_velocity = np.array([0,0,NaCl_10pctRH_gas_flow]),
                                   terminate_on_equilibration=False)
    # Obtain a table giving a history of *all* droplet parameters.
    history_neg_Vx = droplet.complete_trajectory(trajectory)

    #plot_NaCl_sensitivity(history, history_pos_Vx, history_neg_Vx, NaCl_10pctRH_data['all'], 'Horizontal Velocity: 1.232 m/s ± 0.153 m/s')

    ################################

    droplet, trajectory = simulate(NaCl_10pctRH_t,
                                   solution,
                                   T,
                                   RH,
                                   NaCl_10pctRH_R_0[0],
                                   T,
                                   mfs,
                                   initial_position = np.array([0.1/1e3, 0, 1/1e3]),
                                   initial_velocity = np.array([NaCl_10pctRH_v_x0[0],
                                                                0,
                                                                NaCl_10pctRH_v_y0[0] + NaCl_10pctRH_v_y0[1]]),
                                   gravity = np.array([0, 0, 9.80665]),
                                   gas_velocity = np.array([0,0,NaCl_10pctRH_gas_flow]),
                                   terminate_on_equilibration=False)
    # Obtain a table giving a history of *all* droplet parameters.
    history_pos_Vy = droplet.complete_trajectory(trajectory)

    droplet, trajectory = simulate(NaCl_10pctRH_t,
                                   solution,
                                   T,
                                   RH,
                                   NaCl_10pctRH_R_0[0],
                                   T,
                                   mfs,
                                   initial_position = np.array([0.1/1e3, 0, 1/1e3]),
                                   initial_velocity = np.array([NaCl_10pctRH_v_x0[0],
                                                                0,
                                                                NaCl_10pctRH_v_y0[0] - NaCl_10pctRH_v_y0[1]]),
                                   gravity = np.array([0, 0, 9.80665]),
                                   gas_velocity = np.array([0,0,NaCl_10pctRH_gas_flow]),
                                   terminate_on_equilibration=False)
    # Obtain a table giving a history of *all* droplet parameters.
    history_neg_Vy = droplet.complete_trajectory(trajectory)

    #plot_NaCl_sensitivity(history, history_pos_Vy, history_neg_Vy, NaCl_10pctRH_data['all'], 'Vertical Velocity: 0.248 m/s ± 0.089 m/s')

    ################################

    droplet, trajectory = simulate(NaCl_10pctRH_t,
                                   solution,
                                   T,
                                   RH + 0.05,
                                   NaCl_10pctRH_R_0[0],
                                   T,
                                   mfs,
                                   initial_position = np.array([0.1/1e3, 0, 1/1e3]),
                                   initial_velocity = np.array([NaCl_10pctRH_v_x0[0],
                                                                0,
                                                                NaCl_10pctRH_v_y0[0]]),
                                   gravity = np.array([0, 0, 9.80665]),
                                   gas_velocity = np.array([0,0,NaCl_10pctRH_gas_flow]),
                                   terminate_on_equilibration=False)
    # Obtain a table giving a history of *all* droplet parameters.
    history_pos_RH = droplet.complete_trajectory(trajectory)

    droplet, trajectory = simulate(NaCl_10pctRH_t,
                                   solution,
                                   T,
                                   RH - 0.05,
                                   NaCl_10pctRH_R_0[0],
                                   T,
                                   mfs,
                                   initial_position = np.array([0.1/1e3, 0, 1/1e3]),
                                   initial_velocity = np.array([NaCl_10pctRH_v_x0[0],
                                                                0,
                                                                NaCl_10pctRH_v_y0[0]]),
                                   gravity = np.array([0, 0, 9.80665]),
                                   gas_velocity = np.array([0,0,NaCl_10pctRH_gas_flow]),
                                   terminate_on_equilibration=False)
    # Obtain a table giving a history of *all* droplet parameters.
    history_neg_RH = droplet.complete_trajectory(trajectory)

    #plot_NaCl_sensitivity(history, history_pos_RH, history_neg_RH, NaCl_10pctRH_data['all'], 'RH: 10 % ± 5 %')
    ################################

    droplet, trajectory = simulate(NaCl_10pctRH_t,
                                   solution,
                                   T + 1,
                                   RH,
                                   NaCl_10pctRH_R_0[0],
                                   T,
                                   mfs,
                                   initial_position = np.array([0.1/1e3, 0, 1/1e3]),
                                   initial_velocity = np.array([NaCl_10pctRH_v_x0[0],
                                                                0,
                                                                NaCl_10pctRH_v_y0[0]]),
                                   gravity = np.array([0, 0, 9.80665]),
                                   gas_velocity = np.array([0,0,NaCl_10pctRH_gas_flow]),
                                   terminate_on_equilibration=False)
    # Obtain a table giving a history of *all* droplet parameters.
    history_pos_T = droplet.complete_trajectory(trajectory)

    droplet, trajectory = simulate(NaCl_10pctRH_t,
                                   solution,
                                   T - 1,
                                   RH,
                                   NaCl_10pctRH_R_0[0],
                                   T,
                                   mfs,
                                   initial_position = np.array([0.1/1e3, 0, 1/1e3]),
                                   initial_velocity = np.array([NaCl_10pctRH_v_x0[0],
                                                                0,
                                                                NaCl_10pctRH_v_y0[0]]),
                                   gravity = np.array([0, 0, 9.80665]),
                                   gas_velocity = np.array([0,0,NaCl_10pctRH_gas_flow]),
                                   terminate_on_equilibration=False)
    # Obtain a table giving a history of *all* droplet parameters.
    history_neg_T = droplet.complete_trajectory(trajectory)

    #plot_NaCl_sensitivity(history, history_pos_T, history_neg_T, NaCl_10pctRH_data['all'], 'Temperature: 319 K ± 1K')
    ################################

    droplet, trajectory = simulate(NaCl_10pctRH_t,
                                   solution,
                                   T,
                                   RH,
                                   NaCl_10pctRH_R_0[0],
                                   T,
                                   mfs,
                                   initial_position = np.array([0.1/1e3, 0, 1/1e3]),
                                   initial_velocity = np.array([NaCl_10pctRH_v_x0[0],
                                                                0,
                                                                NaCl_10pctRH_v_y0[0]]),
                                   gravity = np.array([0, 0, 9.80665]),
                                   gas_velocity = np.array([0,0,NaCl_10pctRH_gas_flow + 1e-3]),
                                   terminate_on_equilibration=False)
    # Obtain a table giving a history of *all* droplet parameters.
    history_pos_gf = droplet.complete_trajectory(trajectory)

    droplet, trajectory = simulate(NaCl_10pctRH_t,
                                   solution,
                                   T,
                                   RH,
                                   NaCl_10pctRH_R_0[0],
                                   T,
                                   mfs,
                                   initial_position = np.array([0.1/1e3, 0, 1/1e3]),
                                   initial_velocity = np.array([NaCl_10pctRH_v_x0[0],
                                                                0,
                                                                NaCl_10pctRH_v_y0[0]]),
                                   gravity = np.array([0, 0, 9.80665]),
                                   gas_velocity = np.array([0,0,NaCl_10pctRH_gas_flow - 1e-3]),
                                   terminate_on_equilibration=False)
    # Obtain a table giving a history of *all* droplet parameters.
    history_neg_gf = droplet.complete_trajectory(trajectory)

    #plot_NaCl_sensitivity(history, history_pos_gf, history_neg_gf, NaCl_10pctRH_data['all'], 'Gas Flow: 0.0068 m/s ± 0.001 m/s')

if __name__ == '__main__':
    T = 294.15 # Kelvin
    mfs = 0.05
    RH = 0.1

    ### fully parameterised salt
    droplet, trajectory = simulate(NaCl_10pctRH_t,
                                   aqueous_NaCl,
                                   T,
                                   RH,
                                   NaCl_10pctRH_R_0[0],
                                   T,
                                   mfs,
                                   initial_position = np.array([0.1/1e3, 0, 1/1e3]),
                                   initial_velocity = np.array([NaCl_10pctRH_v_x0[0],
                                                                0,
                                                                NaCl_10pctRH_v_y0[0]]),
                                   gravity = np.array([0, 0, 9.80665]),
                                   gas_velocity = np.array([0,0,NaCl_10pctRH_gas_flow]),
                                   terminate_on_equilibration=False)
    # Obtain a table giving a history of *all* droplet parameters.
    history = droplet.complete_trajectory(trajectory)

    ### ideal a_w

    droplet, trajectory = simulate(NaCl_10pctRH_t,
                                   aqueous_NaCl_a_ideal,
                                   T,
                                   RH,
                                   NaCl_10pctRH_R_0[0],
                                   T,
                                   mfs,
                                   initial_position = np.array([0.1/1e3, 0, 1/1e3]),
                                   initial_velocity = np.array([NaCl_10pctRH_v_x0[0],
                                                                0,
                                                                NaCl_10pctRH_v_y0[0]]),
                                   gravity = np.array([0, 0, 9.80665]),
                                   gas_velocity = np.array([0,0,NaCl_10pctRH_gas_flow]),
                                   terminate_on_equilibration=False)
    # Obtain a table giving a history of *all* droplet parameters.
    history_a_ideal = droplet.complete_trajectory(trajectory)

    droplet, trajectory = simulate(NaCl_10pctRH_t,
                                   aqueous_NaCl_d_linear,
                                   T,
                                   RH,
                                   NaCl_10pctRH_R_0[0],
                                   T,
                                   mfs,
                                   initial_position = np.array([0.1/1e3, 0, 1/1e3]),
                                   initial_velocity = np.array([NaCl_10pctRH_v_x0[0],
                                                                0,
                                                                NaCl_10pctRH_v_y0[0]]),
                                   gravity = np.array([0, 0, 9.80665]),
                                   gas_velocity = np.array([0,0,NaCl_10pctRH_gas_flow]),
                                   terminate_on_equilibration=False)
    # Obtain a table giving a history of *all* droplet parameters.
    history_d_linear = droplet.complete_trajectory(trajectory)

    droplet, trajectory = simulate(NaCl_10pctRH_t,
                                   aqueous_NaCl_d_half,
                                   T,
                                   RH,
                                   NaCl_10pctRH_R_0[0],
                                   T,
                                   mfs,
                                   initial_position = np.array([0.1/1e3, 0, 1/1e3]),
                                   initial_velocity = np.array([NaCl_10pctRH_v_x0[0],
                                                                0,
                                                                NaCl_10pctRH_v_y0[0]]),
                                   gravity = np.array([0, 0, 9.80665]),
                                   gas_velocity = np.array([0,0,NaCl_10pctRH_gas_flow]),
                                   terminate_on_equilibration=False)
    # Obtain a table giving a history of *all* droplet parameters.
    history_d_half = droplet.complete_trajectory(trajectory)


# ## R_0, V_x & V_y sensitivity

if __name__ == '__main__':
    def plot_fill_region(ax, x1, y1, x2, y2, colour, label):
        '''
        Useful for plotting regions with different sizes of data
        '''
        ax.fill(np.append(x1, x2[::-1]),
                np.append(y1, y2[::-1]),
                color = colour, alpha = 0.7, label = label)
        return ax

    fig, (ax0, ax1) = plt.subplots(1,2, figsize = figure_size, dpi = resolution)

    #Vx
    ax0 = plot_fill_region(ax0,
                           history_pos_Vx.x/1e-3, history_pos_Vx.z/1e-3,
                           history_neg_Vx.x/1e-3, history_neg_Vx.z/1e-3,
                           '#DC267F', 'V$_x$')
    #Vy
    ax0 = plot_fill_region(ax0,
                           history_pos_Vy.x/1e-3, history_pos_Vy.z/1e-3,
                           history_neg_Vy.x/1e-3, history_neg_Vy.z/1e-3,
                           '#FFB000', 'V$_z$')
    # Radius
    ax0 = plot_fill_region(ax0,
                           history_pos_R0.x/1e-3, history_pos_R0.z/1e-3,
                           history_neg_R0.x/1e-3, history_neg_R0.z/1e-3,
                           '#785EF0', 'R$_0$')

    ax0.plot(history.x/1e-3, history.z/1e-3, 'k', lw = 1, label = 'Simulation')


    #ax0.fill_between(history_pos_R0.x/1e-3, history_pos_R0.z/1e-3, history_neg_R0.z/1e-3, '--', color = 'red', lw = 2, label = '+ve error')

    ax0.errorbar(NaCl_10pctRH_data['all'].x_position_mm_mean,
                 NaCl_10pctRH_data['all'].y_position_mm_mean,
                 NaCl_10pctRH_data['all'].x_position_mm_std,
                 NaCl_10pctRH_data['all'].y_position_mm_std,
                 fmt = 'o', color = 'k', elinewidth = 1, capsize = 2, ms = 3, label = 'FDC data')

    ax0.set_xscale('log')
    ax0.set_yscale('log')
    ax0.set_ylim(100,0.9)
    ax0.set_ylabel('Vertical Position / mm')
    ax0.set_xlabel('Horizontal Position / mm')

    # Radius
    ax1 = plot_fill_region(ax1,
                           history_pos_R0.time, history_pos_R0.radius / 1e-6,
                           history_neg_R0.time, history_neg_R0.radius / 1e-6,
                           '#785EF0', 'R$_0$')
    #Vx
    ax1 = plot_fill_region(ax1,
                           history_pos_Vx.time, history_pos_Vx.radius / 1e-6,
                           history_neg_Vx.time, history_neg_Vx.radius / 1e-6,
                           '#DC267F', 'Vx')
    #Vy
    ax1 = plot_fill_region(ax1,
                           history_pos_Vy.time, history_pos_Vy.radius / 1e-6,
                           history_neg_Vy.time, history_neg_Vy.radius / 1e-6,
                           '#FFB000', 'Vy')

    ax1.plot(history.time, history.radius / 1e-6, 'k', lw = 1, label = 'Simulation')
    ax1.errorbar(NaCl_10pctRH_data['all'].Time_s,
                 NaCl_10pctRH_data['all'].d_e_um_mean / 2,
                 NaCl_10pctRH_data['all'].d_e_um_std / 2,
                 fmt = 'o', color = 'k', elinewidth = 1, capsize = 2, ms = 3, label = 'FDC data')

    ax1.set_ylabel('Radius / µm')
    ax1.set_xlabel('Time / s')
    ax1.legend()

    ax2 = fig.add_axes([0.185, 0.3, 0.23, 0.4])
    ax2.patch.set_alpha(0)
    plt.gca().invert_yaxis()

    # Radius
    ax2 = plot_fill_region(ax2,
                           history_pos_R0.time, history_pos_R0.z/1e-3,
                           history_neg_R0.time, history_neg_R0.z/1e-3,
                           '#785EF0', 'R$_0$')
    #Vx
    ax2 = plot_fill_region(ax2,
                           history_pos_Vx.time, history_pos_Vx.z/1e-3,
                           history_neg_Vx.time, history_neg_Vx.z/1e-3,
                           '#DC267F', 'V$_x$')
    #Vy
    ax2 = plot_fill_region(ax2,
                           history_pos_Vy.time, history_pos_Vy.z/1e-3,
                           history_neg_Vy.time, history_neg_Vy.z/1e-3,
                           '#FFB000', 'V$_z$')

    ax0.plot(history.x/1e-3, history.z/1e-3, 'k', lw = 1, label = 'Simulation')

    ax2.plot(history.time, history.z/1e-3, 'k', lw = 1, label = 'Simulation')

    ax2.errorbar(NaCl_10pctRH_data['all'].Time_s,
                 NaCl_10pctRH_data['all'].y_position_mm_mean,
                 NaCl_10pctRH_data['all'].y_position_mm_std,
                 fmt = 'o', color = 'k', elinewidth = 1, capsize = 2, ms = 3, label = 'FDC data')
    ax2.set_xlabel('Time / s')
    ax2.set_ylabel('Distance Fallen / mm')

    plt.show()

# ## Gas Flow, RH & T sensitivity

if __name__ == '__main__':
    fig, (ax0, ax1) = plt.subplots(1,2, figsize = figure_size, dpi = resolution)

    # Gas Flow
    ax0 = plot_fill_region(ax0,
                           history_pos_gf.x/1e-3, history_pos_gf.z/1e-3,
                           history_neg_gf.x/1e-3, history_neg_gf.z/1e-3,
                           '#5D3A9B', 'Gas Flow')
    #RH
    ax0 = plot_fill_region(ax0,
                           history_pos_RH.x/1e-3, history_pos_RH.z/1e-3,
                           history_neg_RH.x/1e-3, history_neg_RH.z/1e-3,
                           '#40B0A6', 'RH')
    #T
    ax0 = plot_fill_region(ax0,
                           history_pos_T.x/1e-3, history_pos_T.z/1e-3,
                           history_neg_T.x/1e-3, history_neg_T.z/1e-3,
                           '#E66100', 'T')

    ax0.plot(history.x/1e-3, history.z/1e-3, 'k', lw = 1, label = 'Simulation')

    ax0.errorbar(NaCl_10pctRH_data['all'].x_position_mm_mean,
                 NaCl_10pctRH_data['all'].y_position_mm_mean,
                 NaCl_10pctRH_data['all'].x_position_mm_std,
                 NaCl_10pctRH_data['all'].y_position_mm_std,
                 fmt = 'o', color = 'k', elinewidth = 1, capsize = 2, ms = 3, label = 'FDC data')

    ax0.set_xscale('log')
    ax0.set_yscale('log')
    ax0.set_ylim(100,0.9)
    ax0.set_ylabel('Vertical Position / mm')
    ax0.set_xlabel('Horizontal Position / mm')

    # Gas Flow
    ax1 = plot_fill_region(ax1,
                           history_pos_gf.time, history_pos_gf.radius / 1e-6,
                           history_neg_gf.time, history_neg_gf.radius / 1e-6,
                           '#5D3A9B', 'Gas Flow')
    #RH
    ax1 = plot_fill_region(ax1,
                           history_pos_RH.time, history_pos_RH.radius / 1e-6,
                           history_neg_RH.time, history_neg_RH.radius / 1e-6,
                           '#40B0A6', 'RH')
    #T
    ax1 = plot_fill_region(ax1,
                           history_pos_T.time, history_pos_T.radius / 1e-6,
                           history_neg_T.time, history_neg_T.radius / 1e-6,
                           '#E66100', 'T')

    ax1.plot(history.time, history.radius / 1e-6, 'k', lw = 1, label = 'Simulation')
    ax1.errorbar(NaCl_10pctRH_data['all'].Time_s,
                 NaCl_10pctRH_data['all'].d_e_um_mean / 2,
                 NaCl_10pctRH_data['all'].d_e_um_std / 2,
                 fmt = 'o', color = 'k', elinewidth = 1, capsize = 2, ms = 3, label = 'FDC data')

    ax1.set_ylabel('Radius / µm')
    ax1.set_xlabel('Time / s')
    ax1.legend()

    ax2 = fig.add_axes([0.185, 0.3, 0.23, 0.4])
    ax2.patch.set_alpha(0)
    plt.gca().invert_yaxis()

    # Gas Flow
    ax2 = plot_fill_region(ax2,
                           history_pos_gf.time, history_pos_gf.z/1e-3,
                           history_neg_gf.time, history_neg_gf.z/1e-3,
                           '#5D3A9B', 'Gas Flow')
    #RH
    ax2 = plot_fill_region(ax2,
                           history_pos_RH.time, history_pos_RH.z/1e-3,
                           history_neg_RH.time, history_neg_RH.z/1e-3,
                           '#40B0A6', 'RH')
    #T
    ax2 = plot_fill_region(ax2,
                           history_pos_T.time, history_pos_T.z/1e-3,
                           history_neg_T.time, history_neg_T.z/1e-3,
                           '#E66100', 'T')

    ax2.plot(history.time, history.z/1e-3, 'k', lw = 1, label = 'Simulation')

    ax2.errorbar(NaCl_10pctRH_data['all'].Time_s,
                 NaCl_10pctRH_data['all'].y_position_mm_mean,
                 NaCl_10pctRH_data['all'].y_position_mm_std,
                 fmt = 'o', color = 'k', elinewidth = 1, capsize = 2, ms = 3, label = 'FDC data')
    ax2.set_xlabel('Time / s')
    ax2.set_ylabel('Distance Fallen / mm')


    plt.show()

# ## Activity Parameterisation Sensitivity

if __name__ == '__main__':
    fig, (ax0, ax1) = plt.subplots(1,2, figsize = figure_size, dpi = resolution)

    ax0.plot(history.x/1e-3, history.z/1e-3, 'k', lw = 1, label = 'Simulation')

    ax0.errorbar(NaCl_10pctRH_data['all'].x_position_mm_mean,
                 NaCl_10pctRH_data['all'].y_position_mm_mean,
                 NaCl_10pctRH_data['all'].x_position_mm_std,
                 NaCl_10pctRH_data['all'].y_position_mm_std,
                 fmt = 'o', color = 'k', elinewidth = 1, capsize = 2, ms = 3, label = 'FDC data')

    ax0.plot(history_a_ideal.x/1e-3, history_a_ideal.z/1e-3, '--', color = 'r', lw = 2, label = "Raoult's Law")

    ax0.set_xscale('log')
    ax0.set_yscale('log')
    ax0.set_ylim(100,0.9)
    ax0.set_ylabel('Vertical Position / mm')
    ax0.set_xlabel('Horizontal Position / mm')

    ax1.plot(history.time, history.radius / 1e-6, 'k', lw = 1, label = 'E-AIM Activity')

    ax1.plot(history_a_ideal.time, history_a_ideal.radius / 1e-6, '--', color = 'r', lw = 2, label = "Raoult's Law")

    ax1.errorbar(NaCl_10pctRH_data['all'].Time_s,
                 NaCl_10pctRH_data['all'].d_e_um_mean / 2,
                 NaCl_10pctRH_data['all'].d_e_um_std / 2,
                 fmt = 'o', color = 'k', elinewidth = 1, capsize = 2, ms = 3, label = 'FDC data')

    ax1.set_ylabel('Radius / µm')
    ax1.set_xlabel('Time / s')
    ax1.legend()

    ax2 = fig.add_axes([0.185, 0.3, 0.23, 0.4])
    ax2.patch.set_alpha(0)
    plt.gca().invert_yaxis()

    ax2.plot(history.time, history.z/1e-3, 'k', lw = 1, label = 'E-AIM Activity')

    ax2.errorbar(NaCl_10pctRH_data['all'].Time_s,
                 NaCl_10pctRH_data['all'].y_position_mm_mean,
                 NaCl_10pctRH_data['all'].y_position_mm_std,
                 fmt = 'o', color = 'k', elinewidth = 1, capsize = 2, ms = 3, label = 'FDC data')

    ax2.plot(history_a_ideal.time, history_a_ideal.z/1e-3, '--', color = 'r', lw = 2, label = "Raoult's Law")

    ax2.set_xlabel('Time / s')
    ax2.set_ylabel('Distance Fallen / mm')


    plt.show()

# ## Density Paramterisation Sensitivity

if __name__ == '__main__':
    fig, (ax0, ax1) = plt.subplots(1,2, figsize = figure_size, dpi = resolution)

    ax0.plot(history.x/1e-3, history.z/1e-3, 'k', lw = 1, label = 'Simulation')

    ax0.errorbar(NaCl_10pctRH_data['all'].x_position_mm_mean,
                 NaCl_10pctRH_data['all'].y_position_mm_mean,
                 NaCl_10pctRH_data['all'].x_position_mm_std,
                 NaCl_10pctRH_data['all'].y_position_mm_std,
                 fmt = 'o', color = 'k', elinewidth = 1, capsize = 2, ms = 3, label = 'FDC data')


    ax0.plot(history_d_linear.x/1e-3, history_d_linear.z/1e-3, '--', color = 'dodgerblue', lw = 2, label = "Linear Density")
    ax0.plot(history_d_half.x/1e-3, history_d_half.z/1e-3, '-.', color = '#004D40', lw = 2, label = "Half Density")


    ax0.set_xscale('log')
    ax0.set_yscale('log')
    ax0.set_ylim(100,0.9)
    ax0.set_ylabel('Vertical Position / mm')
    ax0.set_xlabel('Horizontal Position / mm')

    ax1.plot(history.time, history.radius / 1e-6, 'k', lw = 1, label = 'E-AIM Density')


    ax1.plot(history_d_linear.time, history_d_linear.radius / 1e-6, '--', color = 'dodgerblue', lw = 2, label = "Linear Density")
    ax1.plot(history_d_half.time, history_d_half.radius / 1e-6, '-.', color = '#004D40', lw = 2, label = "Half Density")

    ax1.errorbar(NaCl_10pctRH_data['all'].Time_s,
                 NaCl_10pctRH_data['all'].d_e_um_mean / 2,
                 NaCl_10pctRH_data['all'].d_e_um_std / 2,
                 fmt = 'o', color = 'k', elinewidth = 1, capsize = 2, ms = 3, label = 'FDC data')

    ax1.set_ylabel('Radius / µm')
    ax1.set_xlabel('Time / s')
    ax1.legend()

    ax2 = fig.add_axes([0.185, 0.3, 0.23, 0.4])
    ax2.patch.set_alpha(0)
    plt.gca().invert_yaxis()

    ax2.plot(history.time, history.z/1e-3, 'k', lw = 1, label = 'E-AIM Density')

    ax2.errorbar(NaCl_10pctRH_data['all'].Time_s,
                 NaCl_10pctRH_data['all'].y_position_mm_mean,
                 NaCl_10pctRH_data['all'].y_position_mm_std,
                 fmt = 'o', color = 'k', elinewidth = 1, capsize = 2, ms = 3, label = 'FDC data')


    ax2.plot(history_d_linear.time, history_d_linear.z/1e-3, '--', color = 'dodgerblue', lw = 2, label = "Linear Density")
    ax2.plot(history_d_half.time, history_d_half.z/1e-3, '-.', color = '#004D40', lw = 2, label = "Half Density")

    ax2.set_xlabel('Time / s')
    ax2.set_ylabel('Distance Fallen / mm')


    plt.show()

# # 4.5.2. Parameter Sensiivty

if __name__ == '__main__':
    def scale_range_for_colorbar(x):
        return (x-x.min())/(x.max() - x.min())

    T = T_freezing + 20
    RH = 0.1
    solution = aqueous_NaCl
    initial_velocity = np.array([1, 0, 0])
    time_limit = 10
    mfs = 0.05
    R0 = 25e-6

    droplet, trajectory = simulate(time_limit,
                                   solution,
                                   T,
                                   RH,
                                   NaCl_10pctRH_R_0[0],
                                   T,
                                   mfs,
                                   initial_position = np.array([0.1/1e3, 0, 1/1e3]),
                                   initial_velocity = initial_velocity,
                                   gravity = np.array([0, 0, 9.80665]),
                                   gas_velocity = np.array([0,0,0]),
                                   terminate_on_equilibration=False)
    # Obtain a table giving a history of *all* droplet parameters.
    history = droplet.complete_trajectory(trajectory)

if __name__ == '__main__':
    def plot_parameter_comparison(label, range_name, history_name, history_w_eff_name, cmap_name = 'plasma_r'):
        fig, (ax0, ax1) = plt.subplots(1,2, figsize = figure_size, dpi = resolution)

        ax0.set_xscale('log')
        ax0.set_yscale('log')
        ax0.set_ylim(200,0.9)
        ax0.set_ylabel('Vertical Position / mm')
        ax0.set_xlabel('Horizontal Position / mm')


        ax1.set_ylabel('Radius / µm')
        ax1.set_xlabel('Time / s')
        #ax1.legend()

        ax2 = fig.add_axes([0.195, 0.3, 0.18, 0.4])
        ax2.patch.set_alpha(0)
        plt.gca().invert_yaxis()

        ax2.set_xlabel('Time / s')
        ax2.set_ylabel('Distance Fallen / mm')

        cmap = mpl.cm.get_cmap(cmap_name)
        norm = mpl.colors.Normalize(vmin=np.min(range_name), vmax=np.max(range_name))

        colors = cmap(scale_range_for_colorbar(range_name))

        for history, history_w_eff, color in zip(history_name, history_w_eff_name, colors):

            #pre efflorescence
            ax0.plot(history[history.time < history_w_eff.tail(1).time.values[0]].x/1e-3,
                     history[history.time < history_w_eff.tail(1).time.values[0]].z/1e-3, color = color, lw = 1, )
            ax1.plot(history[history.time < history_w_eff.tail(1).time.values[0]].time,
                     history[history.time < history_w_eff.tail(1).time.values[0]].radius / 1e-6, color = color, lw = 1, )
            ax2.plot(history[history.time < history_w_eff.tail(1).time.values[0]].time,
                     history[history.time < history_w_eff.tail(1).time.values[0]].z/1e-3, color = color, lw = 1, )

            #post efflorescence
            ax0.plot(history[history.time > history_w_eff.tail(1).time.values[0]].x/1e-3,
                     history[history.time > history_w_eff.tail(1).time.values[0]].z/1e-3, color = color, ls = '--', lw = 1, )
            ax1.plot(history[history.time > history_w_eff.tail(1).time.values[0]].time,
                     history[history.time > history_w_eff.tail(1).time.values[0]].radius / 1e-6, color = color, ls = '--', lw = 1, )
            ax2.plot(history[history.time > history_w_eff.tail(1).time.values[0]].time,
                     history[history.time > history_w_eff.tail(1).time.values[0]].z/1e-3, color = color, ls = '--', lw = 1, )

            #[history.time > history_w_eff.tail(1).time.values[0]]

            if history_w_eff.tail(1).time.max() < history.time.max():
                ax0.scatter(history_w_eff.tail(1).x/1e-3, history_w_eff.tail(1).z/1e-3, marker = '*', s = 50, color = color, lw = 1, )
                ax1.scatter(history_w_eff.tail(1).time, history_w_eff.tail(1).radius / 1e-6, marker = '*', s = 50, color = color, lw = 1, )
                ax2.scatter(history_w_eff.tail(1).time, history_w_eff.tail(1).z/1e-3, marker = '*', s = 50, color = color, lw = 1, )

        plt.colorbar(mpl.cm.ScalarMappable(norm, cmap), ax = [ax0, ax1], orientation = 'vertical', pad = 0.02, aspect = 40, label = label )

        plt.show()

if __name__ == '__main__':
    def generate_parameter_comparison(solution = aqueous_NaCl,
                                      time_limit = 10,
                                      sample_points = 11,
                                      T_min = T_freezing + 20, T_max = T_freezing + 20,
                                      RH_min = 0.3, RH_max = 0.3,
                                      mfs_min = 0.05, mfs_max = 0.05,
                                      R0_min = 25e-6, R0_max = 25e-6,
                                      Vx_min = 1, Vx_max = 1,
                                      Vy_min = 0, Vy_max = 0,
                                      Vz_min = 0, Vz_max = 0):
        """
        For generating and plotting single parameter dependencies.
        Recommended to only vary one parameter at once. Use defaults for all other parameters
        """

        T_range = np.linspace(T_min, T_max, sample_points)
        RH_range = np.linspace(RH_min, RH_max, sample_points)
        mfs_range = np.linspace(mfs_min, mfs_max, sample_points)
        Vx_range = np.linspace(Vx_min, Vx_max, sample_points)
        Vy_range = np.linspace(Vy_min, Vy_max, sample_points)
        Vz_range = np.linspace(Vz_min, Vz_max, sample_points)

        history_list = []
        history_list_with_eff = []

        for T, RH, mfs, Vx, Vy, Vz in zip(T_range, RH_range, mfs_range, Vx_range, Vy_range, Vz_range):

            droplet, trajectory = simulate(time_limit,
                                           solution,
                                           T,
                                           RH,
                                           R0,
                                           T,
                                           mfs,
                                           initial_position = np.array([0.1/1e3, 0, 1/1e3]),
                                           initial_velocity = np.array([Vx, Vy, Vz]),
                                           gravity = np.array([0, 0, 9.80665]),
                                           gas_velocity = np.array([0,0,0]),
                                           rtol = 1e-8,
                                           terminate_on_equilibration=False,
                                           terminate_on_efflorescence=False)
            # Obtain a table giving a history of *all* droplet parameters.
            history = droplet.complete_trajectory(trajectory)
            history_list.append(history)

            droplet, trajectory = simulate(time_limit,
                                           solution,
                                           T,
                                           RH,
                                           R0,
                                           T,
                                           mfs,
                                           initial_position = np.array([0.1/1e3, 0, 1/1e3]),
                                           initial_velocity = np.array([Vx, Vy, Vz]),
                                           gravity = np.array([0, 0, 9.80665]),
                                           gas_velocity = np.array([0,0,0]),
                                           rtol = 1e-8,
                                           terminate_on_equilibration=False,
                                           terminate_on_efflorescence=True)
            # Obtain a table giving a history of *all* droplet parameters.
            history = droplet.complete_trajectory(trajectory)
            history_list_with_eff.append(history)

        return history_list, history_list_with_eff

if __name__ == '__main__':
    plot_parameter_comparison('RH / %', np.linspace(0.25, 0.35,11), cmap_name='plasma', *generate_parameter_comparison(RH_min=0.25, RH_max=0.35))
    plot_parameter_comparison('T / K', np.linspace(273 + 15, 273 +25 ,11), cmap_name='plasma', *generate_parameter_comparison(T_min= 273 + 15, T_max= 273 + 25))
    plot_parameter_comparison('MFS', np.linspace(0, 0.25 ,11), cmap_name='plasma', *generate_parameter_comparison(mfs_min= 0, mfs_max= 0.25))
    plot_parameter_comparison('V$_x$ / ms$^{-1}$', np.linspace(0.9, 1 ,11), cmap_name='plasma', *generate_parameter_comparison(Vx_min= 0.9, Vx_max= 1))

if __name__ == '__main__':


if __name__ == '__main__':


if __name__ == '__main__':


if __name__ == '__main__':


if __name__ == '__main__':


if __name__ == '__main__':


if __name__ == '__main__':


if __name__ == '__main__':


if __name__ == '__main__':

