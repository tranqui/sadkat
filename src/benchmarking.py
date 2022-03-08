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

    for k, RH in zip(Kulmala_k_means, Kulmala_RH_list):
        plt.scatter(RH,
                    k,
                    color = '#004D40', s = 50, marker = 'P')
    plt.plot(np.linspace(0,1,100),
             Kulmala_poly_model(np.linspace(0,1,100)),
             lw = 0.5,
             color = '#004D40')

    for k, RH in zip(Su_k_means, Su_RH_list):
        plt.scatter(RH,
                    k,
                    color = '#1E88E5', s = 50, marker = 'x')
    plt.plot(np.linspace(0,1,100),
             Su_poly_model(np.linspace(0,1,100)),
             lw = 0.5,
             color = '#1E88E5')

    for k, RH in zip(sadkat_k_means, sadkat_RH_list):
        plt.scatter(RH,
                    k,
                    color = '#D81B60', s = 50, marker = 'o')
    plt.plot(np.linspace(0,1,100),
             sadkat_poly_model(np.linspace(0,1,100)),
             lw = 0.5,
             color = '#D81B60')

    plt.xlabel('RH')
    plt.ylabel(r'Evaporation Rate / µm$^2$s$^{-1}$')

    legend_elements = [Line2D([0], [0], marker='P', color='#004D40', label='Kulmala',
                              markerfacecolor='#004D40', markersize=11, lw = 1),
                       Line2D([0], [0], marker='x', color='#1E88E5', label='Su',
                              markerfacecolor='#1E88E5', markersize=11, lw = 1),
                       Line2D([0], [0], marker='o', color='#D81B60', label='SADKAT',
                              markerfacecolor='#D81B60', markersize=11, lw = 1),]

    plt.legend(handles=legend_elements)

    plt.show()

    # comparing temperature suppresssion data
    for k, RH in zip(Kulmala_T_suppression, Kulmala_RH_list):
        plt.scatter(RH,
                    k,
                    color = '#004D40', s = 50, marker = 'P')
    plt.plot(np.linspace(0,1,100),
             Kulmala_T_poly_model(np.linspace(0,1,100)),
             lw = 0.5,
             color = '#004D40')

    for T, RH in zip(Su_T_suppression, Su_RH_list):
        plt.scatter(RH,
                    T,
                    color = '#1E88E5', s = 50, marker = 'x')
    plt.plot(np.linspace(0,1,100),
             Su_T_poly_model(np.linspace(0,1,100)),
             lw = 0.5,
             color = '#1E88E5')
    plt.xlabel('RH')

    for T, RH in zip(sadkat_T_suppression, sadkat_RH_list):
        plt.scatter(RH,
                    T,
                    color = '#D81B60', s = 50)
    plt.plot(np.linspace(0,1,100),
             sadkat_T_poly_model(np.linspace(0,1,100)),
             lw = 0.5,
             color = '#D81B60')
    plt.xlabel('RH')
    plt.ylabel(r'Temperature Suppression / K')
    legend_elements = [Line2D([0], [0], marker='P', color='#004D40', label='Kulmala',
                              markerfacecolor='#004D40', markersize=11, lw = 1),
                       Line2D([0], [0], marker='x', color='#1E88E5', label='Su',
                              markerfacecolor='#1E88E5', markersize=11, lw = 1),
                       Line2D([0], [0], marker='o', color='#D81B60', label='SADKAT',
                              markerfacecolor='#D81B60', markersize=11, lw = 1),]

    plt.legend(handles=legend_elements)
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
    T_range = np.linspace(5 + T_freezing, 25 + T_freezing,11) # Kelvin
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

    def k_surface_plot(T, RH, k_data,
                       fit_k_surface = False,
                       k_ref_value = None,
                       experiment_T = None, experiment_data = None, sadkat_experiment_data = None,
                       delta_T_data = None):


        T_mesh, RH_mesh, k_fit_surface_data, k_surface_parameters, k_surface_parameters_errors = fit_evaporation_surface(T, RH, k_data)

        #angles for viewing 3d plot
        alpha_ang = 30
        beta_ang = 135

        import matplotlib as mpl
        from mpl_toolkits.mplot3d import Axes3D

        # some markers and styles for use in plotting
        marker_list = ['x', 'o', 'v', '^', 's', 'P', 'p', '*', 'D']
        ls_list = [':','-.','--','-']

        colour_ticks = 50
        cmap = mpl.cm.jet
        boundaries = np.arange(k_data.T.min(),k_data.T.max() + 101,colour_ticks)
        norm = mpl.colors.BoundaryNorm(boundaries=boundaries, ncolors=256)

        #3d plot of evaporation rate over temperature and rh
        fig = plt.figure(figsize=(21,13))
        ax = fig.gca(projection='3d')

        #surface of data
        ax.plot_surface(T_mesh, RH_mesh, k_data, alpha = 0.3, cmap = cmap)
        #white contours to show levels
        ax.contour(T_mesh, RH_mesh, k_data, linewidths = 2, colors = 'w', levels = boundaries)

        #scatter data
        for i, Temp in enumerate(T):
            ax.scatter(Temp, RH, k_data.T[i], alpha = 0.9, c = k_data.T[i], cmap = cmap, norm=norm)

        #fitted surface
        if fit_k_surface != False:
            ax.plot_wireframe(T_mesh, RH_mesh, k_fit_surface_data, alpha = 1, color = 'k', lw = 0.5)

        if isinstance(experiment_T, np.ndarray):
                if isinstance(sadkat_experiment_data, np.ndarray):

                    #experiemental temperatures, model data
                    for i, (temp, style) in enumerate(zip(experiment_T, ls_list)):
                        ax.plot(np.full_like(RH,temp), RH, sadkat_experiment_data.T[i], 'r', ls = style, alpha = 0.7, lw = 5, label = str(round(float(temp))) + ' K')

        if isinstance(experiment_T, np.ndarray):
                if isinstance(experiment_data, dict):    

                    #plotting experimental data
                    for i, key in enumerate(experiment_data):
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
                        ax.plot(xs = [experiment_data[key]['exp_props']['T'][0] + T_freezing,
                                      experiment_data[key]['exp_props']['T'][0] + T_freezing],
                                    ys = [experiment_data[key]['exp_props']['RH'].mean() - experiment_data[key]['exp_props']['RH'].std(),
                                          experiment_data[key]['exp_props']['RH'].mean() + experiment_data[key]['exp_props']['RH'].std()],
                                    zs = [(experiment_data[key]['exp_props']['k'].mean()) * 1e12,
                                          (experiment_data[key]['exp_props']['k'].mean()) * 1e12],
                                color = 'k')
                        # k std dev
                        ax.plot(xs = [experiment_data[key]['exp_props']['T'][0] + T_freezing,
                                      experiment_data[key]['exp_props']['T'][0] + T_freezing],
                                    ys = [experiment_data[key]['exp_props']['RH'].mean(),
                                          experiment_data[key]['exp_props']['RH'].mean()],
                                    zs = [(experiment_data[key]['exp_props']['k'].mean() - experiment_data[key]['exp_props']['k'].std()) * 1e12,
                                          (experiment_data[key]['exp_props']['k'].mean() + experiment_data[key]['exp_props']['k'].std()) * 1e12],
                                color = 'k')
                        #plotting data
                        ax.scatter(xs = experiment_data[key]['exp_props']['T'][0] + T_freezing,
                                    ys = experiment_data[key]['exp_props']['RH'].mean(),
                                    zs = experiment_data[key]['exp_props']['k'].mean() * 1e12,
                                    color = 'k' ,marker = marker_list[i], s = 100,
                                    label = str(round(experiment_data[key]['exp_props']['T'][0] + T_freezing,3)) + ' K, ' + str(round(100* experiment_data[key]['exp_props']['RH'].mean())) + ' % RH')

        if k_ref_value != None:

            #k line for given evaporation rate:
            ax.plot3D(*get_constant_k_line(T, k_ref_value, k_surface_parameters), zs = k_ref_value, lw = 5, color = 'k', label = str(round(k_ref_value,3)) + ' µm$^2$s$^{-1}$')

        ax.view_init(alpha_ang, beta_ang)

        ax.set_ylim(0,1)
        ax.set_zlim(0)

        ax.set_xlabel('T / K', labelpad = 25)
        ax.set_ylabel('RH / %', labelpad = 25)
        ax.set_zlabel(r'Evaporation Rate / µm$^2$s$^{-1}$', labelpad = 25)

        plt.legend(bbox_to_anchor = (1.5,1))

        plt.show()

        if isinstance(delta_T_data, np.ndarray):

            cmap_RedBlue = colors.LinearSegmentedColormap.from_list("RedBlue",['r', 'b'])
            temp_colour_ticks = 5
            temp_boundaries = np.arange(delta_T_data.T.min(),delta_T_data.T.max(),temp_colour_ticks)
            temp_norm = mpl.colors.BoundaryNorm(boundaries=temp_boundaries, ncolors=256)

            #3d plot of evaporation rate over temperature and rh
            fig = plt.figure(figsize=(21,13))
            ax = fig.gca(projection='3d')

            #white contours to show levels
            ax.contour(T_mesh, RH_mesh, delta_T_data, linewidths = 2, colors = 'w', levels = temp_boundaries)
            #surface of data
            ax.plot_surface(T_mesh, RH_mesh, delta_T_data, alpha = 0.3, cmap = cmap_RedBlue)

            #scatter data
            for i, Temp in enumerate(T):
                ax.scatter(Temp, RH, delta_T_data.T[i], alpha = 0.9, c = delta_T_data.T[i], cmap = cmap_RedBlue, norm = temp_norm)

            ax.view_init(alpha_ang, beta_ang)

            ax.set_ylim(0,1)

            ax.set_xlabel('T / K', labelpad = 25)
            ax.set_ylabel('RH / %', labelpad = 25)
            ax.set_zlabel(r'Temperature Suppression / K', labelpad = 25)

            plt.legend(bbox_to_anchor = (1.5,1))

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
