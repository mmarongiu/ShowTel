import matplotlib
#matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import matplotlib.colors as colors
from matplotlib.image import imread
import matplotlib as mpl
from matplotlib.ticker import AutoMinorLocator, MultipleLocator
from astropy.constants import *
import numpy as np
from .showtel_distance import *

matplotlib.rc('axes',edgecolor='g')


def plot_radar_sun(dir_out, pos_live, pos_sun, pos_moon, tempo, receiver, dpi_q):

    #pos_sun = (0.3,0.3)
    #pos_moon = (0.1,0.1)
    #pos_live = (0.6,0.4)

    c_sun = plt.Circle(pos_sun, 0.27, color='y')
    c_moon = plt.Circle(pos_moon, 0.27, color='grey')

    fig, ax = plt.subplots()
    plt.style.use('dark_background')
    ax.add_patch(c_sun)
    ax.add_patch(c_moon)
    plt.plot(pos_live[0],pos_live[1],'ro')
    ax.axis('equal')
    ax.grid(color='green')
    ax.set_facecolor('k')
    ax.set_xlabel('Azimuth [deg]', fontsize=18, color='green')
    ax.set_ylabel('Elevation [deg]', fontsize=18, color='green')

    ax.set_title(receiver + ' - ' + tempo, fontsize=18, color='green')
    
    ax.set_xlim((pos_sun[0]-0.3, pos_sun[0]+0.3))
    ax.set_ylim((pos_sun[1]-0.3, pos_sun[1]+0.3))

    ax.tick_params(axis='both', which='major', labelsize=16, colors='green')

    plt.tight_layout()
    plt.savefig(dir_out + 'radar.png', dpi=dpi_q)
    plt.close()

    
def plot_radar(dir_out, time_ref, nomesource, pos_live, pos_comm, res0, res1, altaz_comm, time_source_evo, tmin, tmax,  sun_mode, sun_opt, rfi_tab, rfi_freq, dpi_q, ns):
    # res0 ------- matrix of results (output labelled as 0) of the function "sunmoon_evo" (LIVE POSITION)
    # res1 ------- matrix of results (output labelled as 1) of the function "sunmoon_evo" (GENERIC POSITION)
    # altaz_comm = (array_az, array_alt) [array of floats]
    # pos_live = (az, alt) [float]
    # pos_comm = (az, alt) [float]
    # time_source_evo = array of datetime
    # ns --------- 'N'[north]/'S'[south] visibility

    tempo_min, tempo_max = tmin.datetime, tmax.datetime
    tempo_ref = time_ref.datetime

    if sun_mode == 0:
        fig = plt.figure(figsize=((65.5,11.35)), facecolor='k', edgecolor='g')
        ax1 = plt.subplot(2,1,1)
        ax1.set_facecolor('k')
        ax1.grid(color='g')
        ax1.set_xlabel('Azimuth [deg]', fontsize=20, color='green')
        ax1.set_ylabel('Elevation [deg]', fontsize=20, color='green')

        if ns == 'N':
            ax1.text(145, 83, 'North', fontsize=20, color='green')
            ax1.plot(np.where(pos_comm[0] > 180, pos_comm[0]-360, pos_comm[0]), pos_comm[1], marker='.', color='w', markersize=18, linestyle='none')
            ax1.plot(np.where(res1[:,1].astype('float64') > 180, res1[:,1].astype('float64')-360, res1[:,1].astype('float64')), res1[:,2].astype('float64'), marker='.', color='y', markersize=4, linestyle='none')
            ax1.plot(np.where(res1[:,3].astype('float64') > 180, res1[:,3].astype('float64')-360, res1[:,3].astype('float64')), res1[:,4].astype('float64'), marker='.', color='gray', markersize=4, linestyle='none')
            ax1.plot(np.where(altaz_comm[0].astype('float64') > 180, altaz_comm[0].astype('float64')-360, altaz_comm[0].astype('float64')), altaz_comm[1].astype('float64'), marker='o', color='w', markersize=1, linestyle='none')
            ax1.plot(np.where(res0[1].value > 180, res0[1].value-360, res0[1].value), res0[2].value, marker='o', color='y', markersize=10, linestyle='none')
            ax1.plot(np.where(res0[3].value > 180, res0[3].value-360, res0[3].value), res0[4].value, marker='o', color='gray', markersize=10, linestyle='none')
            ax1.plot(np.where(pos_live[0] > 180, pos_live[0] - 360, pos_live[0]), pos_live[1], 'ro')
        
            ax1.fill_between(np.arange(-200, 200, 10), 6, color='y', alpha=0.2)
            ax1.fill_between(np.arange(-200, 200, 10), 87, 90, color='r', alpha=0.2)
            if isinstance(rfi_tab, str) == False:
                for i in rfi_freq:
                    if i['az_min'] >= 180:
                        az_min = i['az_min']-360
                    else:
                        az_min = i['az_min']
                    if i['az_max'] >= 180:
                        az_max = i['az_max']-360
                    else:
                        az_max = i['az_max']
                    
                    if (i['az_min'] < 180) and (i['az_max'] >= 180):
                        az_arr = [[az_min, 180], [-180, az_max]]
                        el_arr = [i['el_min'], i['el_max']]
                        for ii in az_arr:
                            ax1.fill_between(np.arange(ii[0], ii[1],0.1), el_arr[0], el_arr[1], color='b', alpha=0.5)
                    elif (i['az_min'] == 180) and (i['az_max'] >= 180):
                        az_arr = [-180, az_max]
                        el_arr = [i['el_min'], i['el_max']]
                        ax1.fill_between(np.arange(az_arr[0], az_arr[1], 0.1), el_arr[0], el_arr[1], color='b', alpha=0.5)
                    else:
                        ax1.fill_between(np.arange(az_min, az_max, 0.1), i['el_min'], i['el_max'], color='b', alpha=0.5)

        elif ns == 'S':
            ax1.text(325, 83, 'South', fontsize=20, color='green')
            ax1.plot(pos_comm[0], pos_comm[1], marker='.', color='w', markersize=18, linestyle='none')
            ax1.plot(res1[:,1].astype('float64'), res1[:,2].astype('float64'), marker='.', color='y', markersize=4, linestyle='none')
            ax1.plot(res1[:,3].astype('float64'), res1[:,4].astype('float64'), marker='.', color='gray', markersize=4, linestyle='none')
            ax1.plot(altaz_comm[0].astype('float64'), altaz_comm[1].astype('float64'), marker='o', color='w', markersize=1, linestyle='none')
            ax1.plot(res0[1].value, res0[2].value, marker='o', color='y', markersize=10, linestyle='none')
            ax1.plot(res0[3].value, res0[4].value, marker='o', color='gray', markersize=10, linestyle='none')
            ax1.plot(pos_live[0], pos_live[1], 'ro')
        
            ax1.fill_between(np.arange(0, 400, 10), 6, color='y', alpha=0.2)
            ax1.fill_between(np.arange(0, 400, 10), 87, 90, color='r', alpha=0.2)
            if isinstance(rfi_tab, str) == False:
                for i in rfi_freq:
                    az_min = i['az_min']
                    az_max = i['az_max']
                    ax1.fill_between(np.arange(az_min, az_max, 0.1), i['el_min'], i['el_max'], color='b', alpha=0.5)
                        
        if ns == 'N':
            ax1.set_xlim((-180,180))
            ax1.xaxis.set_ticks([-180, -135, -90, -45, 0, 45, 90, 135, 180], ['180', '225', '270', '315', '0', '45', '90', '135', '180'])    # centrato a 0
        elif ns == 'S':
            ax1.set_xlim((0,360))
            ax1.xaxis.set_major_locator(MultipleLocator(45))
            
        ax1.set_ylim((0,90))
        ax1.tick_params(axis='both', which='both', labelsize=15, colors='g', grid_color='g', grid_alpha=0.8)

    if sun_mode == 1:
        fig = plt.figure(figsize=((15.9,11.35)), facecolor='k', edgecolor='g')
        pos_sun, pos_moon, tempo, receiver = sun_opt[0], sun_opt[1], sun_opt[2], sun_opt[3]
        c_sun = plt.Circle(pos_sun, 0.27, color='y')
        c_moon = plt.Circle(pos_moon, 0.27, color='grey')
        ax1 = plt.subplot(2,1,1)
        plt.style.use('dark_background')
        ax1.add_patch(c_sun)
        ax1.add_patch(c_moon)
        ax1.plot(pos_live[0],pos_live[1],'ro')
        ax1.axis('equal')
        ax1.grid(color='green')
        ax1.set_facecolor('k')
        ax1.set_xlabel('Azimuth [deg]', fontsize=22, color='green')
        ax1.set_ylabel('Elevation [deg]', fontsize=22, color='green')
        ax1.set_title(receiver + ' - ' + tempo, fontsize=25, color='green')
        ax1.set_xlim((pos_sun[0]-0.3, pos_sun[0]+0.3))
        ax1.set_ylim((pos_sun[1]-0.3, pos_sun[1]+0.3))
        if isinstance(rfi_tab, str) == False:
            for i in rfi_tab:
                ax1.fill_between(np.arange(i['az_min'], i['az_max'],0.1), i['el_min'], i['el_max'], color='b', alpha=0.5)
        ax1.tick_params(axis='both', which='major', labelsize=16, colors='green')
        
    from astropy.visualization import time_support
    time_support()

    ax2 = fig.add_subplot(2,1,2)
    ax3 = ax2.twiny()
    ax2.set_facecolor('k')
    ax2.grid(color='g')

    lst = []
    for w in res1[:,0]:
        u0 = w.sidereal_time('mean')
        lst.append(u0)

    new_tick_locations1, new_tick_locations2 = [], []
    for i in res1[:,0][::8]:
        new_tick_locations1 = np.append(new_tick_locations1, i)
        uu = i.sidereal_time('mean')
        new_tick_locations2.append(uu.value)

    time_datetime = []
    for i in res1[:,0]:
        time_datetime.append(i.datetime)

    if sun_mode == 0:
        time_source_datetime = []
        for i in time_source_evo:
            time_source_datetime.append(i.datetime)

        ax2.plot(time_source_datetime , altaz_comm[1].astype('float64'), marker='.', color='w', markersize=2, linestyle='none')
        ax2.plot(time_ref.datetime, pos_comm[1], marker='o', color='w', markersize=10, linestyle='none', label=nomesource)

    ax2.plot(time_datetime, res1[:,2].astype('float64'), marker='.', color='y', markersize=4, linestyle='none')
    ax2.plot(res0[0].datetime, res0[2].value, marker='o', color='y', markersize=10, linestyle='none', label='Sun - now')
    
    ax2.fill_between(time_datetime, 0, 10, color='y', alpha=0.2)
    #ax2.fill_between(time_datetime, 0, -90, color='r', alpha=0.5)
    ax2.fill_between(time_datetime, 80, 90, color='r', alpha=0.3)

    ax2.legend(loc='upper right', numpoints=1, fontsize=14)
    ax2.set_xlabel('UTC Time [hh:mm]', fontsize=22, color='g')
    ax2.set_ylabel('Elevation [deg]', fontsize=22, color='g')

    ax2.tick_params(axis='y', which='major', labelsize=16, colors='g')
    ax2.set_xlim((tempo_min, tempo_max))
    ax2.set_ylim((0,90))

    label1 = []
    for i in new_tick_locations1:
        label1 = np.append(label1,str(i)[11:16])
    
    ax2.set_xticks(time_datetime[::8], labels=label1, fontsize=16, color='g')#, rotation=45)
    ax2.minorticks_off()

    ax2.grid(color='g', ls='solid')

    def tick_function2(x):
        return ["%.2f" % z for z in x]

    ax3.set_xticks(time_datetime[::8], labels=tick_function2(new_tick_locations2), fontsize=16, color='g')#, rotation=45)
    ax3.set_xlim(ax2.get_xlim())
    ax3.minorticks_off()
    ax3.set_xlabel('LST Time [hourangle]', fontsize=22, color='g')
    ax3.grid(False)

    # get the bounds of ax1 and ax2
    x1, y1, w1, h1 = ax1.get_position().bounds
    x2, y2, w2, h2 = ax2.get_position().bounds
    x3, y3, w3, h3 = ax3.get_position().bounds

    # set ax1 width to width of ax2
    if sun_mode == 0:
        ax1.set_position([x1, y1, w2*2, h1])
    
    fig.tight_layout()
    plt.savefig(dir_out + 'radar.png', bbox_inches="tight", dpi=dpi_q)
    plt.close()
