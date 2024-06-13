import numpy as np
from numpy import deg2rad
from astropy import units as u
from astropy.coordinates import Angle, angular_separation
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
from astropy.coordinates import get_sun, get_moon
from astropy.coordinates import solar_system_ephemeris, get_body_barycentric, get_body
from astropy.time import Time
from astropy.constants import au

from datetime import datetime, timedelta
import pytz, time, ephem, sched, threading, socket, os

from astropy.coordinates import solar_system_ephemeris
solar_system_ephemeris.set('jpl')

from tkinter import *
import tkinter as tk
from tkinter import ttk, font
from tkmacosx import Button

from PIL import ImageTk, Image
from pydub import AudioSegment
from pydub.playback import play
from gtts import gTTS

from utilities.showtel_connection import *
from utilities.showtel_discos import *
from utilities.showtel_sound import *
from utilities.showtel_distance import *
#from utilities.showtel_switch import *
from utilities.showtel_forbidden import *
from utilities.showtel_plots import *

import warnings
warnings.filterwarnings("ignore")

update = '2024/06/13'
radiotelescope = 'SRT'                                                                     # SRT or Medicina
delta_time = 4                                                                             # in seconds

tot_eph_h = 24                                                                             # duration of the solar map [in units of hours, float]
step_eph_min = 30                                                                          # step time [in units of minutes, int]

# installare:
# sudo apt-get install ffmpeg libavcodec-extra
# pip install pydub
# conda install gtts


def start_butt():
    log_action('Showtel is ON','black')
    if data_input['app'] == 'DISCOS':
        s, c, ss = on_app_discos()
    elif data_input['app'] == 'Seadas':
        s, c, ss = on_app_seadas()
    
    global running
    running = True

    global stato
    stato = 1
    
    loop(s)

    
def stop_butt():
    log_action('Showtel is OFF','black')
    if data_input['app'] == 'DISCOS':
        s, c, ss = on_app_discos()
    elif data_input['app'] == 'Seadas':
        s, c, ss = on_app_seadas()
    s.close()

    global running
    running = False

    global stato
    stato = 0


def loop(s):
    stile = "Verdana"
    dimensione = 20
    delta_t = int(delta_time*1e3)                                                          # in ms

    file_sys = './utilities/'

    file_alarm = file_sys + 'sound/warning_alarm.wav'
    file_connected = file_sys + 'sound/connected.mp3'
    file_disconnected = file_sys + 'sound/disconnected.mp3'
    file_failure = file_sys + 'sound/failure.mp3'
    file_close = file_sys + 'sound/sunclose.mp3'
    file_rfi = file_sys + 'list_rfi_v1.dat'
    file_rec = file_sys + 'list_receivers.dat'
    
    len_box = 22
    
    if data_input['app'] == 'DISCOS':
        stato = connect_status_discos(s)
    elif data_input['app'] == 'Seadas':
        stato = connect_status_seadas(s)

    if data_input["visibility"] == "SOUTH":
        vis_mode = 'S'
    elif data_input["visibility"] == "NORTH":
        vis_mode = 'N'

    if data_input["quality"] == "LQ":
        q_mode = 1
    elif data_input["quality"] == "HQ":
        q_mode = 2
        
    if stato == 0 or running == False:
        label_text = 'DISCONNECTED'
        data_input.update({"Status": label_text})
        my_canvas.itemconfig(my_oval, fill="red")                                          # Fill the circle with RED

        # Box - left
        var_time.set('')                                                                   # set default path
        var_time_lst.set('')                                                               # set default path
        var_stat.set('')                                                                   # set default path
        txt_stat_box = Entry(tool3_bar, width=len_box-3, textvariable=var_stat, fg="white").place(in_=tool3_bar, relx=0.50, rely=0.40, anchor=W)

        wi2, he2, r2 = 15, 15, 7
        canvas_stat = tk.Canvas(tool3_bar, width=wi2, height=he2, borderwidth=0, bg='SkyBlue1', highlightbackground = 'SkyBlue1')    # Create 200x200 Canvas widget
        canvas_stat.place(relx=0.42, rely=0.33)
        oval_stat = canvas_stat.create_oval(wi2/2-r2, he2/2+r2, wi2/2+r2, he2/2-r2)        # Create a circle on the Canvas

        var_pntg.set('')                                                                   # set default path
        txt_pntg_box = Entry(tool3_bar, width=len_box-3, textvariable=var_pntg, fg="white").place(in_=tool3_bar, relx=0.50, rely=0.60, anchor=W)

        var_rec.set('')                                                                    # set default path

        var_dist.set('')                                                                   # set default path
        txt_dist_box = Entry(tool2l_bar, width=len_box-8, textvariable=var_dist, bg="white", fg="black").place(in_=tool2l_bar, relx=0.70, rely=0.35, anchor=W)

        var_distc.set('')                                                                  # set default path
        var_distm.set('')                                                                  # set default path
        var_distmc.set('')                                                                 # set default path
        var_point_ra.set('')                                                               # set default path
        var_point_dec.set('')                                                              # set default path
        var_point_az.set('')                                                               # set default path
        var_point_el.set('')                                                               # set default path

        canvas_stat2 = tk.Canvas(tool4l_bar, width=wi2, height=he2, borderwidth=0, bg='SkyBlue1', highlightbackground = 'SkyBlue1')  # Create 200x200 Canvas widget
        canvas_stat2.place(relx=0.215, rely=0.502)
        oval_stat2 = canvas_stat2.create_oval(wi2/2-r2, he2/2+r2, wi2/2+r2, he2/2-r2)      # Create a circle on the Canvas

        var_rfi.set('')                                                                    # set default path
        txt_rfi_box = Entry(tool4l_bar, width=len_box+22, textvariable=var_rfi, fg="white").place(in_=tool4l_bar, relx=0.265, rely=0.65, anchor=W)

        # Box - top
        var_name_source.set('')                                                            # set default path
        var_point_az_c.set('')                                                             # set default path
        var_point_el_c.set('')                                                             # set default path
        #var_point_ra_c.set('')                                                             # set default path
        #var_point_dec_c.set('')                                                            # set default path

        radar_frame = radar_f(root)

        data_input.update({"mem_status": (*data_input['mem_status'],data_input['Status'])})
        data_input.update({"mem_status": data_input['mem_status'][-2:]})

        if len(data_input['mem_status']) == 1:
            play(sound_voice(file_disconnected))
            log_action('DISCOS is disconnected!','red')
        else:
            if data_input['mem_status'][-1] != data_input['mem_status'][-2]:
                play(sound_voice(file_disconnected))
                log_action('DISCOS is disconnected!','red')
        
    elif stato == 1 and running == True:
        arr, arr_dict, arr2 = discos_pars(s)
        #print(str(arr))

        #log_action(str(arr),'black')

        radar_frame = radar_f(root)
        
        if len(arr) == 1:
            label_text = 'PLEASE WAIT'
            data_input.update({"Status": label_text})
            my_canvas.itemconfig(my_oval, fill="yellow")                                   # Fill the circle with YELLOW

            data_input.update({"mem_status": (*data_input['mem_status'],data_input['Status'])})
            data_input.update({"mem_status": data_input['mem_status'][-2:]})

            log_action('Waiting for data from DISCOS (runtime?) ...','yellow')
            
        else:
            label_text = 'CONNECTED'
            data_input.update({"Status": label_text})

            my_canvas.itemconfig(my_oval, fill="green")                                    # Fill the circle with GREEN

            data_input.update({"mem_status": (*data_input['mem_status'],data_input['Status'])})
            data_input.update({"mem_status": data_input['mem_status'][-2:]})

            if len(data_input['mem_status']) == 1:
                play(sound_voice(file_connected))
                log_action('DISCOS is connected!','green')
            else:
                if data_input['mem_status'][-1] != data_input['mem_status'][-2]:
                    play(sound_voice(file_connected))
                    log_action('DISCOS is connected!','green')
                
            time_utc = str(arr_dict['SysUTC'])            
            time_arr = loc_time(time_utc, radiotelescope)

            receiver_code = str(arr_dict['ReceiverCode'])
            name_source = str(arr_dict['SourceNameMn'])
            
            #arr_dict['Azimuth'] = '140'    # TEST
            #arr_dict['Elevation'] = '37'   # TEST

            loc_site = time_arr[2]
            altaz = AltAz(obstime=time_arr[3], location=loc_site)

            az0_live, el0_live = arr_dict['Azimuth'], arr_dict['Elevation']                # pointing position (altaz)
            az_live, el_live = float(az0_live), float(el0_live)
            pos_live = (az_live, el_live)
            tool_live = coord_altaz2radec(time_arr[3], (pos_live[1], pos_live[0]), 0, radiotelescope, 0)
            pos_live_radec = (tool_live[5].ra.value, tool_live[5].dec.value)
            tool_live_radec, tool_live_altaz = tool_live[5], tool_live[6]                  # SkyCoord coordinates

            ra0_source, dec0_source = arr_dict['RightAscension'], arr_dict['Declination']  # pointing position (radec)
            radec_string = ra0_source + ' ' + dec0_source
            coord_string = SkyCoord(radec_string, unit=(u.hourangle, u.deg))
            ra_source, dec_source = coord_string.ra.value, coord_string.dec.value
            pos_source = (ra_source, dec_source)
            tool_source_radec = SkyCoord(pos_source[0], pos_source[1], frame='icrs', unit='deg')
            tool_source_altaz = tool_source_radec.transform_to(altaz)

            az_comm, el_comm = arr_dict['CommandedAzimuth'], arr_dict['CommandedElevation']# source position
            pos_comm = (float(az_comm), float(el_comm))
            tool_comm = coord_altaz2radec(time_arr[3], (pos_comm[1], pos_comm[0]), 0, radiotelescope, 0)
            pos_comm_radec = (tool_comm[5].ra.value, tool_comm[5].dec.value)
            tool_comm_radec, tool_comm_altaz = tool_comm[5], tool_comm[6]                  # SkyCoord coordinates

            #dist_sun_pnt, dist_sun_source, dist_moon_pnt, dist_moon_source, tool = calc_angdist_radec(time_utc, radiotelescope, tool_live_radec, tool_source_radec, tot_eph_h, step_eph_min, q_mode)
            #test = calc_angdist(time_utc, 'SRT', az0_live, el0_live, az_comm, el_comm)
            #dist_sun_now, dist_sun_com, dist_moon_now, dist_moon_com = test[0], test[1], test[2], test[3]
            dist_sun_pnt, dist_sun_source, dist_moon_pnt, dist_moon_source, tool = calc_angdist_altaz(time_utc, radiotelescope, tool_live_altaz, tool_comm_altaz, tot_eph_h, step_eph_min, q_mode)

            if q_mode == 2:
                # HQ-mode
                result_ref_sunmoon = tool[6]
                pos_sun, pos_moon = (tool[3].az.value, tool[3].alt.value), (tool[5].az.value, tool[5].alt.value)
                #sun_eph, moon_eph, sun_eph_radec, sun_eph_altaz, moon_eph_radec, moon_eph_altaz, result_sunmoon = horizon_eph(tot_eph_h, time_arr[3], radiotelescope, step_eph_min)
                sun_eph, moon_eph, sun_eph_radec, sun_eph_altaz, moon_eph_radec, moon_eph_altaz, result_sunmoon = skyfield_eph(tot_eph_h, time_arr[3], radiotelescope, step_eph_min)
            elif q_mode == 1:
                # LQ-mode
                result_ref_sunmoon = tool[6]
                pos_sun, pos_moon = (result_ref_sunmoon[1].value, result_ref_sunmoon[2].value), (result_ref_sunmoon[3].value, result_ref_sunmoon[4].value)
                sun_eph, moon_eph = 0, 0
                sun_eph_radec, sun_eph_altaz, moon_eph_radec, moon_eph_altaz, result_sunmoon = tool[0], tool[1], tool[2], tool[3], tool[4]

            try:
                rfi_yes, rfi_source, rfi_tab, rfi_freq = tool_rfi(az_live, el_live, receiver_code, file_rfi, file_rec)
            except:
                rfi_yes, rfi_source, rfi_tab, rfi_freq = 0, 'clean', 'no_tab', 'no_tab'

            try:
                evo_az_source, evo_el_source, source_timevo, timevo_min, timevo_max, evo_radec, evo_altaz = coord_altaz2radec(time_arr[3], (el_comm, az_comm), result_sunmoon[:,0], radiotelescope, 1)
            except:
                evo_az_source, evo_el_source, source_timevo, timevo_min, timevo_max, evo_radec, evo_altaz = 0., 0., 0., 0., 0., 0., 0.
            pos_source_evo = (evo_az_source, evo_el_source)

            if data_input['mode'] == 'NO SOLAR - other':
                plot_radar(file_sys, time_arr[3], name_source, pos_live, pos_comm, result_ref_sunmoon, result_sunmoon, pos_source_evo, source_timevo, timevo_min, timevo_max, 0, 0, rfi_tab, rfi_freq, 70, vis_mode)
                if (dist_sun_pnt.value <= 40) & (dist_sun_pnt.value >= 10):
                    txt_dist_box = Entry(tool2l_bar, width=len_box-8, textvariable=var_dist, bg="yellow", fg="black").place(in_=tool2l_bar, relx=0.70, rely=0.35, anchor=W)
                    log_action('WARNING: The distance SRT-SUN is between 5 and 10 deg!','yellow')
                    play(sound_alarm(file_alarm, 2))
                elif (dist_sun_pnt.value < 10):
                    log_action('ALERT: The distance SRT-SUN is less than 5 deg!','red')
                    txt_dist_box = Entry(tool2l_bar, width=len_box-8, textvariable=var_dist, bg="red", fg="white").place(in_=tool2l_bar, relx=0.70, rely=0.35, anchor=W)
                    play(sound_alarm(file_alarm, 2))
                else:
                    txt_dist_box = Entry(tool2l_bar, width=len_box-8, textvariable=var_dist, bg="white", fg="black").place(in_=tool2l_bar, relx=0.70, rely=0.35, anchor=W)
            elif data_input['mode'] == 'SOLAR           ':
                #plot_radar_sun(file_sys, pos_live, pos_sun, pos_moon, str(time_str), receiver_code, 70)
                plot_radar(file_sys, time_arr[3], name_source, pos_live, 0, result_ref_sunmoon, result_sunmoon, 0, 0, timevo_min, timevo_max, 1, (pos_sun, pos_moon, str(time_arr[0]), receiver_code), rfi_tab, rfi_freq, 70, vis_mode)

            def callback_gen():
                #zoom1, zoom2 = 1.47, 0.832                                                 # 70 dpi
                zoom1, zoom2 = 0.602, 0.567                                                # 100 dpi
                album_data = os.path.join('./utilities/', 'radar.png')
                immagine = Image.open(album_data)                                          # 640 x 380 pixel
                immagine = immagine.resize((int(zoom1*immagine.size[0]), int(zoom2*immagine.size[1])))
                img = ImageTk.PhotoImage(immagine, master=radar_frame)
                panel = Label(radar_frame, image=img)
                panel.place(in_=radar_frame, relx=0.015, rely=0.015)
                panel.configure(image=img)
                panel.image = img

            def callback_sun():
                #zoom = 0.975
                zoom1, zoom2 = 0.562, 0.567                                                # 100 dpi
                album_data = os.path.join('./utilities/', 'radar.png')
                immagine = Image.open(album_data)                                          # 640 x 380 pixel
                immagine = immagine.resize((int(zoom1*immagine.size[0]), int(zoom2*immagine.size[1])))
                img = ImageTk.PhotoImage(immagine, master=radar_frame)
                panel = Label(radar_frame, image=img)
                panel.place(in_=radar_frame, relx=0.015, rely=0.015)
                panel.configure(image=img)
                panel.image = img

            if data_input['mode'] == 'NO SOLAR - other':
                radar_frame.bind("<Return>", callback_gen())
            elif data_input['mode'] == 'SOLAR           ':
                radar_frame.bind("<Return>", callback_sun())
            
            # Box - left
            var_time.set(time_utc)                                                         # set default path
            var_time_lst.set(time_arr[4])                                                  # set default path
            var_stat.set(arr_dict['SystemStatusMn'])                                       # set default path

            wi2, he2, r2 = 15, 15, 7
            canvas_stat = tk.Canvas(tool3_bar, width=wi2, height=he2, borderwidth=0, bg='SkyBlue1', highlightbackground = 'SkyBlue1')   # Create 200x200 Canvas widget
            canvas_stat.place(relx=0.42, rely=0.33)
            oval_stat = canvas_stat.create_oval(wi2/2-r2, he2/2+r2, wi2/2+r2, he2/2-r2)    # Create a circle on the Canvas
            
            if arr_dict['SystemStatusMn'] == 'OK':
                sys_status = 'OK'
                canvas_stat.itemconfig(oval_stat, fill="green")                            # Fill the circle with GREEN
                txt_stat_box = Entry(tool3_bar, width=len_box-3, textvariable=var_stat, bg="green", fg="white").place(in_=tool3_bar, relx=0.50, rely=0.40, anchor=W)
            if arr_dict['SystemStatusMn'] == 'WARNING':
                sys_status = 'WARNING'
                canvas_stat.itemconfig(oval_stat, fill="yellow")                           # Fill the circle with YELLOW
                txt_stat_box = Entry(tool3_bar, width=len_box-3, textvariable=var_stat, bg="yellow", fg="black").place(in_=tool3_bar, relx=0.50, rely=0.40, anchor=W)
            if arr_dict['SystemStatusMn'] == 'FAILURE':
                sys_status = 'FAILURE'
                canvas_stat.itemconfig(oval_stat, fill="red")                              # Fill the circle with RED
                txt_stat_box = Entry(tool3_bar, width=len_box-3, textvariable=var_stat, bg="red", fg="white").place(in_=tool3_bar, relx=0.50, rely=0.40, anchor=W)

            data_input.update({"mem_sys": (*data_input['mem_sys'],sys_status)})
            data_input.update({"mem_sys": data_input['mem_sys'][-2:]})

            if sys_status == 'FAILURE':
                if len(data_input['mem_sys']) == 1:
                    if data_input['mem_sys'][0] == 'FAILURE':
                        play(sound_voice(file_failure))
                        log_action('DISCOS system is in FAILURE status!','red')
                else:
                    if data_input['mem_sys'][-1] != data_input['mem_sys'][-2]:
                        if data_input['mem_sys'][-1] == 'FAILURE':
                            play(sound_voice(file_failure))
                            log_action('DISCOS system is in FAILURE status!','red')

            var_pntg.set(arr_dict['PointingStatusMn'])                                     # set default path
            if arr_dict['PointingStatusMn'] == 'TRACKING':
                txt_pntg_box = Entry(tool3_bar, width=len_box-3, textvariable=var_pntg, bg="green", fg="white").place(in_=tool3_bar, relx=0.50, rely=0.60, anchor=W)
            if arr_dict['PointingStatusMn'] == 'SLEWING':
                txt_pntg_box = Entry(tool3_bar, width=len_box-3, textvariable=var_pntg, bg="yellow", fg="black").place(in_=tool3_bar, relx=0.50, rely=0.60, anchor=W)

            var_rec.set(receiver_code)                                                     # set default path

            if len(receiver_code) == 0:
                txt_rec_box = Entry(tool3_bar, width=len_box-3, textvariable=var_rec, bg="yellow").place(in_=tool3_bar, relx=0.50, rely=0.80, anchor=W)
            else:
                txt_rec_box = Entry(tool3_bar, width=len_box-3, textvariable=var_rec).place(in_=tool3_bar, relx=0.50, rely=0.80, anchor=W)

            var_dist.set(str('%.4f' % dist_sun_pnt.value))                                 # set default path
            var_distc.set(str('%.4f' % dist_sun_source.value))                             # set default path
            var_distm.set(str('%.4f' % dist_moon_pnt.value))                               # set default path
            var_distmc.set(str('%.4f' % dist_moon_source.value))                           # set default path
            var_point_ra.set(ra0_source)                                                   # set default path
            var_point_dec.set(dec0_source)                                                 # set default path
            var_point_az.set(az0_live)                                                     # set default path
            var_point_el.set(el0_live)                                                     # set default path

            #print(ra0_source)                                                   # set default path
            #print(dec0_source)                                                 # set default path
            #print(az0_live)                                                     # set default path
            #print(el0_live)                                                     # set default path

            canvas_stat2 = tk.Canvas(tool4l_bar, width=wi2, height=he2, borderwidth=0, bg='SkyBlue1', highlightbackground = 'SkyBlue1')  # Create 200x200 Canvas widget
            canvas_stat2.place(relx=0.215, rely=0.502)
            oval_stat2 = canvas_stat2.create_oval(wi2/2-r2, he2/2+r2, wi2/2+r2, he2/2-r2)  # Create a circle on the Canvas

            var_rfi.set(rfi_source)                                                        # set default path
            if rfi_source[0] == '':
                txt_rfi_box = Entry(tool4l_bar, width=len_box+22, textvariable='', bg="green").place(in_=tool4l_bar, relx=0.265, rely=0.65, anchor=W)
                canvas_stat2.itemconfig(oval_stat2, fill="green")                          # Fill the circle with GREEN
            else:
                txt_rfi_box = Entry(tool4l_bar, width=len_box+22, textvariable=var_rfi, bg="yellow").place(in_=tool4l_bar, relx=0.265, rely=0.65, anchor=W)
                canvas_stat2.itemconfig(oval_stat2, fill="yellow")                         # Fill the circle with YELLOW
            
            # Box - top
            var_name_source.set(arr_dict['SourceNameMn'])                                  # set default path
            var_point_az_c.set(az_comm)                                                    # set default path
            var_point_el_c.set(el_comm )                                                   # set default path
            #var_point_ra_c.set(ra0_source)                                                 # set default path
            #var_point_dec_c.set(dec0_source)                                               # set default path
    if running:
        root.after(delta_t, loop, s)                                                       # function's name without ()


running = True


# Creation of the widget
root = Tk()                                                                                # Set Tk instance
root.title("Showtel v 0.4")
root.geometry("1128x810")                                                                  # Set the starting size of the window
root.maxsize(1500, 1000)                                                                   # width x height
root.config(bg="DeepSkyBlue3")

data_input = {}
data_input.update({"app": 'DISCOS', "Status": "", "mode": "NO SOLAR - other", "visibility": "SOUTH", "mem_status": [], "mem_sys": [], "quality": 'LQ'})

stile = "Verdana"
dimensione = 10
len_box = 22

img0 = Image.open('./utilities/logo_showtel.png')
zoom = 0.107
img = Image.open('./utilities/logo_showtel.png').resize((int(zoom*img0.size[0]),int(zoom*img0.size[1])))
logo = ImageTk.PhotoImage(img)
panel = Label(root, image = logo)
panel.place(x=9, y=5, relx=0.005, rely=0.01)

copyright_box = Label(root, text="Developed by Dr. Marco Marongiu @ INAF/OAC (Italy) - Last version: " + update, foreground="red", font=(stile, dimensione-4)).place(in_=root, relx=0.007, rely=0.983, anchor=W)
#############################################


# Create DISCOS switch ######################
raggio = 20
wi, he = 130, 90

my_canvas = tk.Canvas(root, width=wi, height=he, bg='DeepSkyBlue1')                        # Create 200x200 Canvas widget
my_canvas.place(x=2, y=52, relx=0.005, rely=0.01)

my_oval = my_canvas.create_oval(wi/1.3-raggio, he/1.85+raggio-20, wi/1.3+raggio, he/1.85-raggio-20)       # Create a circle on the Canvas

button_start = Button(my_canvas, text="START", command=start_butt, fg="green", bg="azure").place(in_=my_canvas, relx=0.3, rely=0.22, anchor=CENTER) # Set a "Start button". The "start_indicators" function is a call-back

button_stop = Button(my_canvas, text="STOP", command=stop_butt, fg="red", bg="azure").place(in_=my_canvas, relx=0.3, rely=0.48, anchor=CENTER)      # Set a "Start button". The "start_indicators" function is a call-back


# Select app to catch the antenna parameters
def change_app(choice):
    choice = select_var.get()
    data_input.update({"app": choice})
    if choice == 'DISCOS':
        opt_menu.config(bg="azure", fg="black", activebackground="azure", activeforeground="black")
        log_action('DISCOS input is selected','brown')
    elif choice == 'Seadas':
        opt_menu.config(bg="turquoise", fg="BLACK", activebackground="turquoise", activeforeground="BLACK")
        log_action('Seadas input is selected','brown')

    return choice

select = ["DISCOS", "Seadas"]
select_var = StringVar()
select_var.set(select[0])

opt_menu = OptionMenu(my_canvas, select_var, *select, command=change_app)
opt_menu.place(in_=my_canvas, relx=0.08, rely=0.63)
opt_menu.config(bg="azure", fg="black", activebackground="azure", activeforeground="black")
#############################################


# Create frames #############################
def left_f1(rr):
    wil, hel = 463, 90
    left_fr = Frame(rr, width=wil, height=hel, bg='DeepSkyBlue1')
    left_fr.place(x=1, y=154, relx=0.005, rely=0.01)
    return left_fr

def left_f2(rr):
    wil, hel = 463, 152
    left_fr = Frame(rr, width=wil, height=hel, bg='DeepSkyBlue1')
    left_fr.place(x=1, y=254, relx=0.005, rely=0.01)
    return left_fr

def left_f3(rr):
    wil, hel = 463, 122
    left_fr = Frame(rr, width=wil, height=hel, bg='DeepSkyBlue1')
    left_fr.place(x=1, y=416, relx=0.005, rely=0.01)
    return left_fr

def left_f4(rr):
    wil, hel = 463, 70
    left_fr = Frame(rr, width=wil, height=hel, bg='DeepSkyBlue1')
    left_fr.place(x=1, y=548, relx=0.005, rely=0.01)
    return left_fr

def rightup_f1(rr):
    wir, her = 320, 140
    rightup_fr = Frame(rr, width=wir, height=her, bg='DeepSkyBlue1')
    rightup_fr.place(x=138, y=5, relx=0.01, rely=0.01)
    return rightup_fr

def rightup_f2(rr):
    wir, her = 640, 140
    rightup_fr = Frame(rr, width=wir, height=her, bg='DeepSkyBlue1')
    rightup_fr.place(x=468, y=5, relx=0.01, rely=0.01)
    return rightup_fr

def radar_f(rr):
    wirad, herad = 640, 464
    radar_fr = Frame(rr, width=wirad, height=herad, bg='DeepSkyBlue1')
    radar_fr.place(x=468, y=154, relx=0.01, rely=0.01)
    return radar_fr

left_frame1 = left_f1(root)
left_frame2 = left_f2(root)
left_frame3 = left_f3(root)
left_frame4 = left_f4(root)
rightup_frame1 = rightup_f1(root)
rightup_frame2 = rightup_f2(root)
radar_frame = radar_f(root)
#############################################


# Create boxes and log-list #################
tool1l_bar = Frame(left_frame1, width=444, height=75, bg='SkyBlue1')
tool1l_bar.place(x=10, rely=0.075)

tool2l_bar = Frame(left_frame2, width=444, height=137, bg='SkyBlue1')
tool2l_bar.place(x=10, rely=0.05)

tool3l_bar = Frame(left_frame3, width=444, height=107, bg='SkyBlue1')
tool3l_bar.place(x=10, rely=0.06)

tool4l_bar = Frame(left_frame4, width=444, height=55, bg='SkyBlue1')
tool4l_bar.place(x=10, rely=0.1)

tool2_bar = Frame(rightup_frame2, width=620, height=125, bg='SkyBlue1')
tool2_bar.place(x=10, rely=0.05)

tool3_bar = Frame(rightup_frame1, width=300, height=125, bg='SkyBlue1')
tool3_bar.place(x=10, rely=0.05)

listbox_loglist = Listbox(root, width=158, height=11, selectmode='extended', bg='grey80')
listbox_loglist.place(x=8, y=634)
scrollbar_loglist_x = Scrollbar(root, orient='horizontal', command=listbox_loglist.xview)
scrollbar_loglist_x.place(x=8, y=768, width=1107, height=13)
scrollbar_loglist_y = Scrollbar(root, orient='vertical', command=listbox_loglist.yview)
scrollbar_loglist_y.place(x=1105, y=633, height=148)
listbox_loglist['xscrollcommand'] = scrollbar_loglist_x.set
listbox_loglist['yscrollcommand'] = scrollbar_loglist_y.set


def get_current_datetime_formatted():
    return datetime.strftime(datetime.now(), "[%Y/%m/%d %H:%M:%S] - ")


def log_action(msg,color):
    listbox_loglist.insert(0, get_current_datetime_formatted() + msg)                      # se al posto di 0 metto tk.END, va come suggerito da Alessandro C., ma poi l'utente deve necessariamente scorrere lo scrollbar per vedere i messaggi più recenti
    listbox_loglist.itemconfig(0, {'fg':color})
#############################################


# Switch to select the solar session mode
def change_solar(choice_sun):
    choice_sun = mode_var.get()
    data_input.update({"mode": choice_sun})
    if choice_sun == 'SOLAR           ':
        mode_menu.config(bg="brown", fg="white", activebackground="brown", activeforeground="white")
        log_action('"SOLAR" session mode is selected','brown')
    elif choice_sun == 'NO SOLAR - other':
        mode_menu.config(bg="turquoise", fg="black", activebackground="turquoise", activeforeground="black")
        log_action('"NO SOLAR - other" session mode is selected','brown')

    return choice_sun

mode = ["SOLAR           ", "NO SOLAR - other"]
mode_var = StringVar()
mode_var.set(mode[1])

mode_menu = OptionMenu(tool2_bar, mode_var, *mode, command=change_solar)
mode_menu.place(in_=tool2_bar, relx=0.73, rely=0.04)
mode_menu.config(bg="turquoise", fg="black", activebackground="turquoise", activeforeground="black")
#############################################

log_action('DISCOS input is selected by default','brown')
log_action('"NO SOLAR - other" session mode is selected by default','brown')


# Switch to select the sky visibility mode
def change_visibility(choice_visibility):
    choice_visibility = visibility_var.get()
    data_input.update({"visibility": choice_visibility})
    if choice_visibility == 'NORTH':
        visibility_menu.config(bg="brown", fg="white", activebackground="brown", activeforeground="white")
        log_action('North sky visibility mode is selected','brown')
    elif choice_visibility == 'SOUTH':
        visibility_menu.config(bg="white", fg="brown", activebackground="white", activeforeground="brown")
        log_action('South sky visibility mode is selected','brown')

    return choice_visibility

visibility = ["NORTH", "SOUTH"]
visibility_var = StringVar()
visibility_var.set(visibility[1])

visibility_menu = OptionMenu(tool2_bar, visibility_var, *visibility, command=change_visibility)
visibility_menu.place(in_=tool2_bar, relx=0.855, rely=0.32)
visibility_menu.config(bg="white", fg="brown", activebackground="white", activeforeground="brown")
#############################################

log_action('South sky visibility mode is selected by default','brown')


# Switch to select the quality of the measurements
def change_quality(choice_quality):
    choice_quality = quality_var.get()
    data_input.update({"quality": choice_quality})
    if choice_quality == 'LQ':
        quality_menu.config(bg="brown", fg="white", activebackground="brown", activeforeground="white")
        log_action('"astropy/JPL-NASA [LQ]" mode is selected','brown')
    elif choice_quality == 'HQ':
        quality_menu.config(bg="white", fg="brown", activebackground="white", activeforeground="brown")
        log_action('"skyfield/JPL-NASA [HQ]" mode is selected','brown')

    return choice_quality

quality = ['LQ', 'HQ']
quality_var = StringVar()
quality_var.set(quality[0])

quality_menu = OptionMenu(tool2_bar, quality_var, *quality, command=change_quality)
quality_menu.place(in_=tool2_bar, relx=0.73, rely=0.32)
quality_menu.config(bg="white", fg="brown", activebackground="white", activeforeground="brown")
#############################################

log_action('"astropy/JPL-NASA [LQ]" mode is selected by default','brown')


# Up - left #################################
source_title_box = Label(tool3_bar, text="Status Panel", font=(stile, dimensione+5), fg='Blue', bg='SkyBlue1').place(in_=tool3_bar, relx=0.04, rely=0.15, anchor=W)

var_stat = StringVar()
var_stat.set('')                                                                           # set default path
stat_box = Label(tool3_bar, text="System Status   ", font=(stile, dimensione), bg='SkyBlue1').place(in_=tool3_bar, relx=0.03, rely=0.40, anchor=W)
txt_stat_box = Entry(tool3_bar, width=len_box-3, textvariable=var_stat).place(in_=tool3_bar, relx=0.50, rely=0.40, anchor=W)

wi2, he2, r2 = 15, 15, 7
canvas_stat = tk.Canvas(tool3_bar, width=wi2, height=he2, borderwidth=0, bg='SkyBlue1', highlightbackground = 'SkyBlue1')   # Create 200x200 Canvas widget
canvas_stat.place(relx=0.42, rely=0.33)
oval_stat = canvas_stat.create_oval(wi2/2-r2, he2/2+r2, wi2/2+r2, he2/2-r2)                # Create a circle on the Canvas

var_pntg = StringVar()
var_pntg.set('')                                                                           # set default path
pntg_box = Label(tool3_bar, text="Pointing Status ", font=(stile, dimensione), bg='SkyBlue1').place(in_=tool3_bar, relx=0.03, rely=0.60, anchor=W)
txt_pntg_box = Entry(tool3_bar, width=len_box-3, textvariable=var_pntg).place(in_=tool3_bar, relx=0.50, rely=0.60, anchor=W)

var_rec = StringVar()
var_rec.set('')                                                                            # set default path
rec_box = Label(tool3_bar, text="Receiver        ", font=(stile, dimensione), bg='SkyBlue1').place(in_=tool3_bar, relx=0.03, rely=0.80, anchor=W)
txt_rec_box = Entry(tool3_bar, width=len_box-3, textvariable=var_rec).place(in_=tool3_bar, relx=0.50, rely=0.80, anchor=W)
#############################################


# Box - left 1 ##############################
source_title_box = Label(tool1l_bar, text="Time Panel", font=(stile, dimensione+5), fg='Blue', bg='SkyBlue1').place(in_=tool1l_bar, relx=0.04, rely=0.15, anchor=W)

var_time = StringVar()
var_time.set('')                                                                           # set default path
time_box = Label(tool1l_bar, text="UTC epoch       ", font=(stile, dimensione), bg='SkyBlue1').place(in_=tool1l_bar, relx=0.03, rely=0.45, anchor=W)
txt_path_box = Entry(tool1l_bar, width=len_box+8, textvariable=var_time).place(in_=tool1l_bar, relx=0.45, rely=0.45, anchor=W)

var_time_lst = StringVar()
var_time_lst.set('')                                                                       # set default path
time_lst_box = Label(tool1l_bar, text="LST time        ", font=(stile, dimensione), bg='SkyBlue1').place(in_=tool1l_bar, relx=0.03, rely=0.77, anchor=W)
txt_time_lst_box = Entry(tool1l_bar, width=len_box+8, textvariable=var_time_lst).place(in_=tool1l_bar, relx=0.45, rely=0.77, anchor=W)
#############################################


# Box - left 2 ##############################
source_title_box = Label(tool2l_bar, text="Distance Panel", font=(stile, dimensione+5), fg='Blue', bg='SkyBlue1').place(in_=tool2l_bar, relx=0.04, rely=0.15, anchor=W)

var_dist = StringVar()
var_dist.set('')                                                                           # set default path
dist_box = Label(tool2l_bar,  text="SUN/telescope (live, deg) ", font=(stile, dimensione), bg='SkyBlue1').place(in_=tool2l_bar, relx=0.03, rely=0.35, anchor=W)
txt_dist_box = Entry(tool2l_bar, width=len_box-8, textvariable=var_dist).place(in_=tool2l_bar, relx=0.70, rely=0.35, anchor=W)

var_distc = StringVar()
var_distc.set('')                                                                          # set default path
distc_box = Label(tool2l_bar, text="SUN/source (live, deg)    ", font=(stile, dimensione), bg='SkyBlue1').place(in_=tool2l_bar, relx=0.03, rely=0.52, anchor=W)
txt_distc_box = Entry(tool2l_bar, width=len_box-8, textvariable=var_distc).place(in_=tool2l_bar, relx=0.70, rely=0.52, anchor=W)

var_distm = StringVar()
var_distm.set('')                                                                          # set default path
distm_box = Label(tool2l_bar,  text="MOON/telescope (live, deg)", font=(stile, dimensione), bg='SkyBlue1').place(in_=tool2l_bar, relx=0.03, rely=0.69, anchor=W)
txt_distm_box = Entry(tool2l_bar, width=len_box-8, textvariable=var_distm).place(in_=tool2l_bar, relx=0.70, rely=0.69, anchor=W)

var_distmc = StringVar()
var_distmc.set('')                                                                         # set default path
distmc_box = Label(tool2l_bar, text="MOON/source (live, deg)   ", font=(stile, dimensione), bg='SkyBlue1').place(in_=tool2l_bar, relx=0.03, rely=0.86, anchor=W)
txt_distmc_box = Entry(tool2l_bar, width=len_box-8, textvariable=var_distmc).place(in_=tool2l_bar, relx=0.70, rely=0.86, anchor=W)
#############################################


# Box - left 3 ##############################
source_title_box = Label(tool3l_bar, text="Pointing Panel", font=(stile, dimensione+5), fg='Blue', bg='SkyBlue1').place(in_=tool3l_bar, relx=0.04, rely=0.15, anchor=W)

var_point_ra = StringVar()
var_point_ra.set('')                                                                       # set default path
point_ra_box = Label(tool3l_bar, text="RA (hms) ", font=(stile, dimensione), bg='SkyBlue1').place(in_=tool3l_bar, relx=0.03, rely=0.50, anchor=W)
txt_point_ra_box = Entry(tool3l_bar, width=len_box-6, textvariable=var_point_ra).place(in_=tool3l_bar, relx=0.21, rely=0.50, anchor=W)

var_point_dec = StringVar()
var_point_dec.set('')                                                                      # set default path
point_dec_box = Label(tool3l_bar, text="DEC (dms) ", font=(stile, dimensione), bg='SkyBlue1').place(in_=tool3l_bar, relx=0.53, rely=0.50, anchor=W)
txt_point_dec_box = Entry(tool3l_bar, width=len_box-6, textvariable=var_point_dec).place(in_=tool3l_bar, relx=0.71, rely=0.50, anchor=W)

var_point_az = StringVar()
var_point_az.set('')                                                                       # set default path
point_az_box = Label(tool3l_bar, text="Az (deg) ", font=(stile, dimensione), bg='SkyBlue1').place(in_=tool3l_bar, relx=0.03, rely=0.75, anchor=W)
txt_point_az_box = Entry(tool3l_bar, width=len_box-6, textvariable=var_point_az).place(in_=tool3l_bar, relx=0.21, rely=0.75, anchor=W)

var_point_el = StringVar()
var_point_el.set('')                                                                       # set default path
point_el_box = Label(tool3l_bar, text="El (deg) ", font=(stile, dimensione), bg='SkyBlue1').place(in_=tool3l_bar, relx=0.53, rely=0.75, anchor=W)
txt_point_el_box = Entry(tool3l_bar, width=len_box-6, textvariable=var_point_el).place(in_=tool3l_bar, relx=0.71, rely=0.75, anchor=W)
#############################################


# Box - left 4 ##############################
source_title_box = Label(tool4l_bar, text="RFI Panel", font=(stile, dimensione+5), fg='Blue', bg='SkyBlue1').place(in_=tool4l_bar, relx=0.04, rely=0.2, anchor=W)

canvas_stat2 = tk.Canvas(tool4l_bar, width=wi2, height=he2, borderwidth=0, bg='SkyBlue1', highlightbackground = 'SkyBlue1')  # Create 200x200 Canvas widget
canvas_stat2.place(relx=0.215, rely=0.502)
oval_stat2 = canvas_stat2.create_oval(wi2/2-r2, he2/2+r2, wi2/2+r2, he2/2-r2)              # Create a circle on the Canvas

var_rfi = StringVar()
var_rfi.set('')                                                                            # set default path
rfi_box = Label(tool4l_bar, text="RFI source", font=(stile, dimensione), bg='SkyBlue1').place(in_=tool4l_bar, relx=0.03, rely=0.65, anchor=W)
txt_rfi_box = Entry(tool4l_bar, width=len_box+22, textvariable=var_rfi).place(in_=tool4l_bar, relx=0.265, rely=0.65, anchor=W)
#############################################


# Box - top #################################
source_title_box = Label(tool2_bar, text="Source Panel", font=(stile, dimensione+5), fg='Blue', bg='SkyBlue1').place(in_=tool2_bar, relx=0.04, rely=0.15, anchor=W)

var_name_source = StringVar()
var_name_source.set('')                                                                    # set default path
name_source_box = Label(tool2_bar, text="Source Name ", font=(stile, dimensione), bg='SkyBlue1').place(in_=tool2_bar, relx=0.02, rely=0.40, anchor=W)
txt_name_source_box = Entry(tool2_bar, width=30, textvariable=var_name_source).place(in_=tool2_bar, relx=0.22, rely=0.40, anchor=W)

var_point_az_c = StringVar()
var_point_az_c.set('')                                                                     # set default path
point_az_c_box = Label(tool2_bar, text="Az (deg) ", font=(stile, dimensione), bg='SkyBlue1').place(in_=tool2_bar, relx=0.02, rely=0.75, anchor=W)
txt_point_az_c_box = Entry(tool2_bar, width=22, textvariable=var_point_az_c).place(in_=tool2_bar, relx=0.22, rely=0.75, anchor=W)

var_point_el_c = StringVar()
var_point_el_c.set('')                                                                     # set default path
point_el_c_box = Label(tool2_bar, text="El (deg) ", font=(stile, dimensione), bg='SkyBlue1').place(in_=tool2_bar, relx=0.52, rely=0.75, anchor=W)
txt_point_el_c_box = Entry(tool2_bar, width=22, textvariable=var_point_el_c).place(in_=tool2_bar, relx=0.72, rely=0.75, anchor=W)

#var_point_ra_c = StringVar()
#var_point_ra_c.set('')                                                                     # set default path
#point_ra_c_box = Label(tool2_bar, text="RA (hms) ", font=(stile, dimensione), bg='SkyBlue1').place(in_=tool2_bar, relx=0.02, rely=0.87, anchor=W)
#txt_point_ra_c_box = Entry(tool2_bar, width=22, textvariable=var_point_ra_c).place(in_=tool2_bar, relx=0.22, rely=0.87, anchor=W)

#var_point_dec_c = StringVar()
#var_point_dec_c.set('')                                                                    # set default path
#point_dec_c_box = Label(tool2_bar, text="DEC (dms) ", font=(stile, dimensione), bg='SkyBlue1').place(in_=tool2_bar, relx=0.52, rely=0.87, anchor=W)
#txt_point_dec_c_box = Entry(tool2_bar, width=22, textvariable=var_point_dec_c).place(in_=tool2_bar, relx=0.72, rely=0.87, anchor=W)

root.mainloop()
