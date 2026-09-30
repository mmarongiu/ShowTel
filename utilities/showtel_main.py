import numpy as np
from numpy import deg2rad
from astropy import units as u
from astropy.coordinates import Angle, angular_separation
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
from astropy.coordinates import solar_system_ephemeris, get_body_barycentric, get_body
from astropy.time import Time
from astropy.constants import au

from datetime import datetime, timedelta
import pytz, time, sched, threading, socket, os, pygame, gc

solar_system_ephemeris.set('jpl')

from tkinter import *
import tkinter as tk
from tkinter import ttk, font
from tkmacosx import Button

from PIL import ImageTk, Image
from gtts import gTTS

from .showtel_connection import *
from .showtel_discos import *
from .showtel_sound import *
from .showtel_distance import *
from .showtel_forbidden import *
from .showtel_plots import *

import warnings
warnings.filterwarnings("ignore")

# Global settings #############
update = '2026/06/04'
radiotelescope = 'SRT'                                                                     # SRT or Medicina
delta_time = 1                                                                             # in seconds

tot_eph_h = 24                                                                             # duration of the solar map [in units of hours, float]
step_eph_min = 30                                                                          # step time [in units of minutes, int]

s0 = connect_port_srtclient()
running = False
stato = 0
timestamps = []

data_input = {}
data_input.update({"app": 'SRTClient', "Status": "", "mode": "NO SOLAR - other", "visibility": "SOUTH", "mem_status": [], "mem_sys": [], "quality": 'AS'})
###############################


def start_butt():
    global running, stato
    log_action('Showtel is ON','black')
    if data_input['app'] == 'DISCOS':
        s, c, ss = on_app_discos()
    elif data_input['app'] == 'Seadas':
        s, c, ss = on_app_seadas()
    elif data_input['app'] == 'SRTClient':
        s = srtclient_pars(s0)

    running = True
    stato = 1
        
    loop(s)

    
def stop_butt():
    global running, stato
    log_action('Showtel is OFF','black')
    my_canvas.itemconfig(my_oval, fill="red")                                              # Fill the circle with RED
    if data_input['app'] == 'DISCOS':
        s, c, ss = on_app_discos()
        s.close()
    elif data_input['app'] == 'Seadas':
        s, c, ss = on_app_seadas()
        s.close()

    running = False
    stato = 0
    timestamps = []
    s0 = None

def loop(s):
    global timestamps, s0, stato
    
    stile = "Verdana"
    dimensione = 22
    dimensione_titolo = 25
    delta_t = int(delta_time*1e3)                                                          # in ms

    #file_sys = './utilities/'
    file_sys = os.path.dirname(os.path.abspath(__file__)) + '/'

    file_alarm = 'sound/warning_alarm.wav'
    file_connected = 'sound/connected.mp3'
    file_disconnected = 'sound/disconnected.mp3'
    file_failure = 'sound/failure.mp3'
    file_close = 'sound/sunclose.mp3'

    rfi_dir = os.path.dirname(os.path.abspath(__file__))
    file_rfi = os.path.join(rfi_dir, 'list_rfi_v1.dat')
    
    rec_dir = os.path.dirname(os.path.abspath(__file__))
    file_rec = os.path.join(rec_dir, 'list_receivers.dat')
    
    len_box = 22
    
    if data_input['app'] == 'DISCOS':
        stato = connect_status_discos(s)
    elif data_input['app'] == 'Seadas':
        stato = connect_status_seadas(s)
    elif data_input['app'] == 'SRTClient':
        s = srtclient_pars(s0)
        arr, arr_dict = s
        arr2 = None
        stato = connect_status_srtclient(arr_dict)

        timestamps.append(arr_dict['SysUTC'])
        timestamps = timestamps[-2:]

        if len(timestamps) == 2:
            if timestamps[-1] == timestamps[-2]:
                stato = 0

                try:
                    s0 = connect_port_srtclient()                                          # Re-open the socket/port
                    s = srtclient_pars(s0)                                                 # Re-read the parameters
                    arr, arr_dict = s
                    stato = connect_status_srtclient(arr_dict)
                except Exception as err_recon:
                    stato = 0

    if data_input["visibility"] == "SOUTH":
        vis_mode = 'S'
    elif data_input["visibility"] == "NORTH":
        vis_mode = 'N'

    if data_input["quality"] == "AS":
        q_mode = 1
    elif data_input["quality"] == "SK":
        q_mode = 2
        
    if stato == 0 or running == False:
        label_text = 'DISCONNECTED'
        data_input.update({"Status": label_text})
        my_canvas.itemconfig(my_oval, fill="red")                                          # Fill the circle with RED

        var_time.set('')
        var_time_lst.set('')
        var_stat.set('')
        var_pntg.set('')
        var_rec.set('')
        var_dist.set('')
        var_distc.set('')
        var_distm.set('')
        var_distmc.set('')
        var_point_ra.set('')
        var_point_dec.set('')
        var_point_az.set('')
        var_point_el.set('')
        var_rfi.set('')
        var_name_source.set('')
        var_point_az_c.set('')
        var_point_el_c.set('')

        txt_stat_box.config(bg="white")
        txt_rec_box.config(bg="white")
        txt_pntg_box.config(bg="white")
        txt_dist_box.config(bg="white")
        txt_rfi_box.config(bg="white")
        
        radar_panel.config(image='')
        radar_panel.image = None

        canvas_stat.itemconfig(oval_stat, fill="SkyBlue1")
        canvas_stat2.itemconfig(oval_stat2, fill="SkyBlue1")

        data_input.update({"mem_status": (*data_input['mem_status'],data_input['Status'])})
        data_input.update({"mem_status": data_input['mem_status'][-2:]})

        if len(data_input['mem_status']) == 1:
            play_audio(file_disconnected)
            log_action('DISCOS is disconnected!','red')
        else:
            if data_input['mem_status'][-1] != data_input['mem_status'][-2]:
                play_audio(file_disconnected)
                log_action('DISCOS is disconnected!','red')
        
    elif stato == 1 and running == True:
        is_srtclient = data_input['app'] == 'SRTClient'
        is_discos = data_input['app'] == 'DISCOS'
        is_seadas = data_input['app'] == 'Seadas'

        if is_discos or is_seadas:
            arr, arr_dict, arr2 = discos_pars(s)
        elif is_srtclient:
            arr_dict.update({'Azimuth': float(arr_dict['Azimuth']),
                             'Elevation': float(arr_dict['Elevation']),
                             'CommandedAzimuth': float(arr_dict['CommandedAzimuth']),
                             'CommandedElevation': float(arr_dict['CommandedElevation']),
                             'AzimuthCorrection': float(arr_dict['AzimuthCorrection']),
                             'AzimuthOffset': float(arr_dict['AzimuthOffset']),
                             'ElevationCorrection': float(arr_dict['ElevationCorrection']),
                             'ElevationOffset': float(arr_dict['ElevationOffset']),
                             'RefractionElevationCorrection': float(arr_dict['RefractionElevationCorrection']),
                             'RaOffset': float(arr_dict['RaOffset']),
                             'DeclOffset': float(arr_dict['DeclOffset']),
                             'GalLongitudeOffset': float(arr_dict['GalLongitudeOffset']),
                             'GalLatitudeOffset': float(arr_dict['GalLatitudeOffset']),
                             })
        
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
                play_audio(file_connected)
                log_action('DISCOS is connected!','green')
            else:
                if data_input['mem_status'][-1] != data_input['mem_status'][-2]:
                    play_audio(file_connected)
                    log_action('DISCOS is connected!','green')

            time_utc_str = str(arr_dict['SysUTC'])
            time_utc = time_utc_str[:-1] if is_srtclient else time_utc_str
            time_arr = loc_time(time_utc, radiotelescope, data_input['app'])

            receiver_code = str(arr_dict['ReceiverCode'])
            name_source = str(arr_dict['SourceName'])

            tab_ric = Table.read(file_rec, format='ascii')
            mask_ric = (receiver_code == tab_ric['nome'])
            tab_ric_select = tab_ric[mask_ric]

            if receiver_code != '':
                thres_critic_ric, thres_standard_ric = float(tab_ric_select['cr_ric'][0]), float(tab_ric_select['st_ric'][0])
            else:
                thres_critic_ric, thres_standard_ric = 0, 0

            #arr_dict['Azimuth'] = 345                  # TEST
            #arr_dict['Elevation'] = 18                 # TEST (mettere 12 per RFI)
            #arr_dict['CommandedAzimuth'] = -14         # TEST

            if arr_dict['CommandedAzimuth'] < 0.:
                arr_dict['CommandedAzimuth'] += 360.
            if arr_dict['Azimuth'] < 0.:
                arr_dict['Azimuth'] += 360.

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

            dist_sun_pnt, dist_sun_source, dist_moon_pnt, dist_moon_source, tool = calc_angdist_altaz(time_utc, radiotelescope, tool_live_altaz, tool_comm_altaz, tot_eph_h, step_eph_min, q_mode, data_input['app'])
            
            result_ref_sunmoon = tool[6]
            if q_mode == 2:
                # HQ-mode
                pos_sun, pos_moon = (tool[3].az.value, tool[3].alt.value), (tool[5].az.value, tool[5].alt.value)
                sun_eph, moon_eph, sun_eph_radec, sun_eph_altaz, moon_eph_radec, moon_eph_altaz, result_sunmoon = skyfield_eph(tot_eph_h, time_arr[3], radiotelescope, step_eph_min)
            elif q_mode == 1:
                # LQ-mode
                pos_sun, pos_moon = (result_ref_sunmoon[1].value, result_ref_sunmoon[2].value), (result_ref_sunmoon[3].value, result_ref_sunmoon[4].value)
                sun_eph, moon_eph = 0, 0
                sun_eph_radec, sun_eph_altaz, moon_eph_radec, moon_eph_altaz, result_sunmoon = tool[:5]

            try:
                rfi_yes, rfi_source, rfi_sp, rfi_tab, rfi_freq = tool_rfi(az_live, el_live, receiver_code, file_rfi, file_rec)
            except:
                rfi_yes, rfi_source, rfi_sp, rfi_tab, rfi_freq = 0, ['clean'], [''], 'no_tab', 'no_tab'

            rfi_sp = 'S' if str(rfi_sp[0]) == 'Y' else ''
            rfi_text = f"{rfi_source[0]} ({rfi_sp})" if rfi_yes == 1 else rfi_source
                
            try:
                evo_az_source, evo_el_source, source_timevo, timevo_min, timevo_max, evo_radec, evo_altaz = coord_altaz2radec(time_arr[3], (el_comm, az_comm), result_sunmoon[:,0], radiotelescope, 1)
            except:
                evo_az_source, evo_el_source, source_timevo, timevo_min, timevo_max, evo_radec, evo_altaz = 0., 0., 0., 0., 0., 0., 0.
            pos_source_evo = (evo_az_source, evo_el_source)

            if data_input['mode'] == 'NO SOLAR - other':
                plot_radar(file_sys, time_arr[3], name_source, pos_live, pos_comm, result_ref_sunmoon, result_sunmoon, pos_source_evo, source_timevo, timevo_min, timevo_max, 0, 0, rfi_tab, rfi_freq, 70, vis_mode)
                if (dist_sun_pnt.value <= thres_standard_ric) & (dist_sun_pnt.value >= thres_critic_ric):
                    txt_dist_box.config(bg="yellow", fg="black")
                    log_action('WARNING: The distance SRT-SUN is between 10 and 40 deg!','yellow')
                elif (dist_sun_pnt.value < thres_critic_ric):
                    log_action('ALERT: The distance SRT-SUN is less than 10 deg!','red')
                    txt_dist_box.config(bg="red", fg="white")
                    play_audio(file_alarm)
                else:
                    txt_dist_box.config(bg="white", fg="black")
            elif data_input['mode'] == 'SOLAR           ':
                plot_radar(file_sys, time_arr[3], name_source, pos_live, 0, result_ref_sunmoon, result_sunmoon, 0, 0, timevo_min, timevo_max, 1, (pos_sun, pos_moon, str(time_arr[0]), receiver_code), rfi_tab, rfi_freq, 70, vis_mode)

            try:
                album_data = os.path.join(file_sys, 'radar.png')
                if os.path.exists(album_data):
                    zoom1, zoom2 = (0.602, 0.567) if data_input['mode'] == 'NO SOLAR - other' else (0.562, 0.567)
            
                    immagine = Image.open(album_data)
                    immagine_res = immagine.resize((int(zoom1*immagine.size[0]), int(zoom2*immagine.size[1])))
                    
                    img = ImageTk.PhotoImage(immagine_res, master=radar_frame)
                    radar_panel.config(image=img)
                    radar_panel.image = img
            except Exception as e:
                print(f"Errore caricamento immagine: {e}")
            
            # Box - left
            var_time.set(time_utc)                                                         # set default path
            var_time_lst.set(time_arr[4])                                                  # set default path
            var_stat.set(arr_dict['SystemStatus'])                                         # set default path
  
            if arr_dict['SystemStatus'] == 'OK':
                sys_status = 'OK'
                canvas_stat.itemconfig(oval_stat, fill="green")                            # Fill the circle with GREEN
                txt_stat_box.config(bg="green", fg="white")
            if arr_dict['SystemStatus'] == 'WARNING':
                sys_status = 'WARNING'
                canvas_stat.itemconfig(oval_stat, fill="yellow")                           # Fill the circle with YELLOW
                txt_stat_box.config(bg="yellow", fg="black")
            if arr_dict['SystemStatus'] == 'FAILURE':
                sys_status = 'FAILURE'
                canvas_stat.itemconfig(oval_stat, fill="red")                              # Fill the circle with RED
                txt_stat_box.config(bg="red", fg="white")

            data_input.update({"mem_sys": (*data_input['mem_sys'],sys_status)})
            data_input.update({"mem_sys": data_input['mem_sys'][-2:]})

            if sys_status == 'FAILURE':
                if len(data_input['mem_sys']) == 1:
                    if data_input['mem_sys'][0] == 'FAILURE':
                        play_audio(file_failure)
                        log_action('DISCOS system is in FAILURE status!','red')
                else:
                    if data_input['mem_sys'][-1] != data_input['mem_sys'][-2]:
                        if data_input['mem_sys'][-1] == 'FAILURE':
                            play_audio(file_failure)
                            log_action('DISCOS system is in FAILURE status!','red')

            var_pntg.set(arr_dict['PointingStatus'])                                       # set default path
            if arr_dict['PointingStatus'] == 'TRACKING':
                txt_pntg_box.config(bg="green", fg="white")
            if arr_dict['PointingStatus'] == 'SLEWING':
                txt_pntg_box.config(bg="yellow", fg="black")

            var_rec.set(receiver_code)                                                     # set default path

            if len(receiver_code) == 0:
                txt_rec_box.config(bg="yellow")
            else:
                txt_rec_box.config(bg="green")

            var_dist.set(str('%.3f' % dist_sun_pnt.value))                                 # set default path
            var_distc.set(str('%.3f' % dist_sun_source.value))                             # set default path
            var_distm.set(str('%.3f' % dist_moon_pnt.value))                               # set default path
            var_distmc.set(str('%.3f' % dist_moon_source.value))                           # set default path
            var_point_ra.set(ra0_source)                                                   # set default path
            var_point_dec.set(dec0_source)                                                 # set default path
            var_point_az.set(f"{az0_live:.3f}" if is_srtclient else az0_live)              # set default path
            var_point_el.set(f"{el0_live:.3f}" if is_srtclient else el0_live)              # set default path

            var_rfi.set(rfi_text)                                                          # set default path
            if rfi_source[0] == '':
                txt_rfi_box.config(bg="green")
                canvas_stat2.itemconfig(oval_stat2, fill="green")                          # Fill the circle with GREEN
            else:
                txt_rfi_box.config(bg="yellow")
                canvas_stat2.itemconfig(oval_stat2, fill="yellow")                         # Fill the circle with YELLOW
            
            # Box - top
            var_name_source.set(arr_dict['SourceName'])                                    # set default path

            var_point_az_c.set(f"{az_comm:.3f}" if is_srtclient else az_comm)              # set default path
            var_point_el_c.set(f"{el_comm:.3f}" if is_srtclient else el_comm)              # set default path

    gc.collect()
    
    if running:
        root.after(delta_t, loop, s)                                                       # function's name without ()
    elif not running:
        return


def get_current_datetime_formatted():
    return datetime.strftime(datetime.now(), "[%Y/%m/%d %H:%M:%S] - ")


def log_action(msg,color):
    listbox_loglist.insert(0, get_current_datetime_formatted() + msg)                      # Using tk.END instead of 0, it works as suggested by Alessandro C., but then the user must scroll down to see the most recent messages
    listbox_loglist.itemconfig(0, {'fg':color})


def radar_f(rr):
    wirad, herad = 640, 464
    radar_fr = Frame(rr, width=wirad, height=herad, bg='DeepSkyBlue1')
    radar_fr.place(x=468, y=154, relx=0.01, rely=0.01)
    return radar_fr

 
def main():
    global root, my_canvas, my_oval, tool3_bar, tool2l_bar, tool4l_bar, tool2_bar, tool1l_bar, tool3l_bar
    global var_stat, var_pntg, var_rec, var_time, var_time_lst, var_dist, var_distc, var_distm, var_distmc
    global var_point_ra, var_point_dec, var_point_az, var_point_el, var_rfi, var_name_source, var_point_az_c, var_point_el_c
    global select_var, opt_menu, mode_var, mode_menu, visibility_var, visibility_menu, quality_var, quality_menu, listbox_loglist
    global txt_dist_box, txt_stat_box, txt_pntg_box, txt_rec_box, txt_rfi_box
    global radar_frame, radar_panel
    global canvas_stat, oval_stat, canvas_stat2, oval_stat2

    global running
    running = True

    # Creation of the widget
    root = Tk()                                                                            # Set Tk instance
    root.title("Showtel v 0.4")
    root.geometry("1128x810")                                                              # Set the starting size of the window
    root.maxsize(1500, 1000)                                                               # width x height
    root.config(bg="DeepSkyBlue3")

    stile = "Verdana"
    dimensione = 12
    dimensione_titolo = 15
    len_box = 22

    current_dir = os.path.dirname(os.path.abspath(__file__))
    logo_path = os.path.join(current_dir, 'logo_showtel.png')
    img0 = Image.open(logo_path)
    zoom = 0.107
    img = Image.open(logo_path).resize((int(zoom*img0.size[0]),int(zoom*img0.size[1])))
    logo = ImageTk.PhotoImage(img)
    panel = Label(root, image = logo)
    panel.place(x=9, y=5, relx=0.005, rely=0.01)

    copyright_box = Label(root, text="Developed by Dr. Marco Marongiu @ INAF/OAC (Italy) - Last version: " + update, foreground="red", font=(stile, dimensione-4)).place(in_=root, relx=0.007, rely=0.983, anchor=W)
    #############################################


    # Create DISCOS switch ######################
    raggio = 20
    wi, he = 130, 90

    my_canvas = tk.Canvas(root, width=wi, height=he, bg='DeepSkyBlue1')                    # Create 200x200 Canvas widget
    my_canvas.place(x=2, y=52, relx=0.005, rely=0.01)

    my_oval = my_canvas.create_oval(wi/1.3-raggio, he/1.85+raggio-20, wi/1.3+raggio, he/1.85-raggio-20)       # Create a circle on the Canvas

    button_start = Button(my_canvas, text="START", command=start_butt, fg="green", bg="azure", width=60).place(in_=my_canvas, relx=0.30, rely=0.18, anchor=CENTER) # Set a "Start button". The "start_indicators" function is a call-back

    button_stop = Button(my_canvas, text="STOP", command=stop_butt, fg="red", bg="azure", width=60).place(in_=my_canvas, relx=0.30, rely=0.44, anchor=CENTER)      # Set a "Start button". The "start_indicators" function is a call-back


    # Select app to catch the antenna parameters
    def change_app(choice):
        choice = select_var.get()
        data_input.update({"app": choice})
        if choice == 'SRTClient':
            opt_menu.config(bg="#007BA7", fg="black", activebackground="#007BA7", activeforeground="white")
            log_action('SRTClient input is selected','brown')
        elif choice == 'DISCOS':
            opt_menu.config(bg="azure", fg="black", activebackground="azure", activeforeground="black")
            log_action('DISCOS input is selected','brown')
        elif choice == 'Seadas':
            opt_menu.config(bg="turquoise", fg="black", activebackground="turquoise", activeforeground="black")
            log_action('Seadas input is selected','brown')

        return choice

    select = ["SRTClient", "DISCOS", "Seadas"]
    select_var = StringVar()
    select_var.set(select[0])

    opt_menu = OptionMenu(my_canvas, select_var, *select, command=change_app)
    opt_menu.place(in_=my_canvas, relx=0.08, rely=0.63)
    opt_menu.config(bg="#007BA7", fg="black", activebackground="#007BA7", activeforeground="white")
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
    radar_panel = Label(radar_frame, bg='DeepSkyBlue1')
    radar_panel.place(in_=radar_frame, relx=0.015, rely=0.015)
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
    scrollbar_loglist_x.place(x=8, y=797, width=1115, height=13)
    scrollbar_loglist_y = Scrollbar(root, orient='vertical', command=listbox_loglist.yview)
    scrollbar_loglist_y.place(x=1111, y=638, height=157)
    listbox_loglist['xscrollcommand'] = scrollbar_loglist_x.set
    listbox_loglist['yscrollcommand'] = scrollbar_loglist_y.set
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

    log_action('SRTClient input is selected by default','brown')
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
    visibility_menu.place(in_=tool2_bar, relx=0.843, rely=0.32)
    visibility_menu.config(bg="white", fg="brown", activebackground="white", activeforeground="brown")
    #############################################

    log_action('South sky visibility mode is selected by default','brown')


    # Switch to select the quality of the measurements
    def change_quality(choice_quality):
        choice_quality = quality_var.get()
        data_input.update({"quality": choice_quality})
        if choice_quality == 'AS':
            quality_menu.config(bg="brown", fg="white", activebackground="brown", activeforeground="white")
            log_action('"astropy/JPL-NASA [AS]" mode is selected','brown')
        elif choice_quality == 'SK':
            quality_menu.config(bg="white", fg="brown", activebackground="white", activeforeground="brown")
            log_action('"skyfield/JPL-NASA [SK]" mode is selected','brown')

        return choice_quality

    quality = ['AS', 'SK']
    quality_var = StringVar()
    quality_var.set(quality[0])

    quality_menu = OptionMenu(tool2_bar, quality_var, *quality, command=change_quality)
    quality_menu.place(in_=tool2_bar, relx=0.73, rely=0.32)
    quality_menu.config(bg="white", fg="brown", activebackground="white", activeforeground="brown")
    #############################################

    log_action('"astropy/JPL-NASA [AS]" mode is selected by default','brown')


    # Up - left #################################
    source_title_box = Label(tool3_bar, text="Status Panel", font=(stile, dimensione_titolo), fg='Blue', bg='SkyBlue1').place(in_=tool3_bar, relx=0.04, rely=0.15, anchor=W)

    var_stat = StringVar()
    var_stat.set('')                                                                       # set default path
    stat_box = Label(tool3_bar, text="System Status   ", font=(stile, dimensione), bg='SkyBlue1').place(in_=tool3_bar, relx=0.03, rely=0.40, anchor=W)
    txt_stat_box = Entry(tool3_bar, width=len_box-7, textvariable=var_stat)
    txt_stat_box.place(in_=tool3_bar, relx=0.50, rely=0.40, anchor=W)

    wi2, he2, r2 = 15, 15, 7
    canvas_stat = tk.Canvas(tool3_bar, width=wi2, height=he2, borderwidth=0, bg='SkyBlue1', highlightbackground = 'SkyBlue1')   # Create 200x200 Canvas widget
    canvas_stat.place(relx=0.43, rely=0.33)
    oval_stat = canvas_stat.create_oval(wi2/2-r2, he2/2+r2, wi2/2+r2, he2/2-r2)            # Create a circle on the Canvas

    var_pntg = StringVar()
    var_pntg.set('')                                                                       # set default path
    pntg_box = Label(tool3_bar, text="Pointing Status ", font=(stile, dimensione), bg='SkyBlue1').place(in_=tool3_bar, relx=0.03, rely=0.60, anchor=W)
    txt_pntg_box = Entry(tool3_bar, width=len_box-7, textvariable=var_pntg)
    txt_pntg_box.place(in_=tool3_bar, relx=0.50, rely=0.60, anchor=W)

    var_rec = StringVar()
    var_rec.set('')                                                                        # set default path
    rec_box = Label(tool3_bar, text="Receiver        ", font=(stile, dimensione), bg='SkyBlue1').place(in_=tool3_bar, relx=0.03, rely=0.80, anchor=W)
    txt_rec_box = Entry(tool3_bar, width=len_box-7, textvariable=var_rec, bg='white')
    txt_rec_box.place(in_=tool3_bar, relx=0.50, rely=0.80, anchor=W)
    #############################################


    # Box - left 1 ##############################
    source_title_box = Label(tool1l_bar, text="Time Panel", font=(stile, dimensione_titolo), fg='Blue', bg='SkyBlue1').place(in_=tool1l_bar, relx=0.04, rely=0.15, anchor=W)

    var_time = StringVar()
    var_time.set('')                                                                       # set default path
    time_box = Label(tool1l_bar, text="UTC epoch       ", font=(stile, dimensione), bg='SkyBlue1').place(in_=tool1l_bar, relx=0.03, rely=0.45, anchor=W)
    txt_path_box = Entry(tool1l_bar, width=len_box+3, textvariable=var_time)
    txt_path_box.place(in_=tool1l_bar, relx=0.46, rely=0.45, anchor=W)

    var_time_lst = StringVar()
    var_time_lst.set('')                                                                   # set default path
    time_lst_box = Label(tool1l_bar, text="LST time        ", font=(stile, dimensione), bg='SkyBlue1').place(in_=tool1l_bar, relx=0.03, rely=0.77, anchor=W)
    txt_time_lst_box = Entry(tool1l_bar, width=len_box+3, textvariable=var_time_lst)
    txt_time_lst_box.place(in_=tool1l_bar, relx=0.46, rely=0.77, anchor=W)
    #############################################


    # Box - left 2 ##############################
    source_title_box = Label(tool2l_bar, text="Distance Panel", font=(stile, dimensione_titolo), fg='Blue', bg='SkyBlue1').place(in_=tool2l_bar, relx=0.04, rely=0.15, anchor=W)

    var_dist = StringVar()
    var_dist.set('')                                                                       # set default path
    dist_box = Label(tool2l_bar,  text="SUN/telescope (live, deg) ", font=(stile, dimensione), bg='SkyBlue1').place(in_=tool2l_bar, relx=0.03, rely=0.35, anchor=W)
    txt_dist_box = Entry(tool2l_bar, width=len_box-8, textvariable=var_dist)
    txt_dist_box.place(in_=tool2l_bar, relx=0.68, rely=0.35, anchor=W)

    var_distc = StringVar()
    var_distc.set('')                                                                      # set default path
    distc_box = Label(tool2l_bar, text="SUN/source (live, deg)    ", font=(stile, dimensione), bg='SkyBlue1').place(in_=tool2l_bar, relx=0.03, rely=0.52, anchor=W)
    txt_distc_box = Entry(tool2l_bar, width=len_box-8, textvariable=var_distc)
    txt_distc_box.place(in_=tool2l_bar, relx=0.68, rely=0.52, anchor=W)

    var_distm = StringVar()
    var_distm.set('')                                                                      # set default path
    distm_box = Label(tool2l_bar,  text="MOON/telescope (live, deg)", font=(stile, dimensione), bg='SkyBlue1').place(in_=tool2l_bar, relx=0.03, rely=0.69, anchor=W)
    txt_distm_box = Entry(tool2l_bar, width=len_box-8, textvariable=var_distm)
    txt_distm_box.place(in_=tool2l_bar, relx=0.68, rely=0.69, anchor=W)

    var_distmc = StringVar()
    var_distmc.set('')                                                                     # set default path
    distmc_box = Label(tool2l_bar, text="MOON/source (live, deg)   ", font=(stile, dimensione), bg='SkyBlue1').place(in_=tool2l_bar, relx=0.03, rely=0.86, anchor=W)
    txt_distmc_box = Entry(tool2l_bar, width=len_box-8, textvariable=var_distmc)
    txt_distmc_box.place(in_=tool2l_bar, relx=0.68, rely=0.86, anchor=W)
    #############################################


    # Box - left 3 ##############################
    source_title_box = Label(tool3l_bar, text="Pointing Panel", font=(stile, dimensione_titolo), fg='Blue', bg='SkyBlue1').place(in_=tool3l_bar, relx=0.04, rely=0.15, anchor=W)

    var_point_ra = StringVar()
    var_point_ra.set('')                                                                   # set default path
    point_ra_box = Label(tool3l_bar, text="RA (hms) ", font=(stile, dimensione), bg='SkyBlue1').place(in_=tool3l_bar, relx=0.03, rely=0.50, anchor=W)
    txt_point_ra_box = Entry(tool3l_bar, width=len_box-9, textvariable=var_point_ra)
    txt_point_ra_box.place(in_=tool3l_bar, relx=0.20, rely=0.50, anchor=W)

    var_point_dec = StringVar()
    var_point_dec.set('')                                                                  # set default path
    point_dec_box = Label(tool3l_bar, text="DEC (dms) ", font=(stile, dimensione), bg='SkyBlue1').place(in_=tool3l_bar, relx=0.50, rely=0.50, anchor=W)
    txt_point_dec_box = Entry(tool3l_bar, width=len_box-9, textvariable=var_point_dec)
    txt_point_dec_box.place(in_=tool3l_bar, relx=0.70, rely=0.50, anchor=W)

    var_point_az = StringVar()
    var_point_az.set('')                                                                   # set default path
    point_az_box = Label(tool3l_bar, text="Az (deg) ", font=(stile, dimensione), bg='SkyBlue1').place(in_=tool3l_bar, relx=0.03, rely=0.75, anchor=W)
    txt_point_az_box = Entry(tool3l_bar, width=len_box-9, textvariable=var_point_az)
    txt_point_az_box.place(in_=tool3l_bar, relx=0.20, rely=0.75, anchor=W)

    var_point_el = StringVar()
    var_point_el.set('')                                                                   # set default path
    point_el_box = Label(tool3l_bar, text="El (deg) ", font=(stile, dimensione), bg='SkyBlue1').place(in_=tool3l_bar, relx=0.50, rely=0.75, anchor=W)
    txt_point_el_box = Entry(tool3l_bar, width=len_box-9, textvariable=var_point_el)
    txt_point_el_box.place(in_=tool3l_bar, relx=0.70, rely=0.75, anchor=W)
    #############################################


    # Box - left 4 ##############################
    source_title_box = Label(tool4l_bar, text="RFI Panel", font=(stile, dimensione_titolo), fg='Blue', bg='SkyBlue1').place(in_=tool4l_bar, relx=0.04, rely=0.2, anchor=W)

    canvas_stat2 = tk.Canvas(tool4l_bar, width=wi2, height=he2, borderwidth=0, bg='SkyBlue1', highlightbackground = 'SkyBlue1')  # Create 200x200 Canvas widget
    canvas_stat2.place(relx=0.234, rely=0.502)
    oval_stat2 = canvas_stat2.create_oval(wi2/2-r2, he2/2+r2, wi2/2+r2, he2/2-r2)          # Create a circle on the Canvas

    var_rfi = StringVar()
    var_rfi.set('')                                                                        # set default path
    rfi_box = Label(tool4l_bar, text="RFI source", font=(stile, dimensione), bg='SkyBlue1').place(in_=tool4l_bar, relx=0.03, rely=0.65, anchor=W)
    txt_rfi_box = Entry(tool4l_bar, width=len_box+11, textvariable=var_rfi)
    txt_rfi_box.place(in_=tool4l_bar, relx=0.292, rely=0.65, anchor=W)
    #############################################


    # Box - top #################################
    source_title_box = Label(tool2_bar, text="Source Panel", font=(stile, dimensione_titolo), fg='Blue', bg='SkyBlue1').place(in_=tool2_bar, relx=0.04, rely=0.15, anchor=W)

    var_name_source = StringVar()
    var_name_source.set('')                                                                # set default path
    name_source_box = Label(tool2_bar, text="Source Name ", font=(stile, dimensione), bg='SkyBlue1').place(in_=tool2_bar, relx=0.02, rely=0.40, anchor=W)
    txt_name_source_box = Entry(tool2_bar, width=30, textvariable=var_name_source)
    txt_name_source_box.place(in_=tool2_bar, relx=0.22, rely=0.40, anchor=W)

    var_point_az_c = StringVar()
    var_point_az_c.set('')                                                                 # set default path
    point_az_c_box = Label(tool2_bar, text="Az (deg) ", font=(stile, dimensione), bg='SkyBlue1').place(in_=tool2_bar, relx=0.02, rely=0.75, anchor=W)
    txt_point_az_c_box = Entry(tool2_bar, width=16, textvariable=var_point_az_c)
    txt_point_az_c_box.place(in_=tool2_bar, relx=0.22, rely=0.75, anchor=W)

    var_point_el_c = StringVar()
    var_point_el_c.set('')                                                                 # set default path
    point_el_c_box = Label(tool2_bar, text="El (deg) ", font=(stile, dimensione), bg='SkyBlue1').place(in_=tool2_bar, relx=0.52, rely=0.75, anchor=W)
    txt_point_el_c_box = Entry(tool2_bar, width=16, textvariable=var_point_el_c)
    txt_point_el_c_box.place(in_=tool2_bar, relx=0.72, rely=0.75, anchor=W)

    gc.disable()

    root.mainloop()


# Questo blocco garantisce che se lanci il file, la funzione main() parta automaticamente
if __name__ == "__main__":
    main()
