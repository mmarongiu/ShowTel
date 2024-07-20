from astropy import units as u
from astropy.coordinates import Angle, angular_separation
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
#from astropy.coordinates import get_sun, get_moon
from astropy.coordinates import solar_system_ephemeris, get_body_barycentric, get_body
from astropy.time import Time, TimeDelta
from astropy.constants import au
from astropy.coordinates.erfa_astrom import erfa_astrom, ErfaAstromInterpolator
#from sunpy.coordinates import frames

import numpy as np
from numpy import deg2rad

from datetime import datetime, timedelta
import pytz, time, sched#ephem, sunpy

#from astropy.coordinates import solar_system_ephemeris
solar_system_ephemeris.set('jpl')

from astroquery.jplhorizons import Horizons
from skyfield import api, almanac


def slr_station(radiot):
    loc = EarthLocation.of_site(radiot).geodetic
    lat_site = 39.492778        #loc.lat.value    # deg
    lon_site = 9.245            #loc.lon.value    # deg
    h_site = 700                #loc.height.value   # meters
    l = EarthLocation.from_geodetic(lon=lon_site,lat=lat_site,height=h_site)
    
    return l, lat_site, lon_site, h_site


def loc_time(datehour, rt):
    # datehour = 'yyyy-ddd-hh:mm:ss.ff' - UTC time - string
    # rt = 'SRT' or 'MED' or 'NOT' - string

    rt_loc = slr_station(rt)[0]              #EarthLocation.of_site(rt)
    local = pytz.timezone('Europe/Rome')
    utc_dt = datetime.strptime(datehour, "%Y-%j-%H:%M:%S.%f")
    loc_time = local.fromutc(utc_dt)
    str_loc_time = str(loc_time)[:-6]
    try:
        loc_dt = datetime.strptime(str_loc_time, "%Y-%m-%d %H:%M:%S.%f")
    except:
        loc_dt = datetime.strptime(str_loc_time, "%Y-%m-%d %H:%M:%S")
    time0_loc = Time(loc_dt, scale='local', location=rt_loc)
    time0_utc = Time(utc_dt, scale='utc', location=rt_loc)
    time0_lst = time0_utc.sidereal_time('mean')

    h_lst = int(time0_lst.hms.h)
    m_lst = int(time0_lst.hms.m)
    s_lst = time0_lst.hms.s

    if h_lst < 10:
        h_lst_str = '0' + str(h_lst)
    else:
        h_lst_str = str(h_lst)

    if m_lst < 10:
        m_lst_str = '0' + str(m_lst)
    else:
        m_lst_str = str(m_lst)

    if s_lst < 10:
        s_lst_str = '0' + '%.3f' % s_lst
    else:
        s_lst_str = '%.3f' % s_lst

    lst_string = h_lst_str + ':' + m_lst_str + ':' + s_lst_str

    return time0_loc, time0_lst, rt_loc, time0_utc, lst_string


def coord_altaz2radec(t_utc, pos, d_t, radiot, multi):
    # t_utc = UTC time [Time object]
    # d_t = array of time elements [astropy Time object]
    # radiot = name of the radio telescope [SRT, Medicina, Noto]
    # pos = (alt, az) [alt and az in units of deg, float]
    # multi = 0 without time range - multi = 1 with time range
    from astropy import units as u

    location = slr_station(radiot)[0]
    #location = EarthLocation.of_site(radiot)                                                                   # geocentric
    frame = AltAz(obstime=t_utc, location=location)

    pos_altaz = SkyCoord(alt=pos[0], az=pos[1], unit=(u.deg, u.deg), frame=frame)
    pos_radec = pos_altaz.transform_to('icrs')

    if multi == 0:
        extra_t = 0 * u.hour
        d_t = 0
        obstime = t_utc + extra_t
        time_min, time_max = t_utc, t_utc
    elif multi == 1:
        obstime = d_t
        time_min, time_max = d_t[0], d_t[-1]

    frame = AltAz(obstime=obstime, location=location)
    evo_altaz = pos_radec.transform_to(frame)

    return evo_altaz.az.value, evo_altaz.alt.value, obstime, time_min, time_max, pos_radec, pos_altaz


def dsun_tool(radiot, utc0_dt):
    d_ref = au                                                                                                 # one astronomical unit in meters
    if (radiot == 'SRT') or (radiot == 'S'):
        l0, l1 = EarthLocation.of_site('SRT').geodetic, EarthLocation.of_site('SRT').geocentric
    elif (radiot == 'MED') or (radiot == 'M'):
        l0, l1 = EarthLocation.of_site('medicina').geodetic, EarthLocation.of_site('medicina').geocentric
    elif (radiot == 'NOTO') or (radiot == 'N'):
        l0, l1 = EarthLocation.of_site('Noto').geodetic, EarthLocation.of_site('Noto').geocentric

    elev = l0.height.value
    h_site = 320
    
    obs = ephem.Observer()
    obs.date = utc0_dt
    obs.lon = deg2rad(l0.lon.value)
    obs.lat = deg2rad(l0.lat.value)
    obs.elevation = elev + h_site
    obs.pressure = 0

    sun = ephem.Sun()
    sun.compute(obs)
    dist_sun = sun.earth_distance*d_ref
    
    return dist_sun


def visibility(ra, dec, time0):
    source = SkyCoord(ra=Angle(ra*u.hourangle), dec=(dec*u.deg)).gcrs
    time0 = Time(time0)

    sun_distance = source.separation(SkyCoord(get_sun(time0).ra,get_sun(time0).dec))
    moon_distance = source.separation(SkyCoord(get_moon(time0).ra,get_moon(time0).dec))
    
    return moon_distance, sun_distance


def _sunmoon_evo(date_t, r_t, delta):
    # date_t -------- 'yyyy-ddd-hh:mm:ss.fff', string (DISCOS output)
    # r_t ----------- 'SRT', string
    # delta --------- range of hour, float

    time0_loc, time_sid, rt_loc, time0_utc, lst_str = loc_time(date_t, r_t)
    
    time_range = np.arange(-delta, delta+1, 0.25)                                                              # I consider the range between -12 hours and +12 hours with respect to the local time of observation

    dt_fast = lambda tt0, tt: tt0 + TimeDelta(tt*u.h)
    dt_loc = dt_fast(time0_loc, time_range)
    dt_utc = dt_fast(time0_utc, time_range)

    altaz = AltAz(obstime=dt_utc, location=rt_loc)

    with solar_system_ephemeris.set('./utilities/de441.bsp'):                                                            # jpl DE441
        sole = get_body('sun', dt_utc, rt_loc)
        luna = get_body('moon', dt_utc, rt_loc)

    sun_altaz = sole.transform_to(altaz)
    moon_altaz = luna.transform_to(altaz)

    result = np.array([dt_loc, dt_utc, sun_altaz.az, sun_altaz.alt, moon_altaz.az, moon_altaz.alt, (dt_utc - time0_utc).value], dtype=object).T
    
    mask_ref = np.where(abs(result[:,-1]) < 0.001)[0]
    ref = result[mask_ref][0]

    return ref, result


def sunmoon_evo(delta, date_t, r_t, step_t):
    # date_t -------- 'yyyy-ddd-hh:mm:ss.fff', string (DISCOS output)
    # r_t ----------- 'SRT', string
    # delta --------- duration of the solar map [in units of hours, float]
    # step_t -------- step time [in units of minutes, int]

    delta_min = 0.5*delta*60                                                                                       # delta in units of minutes
    step_t_min = int(step_t)

    loc_site = slr_station(r_t)[0]

    time0_loc, time_sid, rt_loc, time0_utc, lst_str = loc_time(date_t, r_t)
    
    time_range = np.arange(-delta_min, delta_min+0.1, step_t_min)                                                              # I consider the range between -12 hours and +12 hours with respect to the local time of observation

    dt_fast = lambda tt0, tt: tt0 + TimeDelta(tt*u.minute)
    dt_utc = dt_fast(time0_utc, time_range)

    altaz = AltAz(obstime=dt_utc, location=rt_loc)

    with solar_system_ephemeris.set('./utilities/de441.bsp'):                                                            # jpl DE441
        sun_radec = get_body('sun', dt_utc, rt_loc)
        moon_radec = get_body('moon', dt_utc, rt_loc)

    sun_altaz = sun_radec.transform_to(altaz)
    moon_altaz = moon_radec.transform_to(altaz)

    # angle in decimal mode
    sun_eph_radec = SkyCoord(sun_radec.ra, sun_radec.dec, frame='icrs', unit='deg')
    sun_eph_altaz = SkyCoord(alt = sun_altaz.alt, az = sun_altaz.az, frame='altaz', location=loc_site, unit='deg')

    moon_eph_radec = SkyCoord(moon_radec.ra, moon_radec.dec, frame='icrs', unit='deg')
    moon_eph_altaz = SkyCoord(alt = moon_altaz.alt, az = moon_altaz.az, frame='altaz', location=loc_site, unit='deg')

    result = np.array([dt_utc, sun_eph_altaz.az, sun_eph_altaz.alt, moon_eph_altaz.az, moon_eph_altaz.alt, (dt_utc - time0_utc).value], dtype=object).T

    mask_ref = np.where(abs(result[:,-1]) < 0.001)[0]
    ref = np.array([sun_eph_radec[mask_ref], sun_eph_altaz[mask_ref], moon_eph_radec[mask_ref], moon_eph_altaz[mask_ref]], dtype=object).T
    ref_plot = [dt_utc[mask_ref], sun_eph_altaz[mask_ref].az[0], sun_eph_altaz[mask_ref].alt[0], moon_eph_altaz[mask_ref].az[0], moon_eph_altaz[mask_ref].alt[0]]

    return sun_eph_radec, sun_eph_altaz, moon_eph_radec, moon_eph_altaz, result, ref, ref_plot


def horizon_eph_now(obs_t, radiot):
    # obs_t -------- astropy Time element
    # instrument --- instrument, string

    dt_utc = obs_t
    dt_utc_string = dt_utc.strftime('%Y-%m-%dT%H:%M:%S.%f')

    site_tool = slr_station(radiot)
    loc_site, lat_site, lon_site, h_site = site_tool[0], site_tool[1], site_tool[2], site_tool[3]*1e-3        # altitude is in units of km
    site_dict = {'lat': lat_site, 'lon': lon_site, 'elevation': h_site}

    sole_ref = Horizons(id='Sun', location=site_dict, epochs=dt_utc.jd)
    luna_ref = Horizons(id='301', location=site_dict, epochs=dt_utc.jd)
    sun_ref_eph = sole_ref.ephemerides()
    moon_ref_eph = luna_ref.ephemerides()
    
    try:
        sun_ref_eph_datetime = datetime.strptime(sun_ref_eph['datetime_str'][0], "%Y-%b-%d %H:%M:%S.%f")
    except:
        sun_ref_eph_datetime = datetime.strptime(sun_ref_eph['datetime_str'][0], "%Y-%b-%d %H:%M:%S")
    sun_ref_eph['datetime_str'] = sun_ref_eph_datetime
    sun_ref_eph['datetime_str'].name = 'epoch_datetime'

    try:
        moon_ref_eph_datetime = datetime.strptime(moon_ref_eph['datetime_str'][0], "%Y-%b-%d %H:%M:%S.%f")
    except:
        moon_ref_eph_datetime = datetime.strptime(moon_ref_eph['datetime_str'][0], "%Y-%b-%d %H:%M:%S")
    moon_ref_eph['datetime_str'] = moon_ref_eph_datetime
    moon_ref_eph['datetime_str'].name = 'epoch_datetime'

    # angle in decimal mode
    sun_ref_eph_radec = SkyCoord(sun_ref_eph['RA'][0], sun_ref_eph['DEC'][0], frame='icrs', unit='deg')
    sun_ref_eph_altaz = SkyCoord(alt = sun_ref_eph['EL'][0], az = sun_ref_eph['AZ'][0], frame='altaz', location=loc_site, unit='deg')

    moon_ref_eph_radec = SkyCoord(moon_ref_eph['RA'][0], moon_ref_eph['DEC'][0], frame='icrs', unit='deg')
    moon_ref_eph_altaz = SkyCoord(alt = moon_ref_eph['EL'][0], az = moon_ref_eph['AZ'][0], frame='altaz', location=loc_site, unit='deg')

    result = np.array([dt_utc, sun_ref_eph_altaz.az, sun_ref_eph_altaz.alt, moon_ref_eph_altaz.az, moon_ref_eph_altaz.alt], dtype=object).T
    
    return sun_ref_eph, moon_ref_eph, sun_ref_eph_radec, sun_ref_eph_altaz, moon_ref_eph_radec, moon_ref_eph_altaz, result


def horizon_eph(t_tot, obs_t, radiot, step_t):
    # t_tot -------- duration of the solar map [in units of hours, float]
    # obs_t -------- astropy Time element
    # instrument --- instrument, string
    # step_t ------- step time [in units of minutes, int]
    
    step = str(step_t)
    t_r_down = np.arange(-0.5*t_tot*60, 0, step_t)
    t_r_up = np.arange(0, 0.5*t_tot*60 + 0.1, step_t)
    time_range = np.append(t_r_down, t_r_up)
    
    dt_fast = lambda tt0, tt: tt0 + TimeDelta(tt*u.min)
    dt_utc = dt_fast(obs_t, time_range)
    dt_utc_string = dt_utc.strftime('%Y-%m-%dT%H:%M:%S.%f')

    site_tool = slr_station(radiot)
    loc_site, lat_site, lon_site, h_site = site_tool[0], site_tool[1], site_tool[2], site_tool[3]*1e-3        # altitude is in units of km
    site_dict = {'lat': lat_site, 'lon': lon_site, 'elevation': h_site}

    func_eph = lambda src, sdict, dt_ini, dt_end, st: Horizons(id=src, location=sdict, epochs={'start': dt_ini, 'stop': dt_end, 'step': st + 'm'})
    sole = func_eph('Sun', site_dict, dt_utc_string[0], dt_utc_string[-1], step)
    luna = func_eph('301', site_dict, dt_utc_string[0], dt_utc_string[-1], step)

    sun_eph = sole.ephemerides()
    moon_eph = luna.ephemerides()
    
    sun_eph_datetime, moon_eph_datetime = [], []
    for i in sun_eph['datetime_str']:
        try:
            t = datetime.strptime(i, "%Y-%b-%d %H:%M:%S.%f")
        except:
            t = datetime.strptime(i, "%Y-%b-%d %H:%M:%S")
        sun_eph_datetime = np.append(sun_eph_datetime,t)
    sun_eph['datetime_str'] = sun_eph_datetime
    sun_eph['datetime_str'].name = 'epoch_datetime'

    for i in moon_eph['datetime_str']:
        try:
            t = datetime.strptime(i, "%Y-%b-%d %H:%M:%S.%f")
        except:
            t = datetime.strptime(i, "%Y-%b-%d %H:%M:%S")
        moon_eph_datetime = np.append(moon_eph_datetime,t)
    moon_eph['datetime_str'] = moon_eph_datetime
    moon_eph['datetime_str'].name = 'epoch_datetime'

    # angle in decimal mode
    sun_eph_radec = SkyCoord(sun_eph['RA'], sun_eph['DEC'], frame='icrs', unit='deg')
    sun_eph_altaz = SkyCoord(alt = sun_eph['EL'], az = sun_eph['AZ'], frame='altaz', location=loc_site)

    moon_eph_radec = SkyCoord(moon_eph['RA'], moon_eph['DEC'], frame='icrs', unit='deg')
    moon_eph_altaz = SkyCoord(alt = moon_eph['EL'], az = moon_eph['AZ'], frame='altaz', location=loc_site)

    result = np.array([dt_utc, sun_eph_altaz.az, sun_eph_altaz.alt, moon_eph_altaz.az, moon_eph_altaz.alt, (dt_utc - obs_t).value], dtype=object).T
    
    return sun_eph, moon_eph, sun_eph_radec, sun_eph_altaz, moon_eph_radec, moon_eph_altaz, result


def skyfield_ephem(epoch, ll, de_eph_file, mode):
    # epoch = astropy Time object (isot)
    # ll = astropy EarthLocation object of the observing site
    # de_eph_file = directory path of the .bsp file of the ephemeris [string]
    # mode = 's' [Sun], 'm' [Moon]

    if mode == 's':
        source = 'Sun'
    elif mode == 'm':
        source = 'Moon'
    
    ts = api.load.timescale()
    tempo = ts.utc(int(epoch.value.split('-')[0]), int(epoch.value.split('-')[1]), int(epoch.value.split('-')[2][0:2]), int(epoch.value.split(':')[0][-2:]), int(epoch.value.split(':')[1]), float(epoch.value.split(':')[2]))
    ephem = api.load_file(de_eph_file)
    ephem_source = ephem[source]
    earth = ephem["Earth"]
    location = api.Topos(latitude_degrees=ll.lat.value, longitude_degrees=ll.lon.value, elevation_m=ll.height.value)
    source_pos = (earth + location).at(tempo).observe(ephem_source)                        # astrometric position of the source
    ra, dec, dist = source_pos.radec()
    source_radec = SkyCoord(ra._degrees, dec._degrees, frame='icrs', location=ll, unit='deg')
    altaz = AltAz(obstime=epoch.value, location=ll)
    s0_altaz = source_radec.transform_to(altaz)
    #dist_sun = sunpy.coordinates.sun.earth_distance(epoch)
    #source_altaz = SkyCoord(az=s0_altaz.az, alt=s0_altaz.alt, distance=dist_sun, frame=altaz)
    source_altaz = SkyCoord(az=s0_altaz.az, alt=s0_altaz.alt, frame=altaz)
    
    return source_radec, source_altaz, altaz


def skyfield_eph_now(dt_utc, radiot):
    # obs_t -------- astropy Time element (isot)
    # instrument --- instrument, string

    dt_utc = Time(dt_utc, format='isot')   # astropy Time object (isot)
    dt_utc_string = dt_utc.strftime('%Y-%m-%dT%H:%M:%S.%f')

    loc_site = slr_station(radiot)[0]

    sun_eph_radec, sun_eph_altaz, altaz = skyfield_ephem(dt_utc, loc_site, './utilities/de441.bsp', 's')
    moon_eph_radec, moon_eph_altaz, altaz = skyfield_ephem(dt_utc, loc_site, './utilities/de441.bsp', 'm')

    # angle in decimal mode
    sun_ref_eph_radec = SkyCoord(sun_eph_radec.ra.value, sun_eph_radec.dec.value, frame='icrs', unit='deg')
    sun_ref_eph_altaz = SkyCoord(alt = sun_eph_altaz.alt.value, az = sun_eph_altaz.az.value, frame=altaz, unit='deg')

    moon_ref_eph_radec = SkyCoord(moon_eph_radec.ra.value, moon_eph_radec.dec.value, frame='icrs', unit='deg')
    moon_ref_eph_altaz = SkyCoord(alt = moon_eph_altaz.alt.value, az = moon_eph_altaz.az.value, frame=altaz, unit='deg')

    result = np.array([dt_utc, sun_ref_eph_altaz.az, sun_ref_eph_altaz.alt, moon_ref_eph_altaz.az, moon_ref_eph_altaz.alt], dtype=object).T
    
    return 0., 0., sun_ref_eph_radec, sun_ref_eph_altaz, moon_ref_eph_radec, moon_ref_eph_altaz, result


def skyfield_eph(t_tot, obs_t, radiot, step_t):
    # t_tot -------- duration of the solar map [in units of hours, float]
    # obs_t -------- astropy Time element
    # instrument --- instrument, string
    # step_t ------- step time [in units of minutes, int]

    obs_t = Time(obs_t, format='isot')   # astropy Time object (isot)
    loc_site = slr_station(radiot)[0]
    
    step = str(step_t)
    t_r_down = np.arange(-0.5*t_tot*60, 0, step_t)
    t_r_up = np.arange(0, 0.5*t_tot*60 + 0.1, step_t)
    time_range = np.append(t_r_down, t_r_up)
    
    dt_fast = lambda tt0, tt: tt0 + TimeDelta(tt*u.min)
    dt_utc = dt_fast(obs_t, time_range)
    dt_utc_string = dt_utc.strftime('%Y-%m-%dT%H:%M:%S.%f')

    fine_altaz, fine_sun_ra, fine_sun_dec, fine_sun_alt, fine_sun_az, fine_moon_ra, fine_moon_dec, fine_moon_alt, fine_moon_az = list(), [], [], [], [], [], [], [], []
    for i in dt_utc:
        tool = skyfield_eph_now(i, 'SRT')
        altaz = AltAz(obstime=i, location=loc_site)
        fine_altaz.append(altaz)
        tool_sun_ra, tool_sun_dec, tool_sun_alt, tool_sun_az, tool_moon_ra, tool_moon_dec, tool_moon_alt, tool_moon_az = tool[2].ra.value, tool[2].dec.value, tool[3].alt.value, tool[3].az.value, tool[4].ra.value, tool[4].dec.value, tool[5].alt.value, tool[5].az.value
        fine_sun_ra, fine_sun_dec, fine_sun_alt, fine_sun_az, fine_moon_ra, fine_moon_dec, fine_moon_alt, fine_moon_az = np.append(fine_sun_ra, tool_sun_ra), np.append(fine_sun_dec, tool_sun_dec), np.append(fine_sun_alt, tool_sun_alt), np.append(fine_sun_az, tool_sun_az), np.append(fine_moon_ra, tool_moon_ra), np.append(fine_moon_dec, tool_moon_dec), np.append(fine_moon_alt, tool_moon_alt), np.append(fine_moon_az, tool_moon_az)

    # angle in decimal mode
    sun_eph_radec = SkyCoord(fine_sun_ra, fine_sun_dec, frame='icrs', unit='deg')
    sun_eph_altaz = SkyCoord(alt = fine_sun_alt, az = fine_sun_az, frame='altaz', obstime=dt_utc, location=loc_site, unit='deg')

    moon_eph_radec = SkyCoord(fine_moon_ra, fine_moon_dec, frame='icrs', unit='deg')
    moon_eph_altaz = SkyCoord(alt = fine_moon_alt, az = fine_moon_az, frame='altaz', obstime=dt_utc, location=loc_site, unit='deg')

    result = np.array([dt_utc, sun_eph_altaz.az, sun_eph_altaz.alt, moon_eph_altaz.az, moon_eph_altaz.alt, (dt_utc - obs_t).value], dtype=object).T
    
    return 0., 0., sun_eph_radec, sun_eph_altaz, moon_eph_radec, moon_eph_altaz, result


def calc_angdist(date_t, r_t, az_loc, alt_loc, az_com_loc, alt_com_loc):
    # r_t ----------- 'SRT', string
    # date_t -------- 'yyyy-ddd-hh:mm:ss.fff', string
    # az_loc -------- 'aaa.bbb' (in units of deg), string
    # alt_loc ------- 'aaa.bbb' (in units of deg), string
    # az_com_loc ---- 'aaa.bbb' (in units of deg), string
    # alt_com_loc --- 'aaa.bbb' (in units of deg), string
    
    time0_loc, time_sid, rt_loc, time0_utc, lst_str = loc_time(date_t, r_t)

    ref, result = _sunmoon_evo(date_t, r_t, 12)

    sun_altaz_az, sun_altaz_alt, moon_altaz_az, moon_altaz_alt = ref[2], ref[3], ref[4], ref[5]
    
    az_rt, alt_rt = Angle(az_loc,unit=u.deg), Angle(alt_loc,unit=u.deg)                                        # azimuth/altitude coordinates (radio telescope in real time)
    rt_altaz = SkyCoord(alt=alt_rt,az=az_rt, unit=(u.deg, u.deg), frame='altaz', obstime=time0_utc,location=rt_loc)
    az_rt_com, alt_rt_com = Angle(az_com_loc,unit=u.deg), Angle(alt_com_loc,unit=u.deg)                        # azimuth/altitude coordinates (source in real time)
    rt_altaz_com = SkyCoord(alt=alt_rt_com,az=az_rt_com, unit=(u.deg, u.deg), frame='altaz', obstime=time0_utc,location=rt_loc)

    ang_sep_altaz_now_sun2rt = angular_separation(sun_altaz_az*u.deg,sun_altaz_alt*u.deg,rt_altaz.az.value*u.deg,rt_altaz.alt.value*u.deg).to('deg')         # distance Sun - radio telescope
    ang_sep_altaz_com_sun2rt = angular_separation(sun_altaz_az*u.deg,sun_altaz_alt*u.deg,rt_altaz_com.az.value*u.deg,rt_altaz_com.alt.value*u.deg).to('deg') # distance Sun - source

    ang_sep_altaz_now_moon2rt = angular_separation(moon_altaz_az*u.deg,moon_altaz_alt*u.deg,rt_altaz.az,rt_altaz.alt.value*u.deg).to('deg')                  # distance Moon - radio telescope
    ang_sep_altaz_com_moon2rt = angular_separation(moon_altaz_az*u.deg,moon_altaz_alt*u.deg,rt_altaz_com.az,rt_altaz_com.alt.value*u.deg).to('deg')          # distance Moon - source

    return ang_sep_altaz_now_sun2rt, ang_sep_altaz_com_sun2rt, ang_sep_altaz_now_moon2rt, ang_sep_altaz_com_moon2rt, (sun_altaz_az, sun_altaz_alt), (moon_altaz_az, moon_altaz_alt), time0_loc, time0_utc, ref, result


def calc_angdist_radec(date_t, r_t, radec_pnt, radec_source, delta, step_t, quality):
    # r_t ----------- 'SRT', string
    # date_t -------- 'yyyy-ddd-hh:mm:ss.fff', string
    # radec_pnt ----- SkyCoord element (ra-dec, in units of deg)
    # radec_source -- SkyCoord element (ra-dec, in units of deg)
    # t_tot --------- total time of observation, float
    # t_step -------- time step (ra-dec, in units of minutes), float
    # quality ------- quality of the analysis (1 = astropy/JPL-NASA [LQ-fast] ; 2 = Horizons JPL-NASA [HQ-slow])
    
    time0_loc, time_sid, rt_loc, time0_utc, lst_str = loc_time(date_t, r_t)

    if quality == 1:
        tool = sunmoon_evo(delta, date_t, r_t, step_t)
        sun_eph_radec, sun_eph_altaz, moon_eph_radec, moon_eph_altaz, ref0 = tool[0], tool[1], tool[2], tool[3], tool[5][0]
        sun_ref_eph_radec, sun_ref_eph_altaz, moon_ref_eph_radec, moon_ref_eph_altaz = ref0[0], ref0[1], ref0[2], ref0[3]
    elif quality == 2:
        #tool = horizon_eph_now(time0_utc, r_t)
        tool = skyfield_eph_now(time0_utc, r_t)
        sun_ref_eph, moon_ref_eph, sun_ref_eph_radec, sun_ref_eph_altaz, moon_ref_eph_radec, moon_ref_eph_altaz = tool[0], tool[1], tool[2], tool[3], tool[4], tool[5]

    ang_sep_pnt_sun2rt = angular_separation(sun_ref_eph_radec.ra, sun_ref_eph_radec.dec, radec_pnt.ra, radec_pnt.dec).to('deg')                        # distance Sun - radio telescope
    ang_sep_source_sun2rt = angular_separation(sun_ref_eph_radec.ra, sun_ref_eph_radec.dec, radec_source.ra, radec_source.dec).to('deg')               # distance Sun - source

    ang_sep_pnt_moon2rt = angular_separation(moon_ref_eph_radec.ra, moon_ref_eph_radec.dec, radec_pnt.ra, radec_pnt.dec).to('deg')                     # distance Moon - radio telescope
    ang_sep_source_moon2rt = angular_separation(moon_ref_eph_radec.ra, moon_ref_eph_radec.dec, radec_source.ra, radec_source.dec).to('deg')            # distance Moon - source

    return ang_sep_pnt_sun2rt, ang_sep_source_sun2rt, ang_sep_pnt_moon2rt, ang_sep_source_moon2rt, tool


def calc_angdist_altaz(date_t, r_t, altaz_pnt, altaz_source, delta, step_t, quality):
    # r_t ----------- 'SRT', string
    # date_t -------- 'yyyy-ddd-hh:mm:ss.fff', string
    # altaz_pnt ----- SkyCoord element (alt-az, in units of deg)
    # altaz_source -- SkyCoord element (alt-az, in units of deg)
    # t_tot --------- total time of observation, float
    # t_step -------- time step (ra-dec, in units of minutes), float
    # quality ------- quality of the analysis (1 = astropy/JPL-NASA [LQ-fast] ; 2 = Horizons JPL-NASA [HQ-slow])
    
    time0_loc, time_sid, rt_loc, time0_utc, lst_str = loc_time(date_t, r_t)

    if quality == 1:
        tool = sunmoon_evo(delta, date_t, r_t, step_t)
        sun_eph_radec, sun_eph_altaz, moon_eph_radec, moon_eph_altaz, ref0 = tool[0], tool[1], tool[2], tool[3], tool[5][0]
        sun_ref_eph_radec, sun_ref_eph_altaz, moon_ref_eph_radec, moon_ref_eph_altaz = ref0[0], ref0[1], ref0[2], ref0[3]
    elif quality == 2:
        #tool = horizon_eph_now(time0_utc, r_t)
        tool = skyfield_eph_now(time0_utc, r_t)
        sun_ref_eph, moon_ref_eph, sun_ref_eph_radec, sun_ref_eph_altaz, moon_ref_eph_radec, moon_ref_eph_altaz = tool[0], tool[1], tool[2], tool[3], tool[4], tool[5]

    ang_sep_pnt_sun2rt = angular_separation(sun_ref_eph_altaz.az, sun_ref_eph_altaz.alt, altaz_pnt.az, altaz_pnt.alt).to('deg')                        # distance Sun - radio telescope
    ang_sep_source_sun2rt = angular_separation(sun_ref_eph_altaz.az, sun_ref_eph_altaz.alt, altaz_source.az, altaz_source.alt).to('deg')               # distance Sun - source

    ang_sep_pnt_moon2rt = angular_separation(moon_ref_eph_altaz.az, moon_ref_eph_altaz.alt, altaz_pnt.az, altaz_pnt.alt).to('deg')                     # distance Moon - radio telescope
    ang_sep_source_moon2rt = angular_separation(moon_ref_eph_altaz.az, moon_ref_eph_altaz.alt, altaz_source.az, altaz_source.alt).to('deg')            # distance Moon - source

    return ang_sep_pnt_sun2rt, ang_sep_source_sun2rt, ang_sep_pnt_moon2rt, ang_sep_source_moon2rt, tool
