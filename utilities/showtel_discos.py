import socket
from discos_client import SRTClient


def discos_pars(skt):
    # skt = the output of connect_discos (socket type)
    #skt.settimeout(0.5)
    try:
        skt.sendall(b'antennaParameters')
        arr_pars1 = str(skt.recv(1024))
        arr_pars = arr_pars1[21:-3].split(',')
        #arr_pars2 = str(skt.recv(1024)[18:-1]).split(',')
        #arr_pars2 = str(skt.recv(1024))
    except:
        pass

    if 'arr_pars' in locals():
        arr_pars_dict = {'SysUTC': arr_pars[0], 'SystemStatus': arr_pars[1], 'SourceName': arr_pars[2], 'Azimuth': arr_pars[3], 'Elevation': arr_pars[4], 'RightAscension': arr_pars[5], 'Declination': arr_pars[6], 'GalacticLongitude': arr_pars[7], 'GalacticLatitude': arr_pars[8], 'CommandedAzimuth': arr_pars[9], 'CommandedElevation': arr_pars[10], 'AzimuthError': arr_pars[11], 'AzimuthCorrection': arr_pars[12], 'AzimuthOffset': arr_pars[13], 'ElevationError': arr_pars[14], 'ElevationCorrection': arr_pars[15], 'ElevationOffset': arr_pars[16], 'RefractionElevationCorrection': arr_pars[17], 'RaOffset': arr_pars[18], 'DeclOffset': arr_pars[19], 'GalLongitudeOffset': arr_pars[20], 'GalLatitudeOffset': arr_pars[21], 'ReceiverCode': arr_pars[22], 'LOFrequency': arr_pars[23], 'PointingStatus': arr_pars[24]}
        return arr_pars, arr_pars_dict, arr_pars1
    else:
        return [0], [0], [0]


def srtclient_pars(s):
    mount = s.mount.copy()
    antenna = s.antenna.copy()
    receivers = s.receivers.copy()
    backends = s.backends.copy()
    scheduler = s.scheduler.copy()

    arr_pars = [0, 0]
    pnt_sts = str(scheduler.tracking)
    pnt_txt = 'TRACKING' if pnt_sts == 'True' else 'SLEWING'
    sts_socket = mount.statusSocketConnected
    timestamp = str(antenna.timestamp.iso8601)

    arr_pars_dict = {'Site': str(antenna.site.name), 'SysUTC': timestamp, 'SystemStatus': str(antenna.status), 'SourceName': str(antenna.target.name), 'Azimuth': str(antenna.observedAzimuth), 'Elevation': str(antenna.observedElevation), 'RightAscension': str(antenna.observedRightAscension), 'Declination': str(antenna.observedDeclination), 'GalacticLongitude': str(antenna.observedGalLongitude), 'GalacticLatitude': str(antenna.observedGalLatitude), 'CommandedAzimuth': str(antenna.rawAzimuth), 'CommandedElevation': str(antenna.rawElevation), 'AzimuthCorrection': str(antenna.pointingAzimuthCorrection), 'AzimuthOffset': str(antenna.azimuthOffset), 'ElevationCorrection': str(antenna.pointingElevationCorrection), 'ElevationOffset': str(antenna.elevationOffset), 'RefractionElevationCorrection': str(antenna.refractionCorrection), 'RaOffset': str(antenna.rightAscensionOffset), 'DeclOffset': str(antenna.declinationOffset), 'GalLongitudeOffset': str(antenna.longitudeOffset), 'GalLatitudeOffset': str(antenna.latitudeOffset), 'ReceiverCode': str(receivers.currentSetup), 'Backend': str(backends.currentBackend), 'PointingStatus': pnt_txt, 'ProjectCode': str(scheduler.projectCode), 'Socket': str(sts_socket)}
        
    return arr_pars, arr_pars_dict

        
#'Site': xxx ..................................... SRT.antenna.site.get_value()                        # str
#'SysUTC': arr_pars[0] ........................... SRT.antenna.timestamp.iso8601.get_value()           # str, iso8601 yyyy-mm-ddThh:mm:ss.fffZ UTC
#'SystemStatus': arr_pars[1] ..................... SRT.receivers.status.get_value()                    # str
#'SourceName': arr_pars[2] ....................... SRT.antenna.target.name.get_value()                 # str
#'Azimuth': arr_pars[3] .......................... SRT.antenna.observedAzimuth.get_value()             # float, deg
#'Elevation': arr_pars[4] ........................ SRT.antenna.observedElevation.get_value()           # float, deg
#'RightAscension': arr_pars[5] ................... SRT.antenna.observedRightAscension.get_value()      # str, 
#'Declination': arr_pars[6] ...................... SRT.antenna.observedDeclination.get_value()         # str, 
#'GalacticLongitude': arr_pars[7] ................ SRT.antenna.observedGalLongitude.get_value()        # str, 
#'GalacticLatitude': arr_pars[8] ................. SRT.antenna.observedGalLatitude.get_value()         # str, 
#'CommandedAzimuth': arr_pars[9] ................. SRT.antenna.rawAzimuth.get_value()                  # float, deg
#'CommandedElevation': arr_pars[10] .............. SRT.antenna.rawElevation.get_value()                # float, deg
#'AzimuthError': arr_pars[11] .................... xxx
#'AzimuthCorrection': arr_pars[12] ............... SRT.antenna.pointingAzimuthCorrection.get_value()   # float, deg
#'AzimuthOffset': arr_pars[13] ................... SRT.antenna.azimuthOffset.get_value()               # float, deg
#'ElevationError': arr_pars[14] .................. xxx
#'ElevationCorrection': arr_pars[15] ............. SRT.antenna.pointingElevationCorrection.get_value() # float, deg
#'ElevationOffset': arr_pars[16].................. SRT.antenna.elevationOffset.get_value()             # float, deg
#'RefractionElevationCorrection': arr_pars[17] ... SRT.antenna.refractionCorrection.get_value()        # float, deg
#'RaOffset': arr_pars[18] ........................ SRT.antenna.rightAscensionOffset.get_value()        # float, deg
#'DeclOffset': arr_pars[19] ...................... SRT.antenna.declinationOffset.get_value()           # float, deg
#'GalLongitudeOffset': arr_pars[20] .............. SRT.antenna.longitudeOffset.get_value()             # float, deg
#'GalLatitudeOffset': arr_pars[21] ............... SRT.antenna.latitudeOffset.get_value()              # float, deg
#'ReceiverCode': arr_pars[22] .................... SRT.backends.currentSetup.get_value()               # str
#'Backend': xxx .................................. SRT.backends.currentBackend.get_value()             # str
#'LOFrequency': arr_pars[23] ..................... xxx
#'PointingStatus': arr_pars[24] .................. SRT.antenna.tracking.get_value()                    # Boolean (True = TRACKING; False = SLEWING)
#'ProjectCode': xxx .............................. SRT.scheduler.projectCode.get_value()               # str
