import socket


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
