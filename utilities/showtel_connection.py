import socket

ip_add_discos = '192.168.56.200'
port_discos = 30000

ip_add_seadas = '192.168.56.200'
port_seadas = 30000


def connect_port(ip_add, port):
    # ip_add = str ['192.168.56.200']
    # port = int [30000]
    import socket
    s = socket.socket()
    s.settimeout(2)                                                     # waiting for 2 second the connection answer
    try:
        s.connect((ip_add, port))
        conn = 1
    except:
        conn = 0

    return s, conn


def connect_status_discos(skt):
    # skt = the output of connect_discos (socket type)
    ip_add = ip_add_discos  #'192.168.56.200'
    port = port_discos      #30000
    try:
        test = skt.getpeername()
        if test[0] == ip_add and test[1] == port:
            sts = 1
    except:
        sts = 0
    return sts


def connect_status_seadas(skt):
    # skt = the output of connect_discos (socket type)
    ip_add = ip_add_seadas  #'192.168.56.200'
    port = ip_add_seadas    #30000
    try:
        test = skt.getpeername()
        if test[0] == ip_add and test[1] == port:
            sts = 1
    except:
        sts = 0
    return sts


def on_app_discos():
    s, c = connect_port(ip_add_discos, port_discos)#'192.168.56.200', 30000)                      # giusto
    ss = connect_status_discos(s)

    return s, c, ss


def on_app_seadas():
    s, c = connect_port(ip_add_seadas, port_seadas)#'192.168.56.200', 30000)                      # giusto
    ss = connect_status_seadas(s)

    return s, c, ss
