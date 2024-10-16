#banda             nome   range_freq      freq_cent   freq_cent_SL
#banda K      ---> KKG    18 - 26.5       22.25       23
#banda C-high ---> CCB    5.7 - 7.7       6.7         7.35
#banda C-low  ---> CCG    4.2 - 5.6       4.9         ???
#banda L      ---> LLP    1.3 - 1.8       1.55        1.705
#banda P      ---> PPP    0.305 - 0.410   0.3575      ???

from astropy.table import Table
import numpy as np

def tool_rfi(az, el, name_ric, path_rfi, path_rec):
    tab_rfi = Table.read(path_rfi, format='ascii')
    tab_receivers = Table.read(path_rec, format='ascii')

    mask_receivers = (name_ric == tab_receivers['nome'])
    tab_rec = tab_receivers[mask_receivers]

    freq_min_rec, freq_max_rec = tab_rec['f_min'], tab_rec['f_max']
    freq_min_rfi, freq_max_rfi = tab_rfi['f_min'], tab_rfi['f_max']

    mask_freq = (freq_min_rec <= freq_max_rfi) & (freq_max_rec >= freq_min_rfi)
    tab_freq = tab_rfi[mask_freq]

    mask_single = (az >= tab_rfi['az_min']) & (az <= tab_rfi['az_max']) & (el >= tab_rfi['el_min']) & (el <= tab_rfi['el_max']) & (freq_min_rec <= freq_max_rfi) & (freq_max_rec >= freq_min_rfi)
    tab_single = tab_rfi[mask_single]

    if len(tab_single) == 0:
        rfi = 0
        rfi_name = ['']
        rfi_sp = ['']
    else:
        rfi = 1
        rfi_name = np.array(tab_single['RFI'])
        rfi_sp = np.array(tab_single['Spectrum'])
    
    return rfi, rfi_name, rfi_sp, tab_single, tab_freq
