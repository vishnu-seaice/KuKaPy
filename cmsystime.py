
import time
from datetime import datetime
from math import *
import numpy as np


def cmsystime_xmod(x, m):
    return ((((x%m) + m) + m) % m) 


def cmsystime_mjd2ymd(mjd):
    offset = 2400000.5
    offset_int = floor(offset)

    offset_fra = offset - offset_int
    nn = offset_fra + mjd
    jd_fra = cmsystime_xmod(nn+0.5, 1) - 0.5
    

    nn = nn + offset_int - jd_fra
    nn = nn + (floor(floor((nn - 4479.5)/ 36524.25) * 0.75 + 0.5 ) - 37)

    yr = floor(nn/365.25) - 4712
    dd = floor(cmsystime_xmod(nn-59.25, 365.25))

    mo = floor(cmsystime_xmod(floor((dd+0.5)/30.6) + 2, 12) + 1)
    da = floor(cmsystime_xmod(dd+0.5, 30.6) + 1) + 0.5 + jd_fra

    return yr, mo, da

def cmsystime(arg0):

    cmsystime_months = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', \
                        'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']

    cmsystime_dow = ['Sun', 'Mon', 'Tue', 'Wed', 'Thu', 'Fri', 'Sat']

    MJD_1970 = 40587
    JD_1970 = MJD_1970 + 2400000.5

    gmtime = time.time()
    cltime = datetime.now()
    gmfrac = gmtime % 86400
    gm_mjd = floor(gmtime-gmfrac)/86400 + MJD_1970
    gm_yr, gm_mo, gm_da = cmsystime_mjd2ymd(gm_mjd)
    gm_da = round(gm_da)


    da = cltime.day
    ltimes = np.array([cltime.hour, cltime.minute, cltime.second])
    lfrac = ltimes[2] + 60*(ltimes[1] + 60*ltimes[0])

    tz = lfrac - gmfrac

    if gm_da == da-1:
        tz = tz + 86400
    elif gm_da == da+1:
        tz = tz - 86400
    elif gm_da < da:
        tz = tz - 86400
    elif gm_da > da:
        tz = tz + 86400
        
    cmsystime_timezone = round(tz/60) * 60
    timezone = cmsystime_timezone


    offset = 0
    arg = arg0

    if offset != 0:
        arg = arg + offset
    
    mjd = floor(arg/86400) + MJD_1970
    dsecs = arg-floor(arg/86400)*86400
    hr = floor(dsecs/3600)
    dsecs = dsecs - hr*3600
    mi = floor(dsecs / 60) 
    dsecs = dsecs - mi*60
    se = dsecs
    yr, mo, da = cmsystime_mjd2ymd(mjd)
    dow = cmsystime_xmod((floor(mjd) - 51678), 7)

    try:
        n = len(yr)
    except:
        n = 1
    result = []*n
    
    if len(str(floor(se))) <= 1:

        secc = str(floor(se)) + '0'
    else:
        secc = str(floor(se))

    result = cmsystime_dow[dow] + ' ' + cmsystime_months[mo-1] + ' ' + str(int(da)) + ' ' \
             + str(hr) + ':' + str(mi) + ':' + secc + '  ' + str(yr)

    return result
                        