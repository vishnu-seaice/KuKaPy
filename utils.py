import numpy as np
import math
from math import *
from math import sin, cos, tan, atan, sqrt, atan2
from glob import glob
import re

import matplotlib.pyplot as plt 
from matplotlib.patches import Rectangle


from astropy.convolution import convolve
import scipy.ndimage
import cmath
import scatlib



def earth_radius(lat):
    a_equitorial = 6378137.0          #Earth's radius at the equator
    b_polar = 6356752.3                 #Earth's radius at the poles

    num = (a_equitorial**2*math.cos(lat))**2+(b_polar**2*math.sin(lat))**2
    denom = (a_equitorial*math.cos(lat))**2+(b_polar*math.sin(lat))**2

    R_e = math.sqrt(num/denom)
    
    return R_e

def compute_distance(gps_latitude,gps_longitude,elapsed_time,scatvars,configvars,independent_sample_index):
    lat = gps_latitude
    lon = gps_longitude
    time = elapsed_time

    nvals = len(gps_latitude)

    xoff = np.zeros(nvals)
    yoff = np.zeros(nvals)
    distance = np.zeros(nvals)
    bearing = np.zeros(nvals)
    
    lat1 = math.radians(lat[0])
    lon1 = math.radians(lon[0])
    
    R_e = earth_radius(lat1)

    for i in range(0, nvals):
        lat2 = math.radians(lat[i])
        lon2 = math.radians(lon[i])    
        delta_lat = (lat2-lat1)
        delta_lon = (lon2-lon1)

        afac = math.sin(delta_lat/2)**2+cos(lat1)*cos(lat2)*sin(delta_lon/2.)**2
        
        cfac = 2.0*atan2(sqrt(afac),sqrt(1.-afac))
        
        distance[i] = float(R_e*cfac)
        
        bearing[i] = float(atan2(sin(delta_lon)*cos(lat2),cos(lat1)*sin(lat2)-sin(lat1)*cos(lat2)*(cos(delta_lon))))
                
        if bearing[i] < 0:
            bearing[i] = bearing[i]+2*math.pi
        

        xoff[i]=distance[i]*sin(bearing[i])
        yoff[i]=distance[i]*cos(bearing[i])
    
    smoothfac = 20
    if smoothfac > nvals:
        smoothfac = nvals/6+1
        

    ## NOT ACCURAYE BELOW METHOD
    distance = scipy.ndimage.filters.uniform_filter1d(distance,size=smoothfac, mode="nearest")

    independent_sample_index0 = np.zeros(nvals)
    distance_ref = distance[0]

    # #compute velocity
    vec0 = [i for i in range(nvals-1)]
    vec1 = [i+1 for i in vec0]

    velocity0 = []
    for each in range(len(vec0)):
        velocity0 = (distance[vec1]-distance[vec0])/(elapsed_time[vec1]-elapsed_time[vec0])
    # velocity = [velocity0,velocity0[nvals-2]]
    velocity = np.append(velocity0, velocity0[nvals-2])
    
    Ku_ant_diameter=0.106                   #was incorrectly set to .15   
    Ka_ant_diameter=0.064                   #was incorrectly set to 0.09
    distance_to_independence=Ku_ant_diameter/2.0
    
    if scatvars["chirp_center_frequency_hz"] > 20e09:
        distance_to_independence = Ka_ant_diameter/2.0
    
    min_velocity = 0.5  # test for min velocity to avoid drifting GPS location from appearing as motion

    if configvars["i_proc_ind"] == 0:
        min_velocity =-999.
        distance_to_independence=-999.
    


    count_ind=0
    for i in range(0,nvals):        
        if abs(distance[i]-distance_ref) > distance_to_independence and velocity[i] > min_velocity:            
            independent_sample_index0[count_ind] = i
            count_ind = count_ind+1
            distance_ref = distance[i]

    independent_sample_index = independent_sample_index0[0:count_ind] 
    
    return distance, independent_sample_index

import numpy as np
from scipy.signal import resample
from scipy.interpolate import interp1d

def ResampleLinear1D(original, targetLen):
    original = np.array(original, dtype=np.float)
    index_arr = np.linspace(0, len(original)-1, num=targetLen, dtype=np.float)
    index_floor = np.array(index_arr, dtype=np.int) #Round down
    index_ceil = index_floor + 1
    index_rem = index_arr - index_floor #Remain

    val1 = original[index_floor]
    val2 = original[index_ceil % len(original)]
    interp = val1 * (1.0-index_rem) + val2 * index_rem
    assert(len(interp) == targetLen)
    return interp


def equalize(scatvars,padfac):
    if scatvars["decimation_mode"] == 0:
        Lfac = float(scatvars["decimation"])
        N_predec = float(scatvars["n_gates"]*Lfac) #predecimation number of samples
        N_fft = float(scatvars["n_gates"]*padfac)
        delta_k = 1e-10                     ##added 12/19/19+next line modified
                
        kvec =  []
        for i in range(int(N_predec)):
            kvec.append(float(i) + delta_k)
        kvec = np.asarray(kvec)
        
          
        H_k = []
        ll = []
        for ind, each in enumerate(kvec):
       
            val = ((1.0/Lfac*(1-cmath.exp(complex(0.0, -2*math.pi*kvec[ind]*Lfac/N_predec))) / (1-cmath.exp(complex(0.0,-2*math.pi*kvec[ind]/N_predec))) )) 
            H_k.append(val)
        H_k = np.asarray(H_k)
             
        
        H_k_norm0 = abs(H_k)/max(abs(H_k))
        H_k_norm1 = H_k_norm0[0:int(N_predec/(Lfac))] 
        # h_k_norm = np.resize(H_k_norm1,int(N_fft))
        h_k_norm = ResampleLinear1D(H_k_norm1, int(N_fft))
       
    else:
        h_k_norm = 1.0
        
    return h_k_norm



def to_round(val):
    if (float(val) % 1) >= 0.5:
        x = ceil(val)
    else:
        x = round(val) 
        
    return x
    
    
def read_calfile(calfile, configvars):
    if calfile != "select":
        print(f"******************")
        print(f"calibration file in use: {calfile}")
        print(f"edit name of calibration file in: Ku-(or Ka-)SCAT_IDL_config.txt if necessary")

        scatlib.waitkey(configvars)
        
    if calfile == "select":
        file_names = []
        for ind, name in enumerate(glob('calib_file' + '/*.cal')): 
            file_names.append(name)
            print(ind, " : ", name)
        calfile = file_names[int(input("choose SCAT raw data file by giving the index: "))]
        
    cal_file = open(calfile, 'r')
    

    ainv = np.zeros((2,2),dtype=np.complex_)
    finv = np.zeros((2,2),dtype=np.complex_)
    amat = np.zeros((2,2),dtype=np.complex_)
    fmat = np.zeros((2,2),dtype=np.complex_)

    cal_file_data = cal_file.readlines()
    data_ll = []
    for each in cal_file_data:
        local_ll = []
        p = re.compile('-?\d+.\d+')
        req = p.findall(each)
        for e in req:
            local_ll.append(float(e))
        data_ll.append(local_ll)
        

    ainv[0][0] = data_ll[0][0] + data_ll[0][1]*1j
    ainv[0][1] = data_ll[0][2] + data_ll[0][3]*1j
    ainv[1][0] = data_ll[1][0] + data_ll[1][1]*1j
    ainv[1][1] = data_ll[1][2] + data_ll[1][3]*1j

    finv[0][0] = data_ll[2][0] + data_ll[2][1]*1j
    finv[0][1] = data_ll[2][2] + data_ll[2][3]*1j
    finv[1][0] = data_ll[3][0] + data_ll[3][1]*1j
    finv[1][1] = data_ll[3][2] + data_ll[3][3]*1j


    reference_calibration_loop_power, corner_range,total_corner_power_vv,\
        total_corner_power_hh = data_ll[4][0], data_ll[4][1], data_ll[4][2], data_ll[4][3]
    

    amat[0][0] = data_ll[5][0] + data_ll[5][1]*1j
    amat[0][1] = data_ll[5][2] + data_ll[5][3]*1j
    amat[1][0] = data_ll[6][0] + data_ll[6][1]*1j
    amat[1][1] = data_ll[6][2] + data_ll[6][3]*1j

    fmat[0][0] = data_ll[7][0] + data_ll[7][1]*1j
    fmat[0][1] = data_ll[7][2] + data_ll[7][3]*1j
    fmat[1][0] = data_ll[8][0] + data_ll[8][1]*1j
    fmat[1][1] = data_ll[8][2] + data_ll[8][3]*1j
    
    print('receiver distortion matrix from cal file:')
    print(amat)
    print('cal file transmitter distortion matrix from cal file: ')
    print(fmat)
    print(f'reference cal loop power: {reference_calibration_loop_power}')
    print(f'corner range: {corner_range}')
    print(f'integrated corner reflector power, v: {total_corner_power_vv}')
    print(f'integrated corner reflector power, h: {total_corner_power_hh}')

    cal_file.close() #rcw added

    return ainv,finv,reference_calibration_loop_power,corner_range,total_corner_power_vv,total_corner_power_hh

    

def amat_fmat_compute():
    pass


     


def mueller_compute(c_matrix,c_matrix_vs_bin, L_matrix, L_matrix_vs_bin):
    # find L_matrix averaged over all ranges from individual s-matrix covariances

    # note on notation: s_hhs=s_hh*
    
    s_vv2 = abs(c_matrix[0,0])
    s_vv_s_vhs = c_matrix[1,0]
    s_vv_s_hvs = c_matrix[2,0]
    s_vv_s_hhs = c_matrix[3,0]

    s_vh_s_vvs = c_matrix[0,1]
    s_vh2 =  abs(c_matrix[1,1])
    s_vh_s_hvs = c_matrix[2,1]
    s_vh_s_hhs = c_matrix[3,1]

    s_hv_s_vvs = c_matrix[0,2]
    s_hv_s_vhs =c_matrix[1,2]    
    s_hv2 = abs(c_matrix[2,2])
    s_hv_s_hhs =c_matrix[3,2]

    s_hh_s_vvs = c_matrix[0,3]
    s_hh_s_vhs =c_matrix[1,3]    
    s_hh_s_hvs = c_matrix[2,3]
    s_hh2 = abs(c_matrix[3,3])  


    L_matrix[0,0] = 0.5*(s_vv2+s_hh2+s_hv2+s_vh2)
    L_matrix[1,0] = 0.5*(s_vv2-s_hh2+s_hv2-s_vh2)
    L_matrix[2,0] = s_hv_s_hhs.real + s_vv_s_vhs.real
    L_matrix[3,0] = s_hh_s_hvs.imag + s_vh_s_vvs.imag
    
    L_matrix[0,1] = .5*(s_vv2-s_hh2-s_hv2+s_vh2)
    L_matrix[1,1] = .5*(s_vv2+s_hh2-s_hv2-s_vh2)
    L_matrix[2,1] = s_vv_s_vhs.real - s_hh_s_hvs.real
    L_matrix[3,1] = s_vh_s_vvs.imag + s_hv_s_hhs.imag
    
    L_matrix[0,2] = s_vh_s_hhs.real + s_vv_s_hvs.real
    L_matrix[1,2] = s_vv_s_hvs.real - s_hh_s_vhs.real
    L_matrix[2,2] = s_vv_s_hhs.real + s_vh_s_hvs.real
    L_matrix[3,2] = s_hh_s_vvs.imag + s_vh_s_hvs.imag
    
    L_matrix[0,3] = s_vh_s_hhs.imag + s_vv_s_hvs.imag
    L_matrix[1,3] = s_vv_s_hvs.imag - s_vh_s_hhs.imag
    L_matrix[2,3] = -((s_hh_s_vvs).imag + (s_hv_s_vhs).imag)        #negative sign added 8/1/19 to agree with alt method and keep deg of pol 1 or less    
    L_matrix[3,3] =  (s_vv_s_hhs).real - (s_hv_s_vhs).real
    
    s_vv2 = abs(c_matrix_vs_bin[0,0,:])
    s_vv_s_vhs =c_matrix_vs_bin[1,0,:]
    s_vv_s_hvs =c_matrix_vs_bin[2,0,:]
    s_vv_s_hhs = c_matrix_vs_bin[3,0,:]

    s_vh_s_vvs = c_matrix_vs_bin[0,1,:]
    s_vh2 =  abs(c_matrix_vs_bin[1,1,:])
    s_vh_s_hvs = c_matrix_vs_bin[2,1,:]
    s_vh_s_hhs = c_matrix_vs_bin[3,1,:]

    s_hv_s_vvs = c_matrix_vs_bin[0,2,:]
    s_hv_s_vhs =c_matrix_vs_bin[1,2,:]    
    s_hv2 = abs(c_matrix_vs_bin[2,2,:])
    s_hv_s_hhs =c_matrix_vs_bin[3,2,:]
    
    s_hh_s_vvs = c_matrix_vs_bin[0,3,:]
    s_hh_s_vhs =c_matrix_vs_bin[1,3,:]    
    s_hh_s_hvs = c_matrix_vs_bin[2,3,:]
    s_hh2 = abs(c_matrix_vs_bin[3,3,:])     
    
    

    L_matrix_vs_bin[0,0,:] = .5*(s_vv2+s_hh2+s_hv2+s_vh2)
    L_matrix_vs_bin[1,0,:] = .5*(s_vv2-s_hh2+s_hv2-s_vh2)
    L_matrix_vs_bin[2,0,:] = (s_hv_s_hhs).real + (s_vv_s_vhs).real
    L_matrix_vs_bin[3,0,:] = (s_hh_s_hvs).imag + (s_vh_s_vvs).imag
    
    L_matrix_vs_bin[0,1,:] = .5*(s_vv2-s_hh2-s_hv2+s_vh2)
    L_matrix_vs_bin[1,1,:] = .5*(s_vv2+s_hh2-s_hv2-s_vh2)
    L_matrix_vs_bin[2,1,:] = (s_vv_s_vhs).real - (s_hh_s_hvs).real
    L_matrix_vs_bin[3,1,:] = (s_vh_s_vvs).imag + (s_hv_s_hhs).imag
    
    L_matrix_vs_bin[0,2,:] = (s_vh_s_hhs).real + (s_vv_s_hvs).real
    L_matrix_vs_bin[1,2,:] = (s_vv_s_hvs).real - (s_hh_s_vhs).real
    L_matrix_vs_bin[2,2,:] = (s_vv_s_hhs).real + (s_vh_s_hvs).real
    L_matrix_vs_bin[3,2,:] = (s_hh_s_vvs).imag + (s_vh_s_hvs).imag    
    
    L_matrix_vs_bin[0,3,:] = (s_vh_s_hhs).imag + (s_vv_s_hvs).imag
    L_matrix_vs_bin[1,3,:] = (s_vv_s_hvs).imag - (s_vh_s_hhs).imag
    L_matrix_vs_bin[2,3,:] = -((s_hh_s_vvs).imag + (s_hv_s_vhs).imag)   #neg. sign added 8/1/19 to keep degree of pol <1 and agree with alt method
    L_matrix_vs_bin[3,3,:] =  (s_vv_s_hhs).real - (s_hv_s_vhs).real
    
    
    return L_matrix,L_matrix_vs_bin

def mueller_compute_alt(smat,n_bins,n_blocks_per_group):

    r_mat = np.zeros((4,4), dtype=np.complex_)
                        #in IDL, r_mat(i,j): i=column    #  j=row
    r_mat[0,0] = 1
    r_mat[0,1] = 1
    r_mat[1,0] = 1
    r_mat[1,1] = -1
    r_mat[2,2] = 1
    r_mat[3,2] = 1
    r_mat[2,3] = complex(0,-1)
    r_mat[3,3] = complex(0,1)

    r_mat_inv = np.linalg.inv(r_mat)

    w_mat = np.zeros((4,4), dtype=np.complex_)

                                            #find L_matrix averaged over all ranges
                                            #in IDL, w_mat(i,j): i=column    #  j=row

    
    w_mat[0,0] = sum (sum(np.conj(smat[0,0,:,:]) * smat[0,0,:,:])) /n_blocks_per_group
    w_mat[1,0] = sum (sum(np.conj(smat[1,0,:,:]) * smat[1,0,:,:])) /n_blocks_per_group
    w_mat[2,0] = sum (sum(np.conj(smat[1,0,:,:]) * smat[0,0,:,:])) /n_blocks_per_group
    w_mat[3,0] = sum (sum(np.conj(smat[0,0,:,:]) * smat[1,0,:,:])) /n_blocks_per_group
    
    w_mat[0,1] = sum (sum(np.conj(smat[0,1,:,:]) * smat[0,1,:,:])) /n_blocks_per_group
    w_mat[1,1] = sum (sum(np.conj(smat[1,1,:,:]) * smat[1,1,:,:])) /n_blocks_per_group
    w_mat[2,1] = sum (sum(np.conj(smat[1,1,:,:]) * smat[0,1,:,:])) /n_blocks_per_group
    w_mat[3,1] = sum (sum(np.conj(smat[0,1,:,:]) * smat[1,1,:,:])) /n_blocks_per_group
    
    w_mat[0,2] = sum (sum(np.conj(smat[0,1,:,:]) * smat[0,0,:,:])) /n_blocks_per_group
    w_mat[1,2] = sum (sum(np.conj(smat[1,1,:,:]) * smat[1,0,:,:])) /n_blocks_per_group
    w_mat[2,2] = sum (sum(np.conj(smat[1,1,:,:]) * smat[0,0,:,:])) /n_blocks_per_group
    w_mat[3,2] = sum (sum(np.conj(smat[0,1,:,:]) * smat[1,0,:,:])) /n_blocks_per_group
    
    w_mat[0,3] = sum (sum(np.conj(smat[0,0,:,:]) * smat[0,1,:,:])) /n_blocks_per_group
    w_mat[1,3] = sum (sum(np.conj(smat[1,0,:,:]) * smat[1,1,:,:])) /n_blocks_per_group
    w_mat[2,3] = sum (sum(np.conj(smat[1,0,:,:]) * smat[0,1,:,:])) /n_blocks_per_group
    w_mat[3,3] = sum (sum(np.conj(smat[0,0,:,:]) * smat[1,1,:,:])) /n_blocks_per_group


            #the ## operator multiplies rows of first matrix by columns of second matrix

    L_matrix = (np.dot(r_mat, (np.dot(w_mat,r_mat_inv)))).real
            #find L-matrix vs bin
    L_matrix_vs_bin = np.zeros((4,4,n_bins))
    #     # in IDL, w_mat(i,j): i=column    #  j=row

    for ibin in range(0,n_bins-1):

        w_mat[0,0] = sum(np.conj(smat[0,0,:,ibin]) * smat[0,0,:,ibin]) / n_blocks_per_group
        w_mat[1,0] = sum(np.conj(smat[1,0,:,ibin]) * smat[1,0,:,ibin]) / n_blocks_per_group
        w_mat[2,0] = sum(np.conj(smat[1,0,:,ibin]) * smat[0,0,:,ibin]) / n_blocks_per_group
        w_mat[3,0] = sum(np.conj(smat[0,0,:,ibin]) * smat[1,0,:,ibin]) / n_blocks_per_group
        w_mat[0,1] = sum(np.conj(smat[0,1,:,ibin]) * smat[0,1,:,ibin]) / n_blocks_per_group
        w_mat[1,1] = sum(np.conj(smat[1,1,:,ibin]) * smat[1,1,:,ibin]) / n_blocks_per_group
        w_mat[2,1] = sum(np.conj(smat[1,1,:,ibin]) * smat[0,1,:,ibin]) / n_blocks_per_group
        w_mat[3,1] = sum(np.conj(smat[0,1,:,ibin]) * smat[1,1,:,ibin]) / n_blocks_per_group
        w_mat[0,2] = sum(np.conj(smat[0,1,:,ibin]) * smat[0,0,:,ibin]) / n_blocks_per_group
        w_mat[1,2] = sum(np.conj(smat[1,1,:,ibin]) * smat[1,0,:,ibin]) / n_blocks_per_group
        w_mat[2,2] = sum(np.conj(smat[1,1,:,ibin]) * smat[0,0,:,ibin]) / n_blocks_per_group
        w_mat[3,2] = sum(np.conj(smat[0,1,:,ibin]) * smat[1,0,:,ibin]) / n_blocks_per_group
        w_mat[0,3] = sum(np.conj(smat[0,0,:,ibin]) * smat[0,1,:,ibin]) / n_blocks_per_group
        w_mat[1,3] = sum(np.conj(smat[1,0,:,ibin]) * smat[1,1,:,ibin]) / n_blocks_per_group
        w_mat[2,3] = sum(np.conj(smat[1,0,:,ibin]) * smat[0,1,:,ibin]) / n_blocks_per_group
        w_mat[3,3] = sum(np.conj(smat[0,0,:,ibin]) * smat[1,1,:,ibin]) / n_blocks_per_group

        
        L_matrix_vs_bin[:,:,ibin] = (np.dot(r_mat,np.dot(w_mat,r_mat_inv))).real
    
    return L_matrix, L_matrix_vs_bin


    
        
def process_pol_data_in_footprint(configvars,c_matrix_vs_bin,L_matrix_vs_bin,_range,gate0,gate1,n_bins):
    cor_coe_vs_bin = c_matrix_vs_bin[0,3,:]/sqrt(c_matrix_vs_bin[0,0,:]*c_matrix_vs_bin[3,3,:])   
    cor_coe_xpol_vs_bin = c_matrix_vs_bin[1,2,:]/sqrt(c_matrix_vs_bin[1,1,:]*c_matrix_vs_bin[2,2,:])      
    
    rho_hv_vs_bin = abs(cor_coe_vs_bin)
    rho_hv_xpol_vs_bin = abs(cor_coe_xpol_vs_bin)
    
    phase_hv_vs_bin = atan((cor_coe_vs_bin).imag, (cor_coe_vs_bin).real)/(math.po/180)
    phase_hv_xpol_vs_bin = atan((cor_coe_xpol_vs_bin).imag, (cor_coe_xpol_vs_bin).real)/(math.po/180)
    
    rel_mag_vs_bin = (c_matrix_vs_bin[0,0,:] + c_matrix_vs_bin[3,3,:])/max((c_matrix_vs_bin[0,0,:]+c_matrix_vs_bin[3,3,:]))
    rel_mag_xpol_vs_bin=(c_matrix_vs_bin[1,1,:] + c_matrix_vs_bin[2,2,:])/max((c_matrix_vs_bin[0,0,:]+c_matrix_vs_bin[3,3,:]))
    LDR_vs_bin = (c_matrix_vs_bin[1,1,:] + c_matrix_vs_bin[2,2,:])/(c_matrix_vs_bin[0,0,:] + c_matrix_vs_bin[3,3,:])

    
    
def process_polarimetric_data(scatvars,configvars,calvars, _range,spec,\
    n_blocks_per_group,ngates,smoothfac,igroup, group_index, gate_offset,\
        gate_plot_max,ainv,finv,scan_index,line_el, range_centroid_signal, \
            range_peak_signal, total_power, peak_power, c_matrix_vec, \
                L_matrix_vec, rho_hv_vec, phase_hv_deg_vec, rdepol_dB_vec):
    
    #, base_filename,decon, processed_data_path, current_calibration_loop_power,reference_calibration_loop_power_file ): #rcw added base_filename, decon, 2 x cal loop powers


    c_matrix = np.zeros((4,4), dtype=np.complex_)
    L_matrix = np.zeros((4,4))
    
    avespec_vv = np.zeros(int(ngates))         #average power spectrum (power spectrum=range profile of signal power)
    avespec_hv =  np.zeros(int(ngates))
    avespec_vh =  np.zeros(int(ngates))
    avespec_hh =  np.zeros(int(ngates))
    avespec_cal = np.zeros(int(ngates))
    
    # #average range profiles 

    #print(f"processing group/line number: {igroup}")
        
    for jsweep in range(int(n_blocks_per_group)):

        isweep = int(group_index[igroup]+jsweep)
        sszz = len(avespec_vv)
    
        avespec_vv = avespec_vv + abs(spec[isweep,0,:]**2)[:sszz]
        # avespec_vh = avespec_vh + abs(spec[isweep,1,:]**2)[:sszz]
        # avespec_hv = avespec_hv + abs(spec[isweep,2,:]**2)[:sszz]
        avespec_hv = avespec_hv + abs(spec[isweep,1,:]**2)[:sszz]
        avespec_vh = avespec_vh + abs(spec[isweep,2,:]**2)[:sszz]
        avespec_hh = avespec_hh + abs(spec[isweep,3,:]**2)[:sszz]
        avespec_cal = avespec_cal + abs(spec[isweep,4,:]**2)[:sszz]

      
    avespec_vv = avespec_vv/n_blocks_per_group
    avespec_hv = avespec_hv/n_blocks_per_group
    avespec_vh = avespec_vh/n_blocks_per_group
    avespec_hh = avespec_hh/n_blocks_per_group
    avespec_cal = avespec_cal/n_blocks_per_group
    
    sumspec=(avespec_vv+avespec_hh+avespec_vh+avespec_hv)/4 
    
    
    # rcw save echoes as a text file 
    #name4sav=base_filename[0:7]+'-'+base_filename[7:]+'_igroup_'+str('{:06d}'.format(igroup))+'_line_el_'+str('{:06d}'.format(line_el))
    #import pickle
    # pickle.dump(sumspec, open(processed_data_path+name4sav+'_sumspec.p', 'wb')) #rcw
    # pickle.dump(avespec_cal, open(processed_data_path+name4sav+'_avespec_cal.p', 'wb')) #rcw
    # if decon == 0:
    #     rcw_corr_cal = 1
    # elif decon == 1:
    #     rcw_corr_cal = current_calibration_loop_power/reference_calibration_loop_power_file 
    # openw,wunit,name4sav+'_range_hh_vv_hv_vh.txt',/get_lun
    # printf,wunit,[transpose(range),transpose(avespec_hh),transpose(avespec_vv),transpose(avespec_hv),transpose(avespec_vh)],format='(f10.7,1x,g13.8,1x,g13.8,1x,g13.8,1x,g13.8)'
    # close,wunit
    # free_lun,wunit
    # f = open('/Users/rosie/Documents/mosaic_data/kukapy_test/'+name4sav+'_range_hh_vv_hv_vh.txt', 'w')
    echoes = {'range': _range,
                 'hh': avespec_hh, 
                 'vv': avespec_vv, 
                 'hv': avespec_hv, 
                 'vh': avespec_vh}
    #import pandas as pd
    ## write_df = pd.DataFrame(to_write, columns= ['range', 'hh', 'vv', 'hv', 'vh'])
    ## write_df.to_csv (processed_data_path+name4sav+'_range_hh_vv_hv_vh_decon'+str(decon)+'.txt', index = False, header=False)
    # to_write = str((np.reshape(_range, (_range.shape[0],1)),
    #         np.reshape(avespec_hh, (avespec_hh.shape[0],1)),
    #         np.reshape(avespec_vv, (avespec_vv.shape[0],1)),
    #         np.reshape(avespec_hv, (avespec_hv.shape[0],1)),
    #         np.reshape(avespec_vh, (avespec_vh.shape[0],1))))
            
    # f = open('line_el_'+str(line_el)+'_range_hh_vv_hv_vh.txt', 'w')
    # f.write(to_write)
    # f.close()
    # ;rcw
    
    
    testpeak = scatvars["pedestal_height_m"]
    #print('testpeak',testpeak)
    index = []
    count = 0
    #print('here4')
    #rcw commented out and replaced
    # for each in range(len(_range)):
    #     if _range[each]>testpeak:
    #         index.append(each)
    #         count += 1
    search_margin = int(ngates/50)
    # index_start_search = index[0]-search_margin
    
    index_start_search = min(np.where(_range > testpeak)[0]) - search_margin #rcw new
        
    if index_start_search < ngates/30:
        index_start_search = ngates/30    #avoid DC peak
    
    #print("index_start_search", index_start_search)
    #print(index_start_search)

    proc_thresh_right=10**(configvars["proc_thresh_right_dB"]/10)
    proc_thresh_left=10**(configvars["proc_thresh_left_dB"]/10)    
    
    # max_sumspec = max(sumspec[index_start_search:int(ngates)])
    # index = []
    #rcw commented and replaced
    # for each in range(len(sumspec[index_start_search:int(ngates)])):
    #     if sumspec[each] == max_sumspec: 
    #         index.append(each-index_start_search)
    index = np.where(sumspec[index_start_search:] == max(sumspec[index_start_search:]))[0] #rcw new
    
    #print('index:',index)
    gate_peak = index[0] + index_start_search

    testpeak_cal = -1.0
    #rcw commented and replaced
    # index = []
    # count = 0
    # print('here 5')
    # for each in range(len(_range)):
    #     if _range[each]>testpeak_cal:
    #         index.append(each)
    #         count += 1

    # index_start_search_cal = index[0]
    index_start_search_cal = min(np.where(_range > testpeak_cal)[0]) #rcw new
    
    index_cal = []
    #print('here 6', index_start_search_cal, ngates)
    #rcw commenting and replacing
    # for each in range(index_start_search_cal, int(ngates)):
    #     if avespec_cal[each] == max(avespec_cal[index_start_search_cal:int(ngates)]):
    #         index_cal.append(each-index_start_search_cal)
    # gate_peak_cal = index_cal[0] + index_start_search_cal
    gate_peak_cal = min(np.where(avespec_cal[index_start_search_cal:] == max(avespec_cal[index_start_search_cal:]))) + index_start_search_cal #rcw new
    #print('gate_peak_cal', gate_peak_cal)

    # iloop_left = 0
    # loop0 = 0
    pre_peak_gate = int(gate_peak*.8)   #recoded finding of peak region on 4/27/2020
    
    #index_lo_left = []
    count = 0
    
    #print('here1')   
    
    #rcw commenting and replacing diff for scan and stare modes
    #new
    if max(scan_index) == 1: #scan mode
        iloop_left = 0
        loop0 = 0
        pre_peak_gate = int(gate_peak*.8)   #recoded finding of peak region on 4/27/2020
        
        index_lo_left = []
        count = 0
                
        while count == 0:
            for each in range(pre_peak_gate, gate_peak):
                if sumspec[each]/sumspec[gate_peak] < proc_thresh_left:
                    index_lo_left.append(each - pre_peak_gate)
                    count += 1
        
        
            if count == 0:    
                proc_thresh_left = 2*proc_thresh_left
                if proc_thresh_left > 1 or loop0 > 0:
                    proc_thresh_left = 0.5+loop0/10
                loop0=loop0+1
                if iloop_left > 10:
                    print("cannot find signal region left of peak")
                    exit()
        
                iloop_left = iloop_left+1
        
        gate0=pre_peak_gate+max(index_lo_left)
        
    elif max(scan_index) != 1: #stare mode
        index_lo_left=np.where(sumspec[pre_peak_gate:gate_peak]/sumspec[gate_peak] < proc_thresh_left)[0]
        if len(index_lo_left) > 1: #multiple elements found
            gate0=pre_peak_gate+index_lo_left[-1]
        elif len(index_lo_left) == 1:
            gate0=pre_peak_gate+index_lo_left[0]
        else: #if hasn't found where to cut off on left hand side
            gate0 = 0

    if gate0 > ngates-2:
        gate0 = ngates-2
    if gate0 < 0:
        gate0 = 0
        
    #print('gate0 NEW',gate0)
    
    #old
    # iloop_left = 0
    # loop0 = 0
    # pre_peak_gate = int(gate_peak*.8)   #recoded finding of peak region on 4/27/2020
    
    # index_lo_left = []
    # count = 0
    
    # print('here1')   
    
    # while count == 0:
    #     for each in range(pre_peak_gate, gate_peak):
    #         if sumspec[each]/sumspec[gate_peak] < proc_thresh_left:
    #             index_lo_left.append(each - pre_peak_gate)
    #             count += 1
    
    
    #     if count == 0:    
    #         proc_thresh_left = 2*proc_thresh_left
    #         if proc_thresh_left > 1 or loop0 > 0:
    #             proc_thresh_left = 0.5+loop0/10
    #         loop0=loop0+1
    #         if iloop_left > 10:
    #             print("cannot find signal region left of peak")
    #             exit()
    
    #         iloop_left = iloop_left+1
    
    # gate0=pre_peak_gate+max(index_lo_left)
    
    # if gate0 > ngates-2:
    #     gate0 = ngates-2
    # if gate0 < 0:
    #     gate0 = 0
        
    # print('gate0',gate0)
    

    #new
    index_lo_right = np.where(sumspec[gate_peak:]/max(sumspec[gate_peak:]) < proc_thresh_right)[0]
    if len(index_lo_right) == 0:
        if max(scan_index) == 1: #if scan mode use previous code
            iloop_right = 0
            index_lo_right = []
            count = 0
            #print('here 3i')
            while count == 0:    
                for each in range(gate_peak, int(ngates)):
                    if sumspec[each]/max(sumspec[gate_peak:int(ngates)]) < proc_thresh_right:
                        index_lo_right.append(each-gate_peak)
                        count += 1
            
                if count == 0:
                    proc_thresh_right = 10*proc_thresh_right
                    if proc_thresh_right > 1:
                        proc_thresh_right = 0.5
            
                    if iloop_right > 10:
                        print(f"cannot find signal region right of peak")
                        exit()
            
                iloop_right=iloop_right+1 
        else: #stare mode
            #index_lo_right = [len(_range) - gate_peak] #if stare mode set gate1 to last bin
            index_lo_right = [gate_peak*-1] #rcw December 2022 set it so gate1 = 0 in this case
    #print(gate1,gate0)        
    # print('gate1 NEW', gate1)
    
    #old
    # iloop_right = 0
    # index_lo_right = []
    # count = 0
    # print('here 3i')
    # while count == 0:    
    #     for each in range(gate_peak, int(ngates)):
    #         if sumspec[each]/max(sumspec[gate_peak:int(ngates)]) < proc_thresh_right:
    #             index_lo_right.append(each-gate_peak)
    #             count += 1
    
    #     if count == 0:
    #         proc_thresh_right = 10*proc_thresh_right
    #         if proc_thresh_right > 1:
    #             proc_thresh_right = 0.5
    
    #         if iloop_right > 10:
    #             print(f"cannot find signal region right of peak")
    #             exit()
    
    #     iloop_right=iloop_right+1 
    
    # print('gate1', gate1)

        
    #print('here3a')

    gate1 = gate_peak + index_lo_right[0]
    if gate1-gate0 < 1:
        gate1 = gate0+1
        
    
    if configvars["i_corner_process"] == 1:    
        gate0 = gate_peak - calvars["n_summed_gates"]/2
        gate1 = gate0 + calvars["n_summed_gates"]-1
        
    n_bins = gate1-gate0+1
    
    print('gate0', gate0)
    print('gate1', gate1)


    if n_bins < smoothfac:      # le changed from lt on 4/13/2020
        print(f"n_bins: {n_bins}")
        print(f"n_bins less than smoothfac")
        print(f"reducing smoothing factor to n_bins/2")
        smoothfac = n_bins/2
        
    #print('here4')
    range_peak_signal[igroup] = _range[gate_peak]
    range_centroid_vv = sum(_range[gate0:gate1] * avespec_vv[gate0:gate1])/sum(avespec_vv[gate0:gate1])    
    range_centroid_hh = sum(_range[gate0:gate1] * avespec_hh[gate0:gate1])/sum(avespec_hh[gate0:gate1]) 
    range_centroid_signal[igroup] = (range_centroid_vv+range_centroid_hh)/2.
    
    #print('range_peak_signal, range_centroid_signal', range_peak_signal, range_centroid_signal)
        

    c_matrix_vs_bin = np.zeros((4,4,n_bins), dtype=np.complex_) #covariance matrix at each range within beam footprint
    L_matrix_vs_bin = np.zeros((4,4,n_bins))                        #covariance #matrix at each range within beam footprint
    vmat = np.zeros((2,2,int(n_blocks_per_group),n_bins), dtype=np.complex_)
    smat = np.zeros((2,2,int(n_blocks_per_group),n_bins), dtype=np.complex_)
    
    total_power_vv = sum(avespec_vv[gate0:gate1])
    total_power_hv = sum(avespec_hv[gate0:gate1])
    total_power_vh = sum(avespec_vh[gate0:gate1])
    total_power_hh = sum(avespec_hh[gate0:gate1])
    

    total_power[igroup,0] = total_power_vv
    total_power[igroup,1] = total_power_hv
    total_power[igroup,2] = total_power_vh
    total_power[igroup,3] = total_power_hh
    
    print('total_power[igroup]', total_power[igroup])

    peak_power_vv = max(avespec_vv[gate0:gate1])
    peak_power_hv = max(avespec_hv[gate0:gate1])
    peak_power_vh = max(avespec_vh[gate0:gate1])
    peak_power_hh = max(avespec_hh[gate0:gate1])

    peak_power[igroup,0] = peak_power_vv
    peak_power[igroup,1] = peak_power_hv
    peak_power[igroup,2] = peak_power_vh
    peak_power[igroup,3] = peak_power_hh

    if (gate0 == 0) or (gate1 == 1):
        print('igroup:', igroup, 'gate0 is 0 and gate 1 is 1')
        total_power[igroup,:] = np.nan
        peak_power[igroup,:] = np.nan
    
    
    # ## PLOT    
    
    yrmin = 10 * log10(min(avespec_vv[gate_offset:int(ngates)]))-3
    yrmax = 10*log10(max(np.hstack([avespec_vv[gate_offset:int(ngates)], avespec_hh[gate_offset:int(ngates)]]))) + 3
    gate_plot = int(gate1*1.05)
    if gate_plot > ngates:
        gate_plot=float(.2*ngates)

    if configvars['show_all_plot'] == 1:
        # plotting the points  
        plt.plot(_range[gate_offset: int(gate_plot_max)+1], 10*np.log10(avespec_vv[gate_offset: int(gate_plot_max)+1]), \
                color='black', linestyle='solid', linewidth = 0.8, 
                marker='o', markerfacecolor='blue', markersize=0.1)

        # plotting the points  
        plt.plot(_range[gate_offset: int(gate_plot_max)+1], 10*np.log10(avespec_hh[gate_offset: int(gate_plot_max)+1]), \
                color='red', linestyle='solid', linewidth = 0.8, 
                marker='o', markerfacecolor='blue', markersize=0.1)

        # plotting the points  
        plt.plot(_range[gate_offset: int(gate_plot_max)+1], 10*np.log10(avespec_vh[gate_offset: int(gate_plot_max)+1]), \
                color='blue', linestyle='solid', linewidth = 0.8, 
                marker='o', markerfacecolor='blue', markersize=0.1)
        
        # plotting the points  
        plt.plot(_range[gate_offset: int(gate_plot_max)+1], 10*np.log10(avespec_hv[gate_offset: int(gate_plot_max)+1]), \
                color='green', linestyle='solid', linewidth = 0.8, 
                marker='o', markerfacecolor='blue', markersize=0.1)
        
        plt.ylim(yrmin,yrmax)
        # naming the x axis 
        plt.xlabel('range (m)') 
        # naming the y axis 
        plt.ylabel('relative power (dB)')
        # giving a title to my graph 
        plt.title('ave. power VV blk, HH red, VH blue, HV grn')  


        db_str_vv = str(round((int(1000*np.log10(peak_power_vv)))/100, 1)) + ' dB'
        db_str_hh = str(round((int(1000*np.log10(peak_power_hh)))/100, 1)) + ' dB'
        
        
        rvec1 = [_range[gate0], _range[gate0]]
        rvec2 = [_range[gate1], _range[gate1]]
        yvec = [yrmin, yrmax]
        
        # plotting the points  
        plt.plot(rvec1, yvec, \
                color='cyan', linestyle='solid', linewidth = 0.8)
        plt.plot(rvec2, yvec, \
                color='cyan', linestyle='solid', linewidth = 0.8)
        
        if max(scan_index) == 1:
            el_ang_str = str(round(line_el, 2))

            text1 = 'VV peak: ' + db_str_vv
            text2 = 'HH peak: ' + db_str_hh
            text3 = 'elevation angle (deg): ' + el_ang_str
            extra = Rectangle((0, 0), 1, 1, fc="w", fill=False, edgecolor='none', linewidth=0)
            plt.legend([extra, extra, extra],[text1, text2, text3], loc='upper right', title='')

        else:
            text1 = 'VV peak: ' + db_str_vv
            text2 = 'HH peak: ' + db_str_hh
            extra = Rectangle((0, 0), 1, 1, fc="w", fill=False, edgecolor='none', linewidth=0)
            plt.legend([extra, extra],[text1, text2], loc='upper right', title='')
        plt.show() 

        
    if configvars['show_all_plot'] == 1:
        yrmin=max(10*np.log10(avespec_vv[gate0:(gate1+1)]))+min([configvars["proc_thresh_left_dB"],configvars["proc_thresh_right_dB"]])-5
        yrmax=max(10*np.log10(avespec_vv[gate0:(gate1+1)]))+5

        # plotting the points  
        plt.plot(_range[gate0: int(gate1)+1], 10*np.log10(avespec_vv[gate0:gate1+1]), \
                color='black', linestyle='solid', linewidth = 0.8, 
                marker='o', markerfacecolor='blue', markersize=0.1)
        
        # plotting the points  
        plt.plot(_range[gate0: int(gate1)+1], 10*np.log10(avespec_hh[gate0:gate1+1]), \
                color='red', linestyle='solid', linewidth = 0.8, 
                marker='o', markerfacecolor='blue', markersize=0.1)
        
        # plotting the points  
        plt.plot(_range[gate0: int(gate1)+1], 10*np.log10(avespec_vh[gate0:gate1+1]), \
                color='blue', linestyle='solid', linewidth = 0.8, 
                marker='o', markerfacecolor='blue', markersize=0.1)
        
        # plotting the points  
        plt.plot(_range[gate0: int(gate1)+1], 10*np.log10(avespec_hv[gate0:gate1+1]), \
                color='green', linestyle='solid', linewidth = 0.8, 
                marker='o', markerfacecolor='blue', markersize=0.1)

        gate_plot1=0.4*int(gate0+gate1)

        plt.ylim(yrmin,yrmax)
        # naming the x axis 
        plt.xlabel('range (m)') 
        # naming the y axis 
        plt.ylabel('relative power (dB)')
        # giving a title to my graph 
        plt.title('ave. power VV blk, HH red, VH blue, HV grn')  
        
        text = 'data selected for processing'
        extra = Rectangle((0, 0), 1, 1, fc="w", fill=False, edgecolor='none', linewidth=0)
        plt.legend([extra],[text], loc='lower center', title='')
        plt.show()
        
    
    for jsweep in range(0, int(n_blocks_per_group)):
        isweep = group_index[igroup] + jsweep
        vmat[0,0,jsweep,:] = spec[int(isweep), 0, gate0:(gate1+1)] #vmat(0,0) is the V V return
        vmat[1,0,jsweep,:] = spec[int(isweep), 2, gate0:(gate1+1)] #vmat(1,0) is the V pol return when transmitting H-changed 7/13/19 from vmat(0,1)
        vmat[0,1,jsweep,:] = spec[int(isweep), 1, gate0:(gate1+1)] #vmat(0,1) is the H pol return when transmitting V-changed 7/13/19 from vmat(1,0)
        vmat[1,1,jsweep,:] = spec[int(isweep), 3, gate0:(gate1+1)]       #vmat(1,1) is the H H return

    vmat = np.conj(vmat)

    # # compute complex scattering matrix, Smat# data becomes calibrated by applying ainv and finv
    for j in range(0, int(n_blocks_per_group)):
        for k in range(0, n_bins):
            smat[:,:,j,k] = finv @ vmat[:,:,j,k] @ ainv  #this step computes calibrated 2x2 scattering matrix at each range bin within gate0:gate1       

            
    
    # #compute covariance of S_vv and S_hh
    npoints = int(n_blocks_per_group)*n_bins
    # #display polarization ellipse scatter plot


    if configvars["i_pol_scat"] == 1:
        print("PLOTTING")
    

    smatvec = np.reshape(smat,(4,int(n_blocks_per_group),n_bins), order='F')

    for j in range(0,4):
        for k in range(0,4):                
            c_matrix[j,k] = sum(sum((smatvec[k,:,:]*np.conj(smatvec[j,:,:]))))/int(n_blocks_per_group) #don't divide by n_bins here...we are summing over all range bins# divide by group_index since we're averaging power in time    
    
    for m in range(0,n_bins):
        for j in range(0,4):
            for k in range(0,4):
                c_matrix_vs_bin[j,k,m] = (sum(smatvec[k,:,m]*np.conj(smatvec[j,:,m])))/int(n_blocks_per_group)
                


    
    smooth_vec = np.zeros(3)
    smooth_vec[:] = 1
    smooth_vec[2] = smoothfac

    ## Smoothing
    c_matrix_vs_bin = scipy.ndimage.filters.uniform_filter(abs(c_matrix_vs_bin),smooth_vec) #smoothing over bin (range bin or gate) dimension
       
    
    L_matrix,L_matrix_vs_bin = mueller_compute(c_matrix,c_matrix_vs_bin, L_matrix, L_matrix_vs_bin)

    c_matrix_vec[:,:,igroup] = c_matrix
    L_matrix_vec[:,:,igroup] = L_matrix


    # print(f"Mueller matrix averaged over all range bins: ")
    # print(f"L_matrix: ")
    # print(f"{L_matrix/(max(max(x) for x in L_matrix))}")
    # print(f"L_matrix vs bin for first bin: ")
    # print(L_matrix_vs_bin[:,:,0]/(max(max(x) for x in L_matrix_vs_bin[:,:,0])))

        #the following code was used to check validity of mueller_compute
        #using c_matrix formulation.  These agree when smoothfac=1
    #goto,skip_alt
    L_matrix_alt, L_matrix_vs_bin_alt = mueller_compute_alt(smat,n_bins,n_blocks_per_group)
    # print(f" Mueller matrix averaged over all range bins (alternate form): ")

    # print(f"{L_matrix_alt/(max(max(x) for x in L_matrix_alt))}")
    # print("L_matrix vs bin for first bin: ")
    # print(L_matrix_vs_bin_alt[:,:,0]/(max(max(x) for x in L_matrix_vs_bin_alt[:,:,0])))
    L_matrix = L_matrix_alt
    L_matrix_vs_bin = L_matrix_vs_bin_alt
    
    # # skip_alt:

    cor_coe = (c_matrix[0,3]) / sqrt((c_matrix[0,0]) * (c_matrix[3,3]))   
    cor_coe_xpol = (c_matrix[1,2]) / sqrt((c_matrix[1,1]) * (c_matrix[2,2]))   
    rho_hv = abs(cor_coe)
    
    phase_hv = atan2((cor_coe).imag, (cor_coe).real)

    rdepol = abs(c_matrix[2,1]) / (sqrt(abs(c_matrix[0,0])*(c_matrix[3,3])))
    rco = ((c_matrix[0,0])/(c_matrix[3,3]))
    
    
    # print(f"mag co-polarized correlation coefficient of Smat: {rho_hv}")
    # print(f"phase co-polarized correlation coefficient of Smat (degrees): {(phase_hv)/(math.pi/180)}")
    # print(f"depolarization ratio (dB): {10*log10(rdepol)}")
    # print(f"co-polarized ratio (rco) (dB): {10*log10(rco)}")
    
    rho_hv_vec[igroup] = rho_hv
    phase_hv_deg_vec[igroup] = phase_hv/(math.pi/180)

    rdepol_dB_vec[igroup] = 10*log10(rdepol)
    
    
    if max(scan_index) == 1:
        control_plots = 1
        i_multi_bin = 0
        iline_dum = igroup
    else:
        control_plots = 1
        line_el = 0                         
        i_multi_bin = 0
        iline_dum = 0

    if configvars["i_pol_signature"] == 1:
        polsig(L_matrix,iline_dum,control_plots,line_el,i_multi_bin)
  
    return c_matrix_vec,L_matrix_vec,rho_hv_vec,phase_hv_deg_vec,rdepol_dB_vec,total_power,range_peak_signal,range_centroid_signal, \
        gate_peak,gate_peak_cal, peak_power, echoes




def select_index_az(configvars, azimuth, sweep_count, elapsed_time):
    
    azmin_proc = configvars['azmin_proc']
    azmax_proc = configvars['azmax_proc']

    az_proc_index = []
    count = 0
    for i, each in enumerate(azimuth):
        if each > azmin_proc and each < azmax_proc:
            az_proc_index.append(i)
            count += 1

    if count == 0:
        print('i_az_override = 1')
        print('processing subset of azimuth angles; no azimuth angle fall in specified range')
        #exit() #rcw changed exit to return
        return #rcw changed exit to return

    else:
        sweep_count_override = []
        for each in az_proc_index:
            sweep_count_override.append(sweep_count[each])

    return az_proc_index, sweep_count_override




def plot_range_profiles(i, spec, sweep_count, ngates, gate_offset, gate_plot_max, geometric_peak, elevation, _range, configvars):

    ymin = -90
    ymax = 10
    yoff=-20
    yoff2=-25

    ngates = int(ngates)

    pspec_vv = 20*np.log10(abs(spec[i,0,0:ngates]))
    pspec_hh = 20*np.log10(abs(spec[i,3,0:ngates]))

    # ;check on ordering of subscripts 

    pspec_vh =  20*np.log10(abs(spec[i,2,0:ngates]))
    pspec_hv  = 20*np.log10(abs(spec[i,1,0:ngates]))
    pspec_cal = 20*np.log10(abs(spec[i,4,0:ngates]))



    pspec_noise = 20*np.log10(abs(spec[i,5,0:ngates]))

    phase_vv = []
    phase_hv = []
    phase_vh = []
    phase_hh = []

    for e in spec[i, 0, 0:ngates]:
        val = atan2(e.imag, e)/(math.pi/180) 
        phase_vv.append(val)
    phase_vv = np.array(phase_vv)

    for e in spec[i, 1, 0:ngates]:
        val = atan2(e.imag, e)/(math.pi/180) 
        phase_hv.append(val)
    phase_hv = np.array(phase_hv)
    
    for e in spec[i, 2, 0:ngates]:
        val = atan2(e.imag, e)/(math.pi/180) 
        phase_vh.append(val)
    phase_vh = np.array(phase_vh)
    
    for e in spec[i, 3, 0:ngates]:
        val = atan2(e.imag, e)/(math.pi/180) 
        phase_hh.append(val)
    phase_hh = np.array(phase_hh)
    

    phase_cal = []
    phase_noise = []

    for e in spec[i, 4, 0:ngates]:
        val = atan2(e.imag, e)/(math.pi/180) 
        phase_cal.append(val)
    phase_cal = np.array(phase_cal)

    for e in spec[i, 5, 0:ngates]:
        val = atan2(e.imag, e)/(math.pi/180) 
        phase_noise.append(val)
    phase_noise = np.array(phase_noise)

    phase_vv_hh = phase_vv - phase_hh


    indexvvhh = []
    countvvhh = 0
    for ind, each in enumerate(phase_vv_hh):
        if each < -180.0:
            countvvhh += 1
            indexvvhh.append(ind)

    if countvvhh > 0:
        phase_vv_hh[indexvvhh] = phase_vv_hh[indexvvhh] + 360.0


    indexvvhh = []
    countvvhh = 0
    for ind, each in enumerate(phase_vv_hh):
        if each > 180.0:
            countvvhh += 1
            indexvvhh.append(ind)

    if countvvhh > 0:
        phase_vv_hh[indexvvhh] = phase_vv_hh[indexvvhh] - 360.0


    # plotting the points  
    plt.plot(_range[gate_offset: int(gate_plot_max)+1], pspec_vv[gate_offset:int(gate_plot_max)+1], \
            color='black', linestyle='solid', linewidth = 0.8, 
            marker='o', markerfacecolor='blue', markersize=0.1)

    # plotting the points  
    plt.plot(_range[gate_offset: int(gate_plot_max)+1], pspec_hh[gate_offset:int(gate_plot_max)+1], \
            color='red', linestyle='solid', linewidth = 0.8, 
            marker='o', markerfacecolor='blue', markersize=0.1)

    # plotting the points  
    plt.plot(_range[gate_offset: int(gate_plot_max)+1], pspec_hv[gate_offset:int(gate_plot_max)+1], \
            color='blue', linestyle='solid', linewidth = 0.8, 
            marker='o', markerfacecolor='blue', markersize=0.1)
    
    # plotting the points  
    plt.plot(_range[gate_offset: int(gate_plot_max)+1], pspec_vh[gate_offset:int(gate_plot_max)+1], \
            color='green', linestyle='solid', linewidth = 0.8, 
            marker='o', markerfacecolor='blue', markersize=0.1)
    
    plt.ylim(ymin,ymax)
    # naming the x axis 
    plt.xlabel('range (m)') 
    # naming the y axis 
    plt.ylabel('power (dB)')
    # giving a title to my graph 
    plt.title('Range Profiles  VV blk, HH red, HV blue, VH grn')  

    plt.show() 



    scatlib.waitkey(configvars)
    

    to_plot = []
    for each in  pspec_vv[gate_offset:int(gate_plot_max)+1]:
        val = each - max(pspec_vv[gate_offset:int(gate_plot_max)+1])
        to_plot.append(val)
    to_plot = np.array(to_plot)


    # plotting the points  
    plt.plot(_range[gate_offset: int(gate_plot_max)+1], to_plot, \
            color='black', linestyle='solid', linewidth = 0.8, 
            marker='o', markerfacecolor='blue', markersize=0.1)


    to_plot = []
    for each in  pspec_vv[gate_offset:int(gate_plot_max)+1]:
        val = each - max(pspec_vv[gate_offset:int(gate_plot_max)+1])
        to_plot.append(val)
    to_plot = np.array(to_plot)


    # plotting the points  
    plt.plot(_range[gate_offset: int(gate_plot_max)+1], to_plot, \
            color='blue', linestyle='solid', linewidth = 0.8, 
            marker='o', markerfacecolor='blue', markersize=0.1)

    to_plot = []
    for each in  pspec_vv[5*gate_offset:int(gate_plot_max)+1]:
        val = each - max(pspec_vv[gate_offset:int(gate_plot_max)+1])
        to_plot.append(val)
    to_plot = np.array(to_plot)



    # plotting the points  
    plt.plot(_range[gate_offset: int(gate_plot_max)+1], pspec_cal[gate_offset: int(gate_plot_max)+1], \
            color='black', linestyle='solid', linewidth = 0.8, 
            marker='o', markerfacecolor='blue', markersize=0.1)


    ymin = 0
    ymax = -12
    
    plt.ylim(ymin,ymax)
    # naming the x axis 
    plt.xlabel('range (m)') 
    # naming the y axis 
    plt.ylabel('power (dB)')
    # giving a title to my graph 
    plt.title('internal cal loop')  

    plt.show() 
