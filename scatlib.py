import sys
import file_header
import numpy as np
import math
from math import *
from statistics import median
from scipy import signal
import scipy.ndimage


from pathlib import Path

import matplotlib.pyplot as plt 
from matplotlib.patches import Rectangle

from loading_from_bin import get_file_header_t, get_meta_header_t, data_gathering
import utils as U
from cmsystime import cmsystime


def waitkey(configvars):
    #rcw if configvars["i_batch"] == 1:
    #rcw     pass
    #rcw else:
    #rcw     inp = input("Enter any key to continue, q to quit: ")
    #rcw     if inp in ["q", "Q"]:
    #rcw         exit()
    #rcw     else:
    pass
        

def read_configuration(iscat):
    print("Reading SCAT_IDL_config....")

    configvars = dict()    
    if int(iscat) == 1 : 
        print('Reading Ku-SCAT_IDL_config')
        import Ku_SCAT_IDL_config as file

    if int(iscat) == 2 : 
        print('Reading Ka-SCAT_IDL_config')
        import Ka_SCAT_IDL_config as file  
    
    configvars["instrument_name"] = get_variable(file, 'instrument_name')
    configvars["show_all_plot"] = get_variable(file, 'show_all_plot')

    if iscat == 1 and configvars["instrument_name"] != 'Ku-Scat':
        print('WARNING: instrument name in Ku-SCAT_IDL_config.py must be Ku-Scat')
    elif iscat == 2 and configvars["instrument_name"] != 'Ka-Scat':
        print('WARNING: instrument name in Ka-SCAT_IDL_config.py must be Ka-Scat')


    configvars["i_sub_band"] = get_variable(file, 'i_sub_band')
    configvars["sub_bandwidth"] = get_variable(file, 'sub_bandwidth')
    configvars["sub_center_frequency"] = get_variable(file, 'sub_center_frequency')
    configvars["i_calibrate"] = get_variable(file, 'i_calibrate')
    configvars["calfile"] = get_variable(file, 'calfile')
    configvars["i_batch"] = get_variable(file, 'i_batch')
    configvars["i_raw_plot"] = get_variable(file, 'i_raw_plot')
    configvars["max_display_range"] = get_variable(file, 'max_display_range')
    configvars["i_proc_ind"] = get_variable(file, 'i_proc_ind')
    configvars["n_blocks_ind_per_group"] = get_variable(file, 'n_blocks_ind_per_group')
    configvars["i_az_override"] = get_variable(file, 'i_az_override')
    configvars["azmin_proc"] = get_variable(file, 'azmin_proc')

    configvars["azmax_proc"] = get_variable(file, 'azmax_proc')
    configvars["i_el_override"] = get_variable(file, 'i_el_override')
    configvars["elmin_proc"] = get_variable(file, 'elmin_proc')
    configvars["elmax_proc"] = get_variable(file, 'elmax_proc')

    configvars["i_cal_loop_override"] = get_variable(file, 'i_cal_loop_override')
    configvars["i_corner_cal_override"] = get_variable(file, 'i_corner_cal_override')
    configvars["i_corner_process"] = get_variable(file, 'i_corner_process')
    configvars["i_pol_scat"] = get_variable(file, 'i_pol_scat')
    configvars["i_pol_signature"] = get_variable(file, 'i_pol_signature')
    configvars["i_pol_data_in_footprint"] = get_variable(file, 'i_pol_data_in_footprint')

    configvars["smoothfac"] = get_variable(file, 'smoothfac')
    configvars["i_temp_time_plot"] = get_variable(file, 'i_temp_time_plot')
    configvars["proc_thresh_left_dB"] = get_variable(file, 'proc_thresh_left_dB')
    configvars["proc_thresh_right_dB"] = get_variable(file, 'proc_thresh_right_dB')

    configvars["group_averaging_time"] = get_variable(file, 'group_averaging_time')
    configvars["raw_data_path"] = get_variable(file, 'raw_data_path')
    configvars["processed_data_path"] = get_variable(file, 'processed_data_path')
    return configvars

def get_variable(file,  test_string):
    try:
        return file.extract_config[test_string]
    except Exception as ex:
        print("There is no variable called {} in the selected file", test_string)
    

def read_header(filename, raw_data_path, configvars):
    
    
    file_header_t, sl = get_file_header_t(filename)

    scatvars = file_header_t["config"]
    calvars = file_header_t["calibration"]
    private_config = file_header_t["private_config"]    
    
    print("AD clock frequency (MHz): ",scatvars["ad_clock_frequency_hz"]/1e06)
    print("AD trig delay (ns): ",private_config["ad_trig_delay_ns"])
    print("chirp bandwidth (MHz): ",scatvars["chirp_bandwidth_hz"]/1e06) 
    print("chirp center frequency (MHz): ",scatvars["chirp_center_frequency_hz"]/1e06)
    print("chirp amplitude dbm: ",scatvars["chirp_amplitude_dbm"])
    print("chirp width (ms): ",scatvars["chirp_width_ns"]/1e06)
    print("decimation: ",scatvars["decimation"])
    print("decimation mode: ",scatvars["decimation_mode"])
    print("frame delay (ms): ",scatvars["frame_delay_ns"]/1e06)


    print("number of time samples per sweep (n_gates): ",scatvars["n_gates"])
    print("pulse repetition period (ms): ",scatvars["prp_ns"]/1e06)
    print("nominal range resolution (c_light/2/chirp_bandwidth, m): ",scatvars["range_resolution_m"])
    print("server state, 0=idle# 1=record: ",scatvars["server_state"])
    print("chirp bandwidth (MHz): ",scatvars["chirp_bandwidth_hz"]/1e06)
    print("antenna beamwidth, (deg): ",private_config["antenna_beam_width"])
    print("corner reflector radar cross section (square meters): ",private_config["corner_reflector_sigma"])
    print("pedestal_max_speed (deg/s): ",private_config["pedestal_max_speed"])
    print("pedestal height (m): ",scatvars["pedestal_height_m"])
    
    print("=" * 70)
    
    file_header_size = file_header_t["file_header_size"]    ## find the size of the file
    
    waitkey(configvars)    
    return scatvars ,calvars ,private_config ,file_header_size, sl


def counts_to_voltage(counts):
    voltage=(counts/51.4)/1000.         #volts
    return voltage



def read_raw(configvars, scatvars, filename, file_header_size, sl): #, base_filename, decon, processed_data_path): #rcw added base_filename, decon, processed_data_path
    
    meta_header_size = 146                                  ## find the size of the file
    print("Meta Header Size", meta_header_size)
    
    n_bytes_per_samp=2
    n_pol=6              
    
    block_size = scatvars["n_gates"]*n_pol*n_bytes_per_samp 
    
    file_size = Path(filename).stat().st_size
    
    n_blocks=int((file_size-file_header_size)/(meta_header_size+block_size))

    raw = np.zeros([n_blocks,n_pol,scatvars["n_gates"]])    
    time_sec = np.zeros(n_blocks)
    twelve_volts = np.zeros(n_blocks)   
    input_twelve_volts = np.zeros(n_blocks)   
    five_volts = np.zeros(n_blocks)
    minus_five_volts = np.zeros(n_blocks)
    RF_unit_plate_temp = np.zeros(n_blocks)
    LCD_display_temp = np.zeros(n_blocks)
    power_supply_plate_temp = np.zeros(n_blocks)
    power_supply_plate_temp = np.zeros(n_blocks)
    elevation = np.zeros(n_blocks)
    azimuth = np.zeros(n_blocks)
    sweep_count = np.zeros(n_blocks)
    scan_index = np.zeros(n_blocks)
    transition_flag = np.zeros(n_blocks)
    along_track_tilt = np.zeros(n_blocks)
    cross_track_tilt = np.zeros(n_blocks)
    block_to_block_delta = np.zeros([n_blocks,6])
    gps_latitude = np.zeros(n_blocks)
    gps_longitude = np.zeros(n_blocks)
    plo_20_ghz_lock = np.zeros(n_blocks)
    

    if configvars["i_raw_plot"] > 1:
        print("then window,0,re=2 >>>>>>>>>>>>>>>>>> PLOTTING" )
    if configvars["i_raw_plot"] > 1:
        print("then wshow,0 >>>>>>>>>>>>>>>>>> PLOTTING" )
    if configvars["i_corner_process"] == 1:
        print("i_corner_process=1 assuming you are processing a corner reflector file")
        print("change i_corner_process to 0 in Ku/Ka_SCAT_IDL_config.txt to process regular data")
        waitkey(configvars)

    
    # ## Loop begins
    for i in range(n_blocks):
        meta_header_t, sl = get_meta_header_t(filename, sl)
        
        time_sec[i]= (meta_header_t["timestamp_seconds"])+(meta_header_t["timestamp_nanoseconds"])/1e09    

        
        if i == 0:
            print("start time: ",meta_header_t["timestamp_seconds"], " seconds")
            print("nano-seconds: ",meta_header_t["timestamp_nanoseconds"])

        RF_unit_plate_temp[i] = meta_header_t["project_status"]["rcb_status"]["rf_plate_temp"]
        LCD_display_temp[i]=meta_header_t["project_status"]["rcb_status"]["lcd_display_temp"]
        power_supply_plate_temp[i]=meta_header_t["project_status"]["rcb_status"]["power_supply_plate_temp"]
        azimuth[i] = meta_header_t["project_status"]["ped_status"]["az_pos"]
        elevation[i] = meta_header_t["project_status"]["ped_status"]["el_pos"]
        sweep_count[i]=meta_header_t["project_status"]["ped_status"]["sweep_count"]
        scan_index[i]=meta_header_t["project_status"]["scan_index"]
        transition_flag[i]=meta_header_t["project_status"]["ped_status"]["transition_flag"]
        along_track_tilt[i]=meta_header_t["project_status"]["rcb_status"]["along_track_tilt"]
        cross_track_tilt[i]=meta_header_t["project_status"]["rcb_status"]["cross_track_tilt"]
        gps_longitude[i]=meta_header_t["project_status"]["gps_longitude"]
        gps_latitude[i]=meta_header_t["project_status"]["gps_latitude"]
        plo_20_ghz_lock[i]=meta_header_t["project_status"]["rcb_status"]["plo_20_ghz_lock"]
        five_volts[i]=meta_header_t["project_status"]["rcb_status"]["five_volts"]
        minus_five_volts[i]=meta_header_t["project_status"]["rcb_status"]["minus_five_volts"]
        twelve_volts[i]=meta_header_t["project_status"]["rcb_status"]["twelve_volts"]
        input_twelve_volts[i]=meta_header_t["project_status"]["rcb_status"]["input_twelve_volts"]

        ## cmsystem 
        if (i%100)==0:
            print("reading block number: ",i," time: ", cmsystime(time_sec[i]) )  ## cmsystem_time# will do it later      

        szz = int(block_size/2/n_pol)

        gathered_data, sl = data_gathering(filename, sl, szz)

        hh = gathered_data["hh"]
        vv = gathered_data["vv"]
        hv = gathered_data["hv"]
        vh = gathered_data["vh"]
        cal = gathered_data["cal"]
        noise = gathered_data["noise"]
        
        #rcw adding code to save out the raw data
        # save_folder = '/Users/rosie/Documents/mosaic/mosaic_data/kukapy_output/'
        # name4sav=base_filename[0:7]+'-'+base_filename[7:]+'_block'+str(i)
        # np.savetxt(name4sav+"_hh.csv", hh, delimiter=" ")
        # np.savetxt(name4sav+"_vv.csv", vv, delimiter=" ")
        # np.savetxt(name4sav+"_hv.csv", hv, delimiter=" ")
        # np.savetxt(name4sav+"_vh.csv", vh, delimiter=" ")
        #end rcw
                
        raw[i,0,:] = counts_to_voltage(vv)
        raw[i,1,:] = counts_to_voltage(hv)
        raw[i,2,:] = counts_to_voltage(vh)
        raw[i,3,:] = counts_to_voltage(hh)
        raw[i,4,:] = counts_to_voltage(cal)
        raw[i,5,:] = counts_to_voltage(noise)
        
        pol_label = list()
        pol_label.append('raw data, VV')
        pol_label.append('raw data, HV')
        pol_label.append('raw data, VH')
        pol_label.append('raw data, HH')
        pol_label.append('raw data, CAL')
        pol_label.append('raw data, noise')
        
        if configvars["i_raw_plot"] > 1 and i == 0:
            ymax = max(raw)
            ymin = min(raw)
            
            print(f"block: {i}")
            
#     #rcw adding code to save out the raw arrays and deconvolve if wanted
# #    save_folder = '/Users/rosie/Documents/mosaic/mosaic_data/kukapy_output/'
#     save_folder = processed_data_path

#     if decon == 0: 
#         for raw_layer_number in range(6):
#             #save and plot raw data
#             name4sav=base_filename[0:7]+'-'+base_filename[7:]+'_layer'+str(raw_layer_number)+'_raw_decon0_t0py'
#             np.savetxt(save_folder+name4sav+".csv", raw[:,raw_layer_number,:], delimiter=" ")
#             plt.imshow(raw[:,raw_layer_number,:])
#             plt.title('scatlib plot of raw'+str(raw_layer_number))
#             plt.show()    
        
#     if decon == 1:
#         #deconvolve raw data
#         import deconvolveTomModRosie as deconvolve
#         raw = deconvolve.runDeconvolution(raw, base_filename[0:2])        
#         #save and plot deconvolved raw data
#         for raw_layer_number in range(6):
#             name4sav=base_filename[0:7]+'-'+base_filename[7:]+'_layer'+str(raw_layer_number)+'_raw_decon1_t0py'
#             np.savetxt(save_folder+name4sav+".csv", raw[:,raw_layer_number,:], delimiter=" ")       
#             plt.imshow(raw[:,raw_layer_number,:])
#             plt.title('scatlib plot of raw_decon'+str(raw_layer_number))
#             plt.show()          
#     #end rcw
            
    index = []
    for ind, each in enumerate(time_sec):
        if each != 0: 
            index.append(ind)

            
    elapsed_time = []
    for ind, each in enumerate(index):
        elapsed_time.append(time_sec[each] - time_sec[0])
    elapsed_time = np.asarray(elapsed_time)
    
    
    if configvars["i_temp_time_plot"] == 1:
        print("Some Plotting.....")

        waitkey(configvars)
        
    independent_sample_index = np.asarray([i for i in range(n_blocks)])
    

    if max(scan_index) != 1:
        distance, independent_sample_index = U.compute_distance(gps_latitude,gps_longitude,elapsed_time,scatvars,configvars,independent_sample_index)
    else:
        distance = 0
        
    if configvars["i_az_override"] == 1:
        az_proc_index, sweep_count_override = U.select_index_az(configvars, azimuth, sweep_count, elapsed_time)
        sweep_count_override = [int(i) for i in sweep_count_override]
    else:
        az_proc_index, sweep_count_override = 0, 0


    return raw, scan_index, sweep_count, transition_flag, elevation, \
          n_blocks, n_pol, elapsed_time,time_sec,gps_latitude,gps_longitude,\
            along_track_tilt,cross_track_tilt, independent_sample_index,distance,\
                az_proc_index, sweep_count_override
    
        

def compute_range_profiles(n_blocks,n_pol,elevation,scan_index,\
                            sweep_count,raw,scatvars,private_config, \
                            configvars):

    if configvars["i_raw_plot"] > 1:
        print("Some plotting")
        
    max_display_range=configvars["max_display_range"]
    
    raw_length = len(raw[0,0,:])
    fac = log10(raw_length)/log10(2)
    
    print(f"raw length {raw_length}")
    print(f"fac: {fac}")
    
    upfac = 2
    if raw_length > 4096:
        upfac = 1
    ipad=1                  # set to 1 to use zero padding
    
    if ipad == 1:
        fftlen = (2.**(floor(fac) + upfac)) 
        print(f"fft length: {fftlen}")
        sizepad = fftlen-raw_length
        zeros = np.zeros(int(sizepad))
        padfac=float(fftlen)/raw_length     #used for digital rx block-averaging equalization 
    else:
        fftlen = raw_length
        padfac = 1
        print(f"NOTICE: no zero padding.  Just use for testing")
        inp = input(f"enter c to continue")
        if inp != 'c' or inp != 'C':
            exit()
        

    if configvars["instrument_name"] == 'Ku-Scat':
        range_offset_hh = -.01
        range_offset_vh = 0
        range_offset_hv = 0
        ped_el_offset_inches = 12.6065
        hyp_inches = 46
        hyp_phi_deg = 0
    
    if configvars["instrument_name"] == 'Ka-Scat':
        range_offset_hh = 0
        range_offset_vh = 0
        range_offset_hv = 0
        ped_el_offset_inches = 12.6065
        hyp_inches = 46
        hyp_phi_deg = 0


    range_gate_spacing = scatvars["range_resolution_m"]*raw_length/fftlen
    radar_range_offset = private_config["range_to_antenna_m"]
    gate_offset = int(radar_range_offset/range_gate_spacing)
    
    gate_offset_hv = round(range_offset_hv/range_gate_spacing)
    gate_offset_vh = round(range_offset_vh/range_gate_spacing)
    gate_offset_hh = round(range_offset_hh/range_gate_spacing)   

    print(f"gate_offset_hv: {gate_offset_hv} ")
    print(f"gate_offset_vh: {gate_offset_vh} ")
    print(f"gate_offset_hh: {gate_offset_hh} ")

    ngates = fftlen/2
    gate_plot_max = ngates-1
    
    print(f"raw range resolution: {scatvars['range_resolution_m']}, range_gate_spacing with zero pad: {range_gate_spacing}")

    spec = np.zeros((n_blocks,n_pol,int(fftlen)), dtype=np.complex_)

    _range = []
    for i in range(int(ngates)):
        _range.append((i * range_gate_spacing)-radar_range_offset)
    _range = np.asarray(_range)
    print(f"maximum range (m): {max(_range)}")
    
    if max_display_range > max(_range):
        max_display_range = max(_range)

    itest = [] 
    countr = 0
    for ind, each in enumerate(range(int(max_display_range))):
        if each < _range[ind]:
            itest.append(ind)
            countr += 1

    if countr > 0:
        gate_plot_max=itest[0]
    print(f"maximum displayed range (m): {max_display_range}")

    weights = np.hanning(scatvars["n_gates"])
        
    if configvars["i_sub_band"] == 1:
        sub_band_fac = configvars["sub_bandwidth"]/(scatvars["chirp_bandwidth_hz"]/1e06)
        weights_sub = np.hanning(scatvars["n_gates"]*sub_band_fac)
        sweep_start_freq = (scatvars["chirp_center_frequency_hz"]-scatvars["chirp_bandwidth_hz"]/2)/1e6
        sweep_end_freq = sweep_start_freq + scatvars["chirp_bandwidth_hz"]/1e06
        sweep_bandwidth = scatvars["chirp_bandwidth_hz"]/1e06
        sub_band_start_freq = configvars["sub_center_frequency"]-configvars["sub_bandwidth"]/2.
        weights_start_index = int(scatvars["n_gates"]*(sub_band_start_freq-sweep_start_freq)/(sweep_bandwidth))
        weights_end_index = int(weights_start_index+len(weights_sub))
        weights = np.zeros(scatvars["n_gates"])
        print(len(weights), len(weights_sub), weights_start_index, weights_end_index)
        weights[weights_start_index:weights_end_index] = weights_sub

    range_index = np.zeros(n_blocks)

    # # #set up variables to compute antenna height

    ped_el_offset_m = ped_el_offset_inches*.0254  
    hyp_m = hyp_inches*.0254
    hyp_phi = radians(hyp_phi_deg)

    spec_coh_ave = np.zeros((n_blocks,n_pol,int(fftlen)), dtype=np.complex_)
    spec_coh_ave0 = np.zeros((n_pol,int(fftlen)), dtype=np.complex_)


    print(f"number of records: {n_blocks}")
    
    
    if scatvars["decimation_mode"] == 0:
        H_k_norm = U.equalize(scatvars,padfac)
    else:
        H_k_norm=1.0

                                        ##compute antenna phase center height
    ped_elevation_deg = np.asarray(elevation)
    pos_height = scatvars["pedestal_height_m"]

       
    height = []
    for i in ped_elevation_deg:
        v = pos_height + ped_el_offset_m + hyp_m*sin(radians(i)-hyp_phi)
        height.append(v)
    height = np.asarray(height)
    
    geometric_peak = []
    for i, each in enumerate(height):
        geometric_peak.append(height[i]/cos(radians(elevation[i])))
    geometric_peak = np.asarray(geometric_peak)
    
    
    if gate_plot_max > ngates-1:
        gate_plot_max = ngates-1
    i_suppress = 0
    sweep_count_old = 0
    print(f"processing data into range profiles: ")
    if configvars["i_raw_plot"] != 0:
        print("Printing record numbers....")
        

    z_ohms = 50.0                     
                                    ## #pscale preserves peak power
    Watts_to_mW = 1000.
    pscale = Watts_to_mW*2.0*(fftlen/U.to_round(np.sum(weights)))**2/z_ohms
    vscale = sqrt(pscale)


    for i in range(n_blocks):
        if (i%100) == 0 or i==0:
            print(f"processing record number: {i}")
        for j in range(4):
            if ipad == 1:
                rd = weights*np.squeeze(raw[i,j,:]) # zeros
                rawpad = np.zeros(len(rd) + len(zeros))
                rawpad[:len(rd)] = rd
            else:
                rawpad = weights*np.squeeze(raw[i,j,:])
            rawpad = np.array(rawpad)      

            spec[i,j,:] = (np.fft.ifft(rawpad)-spec_coh_ave[i,j,:])*vscale/H_k_norm       #subtract sky noise from ground data and equalize result to account for block averaging in digital rx
                       #correct for fine offset in range 
        spec[i,1,:] = np.roll(spec[i,1,:],gate_offset_hv)
        spec[i,2,:] = np.roll(spec[i,2,:],gate_offset_vh)
        spec[i,3,:] = np.roll(spec[i,3,:],gate_offset_hh) 

      
        for j in range(4, 5):
            if ipad == 1:
                rd = weights*np.squeeze(raw[i,j,:]) # zeros
                rawpad = np.zeros(len(rd) + len(zeros))
                rawpad[:len(rd)] = rd
            else:
                rawpad = weights*np.squeeze(raw[i,j,:])
            rawpad = np.array(rawpad)
            spec[i,j,:] = np.fft.ifft(rawpad)*vscale/H_k_norm           

        
        if configvars["i_raw_plot"] == 1 and scan_index[i] == 1: 
            if sweep_count[i] != sweep_count_old:
                U.plot_range_profiles(i, spec, sweep_count, ngates, gate_offset, gate_plot_max, geometric_peak, elevation, _range, configvars)
        
                
                
        if (configvars["i_raw_plot"] == 2 and abs(scan_index[i]) == 1) or configvars["i_raw_plot"] == 3:            
            U.plot_range_profiles(i, spec, sweep_count, ngates, gate_offset, gate_plot_max, geometric_peak, elevation, _range, configvars)

        sweep_count_old=sweep_count[i]
        
    
    pll = abs(spec[0,4,gate_offset:int(ngates-1)])

    pspec_cal = []
    for each in pll:
        pspec_cal.append(20*log10(each))
    pspec_cal = np.asarray(pspec_cal)
    
    pspec_cal1 = 20*np.log10(pll)
    
    print('max(pspec_cal1 - pspec_cal)', max(pspec_cal1 - pspec_cal))
    
    gate_max_cal = []
    for ind, each in enumerate(pspec_cal):
        if each == max(pspec_cal):
            gate_max_cal.append(gate_offset + ind)
    gate_max_cal = np.asarray(gate_max_cal)[0]
    
    gate_max_cal1 = np.where(pspec_cal == max(pspec_cal))[0]
    if len(gate_max_cal1) > 1:
        gate_max_cal1[0]
    gate_max_cal1 = np.array(gate_max_cal1 + gate_offset)
    
    print('gate_max_cal1 - gate_max_cal', gate_max_cal1 - gate_max_cal)
    
    index_surface_vec = []
    count = 0
    for ind, each in enumerate(scan_index):
        if each == 1:
            index_surface_vec.append(ind)
            count += 1
    index_surface_vec = np.asarray(index_surface_vec)
    # index_surface = index_surface_vec[0]
     
    return range_gate_spacing,ngates,gate_offset,gate_plot_max,pos_height,height, \
        _range,spec,gate_max_cal,geometric_peak 
        
        
        
def calibrate(calvars,\
                configvars,gate_offset,gate_plot_max,spec,ngates,n_blocks,n_pol,\
                    _range,gate_max_cal, range_gate_spacing,scatvars,scan_index,sweep_count,\
                        elevation,base_filename,calfile):

    calpwr_i = abs(spec[0:n_blocks-1,4,gate_max_cal])**2 
    median_cal_pwr = median(calpwr_i)
    calpwr_i_db = []
    for ind, each in enumerate(calpwr_i):
        calpwr_i_db.append(10*log10(calpwr_i[ind]))
    median_cal_pwr_db = 10*log10(median_cal_pwr)
    
    current_calibration_loop_power = sum(calpwr_i)/n_blocks    
    print(f"mean calibration power {10*log10(current_calibration_loop_power)} dB")
    
    
    if configvars["i_calibrate"] == 0:
        ainv,finv,reference_calibration_loop_power,corner_range,total_corner_power_vv,total_corner_power_hh = U.read_calfile(calfile, configvars)
        print(f"mean calibration power stored in calibration file {10*log10(reference_calibration_loop_power)} dB")
        


    if configvars["i_calibrate"] == 1:
        print("Did not completed the part when configvars is => 1")
        print("Please exit the code")
        waitkey(configvars)
        
    return reference_calibration_loop_power, current_calibration_loop_power,ainv,finv, corner_range, \
            total_corner_power_vv,total_corner_power_hh


def polarimetric_processing_stare_by_independent_samples(scatvars,\
            configvars,calvars,ngates,range_gate_spacing,_range,gate_offset,gate_plot_max,spec,ainv,finv,\
            n_blocks, elapsed_time, scan_index, independent_sample_index, distance, base_filename):
    
    n_blocks_ind = len(independent_sample_index)
    n_bins_max = len(_range)

    n_groups = round(n_blocks_ind/configvars['n_blocks_ind_per_group'])   #eventually, couple this to time and/or distance tranvelled
    n_blocks_per_group = configvars['n_blocks_ind_per_group']

    c_matrix_vec = np.zeros((4,4,int(n_groups)), dtype=np.complex_)
    L_matrix_vec = np.zeros((4,4,n_groups))   #Mueller matrix for average of all gates
    total_power = np.zeros((n_groups,4))      #stores VV/HV/VH/HH power from ground target
    peak_power = np.zeros((n_groups,4))       #stores VV/HV/VH/HH peak power from ground target
    rho_hv_vec = np.zeros((n_groups))
    phase_hv_deg_vec = np.zeros((n_groups))
    rdepol_dB_vec = np.zeros((n_groups))
    range_peak_signal = np.zeros((n_groups))
    range_centroid_signal = np.zeros((n_groups))
    

    ngates = int(ngates)

    hh_image0 = np.zeros((n_blocks_ind,ngates))
    vh_image0 = np.zeros((n_blocks_ind,ngates))
    hv_image0 = np.zeros((n_blocks_ind,ngates))
    vv_image0 = np.zeros((n_blocks_ind,ngates))
    
    group_index = []
    for each in range(n_groups):
        group_index.append(each*n_blocks_per_group) 
        
    smoothfac=configvars["smoothfac"]


    if smoothfac < 1:
        smoothfac=1
    if configvars["i_corner_process"] == 1 and smoothfac != 1:
        smoothfac=1
        print('setting smoothing factor to 1 when processing corner reflector data')
        waitkey(configvars)
    

    independent_sample_index = [int(i) for i in independent_sample_index]
    # print(independent_sample_index)
    # exit()

    spec_independent = spec[independent_sample_index,:,:]
    line_el = 0


    # print(group_index)
    # print(n_blocks_per_group)
    # print(n_groups)
    # print(spec_independent.shape)
    # exit()
    
    #rcw create dict to contain echoes
    echoes_all = {'range': _range,
                'hh': np.full((len(_range), n_groups), np.nan),
                'vh': np.full((len(_range), n_groups), np.nan),
                'hv': np.full((len(_range), n_groups), np.nan), 
                'vv': np.full((len(_range), n_groups), np.nan)}  
    #end rcw                             

    for igroup in range(0,n_groups):
        print('igroup: ', igroup)

        try:
            c_matrix_vec,L_matrix_vec,rho_hv_vec,phase_hv_deg_vec,rdepol_dB_vec,total_power,range_peak_signal,range_centroid_signal, \
            gate_peak,gate_peak_cal, peak_power, echoes = U.process_polarimetric_data(scatvars,configvars,calvars,_range,spec_independent,n_blocks_per_group,ngates,smoothfac, \
            igroup,group_index,gate_offset,gate_plot_max, ainv,finv, scan_index,line_el, range_centroid_signal, range_peak_signal,\
            total_power,peak_power, c_matrix_vec,L_matrix_vec,rho_hv_vec,phase_hv_deg_vec,rdepol_dB_vec) 
                
            hh_image0[group_index[igroup]:group_index[igroup]+n_blocks_per_group,:] = 20*np.log10(abs(spec_independent[group_index[igroup]:group_index[igroup] + n_blocks_per_group,3,0:ngates]))
            vv_image0[group_index[igroup]:group_index[igroup]+n_blocks_per_group,:] = 20*np.log10(abs(spec_independent[group_index[igroup]:group_index[igroup] + n_blocks_per_group,0,0:ngates]))
            hv_image0[group_index[igroup]:group_index[igroup]+n_blocks_per_group,:] = 20*np.log10(abs(spec_independent[group_index[igroup]:group_index[igroup] + n_blocks_per_group,1,0:ngates]))
            vh_image0[group_index[igroup]:group_index[igroup]+n_blocks_per_group,:] = 20*np.log10(abs(spec_independent[group_index[igroup]:group_index[igroup] + n_blocks_per_group,2,0:ngates]))
            
            #rcw fill echoes_all
            echoes_all['hh'][:,igroup] = echoes['hh']
            echoes_all['vh'][:,igroup] = echoes['vh']
            echoes_all['hv'][:,igroup] = echoes['hv']
            echoes_all['vv'][:,igroup] = echoes['vv']
            #end rcw
        except:
            pass
    # exit() 

    HH_image = hh_image0[0:max(group_index)+n_blocks_per_group-1,:]
    VV_image = vv_image0[0:max(group_index)+n_blocks_per_group-1,:]
    HV_image = hv_image0[0:max(group_index)+n_blocks_per_group-1,:]
    VH_image = vh_image0[0:max(group_index)+n_blocks_per_group-1,:]
    LDR_image = (VH_image+HV_image)/2-(HH_image+VV_image)/2
    HHVV_ave_image = (HH_image+VV_image)/2    
    
    min_range = 1.0
    max_range = 3.0
    index_range_min = []
    for ind, each in enumerate(_range):
        if each > min_range:
            index_range_min.append(ind)

    index_range_max = []
    for ind, each in enumerate(_range):
        if each > max_range:
            index_range_max.append(ind)

    distance_r = []
    for e in independent_sample_index:
        distance_r.append(distance[e])
    
    x_range = [min(distance_r), max(distance_r)]
    y_range = [_range[index_range_max[0]], _range[index_range_min[0]]]

    xlabel = "approximate distance along transect (m)"
    ylabel = "radar range for antenna (m)"
    smoothfac = 1

    # flag = True
    # while flag == True:
    #     flag = False

    hh_image_s = scipy.ndimage.filters.uniform_filter(HH_image[:, index_range_min[0]:index_range_max[0]], size=[smoothfac,1])
    hh_image_smooth = np.resize(hh_image_s, (200, 200))


    hv_image_s = scipy.ndimage.filters.uniform_filter(HV_image[:, index_range_min[0]:index_range_max[0]], size=[smoothfac,1])
    hv_image_smooth = np.resize(hv_image_s, (200, 200))


    vh_image_s = scipy.ndimage.filters.uniform_filter(VH_image[:, index_range_min[0]:index_range_max[0]], size=[smoothfac,1])
    vh_image_smooth = np.resize(vh_image_s, (200, 200))


    vv_image_s = scipy.ndimage.filters.uniform_filter(VV_image[:, index_range_min[0]:index_range_max[0]], size=[smoothfac,1])
    vv_image_smooth = np.resize(vv_image_s, (200, 200))

    hhvv_ave_image_s = scipy.ndimage.filters.uniform_filter(HHVV_ave_image[:, index_range_min[0]:index_range_max[0]], size=[smoothfac,1])
    hhvv_ave_image_smooth = np.resize(hhvv_ave_image_s, (200, 200))

    LDR_image_s = scipy.ndimage.filters.uniform_filter(LDR_image[:, index_range_min[0]:index_range_max[0]], size=[smoothfac,1])
    LDR_image_smooth = np.resize(LDR_image_s, (200, 200))



    data_max = np.amax(HH_image)
    data_min = data_max - 50.0
    ct_string = 'relative power (db)'

    vplot = hh_image_smooth[::-1]
    title_str = 'HH power along transect'
    print("Graph is not loading...Enter 99 to skip")
    ## Image plot

    vplot = vv_image_smooth[::-1]
    title_str = 'VV power along transect'
    ## Image plot

    ct_strin = 'LDR (dB)'
    data_max = 5
    data_min = -40
    vplot = LDR_image_smooth[::-1]
    title_str = 'Linear Deplorization Ratio along transect'
    ## Image plot

    if configvars["i_batch"] == 0:
        smoothfac = int(input("Enter smoothing factor for 2D imaage, 99 to exit: "))
        smoothfac = 99
    else:
        smoothfac = 99

    # if smoothfac != 99:
    #     flag = True

    if configvars['show_all_plot'] == 1:
        
        if n_groups > 1:
            # plotting the points  
            plt.plot(rho_hv_vec)
            plt.xlabel('group number') 
            plt.ylabel('correlation coeff. mag.')
            plt.title('VV/HH cor. coefficient vs. group')          
            plt.show()
            
            # plotting the points  
            plt.plot(phase_hv_deg_vec)
            plt.xlabel('group number') 
            plt.ylabel('correlation coeff. phase (deg)')
            plt.title('phase VV/HH cor. coe. vs group')          
            plt.show()
            
            
            # plotting the points  
            plt.plot(rdepol_dB_vec)
            plt.xlabel('group number') 
            plt.ylabel('depolarization ratio (dB)')
            plt.title('depolarization vs elevation group')          
            plt.show()

    return range_peak_signal,line_el,n_groups,rho_hv_vec,\
            phase_hv_deg_vec,total_power,range_centroid_signal, L_matrix_vec,c_matrix_vec, \
                group_index, n_blocks_per_group, echoes_all #rcw added echoes_all



def polarimetric_processing_stare(scatvars,configvars,calvars,ngates,range_gate_spacing,\
                    _range,gate_offset,gate_plot_max,spec,ainv,finv,n_blocks, elapsed_time,\
                    scan_index):#, base_filename, decon, processed_data_path,\
                        #current_calibration_loop_power,reference_calibration_loop_power_file): #rcw addded base_filename, decon, 2x cal powers


    n_bins_max = len(_range)
    
    # n_groups = int(round(max(elapsed_time)/configvars["group_averaging_time"]))   #eventually, couple this to time and/or distance tranvelled #rcw
    # n_blocks_per_group = round(n_blocks/n_groups) #rcw
    n_groups = n_blocks   #rcw changed from above two lines 
    n_blocks_per_group = 1 #rcw changed from above two lines
   
    n_groups = int(n_groups)
    c_matrix_vec = np.zeros((4,4,n_groups), dtype=np.complex_)
    
    L_matrix_vec = np.zeros((4,4,n_groups))   #Mueller matrix for average of all gates
    total_power = np.zeros((n_groups,4))      #stores VV/HV/VH/HH power from ground target
    peak_power = np.zeros((n_groups,4))       #stores VV/HV/VH/HH peak power from ground target
    rho_hv_vec = np.zeros((n_groups))
    phase_hv_deg_vec = np.zeros((n_groups))
    rdepol_dB_vec = np.zeros((n_groups))
    range_peak_signal = np.zeros((n_groups))
    range_centroid_signal = np.zeros((n_groups))


    group_index = []
    for each in range(n_groups):
        group_index.append(each*n_blocks_per_group) 
        
        
    smoothfac = configvars["smoothfac"]

    if smoothfac < 1:
        smoothfac=1
    if configvars["i_corner_process"] == 1 and smoothfac != 1:
        smoothfac=1
        print('setting smoothing factor to 1 when processing corner reflector data')
        waitkey(configvars)

    line_el = 0

    
    # print(n_bins_max)
    # print(n_groups)
    # print(n_blocks_per_group)
    # exit()
    
        
    #rcw create dict to contain echoes
    echoes_all = {'range': _range,
                'hh': np.full((len(_range), n_groups), np.nan),
                'vh': np.full((len(_range), n_groups), np.nan),
                'hv': np.full((len(_range), n_groups), np.nan), 
                'vv': np.full((len(_range), n_groups), np.nan)}  
    #end rcw                             

    for  igroup  in range(n_groups):
        c_matrix_vec,L_matrix_vec,rho_hv_vec,phase_hv_deg_vec,rdepol_dB_vec,total_power,range_peak_signal,range_centroid_signal, \
        gate_peak,gate_peak_cal, peak_power, echoes = U.process_polarimetric_data(scatvars,configvars,calvars,_range,spec,n_blocks_per_group,ngates,smoothfac, igroup,group_index,gate_offset, gate_plot_max, ainv,finv,scan_index,\
            line_el, range_centroid_signal,range_peak_signal, total_power,\
            peak_power, c_matrix_vec,L_matrix_vec,rho_hv_vec,phase_hv_deg_vec,rdepol_dB_vec)#,
            #base_filename, decon, processed_data_path, current_calibration_loop_power,reference_calibration_loop_power_file) #rcw added echoes and base_filename onwards
       
        #rcw fill echoes_all
        echoes_all['hh'][:,igroup] = echoes['hh']
        echoes_all['vh'][:,igroup] = echoes['vh']
        echoes_all['hv'][:,igroup] = echoes['hv']
        echoes_all['vv'][:,igroup] = echoes['vv']
        #end rcw
    ## Plotting

    if configvars['show_all_plot'] == 1:
        
        if n_groups > 1:
            # plotting the points  
            plt.plot(rho_hv_vec)
            plt.xlabel('elevation angle (deg)') 
            plt.ylabel('correlation coeff. mag.')
            plt.title('VV/HH cor. coefficient vs. elevation angle')          
            plt.show()
            
            # plotting the points  
            plt.plot(phase_hv_deg_vec)
            plt.xlabel('elevation angle (deg)') 
            plt.ylabel('correlation coeff. phase (deg)')
            plt.title('phase VV/HH cor. coe. vs elevation')          
            plt.show()
            
            
            # plotting the points  
            plt.plot(rdepol_dB_vec)
            plt.xlabel('elevation angle (deg)') 
            plt.ylabel('depolarization ratio (dB)')
            plt.title('depolarization vs elevation angle')          
            plt.show()
    

    return  range_peak_signal,line_el,n_groups,rho_hv_vec,\
            phase_hv_deg_vec,total_power,range_centroid_signal, L_matrix_vec,c_matrix_vec, \
                group_index, n_blocks_per_group, echoes_all #rcw added echoes_all
            
def polarimetric_processing_scan(scatvars,configvars,\
            calvars,geometric_peak,scan_index,\
            sweep_count,transition_flag,elevation,ngates,\
            range_gate_spacing,_range,height,gate_offset,\
            gate_plot_max,spec,ainv,finv, az_proc_index,\
                sweep_count_override):


    nlines = max(sweep_count)+1         # number of surface scan lines
    if nlines == 0:
        nlines=1
    start_line=0
    
    
    if configvars["i_el_override"] == 1:
        sweep_count_max = int(max(sweep_count))
                
        elvals = np.zeros(sweep_count_max+1)
        

        for i in range(0,sweep_count_max+1):
            index, count = [], 0
            for each in range(len(sweep_count)):

                if sweep_count[each] == i and transition_flag[each] == 0 and scan_index[each] == 1:
                    index.append(each)
                    count += 1

            each_elevation = []
            for ii in index:
                each_elevation.append(elevation[ii])
            
            elvals[i] = U.to_round(median(each_elevation))

        index0 = []
        count0 = 0
        for i in range(len(elvals)):
            if elvals[i] >= configvars["elmin_proc"]:
                index0.append(i)
                count0 += 1

        if count0 > 0:
            start_line=index0[0]
        else:
            print("elmin_proc is greater than all elevation angles")
            exit()
        
        index1 = []
        count1 = 0
        for i in range(len(elvals)):
            if elvals[i] > configvars["elmax_proc"]:
                index1.append(i)
                count1 += 1

        if count1 > 0:
            stop_line = index1[0]-1
        else:
            print("elmax_proc is greater than all elevation angles# processing up to max elevation")
            stop_line = sweep_count_max
        nlines = stop_line-start_line+1

    if configvars["i_corner_process"] == 1:
        nlines=1                # only look at corner reflector data if i_corner_process=1
     
    c_matrix_vec = np.zeros((4,4,int(nlines)), dtype=np.complex_) #covariance matrix
    L_matrix_vec = np.zeros((4,4,int(nlines)))     #Mueller matrix
    total_power=np.zeros((int(nlines),4))          #stores VV/HV/VH/HH power from ground target
    rho_hv_vec=np.zeros(int(nlines))
    phase_hv_deg_vec=np.zeros(int(nlines))
    rdepol_dB_vec=np.zeros(int(nlines))
    line_index = np.zeros(int(nlines))   #sweep number at the beginning of lines having data (even numbered lines have data)
    block_count = np.zeros(int(nlines))  # number of VV/HV/VH/HH/cal/noise data blocks per sweep 
    line_elevation = np.zeros(int(nlines)) #antenna elevation angle
    line_height = np.zeros(int(nlines))
    range_peak_signal = np.zeros(int(nlines))
    range_centroid_signal = np.zeros(int(nlines))
    total_power = np.zeros((int(nlines),4))  #stores VV/HV/VH/HH power from ground target
    peak_power = np.zeros((int(nlines),4))   #stores VV/HV/VH/HH peak power from ground target



    if nlines > 1 or configvars["i_corner_process"] == 0:
        for iline0 in range(0, int(nlines)):
            iline = iline0 + start_line
            index = []
            count = 0
            for each in range(len(sweep_count)):
                if sweep_count[each] == iline and transition_flag[each] == 0 and scan_index[each] == 1:
                    index.append(each)
                    count += 1
            

            if configvars['i_az_override'] == 1:
                count = 0
                index0 = []
                for ind, each in enumerate(sweep_count_override):
                    if each == iline:
                        index0.append(ind)
                        count += 1

                index = []
                for each in index0:
                    index.append(az_proc_index[each])

            if count == 0:
                count=1
                index[0] = len(elevation)-1
                
            line_index[iline0] = index[0]
            block_count[iline0] = count
            line_elevation[iline0] = elevation[index[0]]
            line_height[iline0] = height[index[0]]          
    else:
        ## Not verified as this is not the condition
        index = []
        count = 0
        for each in range(scan_index):
            if scan_index[each] == -1:
                index.append(each)
                count += 1

        if count == 0:
            count = 1
            index[0] = 0
            
        line_index[0] = index[0]
        block_count[0] = 1         ##count
        line_elevation[0] = elevation[index[0]]


    smoothfac = configvars["smoothfac"]
    
    if smoothfac < 1:
        smoothfac = 1
    if configvars["i_corner_process"] == 1 and smoothfac != 1:
        smoothfac=1
        print('setting smoothing factor to 1 when processing corner reflector data')
        waitkey(configvars)
        
    #rcw create dict to contain echoes
    echoes_all = {'range': _range,
                'hh': np.full((len(_range), int(nlines)), np.nan),
                'vh': np.full((len(_range), int(nlines)), np.nan),
                'hv': np.full((len(_range), int(nlines)), np.nan), 
                'vv': np.full((len(_range), int(nlines)), np.nan)}  
    #end rcw                             


    for iline in range(int(nlines)):
        n_blocks_per_line = block_count[iline]
        line_el = line_elevation[iline]

        c_matrix_vec,L_matrix_vec,rho_hv_vec,phase_hv_deg_vec,rdepol_dB_vec,total_power,range_peak_signal,range_centroid_signal, \
        gate_peak,gate_peak_cal, peak_power, echoes = U.process_polarimetric_data(scatvars,configvars,calvars,\
                        _range,spec,n_blocks_per_line,ngates,smoothfac,iline, line_index, gate_offset,gate_plot_max,ainv,finv,scan_index,line_el, \
                            range_centroid_signal, range_peak_signal, total_power, peak_power, c_matrix_vec, L_matrix_vec, rho_hv_vec, phase_hv_deg_vec, rdepol_dB_vec)#, current_calibration_loop_power,reference_calibration_loop_power_file) 
    
        #rcw fill echoes_all
        echoes_all['hh'][:,iline] = echoes['hh']
        echoes_all['vh'][:,iline] = echoes['vh']
        echoes_all['hv'][:,iline] = echoes['hv']
        echoes_all['vv'][:,iline] = echoes['vv']
        #end rcw
    
    waitkey(configvars)
    
    
    if configvars['show_all_plot'] == 1:
    
        if nlines > 1:
            # plotting the points  
            plt.plot(line_elevation, rho_hv_vec)
            plt.xlabel('elevation angle (deg)') 
            plt.ylabel('correlation coeff. mag.')
            plt.title('VV/HH cor. coefficient vs. elevation angle')          
            plt.show()
            
            # plotting the points  
            plt.plot(line_elevation, phase_hv_deg_vec)
            plt.xlabel('elevation angle (deg)') 
            plt.ylabel('correlation coeff. phase (deg)')
            plt.title('phase VV/HH cor. coe. vs elevation')          
            plt.show()
            
            
            # plotting the points  
            plt.plot(line_elevation, rdepol_dB_vec)
            plt.xlabel('elevation angle (deg)') 
            plt.ylabel('depolarization ratio (dB)')
            plt.title('depolarization vs elevation angle')          
            plt.show()
        
    
    waitkey(configvars)
    
    return range_peak_signal, line_elevation,line_height, nlines,rho_hv_vec,phase_hv_deg_vec,\
            total_power,range_centroid_signal, L_matrix_vec, c_matrix_vec, echoes_all #rcw addedd echoes_all
            
            
def nrcs_compute_stare(scatvars,calvars,configvars,private_config,c_matrix,range_peak_signal,\
        range_centroid_signal,corner_range,pos_height,n_groups,reference_calibration_loop_power,current_calibration_loop_power,\
        L_matrix,processed_data_filename,rho_hv_vec,phase_hv_deg_vec,total_power,total_corner_power_vv,total_corner_power_hh,\
        calfile,time_sec,gps_latitude,gps_longitude,group_index,along_track_tilt,cross_track_tilt,n_blocks_per_group,independent_sample_index):
    
    
    total_corner_power_vv_file = 10**(calvars["corner_reflector_vv_power_dbm"]/10.)
    total_corner_power_hh_file = 10**(calvars["corner_reflector_hh_power_dbm"]/10.)
    corner_range_file = calvars["corner_reflector_range_m"]
    reference_calibration_loop_power_file = 10**(calvars["cal_peak_dbm"]/10.) 
    
    jump_to = False

    while jump_to != False:
        if abs(calvars["corner_reflector_vv_power_dbm"] - 10*np.log10(total_corner_power_vv)) > 0.1 or \
        abs(calvars["corner_reflector_hh_power_dbm"] - 10*np.log10(total_corner_power_hh)) > 0.1:

            
            print('corner powers differ between current calibration file and values in data file: ')
            print('data file values: ')
            print('corner reflector power vv: ',total_corner_power_vv_file )
            print('corner reflector power hh: ',total_corner_power_hh_file )
            print('corner reflector range: ',corner_range_file)
            print('calibration file values: ')
            print('corner reflector power vv: ',total_corner_power_vv) 
            print('corner reflector power hh: ',total_corner_power_hh )
            print('corner reflector range: ',corner_range)
            print('this is just a heads up.')
            
            
            waitkey(configvars)
    
    pi = math.pi
    ln_2 = np.log10(2)/np.log10(exp(1))

    nrcs = np.zeros((4,n_groups))
    antenna_beamwidth_rad = private_config["antenna_beam_width"]*(pi/180)
    

    corr_cal = current_calibration_loop_power/reference_calibration_loop_power_file
    print('******************************')
    print('reference calibration power from data file (dB):',10*np.log10(reference_calibration_loop_power_file))
    print('reference calibration power from cal file(dB):',10*np.log10(reference_calibration_loop_power))
    print('current calibration power (dB): ',10*np.log10(current_calibration_loop_power))
    print('drift correction factor (linear): ',corr_cal)
    print('*****************************')
    
     
    corner_sigma = private_config["corner_reflector_sigma"]
    scale_factor = np.zeros((n_groups))



    psi = []
    for each in range(len(along_track_tilt)):
        psi.append( atan(sqrt(((tan(along_track_tilt[each]*(pi/180)))**2) + tan(cross_track_tilt[each]*(pi/180))**2)) )   


    for igroup in range(0,n_groups):

        ss = [i for i in range(int(n_blocks_per_group))] 
        vec = [i+group_index[igroup] for i in ss]

        find_mean = []
        for each in vec:
            find_mean.append(psi[int(each)])
            
        mean_tilt = np.mean(find_mean)
        
        scale_factor[igroup] = 8*ln_2*range_centroid_signal[igroup]**2*corner_sigma*cos(mean_tilt)/(pi*corner_range_file**4*antenna_beamwidth_rad**2)/corr_cal 
        
        if configvars["i_corner_process"] == 1:
            ## Have to check
            scale_factor[igroup] = range_peak_signal[igroup]**4*corner_sigma/corner_range_file**4
                  
        for j in range(0,4):
            if j == 0:
                cr_power = total_corner_power_vv_file
            if j == 1 or j == 2:
                cr_power = (total_corner_power_vv_file+total_corner_power_hh_file)/2
            if j == 3:
                cr_power = total_corner_power_hh_file
        
            pratio = total_power[igroup,j]/cr_power
            nrcs[j,igroup] = scale_factor[igroup]*pratio



    nrcs_db = 10*np.log10(nrcs)
    if n_groups > 1:
        pass
    
        # wset,0
        # plot,nrcs_dB(0,*),yrange=[-60,10],xtitle='group number',ytitle='NRCS, dB m2/m2',charsize=1.6,title='NRCS VV blk# HH red# VH blue# HV grn#',xs=1
        # oplot,nrcs_dB(1,*),color=80,line=2
        # oplot,nrcs_dB(2,*),color=40,line=2
        # oplot,nrcs_dB(3,*),color=160
        # wset,1
        # plot,rho_hv_vec,xtitle='group number',ytitle='correlation coeff. mag.',title='VV/HH cor. coefficient vs. group',charsize=1.4  
        # wset,2
        # plot,phase_hv_deg_vec,xtitle='group number',ytitle='correlation coeff. phase (deg)',title='phase VV/HH cor. coe. vs group',charsize=1.4,yr=[-180,180],ys=1
    else:
        print('NRCS (dB): ')
        print(nrcs_db)
        print('scale factor: ',scale_factor)


    #store processed data
    extension='.stare'
    if configvars["i_corner_process"] == 1:
        extension='.stare.corner'


    file = open(processed_data_filename+extension, 'w')
    
    print('summary data stored in: ',processed_data_filename+extension)
    
    file.write('number of data groups: ')
    file.write('\n')
    file.write(str(n_groups))
    file.write('\n')
    file.write('\n')


    file.write('n_blocks_per_group')
    file.write('\n')
    file.write(str(n_blocks_per_group))
    file.write('\n')
    file.write('\n')


    file.write('range peak signal: ')
    file.write('\n')
    ii = 1
    for each in range_peak_signal:
        file.write(str(round(each, 4)) + ' ') 
        if ii == 6:
            file.write('\n')
            ii = 1
        else:
            ii += 1
    # file.write(str(range_peak_signal))
    file.write('\n')
    file.write('\n')


    file.write('range centroid signal: ')
    file.write('\n')
    ii = 1
    for each in range_centroid_signal:
        file.write(str(round(each, 4)) + ' ') 
        if ii == 6:
            file.write('\n')
            ii = 1
        else:
            ii += 1
    # file.write(str(range_centroid_signal))
    file.write('\n')
    file.write('\n')


    file.write('corner reflector sigma (square meters): ')
    file.write('\n')
    file.write(str(corner_sigma))
    file.write('\n')
    file.write('\n')


    file.write('6 dB two-way antenna beamwidth (deg): ')
    file.write('\n')
    file.write(str(private_config["antenna_beam_width"]))
    file.write('\n')
    file.write('\n')


    file.write('current cal loop power (dB): ')
    file.write('\n')
    file.write(str(round(10*np.log10(current_calibration_loop_power), 4)))
    file.write('\n')
    file.write('\n')


    file.write('reference cal loop power (dB): ')
    file.write('\n')
    file.write(str(round(10*np.log10(reference_calibration_loop_power), 4)))
    file.write('\n')
    file.write('\n')


    file.write('corner range (m): ')
    file.write('\n')
    file.write(str(corner_range))
    file.write('\n')
    file.write('\n')


    for igroup in range(0,n_groups):

        try:

            ss = [i for i in range(n_blocks_per_group)] 
            vec0 = [i+group_index[igroup] for i in ss]
            vec=independent_sample_index[vec0]
                    
            file.write('group number: ')
            file.write('\n')
            file.write(str(igroup))
            file.write('\n')
            file.write('\n')
            file.write('group index: ')
            file.write('\n')
            file.write(str(group_index[igroup]))
            file.write('\n')
            file.write('\n')
            file.write('start time of group: ')
            file.write('\n')
            file.write(str(round(time_sec[group_index[igroup]], 4)))
            file.write('\n')
            file.write('\n')

            file.write('end time of group: ')
            file.write('\n')
            file.write(str(round(time_sec[group_index[igroup]]+n_blocks_per_group-1, 4)))
            file.write('\n')
            file.write('\n')


            file.write('mean latitude: ')
            file.write('\n')
            # file.write(str(np.mean(gps_latitude[vec])))
            required = []
            for each in vec:
                required.append(gps_latitude[int(each)])
            file.write(str(round(np.mean(required), 7))) #rcw increased precision from 4 to 7
            file.write('\n')
            file.write('\n')


            file.write('mean longitude: ')
            file.write('\n')
            # file.write(str(np.mean(gps_longitude[vec])))
            required = []
            for each in vec:
                required.append(gps_longitude[int(each)])
            file.write(str(round(np.mean(required), 7))) #rcw increased precision from 4 to 7
            file.write('\n')
            file.write('\n')


            file.write('mean along track tilt: ')
            file.write('\n')
            # file.write(str(np.mean(along_track_tilt[vec])))
            required = []
            for each in vec:
                required.append(along_track_tilt[int(each)])
            file.write(str(round(np.mean(required), 4)))
            file.write('\n')
            file.write('\n')



            file.write('standard deviation along track tilt: ')
            file.write('\n')
            # file.write(str(np.std(along_track_tilt[vec])))
            required = []
            for each in vec:
                required.append(along_track_tilt[int(each)])
            file.write(str(round(np.std(required), 4)))
            file.write('\n')
            file.write('\n')


            file.write('mean cross track tilt: ')
            file.write('\n')
            # file.write(str(np.mean(cross_track_tilt(vec))))
            required = []
            for each in vec:
                required.append(cross_track_tilt[int(each)])
            file.write(str(round(np.mean(required), 4)))
            file.write('\n')
            file.write('\n')


            file.write('standard deviation cross track tilt: ')
            file.write('\n')
            # file.write(str(np.std[cross_track_tilt[vec]]))
            required = []
            for each in vec:
                required.append(cross_track_tilt[int(each)])
            file.write(str(round(np.std(required), 4)))
            file.write('\n')
            file.write('\n')

            file.write('mean tilt relative to vertical axis: ')
            file.write('\n')
            # file.write(str(np.mean(psi(vec)/(pi/180))))
            required = []
            for each in vec:
                required.append(psi[int(each)]/(pi/180))
            file.write(str(round(np.mean(required), 4))) 

            file.write('\n')
            file.write('\n')
            file.write('standard deviation of tilt relative to vertical axis: ')
            file.write('\n')
            # file.write(str(np.std(psi(vec)/(pi/180))))
            required = []
            for each in vec:
                required.append(psi[int(each)]/(pi/180))
            file.write(str(round(np.std(required), 4))) 

            file.write('\n')
            file.write('\n')
            file.write('Mueller matrix ')
            file.write('\n')

            for i in range(0,4):  
                # L_matrix_norm = L_matrix[:,:,igroup]/max(L_matrix[:,:,igroup])
                L_matrix_norm = L_matrix[:,:,igroup]/max(max(l_matr) for l_matr in L_matrix[:,:,igroup])
                file.write(str(round(L_matrix_norm[0, i], 6)) + ' ')
                file.write(str(round(L_matrix_norm[1, i], 6)) + ' ')
                file.write(str(round(L_matrix_norm[2, i], 6)) + ' ')
                file.write(str(round(L_matrix_norm[3, i], 6)) + ' ')
                file.write('\n')
            
            
            ## Covariance
            file.write('\n')
            file.write('\n')
            file.write('Normalized Covariance matrix ')
            file.write('\n')
            for i in range(0,4):      
                c_matrix_norm = c_matrix[:,:,igroup]/max(max(abs(c_matr)) for c_matr in c_matrix[:,:,igroup])

                real0 = round((c_matrix_norm[0,i]).real, 5)
                real1 = round((c_matrix_norm[1,i]).real, 5)
                real2 = round((c_matrix_norm[2,i]).real, 5)
                real3 = round((c_matrix_norm[3,i]).real, 5)
                imag0 = round((c_matrix_norm[0,i]).imag, 5)
                imag1 = round((c_matrix_norm[1,i]).imag, 5)
                imag2 = round((c_matrix_norm[2,i]).imag, 5)
                imag3 = round((c_matrix_norm[3,i]).imag, 5)

                print_0 = '(' + str(real0) + ', ' + str(imag0) + ')' + ' ' + \
                        '(' + str(real1) + ', ' + str(imag1) + ')' + ' ' + \
                        '(' + str(real2) + ', ' + str(imag2) + ')' + ' ' + \
                        '(' + str(real3) + ', ' + str(imag3) + ')' 
                
                
                file.writelines(print_0)
                file.writelines('\n')

                    
            file.writelines('\n')
            file.writelines('\n')

            file.write('magnitude of HH/VV correlation coefficient: ')
            file.writelines('\n')
            file.write(str(round(rho_hv_vec[igroup], 5)))
            file.writelines('\n')
            file.writelines('\n')

            file.write('phase of HH/VV correlation coefficient (deg): ')
            file.writelines('\n')
            file.write(str(round(phase_hv_deg_vec[igroup], 5)))
            file.writelines('\n')
            file.writelines('\n')

        
            if configvars["i_corner_process"] == 0:
                file.write('Normalized Radar Cross Section (dBm2/m2) columns: VV, HV, VH, HH# rows: elevations angles : ')
            if configvars["i_corner_process"] == 1:
                file.write('Radar Cross Section (dBm2): VV, HV, VH, HH: ')
            

            file.writelines('\n')

            for each in nrcs_db[:,igroup]:
                file.write(str(round(each, 5)) + '  ')  

            file.write('\n')
            file.write('\n')
        
        except:
            pass
    
    file.writelines('\n')
    file.writelines('\n')
    file.write('calibration filename')
    file.writelines('\n')
    file.write(calfile)
    file.close()

    print('done')
    
    return nrcs_db

            
            
def nrcs_compute_scan(scatvars,calvars,configvars,private_config,c_matrix,range_peak_signal,\
                    range_centroid_signal,corner_range,pos_height,line_elevation,nlines,reference_calibration_loop_power,\
                    current_calibration_loop_power,L_matrix,processed_data_filename,rho_hv_vec,\
                    phase_hv_deg_vec,total_power,line_height,total_corner_power_vv,total_corner_power_hh,gps_latitude,gps_longitude,calfile):
    
    
    total_corner_power_vv_file = 10**(calvars['corner_reflector_vv_power_dbm']/10.)
    total_corner_power_hh_file= 10**(calvars['corner_reflector_hh_power_dbm']/10.)
    reference_calibration_loop_power_file = 10**(calvars['cal_peak_dbm']/10.) 
    corner_range_file = calvars['corner_reflector_range_m']
    
    if (abs(calvars['corner_reflector_vv_power_dbm'] - 10*log10(total_corner_power_vv))) > 0.1 or \
        (abs(calvars['corner_reflector_hh_power_dbm'] - 10*log10(total_corner_power_hh))) > 0.1:    
            
        print(f"corner powers differ between current calibration file and values in data file:")
        print(f"data file values:")
        print(f"corner reflector power vv: {total_corner_power_vv_file}") 
        print(f"corner reflector power hh: {total_corner_power_hh_file}") 
        print(f"corner reflector range: {corner_range_file}")
        print(f"calibration file values: ")
        print(f"corner reflector power vv: {total_corner_power_vv}") 
        print(f"corner reflector power hh: {total_corner_power_hh}") 
        print(f"corner reflector range: {corner_range}")
        print(f"this is just a heads up.")
    
        waitkey(configvars)
        
        
    ## PLot
    
    pi = math.pi
    ln_2 = log10(2)/log10(exp(1))
    nlines = int(nlines)
    print(line_elevation)
    print(type(line_elevation))
    print(line_elevation.shape)
    theta_i = [radians(i) for i in line_elevation]
    nrcs = np.zeros((4,nlines))
    nrcs_alt = np.zeros((4,nlines))
    antenna_beamwidth_rad = radians(private_config["antenna_beam_width"])

    corr_cal = current_calibration_loop_power/reference_calibration_loop_power_file

    if configvars['i_cal_loop_override'] == 0:
        corr_cal=1.0
    
    print(f"******************************")
    print(f"reference calibration power from data file (dB): {10*log10(reference_calibration_loop_power_file)}")
    print(f"reference calibration power from cal file(dB): {10*log10(reference_calibration_loop_power)}")  
    print(f"current calibration power (dB):  {10*log10(current_calibration_loop_power)}")
    print(f"drift correction factor (linear): {corr_cal}")
    
    scale_factor = np.zeros(nlines)
    scale_factor_alt = np.zeros(nlines)
    corner_sigma = private_config["corner_reflector_sigma"]



    for iline in range(nlines):


        scale_factor[iline] = 8*ln_2*range_centroid_signal[iline]**2*corner_sigma*cos(theta_i[iline])/(pi*corner_range_file**4*antenna_beamwidth_rad**2)/corr_cal # 

        if configvars['i_corner_process'] == 1:
            scale_factor[iline] = range_peak_signal[iline]**4 * corner_sigma/corner_range_file**4
        
        # print("scale_factor")
        # print(scale_factor)
        # print(scale_factor.shape)
        # exit()

        for j in range(0, 4):
            if j == 0:
                cr_power = total_corner_power_vv_file
            if j == 1 or j == 2:
                cr_power = (total_corner_power_vv_file+total_corner_power_hh_file)/2
            if j == 3:
                cr_power = total_corner_power_hh_file

            pratio = total_power[iline,j]/cr_power        #float(c_matrix(j,j,iline))            
            nrcs[j,iline] = scale_factor[iline]*pratio

            # print("pratio")
            # print(pratio)
            # print("cr_power")
            # print(cr_power)
            # print("scale_factor")
            # print(scale_factor)
            # exit()



    # print(nrcs[0, 0])
    # print(nrcs[1, 2])
    # print(nrcs[2, 4])
    # print(nrcs[3, 11])

    # print(nrcs.shape)
    # exit()

    nrcs_cross_pol = (nrcs[1,:] + nrcs[2,:])/2.0
        
    nrcs_db = []
    for nn in nrcs:
        nndb = []
        for n in nn:
            nndb.append(10*log10(n))
        nrcs_db.append(nndb)      
       
    
    # #store processed data
    extension='.scan'
    if configvars['i_corner_process'] == 1:
        extension='.scan.corner'
    
    print(f"summary data stored in: {processed_data_filename+extension}")
    
    
    out_file = open(processed_data_filename+extension, 'w')
    
    ## Number of elevation angles
    out_file.writelines('number of elevation angles: ')
    out_file.writelines('\n')
    out_file.writelines(str(nlines))
    out_file.writelines('\n')
    out_file.writelines('\n')
    
    ## Elevation angles
    out_file.writelines('elevation angles (deg)')
    out_file.writelines('\n')
    to_print = ''
    for le in line_elevation:
        to_print += str(round(le, 5)) + '\t'        
    out_file.writelines(to_print)
    out_file.writelines('\n')
    out_file.writelines('\n')
    
    
    
    out_file.writelines('range peak signal (m): ')
    out_file.writelines('\n')
    to_print = ''
    for le in range_peak_signal:
        to_print += str(round(le, 5)) +'\t'  
    out_file.writelines(to_print)
    out_file.writelines('\n')
    out_file.writelines('\n')
    
    
    out_file.writelines('corner reflector sigma (square meters): ')
    out_file.writelines('\n')
    out_file.writelines(str(round(corner_sigma, 5)))
    out_file.writelines('\n')
    out_file.writelines('\n')
    
    out_file.writelines('6 dB two-way antenna beamwidth (deg): ')
    out_file.writelines('\n')
    out_file.writelines(str(round(private_config['antenna_beam_width'], 5)))
    out_file.writelines('\n')
    out_file.writelines('\n')


    out_file.writelines('current cal loop power (dB): ')
    out_file.writelines('\n')
    out_file.writelines(str(round(10*log10(current_calibration_loop_power), 5)))
    out_file.writelines('\n')
    out_file.writelines('\n')


    out_file.writelines('reference cal loop power (dB): ')
    out_file.writelines('\n')
    out_file.writelines(str(round(10*log10(reference_calibration_loop_power), 5)))
    out_file.writelines('\n')
    out_file.writelines('\n')
    
    out_file.writelines('corner range (m): ')
    out_file.writelines('\n')
    out_file.writelines(str(round(corner_range, 5)))
    out_file.writelines('\n')
    out_file.writelines('\n')
    
    out_file.writelines('latitude: ')
    out_file.writelines('\n')
    out_file.writelines(str(round(np.mean(gps_latitude), 5)))
    out_file.writelines('\n')
    out_file.writelines('\n')
    
    out_file.writelines('longitude: ')
    out_file.writelines('\n')
    out_file.writelines(str(round(np.mean(gps_longitude), 5)))
    out_file.writelines('\n')
    out_file.writelines('\n')
   
    
    
    
    out_file.writelines('Mueller matrix for various elevation angles: ')
    out_file.writelines('\n')
    
        
    for iline in range(0,nlines):
        for i in range(0,4):      
            L_matrix_norm = L_matrix[:,:,iline]/max(max(l_matr) for l_matr in L_matrix[:,:,iline])
            
            if i == 0:
                out_file.write(str(round(L_matrix_norm[0,i], 5)))
                out_file.write('\t')
                out_file.write(str(round(L_matrix_norm[1,i], 5)))
                out_file.write('\t')
                out_file.write(str(round(L_matrix_norm[2,i], 5)))
                out_file.write('\t')
                out_file.write(str(round(L_matrix_norm[3,i], 5)))

            if i == 1:
                out_file.write(str(round(L_matrix_norm[0,i], 5)))
                out_file.write('\t')
                out_file.write(str(round(L_matrix_norm[1,i], 5)))
                out_file.write('\t')
                out_file.write(str(round(L_matrix_norm[2,i], 5)))
                out_file.write('\t')
                out_file.write(str(round(L_matrix_norm[3,i], 5)))
                
            if i == 2:
                out_file.write(str(round(L_matrix_norm[0,i], 5)))
                out_file.write('\t')
                out_file.write(str(round(L_matrix_norm[1,i], 5)))
                out_file.write('\t')
                out_file.write(str(round(L_matrix_norm[2,i], 5)))
                out_file.write('\t')
                out_file.write(str(round(L_matrix_norm[3,i], 5)))

            if i == 3:
                out_file.write(str(round(L_matrix_norm[0,i], 5)))
                out_file.write('\t')
                out_file.write(str(round(L_matrix_norm[1,i], 5)))
                out_file.write('\t')
                out_file.write(str(round(L_matrix_norm[2,i], 5)))
                out_file.write('\t')
                out_file.write(str(round(L_matrix_norm[3,i], 5)))

            out_file.writelines('\n')    
        out_file.writelines('\n')
        out_file.writelines('\n')

    out_file.writelines('\n')
    out_file.writelines('\n')
    out_file.writelines('Covariance matrix for various elevation angles: ')
    out_file.writelines('\n')
    
    
    for iline in range(0,nlines):
        for i in range(0,4):      
            c_matrix_norm = c_matrix[:,:,iline]/max(max(abs(c_matr)) for c_matr in c_matrix[:,:,iline])

            real0 = round((c_matrix_norm[0,i]).real, 5)
            real1 = round((c_matrix_norm[1,i]).real, 5)
            real2 = round((c_matrix_norm[2,i]).real, 5)
            real3 = round((c_matrix_norm[3,i]).real, 5)
            imag0 = round((c_matrix_norm[0,i]).imag, 5)
            imag1 = round((c_matrix_norm[1,i]).imag, 5)
            imag2 = round((c_matrix_norm[2,i]).imag, 5)
            imag3 = round((c_matrix_norm[3,i]).imag, 5)

            print_0 = '(' + str(real0) + ', ' + str(imag0) + ')' + ' ' + \
                      '(' + str(real1) + ', ' + str(imag1) + ')' + ' ' + \
                      '(' + str(real2) + ', ' + str(imag2) + ')' + ' ' + \
                      '(' + str(real3) + ', ' + str(imag3) + ')' 
            
            
            out_file.writelines(print_0)
            out_file.writelines('\n')

            
        out_file.writelines('\n')
            
        
    out_file.writelines('\n')
 

    out_file.writelines('magnitude of HH/VV correlation coefficient: ')
    out_file.writelines('\n')
    to_print = ''
    for ind, le in enumerate(rho_hv_vec):
        to_print += str(round(le, 5)) +'\t'
          
    out_file.writelines(to_print)
    out_file.writelines('\n')
    out_file.writelines('\n')

    
    out_file.writelines('phase of HH/VV correlation coefficient (deg): ')
    out_file.writelines('\n')
    to_print = ''
    for ind, le in enumerate(phase_hv_deg_vec):
        to_print += str(round(le, 5)) +'\t'
  
    out_file.writelines(to_print)
    out_file.writelines('\n')
    out_file.writelines('\n')

    
    if configvars['i_corner_process'] == 0:
        out_file.writelines('Normalized Radar Cross Section (dbm2/m2) columns: VV, HV, VH, HH# rows: elevations angles : ')
    if configvars['i_corner_process'] == 1:
        out_file.writelines('Radar Cross Section (dbm2): VV, HV, VH, HH: ')

    out_file.writelines('\n')
    to_print = ''
    
    nrcs_db = np.asarray(nrcs_db).transpose()

    # print(nrcs_db)
    # exit()
    
    for m in range(nrcs_db.shape[0]):
        for n in range(nrcs_db.shape[1]):
            to_print += str(round(nrcs_db[m][n], 5)) +'\t'
        to_print += '\n'            
        
          
    out_file.writelines(to_print)
    out_file.writelines('\n')
    out_file.writelines('\n')
               
    out_file.writelines('calibration filename ')
    out_file.writelines('\n')
    out_file.writelines(str(calfile))
    out_file.writelines('\n')
    out_file.close() #rcw added
    
    print('___'*30)
    print('Done!')
    print('___'*30)
    
    return nrcs_db