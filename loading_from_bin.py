from struct import unpack
import struct

import numpy as np


def get_file_header_t(file_path):
    with open(file_path, 'rb') as f:

        struct_fmt = '<6Q4dq3L2q2Lbdq2d2L256bl256bl2b128b2db5dqb7d2q2db3dL9dq3L'
        struct_len = struct.calcsize(struct_fmt)
        
        data = unpack(struct_fmt, f.read(struct_len)) 
        
        config_t = {
                    "ad_clock_frequency_hz": data[6], 
                    "chirp_amplitude_dbm": data[7], 
                    "chirp_bandwidth_hz": data[8], 
                    "chirp_center_frequency_hz": data[9], 
                    "chirp_width_ns": data[10],
                    "decimation": data[11],
                    "decimation_mode": data[12], 
                    "decimation_setting": data[13],
                    "file_roll_interval_ns": data[14],
                    "frame_delay_ns": data[15], 
                    "n_gates": data[16],
                    "n_gates_override": data[17], 
                    "n_gates_use_override": data[18],
                    "pedestal_height_m": data[19],
                    "prp_ns": data[20],
                    "range_resolution_m": data[21], 
                    "range_to_corner_reflector_m": data[22], 
                    "server_state": data[23], 
                    "switch_mode": data[24],
        }
        
        private_config_t = {
                            "pedestal_device_file_path" : np.array(list(data[25:153])),                            
                            "rcb_host": np.array(list(data[153:281])), 
                            "rcb_port": data[281], 
                            "rcb_scope_data_device_file_path": np.array(list(data[282:410])), 
                            "gps_device_file_path": np.array(list(data[410:538])), 
                            "scat_id": data[538],
                            "always_master": data[539], 
                            "peer_system_enabled": data[540], 
                            "peer_system_host": np.array(list(data[541:669])),
                            "antenna_beam_width": data[669], 
                            "corner_reflector_sigma": data[670], 
                            "heater_control_mode": data[671],
                            "heater_control_set_point": data[672], 
                            "heater_control_lcd_panel_set_point": data[673], 
                            "heater_control_hysteresis": data[674], 
                            "heater_control_proportional_gain": data[675],
                            "heater_control_integral_gain": data[676], 
                            "heater_control_minimum_pulse_width_ns": data[677],  
                            "heater_control_pulse_period_multiplier": data[678],  
                            "cooler_set_point": data[679],  
                            "cooler_hysteresis": data[680],  
                            "pedestal_max_speed": data[681],  
                            "pedestal_azimuth_offset": data[682],  
                            "pedestal_elevation_offset": data[683],  
                            "range_to_antenna_m": data[684],  
                            "range_to_antenna_offset_m": data[685],   
                            "ad_trig_delay_ns": data[686],   
                            "synthesizer_trig_delay_ns": data[687],   
                            "along_track_tilt_offset": data[688],  
                            "cross_track_tilt_offset": data[689],  
                            "calculate_ground_indices": data[690],  
                            "peak_detector_copol_threshold_db": data[691],  
                            "peak_detector_crosspol_threshold_db": data[692],  
                            "peak_detector_range_window_width_m": data[693],
            
        }
        
        calibration_t = {
                        "timestamp_seconds": data[694],
                        "corner_reflector_vv_power_dbm": data[695],
                        "corner_reflector_hh_power_dbm": data[696],
                        "cal_peak_dbm": data[697],
                        "corner_reflector_range_m": data[698],
                        "corner_reflector_az": data[699],
                        "corner_reflector_el": data[700],
                        "chirp_amplitude_dbm": data[701],
                        "chirp_bandwidth_hz": data[702],
                        "chirp_center_frequency_hz": data[703],
                        "chirp_width_ns": data[704],
                        "decimation": data[705],
                        "decimation_mode": data[706],
                        "n_summed_gates": data[707],   
        }
        
        file_header_t = {   "file_header_size": data[0], 
                            "meta_header_size": data[1], 
                            "data_size": data[2], 
                            "version_major": data[3], 
                            "version_minor": data[4], 
                            "version_patch": data[5], 
                            "config" : config_t,
                            "private_config": private_config_t,  
                            "calibration": calibration_t ,
                            }
                  
    return file_header_t, struct_len


def get_meta_header_t(file_path, sl):
    
    with open(file_path, 'rb') as f:
        f.seek(sl)
        struct_fmt = '<Q2Lbb9fL2b2l6f8blBlb3l3d'
        struct_len = struct.calcsize(struct_fmt)
        data = unpack(struct_fmt, f.read(struct_len))  

        
        
        rcb_status_t = {
                        "plo_20_ghz_lock": data[4],
                        "five_volts": data[5],
                        "minus_five_volts": data[6], 
                        "twelve_volts": data[7],
                        "input_twelve_volts": data[8],
                        "power_supply_plate_temp": data[9], 
                        "lcd_display_temp": data[10],
                        "rf_plate_temp": data[11],
                        "cross_track_tilt": data[12],
                        "along_track_tilt": data[13],
                        "error_bits": data[14],
                        "status_bits": data[15],
                        }
        ped_status_t = {
                        "az_mode": data[17], 
                        "el_mode": data[18], 
                        "az_pos": data[19], 
                        "el_pos": data[20], 
                        "az_vel": data[21], 
                        "el_vel": data[22], 
                        "az_current": data[23], 
                        "el_current": data[24], 
                        "az_at_ccw_hardware_limit":data[25],  
                        "az_at_cw_hardware_limit": data[26],  
                        "az_at_ccw_software_limit":data[27],  
                        "az_at_cw_software_limit": data[28],  
                        "el_at_ccw_hardware_limit":data[29],  
                        "el_at_cw_hardware_limit": data[30],  
                        "el_at_ccw_software_limit":data[31],  
                        "el_at_cw_software_limit": data[32],  
                        "sweep_count": data[33], 
                        "transition_flag": data[34], 
                        }
        
        project_status_t = {
                        "rcb_status_valid": data[3], 
                        "rcb_status": rcb_status_t, 
                        "ped_status_valid": data[16], 
                        "ped_status": ped_status_t, 
                        "scan_index": data[35], 
                        "gps_status_valid": data[36], 
                        "gps_solution_status": data[37], 
                        "gps_week": data[38], 
                        "gps_milliseconds": data[39],  
                        "gps_latitude": data[40], 
                        "gps_longitude": data[41], 
                        "gps_height": data[42], 
                        }   
        
        meta_header_t = { 
                        "index": data[0], 
                        "timestamp_seconds": data[1],  
                        "timestamp_nanoseconds": data[2],  
                        "project_status": project_status_t 
                        }
        
        new_sl = sl + struct_len   
    return meta_header_t, new_sl
    

def data_gathering(file_path, sl, szz):
    with open(file_path, 'rb') as f:
        f.seek(sl)
        
        struct_fmt = '<' + str(szz*6) + 'h'
        struct_len = struct.calcsize(struct_fmt)
        data = unpack(struct_fmt, f.read(struct_len))  
        
        gathered_data ={
            "hh": np.array(list(data[0:szz])),
            "vv": np.array(list(data[szz:(szz*2)])),
            "hv": np.array(list(data[(szz*2):(szz*3)])),
            "vh": np.array(list(data[(szz*3):(szz*4)])),
            "cal": np.array(list(data[(szz*4):(szz*5)])),
            "noise": np.array(list(data[(szz*5):(szz*6)])),
        }
        new_sl = sl + struct_len    
        return gathered_data, new_sl

