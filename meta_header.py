rc_status_t = {   
'plo_20_ghz_lock' : 0   ,   
'five_volts' : 0.0,   
'minus_five_volts' : 0.0,   
'twelve_volts' : 0.0,   
'input_twelve_volts' : 0.0,   
'power_supply_plate_temp' : 0.0,   
'lcd_display_temp' : 0.0,   
'rf_plate_temp' : 0.0,   
'cross_track_tilt' : 0.0,   
'along_track_tilt' : 0.0,   
'error_its' : 0,   
'status_its' : 0      
}

ped_status_t = {   
'az_mode': 0 ,   
'el_mode': 0 ,   
'az_pos': 0.0,   
'el_pos': 0.0,   
'az_vel': 0.0,   
'el_vel': 0.0,   
'az_current': 0.0,   
'el_current': 0.0,   
'az_at_ccw_hardware_limit': 0   ,   
'az_at_cw_hardware_limit': 0   ,   
'az_at_ccw_software_limit': 0   ,   
'az_at_cw_software_limit': 0   ,   
'el_at_ccw_hardware_limit': 0   ,   
'el_at_cw_hardware_limit': 0   ,   
'el_at_ccw_software_limit': 0   ,   
'el_at_cw_software_limit': 0   ,   
'sweep_count': 0 ,   
'transition_flag': 0      
}

project_status_t = {   
'rc_status_valid': 0   ,   
'rc_status': rc_status_t,   
'ped_status_valid': 0,   
'ped_status': ped_status_t,   
'scan_index': 0 ,   
'gps_status_valid': 0   ,   
'gps_solution_status': 0 ,   
'gps_week': 0 ,   
'gps_milliseconds': 0 ,   
'gps_latitude': 0.0  ,   
'gps_longitude': 0.0  ,   
'gps_height': 0.0     
}

meta_header_t = {   
'index': 0     ,   
'timestamp_seconds': 0    ,   
'timestamp_nanoseconds': 0    ,   
'project_status': project_status_t   
}