extract_config = {
    "instrument_name" :"Ku-Scat", 
    "i_sub_band" : 1, #set to 1 to process sub-band
    "sub_bandwidth":320, #bandwidth of sub-band in MHz
    "sub_center_frequency":13575, #center frequency of sub-band in MHz
    "i_calibrate":0, #0 to read existing cal file;1: find tx/rx distortion matrix using isotropic medium
    "calfile":"calib_file/Ku-Scat-20190731-145910.cal", #filename of calibration file. write "select" if you want to pick a file manually
    "i_corner_cal_override" : 0,
    "i_cal_loop_override" : 0,
    "i_batch":1, #0 for single file; 1 for batch processing
    "i_raw_plot":0, #0 skip plotting; 1 plot one range profile for each elevation; 2 plot all range profiles; 3 to plot non-scanning data
    "max_display_range":20, #maximum range in meters to display
    "i_proc_ind": 0,    #set to 1 to only process independent samples due to motion in stare data
    "n_blocks_ind_per_group":40,      #number of ind. samps per group when i_proc_ind:1
    "i_az_override":0,     # set to 1 to process a limited azimuth segment
    "azmin_proc":-5, #minimum azimuth angle to process when i_az_override:1
    "azmax_proc":45,            #maximum azimuth angle to process when i_az_override:1
    "i_el_override":0,         #set to 1 to process a limited range of elevation angles 
    "elmin_proc": 0,        #minimum elevation angle to process when i_el_override:1
    "elmax_proc":50,        #maximum elevation angle to process when i_el_override:1
    "i_corner_process": 0, #set to 1 to process Mueller matrix, etc. for corner reflector data
    "i_pol_scat":0, #1 to show polarization scatter plots
    "i_pol_signature": 0, #1 to show polarization signatures
    "i_pol_data_in_footprint": 0,    #1 to show polarimetric parameters vs range bin within footprint
    "smoothfac": 25,  #boxcar average length to smooth polarimetric parameters vs range bin 
    "i_temp_time_plot": 0,    #1 to plot RF unit temp vs time and data block time
    "proc_thresh_left_dB" :-20.0,    # i#nclude all data from this threshold left of peak to peak 
    "proc_thresh_right_dB" :-40.0,     #include all data from right of peak to this threshold right of peak 
    "group_averaging_time":10.0,         # s:averaging time in stare mode
    "raw_data_path":"/Volumes/bigLacie/all_stare/", #data file path for location of raw data
    "processed_data_path": "/Volumes/bigLacie/processed_stare_cs2_alt_new_nrcs/",     #data file path for location of processed data   
    "show_all_plot" : 0, 
}
