import os
from glob import glob

from scatlib import read_configuration, read_header, read_raw, compute_range_profiles, calibrate, \
    polarimetric_processing_scan, nrcs_compute_scan, polarimetric_processing_stare, nrcs_compute_stare , \
        polarimetric_processing_stare_by_independent_samples


def main():
    iscat = input("enter 1 for Ku-Scat; 2 for Ka-Scat: ")
    
    # Reading configurations from configvars in SCATlib.pro
    configvars = read_configuration(iscat) 

    instrument_name = configvars["instrument_name"]
    calfile = configvars["calfile"]    #pass this as plain variable since it may get modified if creating a new calfile
    raw_data_path = configvars["raw_data_path"]
    processed_data_path = configvars["processed_data_path"]

    if iscat == "1":
        extension_for_files = '/KU*.dat'
    else:
        extension_for_files = '/KA*.dat'
    
    
    if configvars["i_batch"] == 0:
        file_names= list()
        for ind, name in enumerate(glob(raw_data_path + extension_for_files)): 
            file_names.append(name)
            print(ind, " : ", name.split('/')[-1])
        filename = file_names[int(input("choose SCAT raw data file by giving the index: "))]
        files = [filename]
        filecount = len(files)
    else:
        files, filecount = list(), 0
        for name in glob(raw_data_path + extension_for_files): 
            files.append(name)
            filecount += 1
        files.sort() #rcw added to put in order
    
    for ifile in range(0,filecount):
        filename = files[ifile]
        print(filename)
        print("=" * 70)
        
            
        base_filename = instrument_name + filename.split('/')[-1][-19:-4] #rcw
        
        ## READ_HEADER
        print("READ_HEADER")
        scatvars ,calvars ,private_config ,file_header_size, sl = read_header(filename, raw_data_path, configvars)
        
        #rcw pickle scatvars ,calvars ,private_config
        #pickle.dump(scatvars, open(processed_data_path+base_filename+'_scatvars.p', 'wb')) 
        #pickle.dump(calvars, open(processed_data_path+base_filename+'_calvars.p', 'wb')) 
        #pickle.dump(configvars, open(processed_data_path+base_filename+'_configvars.p', 'wb')) 
        #rcw pickling done

        # rcw get the reference_calibration_loop_power_file         
        reference_calibration_loop_power_file = 10**(calvars["cal_peak_dbm"]/10.)  

        ## READ_RAW
        print("READ_RAW")
    
        raw, scan_index, sweep_count, transition_flag, elevation, \
          n_blocks, n_pol, elapsed_time,time_sec,gps_latitude,gps_longitude,\
            along_track_tilt,cross_track_tilt, independent_sample_index,distance,\
                az_proc_index, sweep_count_override =  read_raw(configvars, scatvars,\
                    filename, file_header_size, sl)#, base_filename, decon, processed_data_path) #rcw added base_filename, decon, processed_data_path
        
        
        print('RRR',n_blocks, len(gps_latitude), len(elevation), len(time_sec))
        #rcw adding code to save out the raw arrays and deconvolve if wanted
        import numpy as np
        
        #for decon in [0]: #rcw switch off use of deconvolved data for sub-banding because we don't have the deconvolution waveforms
        for decon in [0,1]: #rcw switch on/off use of deconvolved data

            # if decon == 0: 
            #     for raw_layer_number in range(6):
                    #save and plot raw data
                    # name4sav=base_filename[0:7]+'-'+base_filename[7:]+'_layer'+str(raw_layer_number)+'_raw_decon0_t0py'
                    # np.savetxt(processed_data_path+name4sav+".csv", raw[:,raw_layer_number,:], delimiter=" ")
                    # plt.imshow(raw[:,raw_layer_number,:])
                    # plt.title('scatproc plot of raw'+str(raw_layer_number))
                    # plt.show()    
                
            if decon == 1:
                #deconvolve raw data
                import deconvolveTomModRosie as deconvolve
                raw = deconvolve.runDeconvolution(raw, base_filename[0:2])        
                #save and plot deconvolved raw data
                # for raw_layer_number in range(6):
                #     name4sav=base_filename[0:7]+'-'+base_filename[7:]+'_layer'+str(raw_layer_number)+'_raw_decon1_t0py'
                #     np.savetxt(processed_data_path+name4sav+".csv", raw[:,raw_layer_number,:], delimiter=" ")       
                    # plt.imshow(raw[:,raw_layer_number,:])
                    # plt.title('scatproc plot of raw_decon'+str(raw_layer_number))
                    # plt.show()          
            #end rcw save/plot raw data and deconvolve
    
            #rcw        base_filename = instrument_name + filename.split('/')[-1][-19:-4]
            processed_data_filename = processed_data_path+'/'+base_filename+'_decon'+str(decon)+'.nrcs'
            print(f"processed_data_filename: {processed_data_filename}")

    
            print("Compute Range Profiles")
            ## compute range profiles of reflectivity from raw data
            range_gate_spacing,ngates,gate_offset,gate_plot_max,pos_height,height, \
            _range,spec,gate_max_cal,geometric_peak = compute_range_profiles(n_blocks,\
                                                        n_pol,elevation,scan_index,\
                                                        sweep_count,raw,scatvars,private_config\
                                                        ,configvars)
            
            print("Calibrate")
            # ## the calibrate procedure is used to generate the transmit distortion matrix 
            reference_calibration_loop_power, current_calibration_loop_power,ainv,finv, corner_range, \
                total_corner_power_vv,total_corner_power_hh = calibrate(calvars,\
                    configvars,gate_offset,gate_plot_max,spec,ngates,n_blocks,n_pol,\
                        _range,gate_max_cal, range_gate_spacing,scatvars,scan_index,sweep_count,\
                            elevation,base_filename,calfile)
                    
            
            #pickle.dump(current_calibration_loop_power, open(processed_data_path+base_filename+'_current_calibration_loop_power.p', 'wb'))  #rcw
            
            #compute calibrated covariance matrix and Mueller matrix containing average polarimetric
            #scattering properties of target
            
            if max(scan_index) == 1:
                print("Polarimetric Scan")
                range_peak_signal,line_elevation,line_height,nlines,rho_hv_vec,\
                phase_hv_deg_vec,total_power,range_centroid_signal, L_matrix,c_matrix, echoes_all = polarimetric_processing_scan(scatvars,configvars,\
                            calvars,geometric_peak,scan_index,\
                            sweep_count,transition_flag,elevation,ngates,\
                            range_gate_spacing,_range,height,gate_offset,\
                            gate_plot_max,spec,ainv,finv, az_proc_index,\
                                sweep_count_override)
            
    
    
            elif max(scan_index) != 1 and configvars["i_proc_ind"] == 0:
                print("Polarimetric Scan Stare")
                
                range_peak_signal,line_elevation,n_groups,rho_hv_vec,\
                phase_hv_deg_vec,total_power,range_centroid_signal, L_matrix,c_matrix, \
                    group_index, n_blocks_per_group, echoes_all \
                    = polarimetric_processing_stare(scatvars,configvars,calvars,ngates,range_gate_spacing,\
                        _range,gate_offset,gate_plot_max,spec,ainv,finv,n_blocks, elapsed_time,\
                        scan_index)#, base_filename, decon, processed_data_path,\
                            #current_calibration_loop_power,reference_calibration_loop_power_file)   #rcw added base_filename, decon
    
    
            elif max(scan_index) != 1 and configvars["i_proc_ind"] == 1:
                print("Polarimetric Scan Stare by Independant Samples")
    
                range_peak_signal,line_elevation,n_groups,rho_hv_vec,\
                phase_hv_deg_vec,total_power,range_centroid_signal, L_matrix,c_matrix, \
                group_index, n_blocks_per_group, echoes_all = polarimetric_processing_stare_by_independent_samples(scatvars,\
                configvars,calvars,ngates,range_gate_spacing,_range,gate_offset,gate_plot_max,spec,ainv,finv,\
                n_blocks, elapsed_time, scan_index, independent_sample_index, distance)            
            
            
            if max(scan_index) == 1:
                print("nrcs_compute_scan")
                nrcs_db = nrcs_compute_scan(scatvars,calvars,configvars,private_config,c_matrix,range_peak_signal,\
                        range_centroid_signal,corner_range,pos_height,line_elevation,nlines,reference_calibration_loop_power,\
                        current_calibration_loop_power,L_matrix,processed_data_filename,rho_hv_vec,\
                        phase_hv_deg_vec,total_power,line_height,total_corner_power_vv,total_corner_power_hh,gps_latitude,gps_longitude,calfile)            
    
            elif max(scan_index) != 1:
                print("nrcs_compute_stare")
                nrcs_db = nrcs_compute_stare(scatvars,calvars,configvars,private_config,c_matrix,range_peak_signal,\
                    range_centroid_signal,corner_range,pos_height,n_groups,reference_calibration_loop_power,current_calibration_loop_power,\
                    L_matrix,processed_data_filename,rho_hv_vec,phase_hv_deg_vec,total_power,total_corner_power_vv,total_corner_power_hh,\
                    calfile,time_sec,gps_latitude,gps_longitude,group_index,along_track_tilt,cross_track_tilt,n_blocks_per_group,independent_sample_index)
                
            # print('echoes_all')
            # print(echoes_all['vv'])
            # rcw adding netcdf file storage
            if decon == 0:

                from netCDF4 import Dataset
                
                #create file
                if max(scan_index) == 1:
                    nfilename = processed_data_path + 'kuka_scan_decon_'+base_filename+'.nc'
                else:
                    nfilename = processed_data_path + 'kuka_stare_decon_'+base_filename+'.nc'
                print('creating netcdf file: ', nfilename)
                ncfile = Dataset(nfilename, mode='w', format='NETCDF4')
                
                #general info
                ncfile.title='KuKa combined echo and summary data'
                ncfile.use_permission = 'Permission of project PI (Prof. Julienne Stroeve j.stroeve@ucl.ac.uk) is required to access, analyse and publish data before 1st January 2023. See MOSAiC Data Policy https://mosaic-expedition.org/wp-content/uploads/2020/12/mosaic_datapolicy.pdf.'
                ncfile.assistance = 'For data processing assistance, contact Dr. Rosemary Willatt (r.willatt@ucl.ac.uk) (Stare Mode), Dr. Vishnu Nandan (vishnu.nandan@umanitoba.ca) (Scan Mode), Dr Thomas Newman (t.newman@ucl.ac.uk) (Deconvolution)'
                ncfile.file_processed_with = 'KuKaPy translated from ProSensing IDL code by Vishnu Nandan with additions by Rosemary Willatt and Thomas Newman'
                ncfile.file_created_by = 'Willatt UCL' 
                import datetime as dt
                ncfile.date_created = dt.date.today().strftime('%Y%m%d')
                ncfile.data_type = 'KuKaPy output'
                ncfile.kuka_operators = 'Stefan Hendricks, Gunnar Spreen and Oguz Demir (leg 1), Julienne Stroeve, Vishnu Nandan, Rasmus Tonboe and Marcus Huntemann (leg 2) \
                    Aikaterini Tavri and Mallik Mahmud (leg 4), Gunnar Spreen (leg 5)'
                    

                ncfile.current_calibration_loop_power = current_calibration_loop_power
                
                #groups to contain the vars for setup
                sv = ncfile.createGroup('scatvars')
                for k,v in scatvars.items():
                    setattr(sv, k, v)
                cv = ncfile.createGroup('calvars')
                for k,v in calvars.items():
                    setattr(cv, k, v)
                co = ncfile.createGroup('configvars')
                for k,v in configvars.items():
                    setattr(co, k, v)
                
                #create dimensions
                if max(scan_index) ==1: #elevation angles for scan data
                    sample_dim = ncfile.createDimension('sample', nlines)
                elif max(scan_index) !=1: #n_groups for stare data
                    sample_dim = ncfile.createDimension('sample', n_groups) 
                
                range_dim = ncfile.createDimension('range', len(_range)) 
                pol_dim = ncfile.createDimension('pol', 4)
                
                #fill values
                polarisation = ncfile.createVariable('polarisation', np.str, ('pol',))
                polarisation[:] = np.array(['VV', 'HV', 'VH', 'HH'])
                                


                if max(scan_index) ==1: #elevation angles for scan data
                    ncfile.nlines = nlines
                    
                    nelevation_angles = ncfile.createVariable('elevation_angles', np.float64, ('sample'))
                    nelevation_angles.units = 'degrees'
                    # print('elevation_angles', elevation_angles)
                    nelevation_angles[:] = line_elevation
                    
                    ncfile.ave_lat = np.mean(gps_latitude)
                    ncfile.ave_lon = np.mean(gps_longitude)
                    
                elif max(scan_index) !=1: #tilts for stare data
                    ncfile.n_groups = n_groups
                    ncfile.n_blocks_per_group = n_blocks_per_group
                    
                    along_tilt = ncfile.createVariable('along_tilt', np.float64, ('sample',))
                    along_tilt.units = 'degrees'
                    along_tilt.long_name = 'along track tilt'
                    along_tilt[:] = along_track_tilt
                    
                    cross_tilt = ncfile.createVariable('cross_tilt', np.float64, ('sample',))
                    cross_tilt.units = 'degrees'
                    cross_tilt.long_name = 'cross track tilt'
                    cross_tilt[:] = cross_track_tilt  
                    
                    start_time = ncfile.createVariable('start_time', np.float64, ('sample',))
                    start_time.units = 'group start time in seconds since 1970-1-1'
                    start_time.long_name = 'start times'
                    start_time[:] = time_sec
                    
                    lat = ncfile.createVariable('lat', np.float64, ('sample',))
                    lat.units = 'degrees_north'
                    lat.long_name = 'latitudes'
                    lat[:] = gps_latitude
                    
                    lon = ncfile.createVariable('lon', np.float64, ('sample',))
                    lon.units = 'degrees_east'
                    lon.long_name = 'longitudes'
                    lon[:] = gps_longitude
                    
                nrange = ncfile.createVariable('range', np.float64, ('range',))
                nrange.units = 'metres'
                nrange.long_name = 'range'
                nrange[:] = _range
        
                vv_power_decon0 = ncfile.createVariable('vv_power_decon0',np.float64,('range','sample')) 
                vv_power_decon0.units = 'linear power non deconvolved' 
                vv_power_decon0.standard_name = 'VV power non deconvolved' 
                vv_power_decon0[:] = echoes_all['vv']
                
                hv_power_decon0 = ncfile.createVariable('hv_power_decon0',np.float64,('range','sample')) 
                hv_power_decon0.units = 'linear power non deconvolved' 
                hv_power_decon0.standard_name = 'HV power non deconvolved'    
                hv_power_decon0[:] = echoes_all['hv']

                vh_power_decon0 = ncfile.createVariable('vh_power_decon0',np.float64,('range','sample')) 
                vh_power_decon0.units = 'linear power non deconvolved' 
                vh_power_decon0.standard_name = 'VH power non deconvolved'
                vh_power_decon0[:] = echoes_all['vh']
                
                hh_power_decon0 = ncfile.createVariable('hh_power_decon0',np.float64,('range','sample')) 
                hh_power_decon0.units = 'linear power non deconvolved' 
                hh_power_decon0.standard_name = 'HH power non deconvolved' 
                hh_power_decon0[:] = echoes_all['hh']

                range_peak_signal_decon0 = ncfile.createVariable('range_peak_signal_decon0', np.float64, ('sample',))
                range_peak_signal_decon0.units = 'metres'
                range_peak_signal_decon0[:] = range_peak_signal

                #'Normalized Radar Cross Section (dBm2/m2) columns: VV, HV, VH, HH#'
                nrcs_db_decon0 = ncfile.createVariable('nrcs_decon0',np.float64,( 'pol','sample'))
                nrcs_db_decon0.standard_name = 'Normalized Radar Cross Section (dBm2/m2) columns: VV, HV, VH, HH non-deconvolved data'
                if max(scan_index) == 1:
                    nrcs_db_decon0[:] = nrcs_db.T
                elif max(scan_index) != 1:
                    nrcs_db_decon0[:] = nrcs_db
                
                #HH/VV correlation coefficient - magnitude
                rho_hv_vec_decon0 = ncfile.createVariable('rho_hv_vec_decon0',np.float64,('sample'))
                rho_hv_vec_decon0.standard_name = 'magnitude of HH/VV correlation coefficient non-deconvolved data'
                rho_hv_vec_decon0[:] = rho_hv_vec
                
                #HH/VV correlation coefficient - phase
                phase_hv_deg_vec_decon0 = ncfile.createVariable('phase_hv_deg_vec_decon0',np.float64,( 'sample'))
                phase_hv_deg_vec_decon0.standard_name = 'phase of HH/VV correlation coefficient (deg) non-deconvolved data'
                phase_hv_deg_vec_decon0[:] = phase_hv_deg_vec
                
                # print('rho_hv_vec', rho_hv_vec.shape, rho_hv_vec)
                # print('phase_hv_deg_vec', phase_hv_deg_vec.shape, phase_hv_deg_vec)
                
                print('echoes shape', echoes_all['vh'].shape)
                
                # np.savetxt(processed_data_path+base_filename+'echoes_all_vh_decon0.csv', echoes_all['vh'], delimiter=" ") 
                # np.savetxt(processed_data_path+base_filename+'echoes_all_hv_decon0.csv', echoes_all['hv'], delimiter=" ") 
                # np.savetxt(processed_data_path+base_filename+'echoes_all_vv_decon0.csv', echoes_all['vv'], delimiter=" ") 
                # np.savetxt(processed_data_path+base_filename+'echoes_all_hh_decon0.csv', echoes_all['hh'], delimiter=" ") 
            

            elif decon == 1:
                
                rcw_corr_cal = current_calibration_loop_power/reference_calibration_loop_power_file
                
                vv_power = ncfile.createVariable('vv_power',np.float64,('range','sample')) 
                vv_power.units = 'linear power deconvolved' 
                vv_power.standard_name = 'VV power deconvolved' 
                vv_power[:] = echoes_all['vv']/rcw_corr_cal
                
                hv_power = ncfile.createVariable('hv_power',np.float64,('range','sample')) 
                hv_power.units = 'linear power deconvolved' 
                hv_power.standard_name = 'HV power deconvolved'    
                hv_power[:] = echoes_all['hv']/rcw_corr_cal

                vh_power = ncfile.createVariable('vh_power',np.float64,('range','sample')) 
                vh_power.units = 'linear power deconvolved' 
                vh_power.standard_name = 'VH power deconvolved'
                vh_power[:] = echoes_all['vh']/rcw_corr_cal
                
                hh_power = ncfile.createVariable('hh_power',np.float64,('range','sample')) 
                hh_power.units = 'linear power deconvolved' 
                hh_power.standard_name = 'HH power deconvolved' 
                hh_power[:] = echoes_all['hh']/rcw_corr_cal

                range_peak_signal1 = ncfile.createVariable('range_peak_signal', np.float64, ('sample',))
                range_peak_signal1.units = 'metres'
                range_peak_signal1[:] = range_peak_signal
                
                #'Normalized Radar Cross Section (dBm2/m2) columns: VV, HV, VH, HH#'
                nrcs_db1 = ncfile.createVariable('nrcs_db',np.float64,( 'pol','sample'))
                nrcs_db1.standard_name = 'Normalized Radar Cross Section (dBm2/m2) columns: VV, HV, VH, HH DECONVOLVED data'
                if max(scan_index) == 1:
                    nrcs_db1[:] = nrcs_db.T
                elif max(scan_index) != 1:
                    nrcs_db1[:] = nrcs_db                
                #HH/VV correlation coefficient - magnitude
                rho_hv_vec1 = ncfile.createVariable('rho_hv_vec',np.float64,('sample'))
                rho_hv_vec1.standard_name = 'magnitude of HH/VV correlation coefficient DECONVOLVED data'
                rho_hv_vec1[:] = rho_hv_vec
                
                #HH/VV correlation coefficient - phase
                phase_hv_deg_vec1 = ncfile.createVariable('phase_hv_deg_vec',np.float64,('sample'))
                phase_hv_deg_vec1.standard_name = 'phase of HH/VV correlation coefficient DECONVOLVED data'
                phase_hv_deg_vec1[:] = phase_hv_deg_vec

                print(ncfile)
                ncfile.close(); print('netcdf file save complete')     
                
                # np.savetxt(processed_data_path+base_filename+'echoes_all_vh_decon1.csv', echoes_all['vh'], delimiter=" ") 
                # np.savetxt(processed_data_path+base_filename+'echoes_all_hv_decon1.csv', echoes_all['hv'], delimiter=" ") 
                # np.savetxt(processed_data_path+base_filename+'echoes_all_vv_decon1.csv', echoes_all['vv'], delimiter=" ") 
                # np.savetxt(processed_data_path+base_filename+'echoes_all_hh_decon1.csv', echoes_all['hh'], delimiter=" ") 
                
                #end rcw netcdf additions
        
if __name__ =="__main__":
    main()
    