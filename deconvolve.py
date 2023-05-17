#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May  6 11:27:00 2021

@author: Rosie Willatt coding Thomas Newman's deconvolution and coherent noise subtraction processing, translating from Matlab
"""

# for reference: raw structure from scatlib.py:
       #  raw[i,0,:] = counts_to_voltage(vv)
       #  raw[i,1,:] = counts_to_voltage(hv)
       #  raw[i,2,:] = counts_to_voltage(vh)
       #  raw[i,3,:] = counts_to_voltage(hh)
       #  raw[i,4,:] = counts_to_voltage(cal)
       #  raw[i,5,:] = counts_to_voltage(noise)



def coh_noise_sub(raw_layer):
    
    waveletDecompositionLevel = 4 #Set wavelet decompostion level
    
    thisDataTranpose = raw_layer.T #transpose array
    
    thisDataTranposeHill = hilbert(thisDataTranpose) #take Hilbert transform
    
    thisDataTranposeHillAbs = np.abs(thisDataTranposeHill) #take magnitudes
    
    thisDataTranposeHillAbsMean = np.mean(thisDataTranposeHillAbs, axis = 1) #double check correct axis in python!
    
    
    
    ######HERE TO ADD FURTHER CODE!
    
    
    
    return coh_subd



        
def deconvolve(raw, band):  
    
    import pandas as pd
    import numpy as np
    from scipy.signal import hilbert
    import matplotlib.pyplot as plt
    
    print('raw.shape',raw.shape)
    
    decon = np.full(raw.shape, np.nan) #this will contain the deconvolved data
    coh_subd = np.full(raw.shape, np.nan)  #this will contain the raw data after coherent noise subtraction
    
    for j in range(4): #loop over 4 polarisations, don't do cal or noise
    
        # find correct deconvolution waveform
        file_path = '/Users/rosie/Documents/mosaic/mosaic_data/New_deconvolution_waveforms/'
        if j ==  0:
            decon_wf_df = pd.read_csv(file_path+'kukaDeconvolutionWaveform_'+band.lower()+'_vv.csv', header = None, delimiter = ',')
        elif j == 1:
            decon_wf_df = pd.read_csv(file_path+'kukaDeconvolutionWaveform_'+band.lower()+'_hv.csv', header = None, delimiter = ',')
        elif j == 2:
            decon_wf_df = pd.read_csv(file_path+'kukaDeconvolutionWaveform_'+band.lower()+'_vh.csv', header = None, delimiter = ',')
        elif j == 3:
            decon_wf_df = pd.read_csv(file_path+'kukaDeconvolutionWaveform_'+band.lower()+'_hh.csv', header = None, delimiter = ',')
        
        # remove coherent noise
        coh_subd[:,j,:] = coh_noise_sub(raw[:,j,:]) #process correct layer of data to remove coherent noise


        # convert to complex numpy array
        decon_wf = pd.DataFrame(decon_wf_df).to_numpy() #first convert dataframe to numpy
        dwf_comp = decon_wf[:,0] + 1j * decon_wf[:,1] #convert deconvolution waveform to complex

        # print('decon.shape',decon.shape)
        
        for i in range(raw.shape[0]): #loop through the samples/echoes

            #next two lines removed as were pre-coherent noise subtraction, left for reference only
            # raw_comp = hilbert(raw[i,j,:]) #convert to complex
            # decon[i,j,:] = raw_comp * dwf_comp #multiply elementwise to get the decovolved raw data

           decon[i,j,:] = coh_subd[i,:] * dwf_comp #multiply elementwise to get the decovolved raw data
            
    decon[:,4,:] = raw[:,4,:]
    decon[:,5,:] = raw[:,5,:]
    
    for j in range (6):
        plt.imshow(decon[:,j,:])
        plt.title('decon plot:'+str(j))
        plt.show()
        
        plt.imshow(raw[:,j,:])
        plt.title('raw plot:'+str(j))
        plt.show()

    return decon #Tom here this is returning the complex array but we probably need to explicitly convert to real as discussed!
            


