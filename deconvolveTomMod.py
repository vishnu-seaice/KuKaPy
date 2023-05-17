#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May  6 11:27:00 2021

@author: Tom modifing Rosie Willatt coding Thomas Newman's deconvolution and
coherent noise subtraction processing, translating from Matlab
"""

# #############################################################################
# LOAD PACKAGES
# #############################################################################
import pandas as pd
import numpy as np
import pywt
from scipy.signal import hilbert
import matplotlib.pyplot as plt 

# #############################################################################
# DEFINE FUNCTION: COHERENT NOISE REMOVAL
# #############################################################################


def runDeconvolution(raw, kukaBand):  # raw array and frequency band (Ka or Ku)

    pols = ['vv', 'hv', 'vh', 'hh']
    
    decon = np.full(raw.shape, np.nan) #this will contain the deconvolved data

    for j in range(4):  # loop over 4 polarisations, don't do cal or noise

        thisTestData = raw[:, j, :]

        # Remove coherent noise from the data
        thisTestDataCoh = removeCoherentNoise(thisTestData)
        np.savetxt('/Users/rosie/Documents/mosaic/mosaic_data/kukapy_output/'+'thisTestDataCoh_'+str(j)+".csv", thisTestDataCoh, delimiter=" ")

        # load decon waveforms for current band (KA or KU)
        thisBandDeconWaveforms = loadDeconWaveforms(kukaBand)

        # Apply deconvolution waveforms to data
        kukaPol = pols[j]
        thisTestDataCohDecon = applyDeconWaveforms(
            thisTestDataCoh, thisBandDeconWaveforms, kukaPol)

        decon[:, j, :] = np.transpose(thisTestDataCohDecon)

    decon[:, 4, :] = raw[:, 4, :]
    decon[:, 5, :] = raw[:, 5, :]
    
    

    for j in range (6):
        plt.imshow(decon[:,j,:])
        plt.title('decon plot:'+str(j))
        plt.show()
        
        plt.imshow(raw[:,j,:])
        plt.title('raw plot:'+str(j))
        plt.show()

    return decon 

        
# #############################################################################
# DEFINE FUNCTION: COHERENT NOISE REMOVAL
# #############################################################################
def removeCoherentNoise(inputData):
    
    #print('inputData.shape',inputData.shape)
    
    inputData = np.transpose(inputData)

    # Set the wavelet decompostion level
    waveletDecompositionLevel = 3  # (!! THIS SHOULD BE 3 FOR KUKA DATA !!)

    # The number of input data columns
    nCols = np.shape(inputData)[1]
    # print('nCols',nCols)
    # print('inputData.shape',inputData.shape)
    
    # .........................................................................
    # Preallocation first create empty
    cohNoiseArray = np.empty_like(inputData)
    
    # Then convert empty to nan
    cohNoiseArray[:] = np.nan
    # .........................................................................
    
    # =========================================================================
    # Remove coherent noise from each column
    # =========================================================================
    for iCol in range(nCols):
        
        # The current column
        thisCol = inputData[:,iCol]
        
        # ---------------------------------------------------------------------
        # Perform the multilevel wavelet decomposition using wavedec
        # ---------------------------------------------------------------------
        # Perform the multilevel wavelet decomposition using wavedec
        coeffs = pywt.wavedec(thisCol,'sym4', mode='sym', level=waveletDecompositionLevel)

        # Zero out all the detail coefficients from wavedec
        coeffs[1] = np.zeros_like(coeffs[1])
        coeffs[2] = np.zeros_like(coeffs[2])
        coeffs[3] = np.zeros_like(coeffs[3])

        # Perform multilevel reconstruction using waverec to obtain approximations
        a4Dwt = pywt.waverec(coeffs, 'sym4')
        
        # Assign to array
        # print('a4Dwt.shape',a4Dwt.shape)
        # print(iCol)
        cohNoiseArray[:,iCol] = a4Dwt
        # ---------------------------------------------------------------------
        
    # =========================================================================

    # Remove the coherent noise array from original data
    cohCorrectedArray = inputData - cohNoiseArray

    # Specify the output data
    outputData = cohCorrectedArray

    return outputData


# #############################################################################
# DEFINE FUNCTION: LOAD DECON WAVEFORMS
# #############################################################################
def loadDeconWaveforms(kukaBand):

    # Define the deconvolution file path
    deconFilePath = '/Users/rosie/Documents/mosaic/mosaic_data/waveforms/' #rcw

    # Define the pol vector
    polVect = ["vv", "hv", "vh", "hh"]

    # Setup dictionary for deconvolution waveforms
    deconWaveforms = {}

    # Loop over all polarisations
    for thisPol in polVect:

        # Define the current pol file
        thisDeconPolFilename = 'kukaDeconvolutionWaveform_' + kukaBand.lower() + '_' + \
            thisPol + '.csv'

        # Get pol data
        thisDeconPolData = pd.read_csv(
            deconFilePath + thisDeconPolFilename, header=None, delimiter=',')

        # Convert pol data to complex
        thisDeconPolDataComplex = thisDeconPolData[0] + 1j*thisDeconPolData[1]

        # Add to deconWaveforms dictionary
        deconWaveforms[thisPol] = thisDeconPolDataComplex

    # Specify output data
    outputData = deconWaveforms

    return outputData


# #############################################################################
# DEFINE FUNCTION: APPLY DECON WAVEFORMS
# #############################################################################
def applyDeconWaveforms(inputData, deconWaveforms, kukaPol):

    # Find the number of columns (traces) in input array
    nColsInputData = np.shape(inputData)[1]

    # Extract the appropriate deconvolution waveform
    thisDeconWaveform = deconWaveforms[kukaPol]

    # Turn deconvolution waveform in an array by tiling it
    thisDeconWaveformArray = np.transpose(
        np.tile(thisDeconWaveform, (nColsInputData, 1)))

    # Take the Hilbert transform of the input array
    thisDataHill = hilbert(inputData, axis=0)

    # Multiply hilbert transformed input array by deconvoltion waveform array
    thisDataHillDecon = thisDataHill * thisDeconWaveformArray

    # Convert back to real
    thisDataDecon = np.real(thisDataHillDecon)

    # Specify output data
    outputData = thisDataDecon

    return outputData

# #############################################################################


# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# BEGIN CODE
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# =============================================================================
# Load test data
# =============================================================================
# for reference: raw structure from scatlib.py:
#  raw[i,0,:] = counts_to_voltage(vv)
#  raw[i,1,:] = counts_to_voltage(hv)
#  raw[i,2,:] = counts_to_voltage(vh)
#  raw[i,3,:] = counts_to_voltage(hh)
#  raw[i,4,:] = counts_to_voltage(cal)
#  raw[i,5,:] = counts_to_voltage(noise)

# .............................................................................
# Test data
# .............................................................................
# KA-SCAT-20200116-120702-vv
# KA-SCAT-20200116-120702-hv
# KA-SCAT-20200116-120702-vh
# KA-SCAT-20200116-120702-hh
# ............................................................................
# KU-SCAT-20200116-120700-vv
# KU-SCAT-20200116-120700-hv
# KU-SCAT-20200116-120700-vh
# KU-SCAT-20200116-120700-hh
# .............................................................................

# # Define test data path
# testDataPath = '/Volumes/uclDataDiskB/Mosaic/KuKaData/TEST/STARE/20200116/'

# # Specify test data file
# testDataFile = 'KA-SCAT-20200116-120702-vh'

# # Extract kuka band
# kukaBand = testDataFile[0:2]

# # Extract kuka pol
# kukaPol = testDataFile[-2:]

# # Define test data file path
# testDataFilePath = testDataPath + testDataFile + '.csv'

# # Load test data
# thisTestDataPd = pd.read_csv(testDataFilePath, header=None, delimiter=',')

# # first convert dataframe to numpy
# thisTestData = pd.DataFrame(thisTestDataPd).to_numpy()


# # =============================================================================
# # Process data
# # =============================================================================
# # Remove coherent noise from the data
# thisTestDataCoh = removeCoherentNoise(thisTestData)

# # load decon waveforms for current band (KA or KU)
# thisBandDeconWaveforms = loadDeconWaveforms(kukaBand)

# # Apply deconvolution waveforms to data
# thisTestDataCohDecon = applyDeconWaveforms(
#     thisTestDataCoh, thisBandDeconWaveforms, kukaPol)

# # Save to compare with matlab
# np.savetxt(testDataFile+'-decon.csv', thisTestDataCohDecon, delimiter=",")

# =============================================================================


# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# END OF CODE
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
