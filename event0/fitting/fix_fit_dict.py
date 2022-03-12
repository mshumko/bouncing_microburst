# -*- coding: utf-8 -*-
"""
Created on Tue Apr  4 13:55:45 2017

@author: mike
"""
import pickle
import numpy as np

# Fix the fit paramaters dictionary to remove the background Gaussians.
for sc_id in [3,4]:
    
    saveDir = '/home/mike/FIREBIRD/bouncing_packet_microburst/'
    saveName = 'FU' + str(sc_id) + '_fit_dictionary.pickle'
    with open(saveDir + saveName, 'rb') as handle:
        fitParams = pickle.load(handle)
        
    # Remove wide or incorrect Gaussians
    if sc_id is 3:
        for key in fitParams['CH0'].keys():
           fitParams['CH0'][key] = np.delete(
           fitParams['CH0'][key], (5, 6), axis = 0)
           fitParams['CH1'][key] = np.delete(fitParams['CH1'][key], 
            (len(fitParams['CH1'][key])-1, len(fitParams['CH1'][key])-2), axis = 0)
           fitParams['CH2'][key] = np.delete(
           fitParams['CH2'][key], (len(fitParams['CH2'][key])-1), axis = 0)
           
    elif sc_id is 4:
        for key in fitParams['CH0'].keys():   
            fitParams['CH0'][key] = np.delete(
            fitParams['CH0'][key], (len(fitParams['CH0'][key])-1), axis = 0)
            fitParams['CH1'][key] = np.delete(
            fitParams['CH1'][key], (len(fitParams['CH1'][key])-2), axis = 0)
            
    with open(saveDir + saveName, 'wb') as toFile:
        pickle.dump(fitParams, toFile, protocol=pickle.HIGHEST_PROTOCOL)