#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 11 12:03:36 2026

@author: zjpeters
"""
import os
import json
from glob import glob
import pandas as pd
import numpy as np
import scipy.sparse as sp_sparse
currentLocation = os.path.dirname(os.path.abspath(__file__))
rawdata = os.path.join('/','media','zjpeters','Expansion','eklein','rawdata')
derivatives = os.path.join('/','media','zjpeters','Expansion','eklein','derivatives')
#%% function to load data for regional samples

def load3dRegionalSample(save_path):
    sample = {}
    jsonLoc = glob(os.path.join(save_path, "*_processing_information.json"))[0]
    with open(jsonLoc) as jsonData:
        jsonDict = json.loads(jsonData.read())
        sample['sampleID'] = jsonDict['sampleID']
        sample['save_path'] = save_path
        sample['age'] = jsonDict['age']
        sample['sex'] = jsonDict['sex']
        sample['original_filename'] = jsonDict['original_filename']
        sample['geneList'] = jsonDict['geneList']
        sample['cellTypeClassNames'] = jsonDict['cellTypeClassNames']
        sample['ccfRegionsIncluded'] = jsonDict['ccfRegionsIncluded']
    ccfCoorCsvLoc = glob(os.path.join(save_path, "*ccf_coordinates.csv"))[0]
    ccfCsvData = pd.read_csv(ccfCoorCsvLoc)
    sample['ccfCoordinates'] = np.array([ccfCsvData['x_CCF'], ccfCsvData['y_CCF'], ccfCsvData['z_CCF']]).T
    sample['cellTypeClassInt'] = np.array(ccfCsvData['cellTypeClassInt'], dtype='int32')
    sample['cellTypeClassColorHex'] = np.array(ccfCsvData['cellTypeClassColorHex'], dtype=np.dtype('U17'))
    sample['cell_label'] = np.array(ccfCsvData['cell_label'], dtype=np.dtype('U7'))
    sample['ccfID'] = np.array(ccfCsvData['ccfID'], dtype='int64')
    sample['z_coord'] = np.array(ccfCsvData['z_coord'], dtype='int64')
    geneMatrixLoc = glob(os.path.join(save_path, "*_gene_matrix.npz"))[0]
    sample['geneMatrix'] = sp_sparse.load_npz(geneMatrixLoc)
    return sample

#%% load samples

filenames = ['mouse_609882_registered_082725_amygdalar', 'mouse_609889_registered_082725_amygdalar',
             'mouse_638850_registered_082725_amygdalar', 'mouse_687997_registered_082725_amygdalar']

allSamples = {}
for i, filename in enumerate(filenames):
    allSamples[i] = load3dRegionalSample(os.path.join(rawdata,filename))
