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
rawdata = os.path.join('C:',os.sep,'Users','Rena','Downloads','Uni','Spatial Transcriptomics',"amygdalaSpatialInvestigationElliKlein","rawdata")
derivatives = os.path.join('C:',os.sep,'Users','Rena','Downloads','Uni','Spatial Transcriptomics',"amygdalaSpatialInvestigationElliKlein","derivatives")

#%% function to load data for regional samples

def load3dRegionalSample(save_path):
    """
    Load 3D regional sample created by using the saveRegionalCellsAsSample
    function. Sample includes: sampleID, save_path, age, sex, original_filename,
    geneList, cellTypeClassNames, ccfRegionsIncluded, ccfCoordinates,
    cellTypeClassInt, cellTypeClassColorHex, cell_label, ccfID, and geneMatrix

    Parameters
    ----------
    save_path : path
        Location of sample to be loaded to workspace.

    Returns
    -------
    dict
        Dictionary containing sample information used for later analyses.
    """
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

foldernames = ['mouse_609882_registered_082725_amygdalar', 'mouse_609889_registered_082725_amygdalar',
             'mouse_638850_registered_082725_amygdalar', 'mouse_687997_registered_082725_amygdalar']

allSamples = {}
for i, foldername in enumerate(foldernames):
    print(i, foldername)
    allSamples[foldername] = load3dRegionalSample(os.path.join(rawdata, foldername))

#%% count cells of one specific cell type in one sample
cellTypeIdx = allSamples['mouse_609882_registered_082725_amygdalar']["cellTypeClassNames"].index("01 IT-ET Glut")
nOfCellType = len(np.where(allSamples['mouse_609882_registered_082725_amygdalar']['cellTypeClassInt'] == cellTypeIdx)[0])
#%% write a function to find all cells of a cell Type
cellTypeIdxList = []
for cellTypeName in allSamples['mouse_609882_registered_082725_amygdalar']["cellTypeClassNames"]:
    cellTypeSplit = cellTypeName.split(" ")
    if "Glut" in cellTypeSplit:
        print(cellTypeName)
        cellTypeIdx = allSamples['mouse_609882_registered_082725_amygdalar']["cellTypeClassNames"].index(cellTypeName)
        cellTypeIdxList.append(cellTypeIdx)
# convert list into np.array
cellTypeIdxList = np.array(cellTypeIdxList, dtype = "int32")
       
nOfCellType = 0
for cellTypeIdx in cellTypeIdxList:
    nOfCellType += len(np.where(allSamples['mouse_609882_registered_082725_amygdalar']['cellTypeClassInt'] == cellTypeIdx)[0]) 
    print(nOfCellType)



#%%

def countAllCellsOfCellType(sample, cellType):
    cellTypeIdxList = []
    for cellTypeName in sample["cellTypeClassNames"]:
        cellTypeSplit = cellTypeName.split(" ")
        if cellType in cellTypeSplit:
            #print(cellTypeName)
            cellTypeIdx = sample["cellTypeClassNames"].index(cellTypeName)
            cellTypeIdxList.append(cellTypeIdx)
        elif cellType == cellTypeName:
            cellTypeIdx = sample["cellTypeClassNames"].index(cellTypeName)
            cellTypeIdxList.append(cellTypeIdx)
            
    # convert list into np.array
    cellTypeIdxList = np.array(cellTypeIdxList, dtype = "int32")
           
    nOfCellType = 0
    for cellTypeIdx in cellTypeIdxList:
        nOfCellType += len(np.where(sample['cellTypeClassInt'] == cellTypeIdx)[0]) 
    return nOfCellType
x = countAllCellsOfCellType(allSamples['mouse_609882_registered_082725_amygdalar'], "GABA")

# recreate for Loop for each sample to count Glut cells, rename x and append it to a list 
#list
for sampleID in foldernames:
    allSamples[sampleID]
        
    












