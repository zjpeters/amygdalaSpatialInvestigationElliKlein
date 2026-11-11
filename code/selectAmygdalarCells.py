#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 27 07:35:40 2026

@author: zjpeters
"""
import os
import pandas as pd
import numpy as np
import h5py
import scipy.sparse as sp_sparse
from scipy.sparse import csr_matrix
from matplotlib import pyplot as plt
import json
from glob import glob
# next line identifies the path for the directory where this script is located,
# since it will be different between computers, i.e. my computer vs yours
currentLocation = os.path.dirname(os.path.abspath(__file__))
# this line combines the above location with where directory where the 
# "allen_ccf_annotation.csv" file is located
ccfLocation = os.path.join(currentLocation,'data','ccf')
overlappingGenes = list(pd.read_csv(os.path.join(currentLocation, 'data','list_of_overlapping_genes.csv'))['gene_name'])
#%% code we worked on
# load the CCF information from the csv file
ccfInfo = pd.read_csv(os.path.join(ccfLocation, 'allen_ccf_annotation.csv'))
# convert names from pandas dataframe into a python list
ccfList = list(ccfInfo["name"])

# assigns the index where "Cortical amygdalar area" is stored in the ccfList
regionIdx = ccfList.index("Cortical amygdalar area")

# the next few lines use the try and except functionality to look for a name that 
# isn't present in the list of regions and prints that the region isn't present
# if not found 
regionName = 'Fos'
try:
    regionIdx = ccfList.index(regionName)
    print(regionIdx)
except ValueError:
    # f-strings (formatted as print(f'{someVariable} and some string')) allow 
    # the addition of variables to strings, where the variable is added in 
    # curly brackets {}
    print(f'{regionName} is not present')

# the childNode lines are only for display, and won't be used in next section
childNode = ccfInfo["structure_id_path"][regionIdx]
childNodeSplit = childNode.split("/")
# use a single region name to practice using the split function for strings
nameNode = ccfInfo["name"][160]
# split a string using spaces (" ")
nameNodeSplit = nameNode.split(" ")
# create an empty list to store individual words
wordList = []
for splitName in nameNodeSplit: 
    # casefold converts all letters to lowercase
    print(splitName.casefold())
    commaSplit = splitName.split(",")
    if len(commaSplit) == 1:
        wordList.append(commaSplit[0].casefold())
    elif len(commaSplit) == 2:
        wordList.append(commaSplit[0].casefold())
        

#%% create loop to identify amygdalar areas in list
"""
In this section I would like you to:
    1. Write a for loop that loops over each of the regions contained in ccfInfo['name']
    2. Use the code above to create a list of all words in that make up the 
    region name that you are looping over
    3. Identify whether that region name contains the word "amygdalar"
    4. If the name contains "amygdalar", print the name from ccfList to the screen/console
    5. If the name does NOT contain "amygdalar", print "Does not contain Amygdalar"
"""
ccfIdList = []
for nameNode in ccfList:
    nameNodeSplit = nameNode.split(" ")
    # create an empty list to store individual words
    wordList = []
    for splitName in nameNodeSplit: 
        commaSplit = splitName.split(",")
        if len(commaSplit) == 1:
            wordList.append(commaSplit[0].casefold())
        elif len(commaSplit) == 2:
            wordList.append(commaSplit[0].casefold())
    if 'amygdalar' in wordList:
        # amygIdx = wordList.index("amygdalar")
        # print(nameNode)
        regionIdx = ccfList.index(nameNode)
        ccfIdList.append(ccfInfo['id'][regionIdx])
ccfIdList = np.array(ccfIdList, dtype='int64')



#%% 
def loadRegionalCellsFromSample(h5adLocation, regionName):
    
    # prepare CCF info to look for specific region
    ccfInfo = pd.read_csv(os.path.join(ccfLocation, 'allen_ccf_annotation.csv'))
    
    # regionIdx = np.where(ccfInfo['name'] == regionName)[0]
    # structureID = int(ccfInfo['id'][regionIdx].iloc[0])
    # combStructID = []
    
    # for StructureIDpath in enumerate(ccfInfo['structure_id_path']):
    #     splitStructureID = StructureIDpath[1].split('/')

    #     if np.any(np.array(splitStructureID) == f'{structureID}'):
    #         combStructID.append(ccfInfo['id'][StructureIDpath[0]])
    
    # combStructID = np.array(combStructID)
    ccfList = list(ccfInfo["name"])
    ccfIdList = []
    regionNameList = []
    for nameNode in ccfList:
        nameNodeSplit = nameNode.split(" ")
        # create an empty list to store individual words
        wordList = []
        for splitName in nameNodeSplit: 
            commaSplit = splitName.split(",")
            if len(commaSplit) == 1:
                wordList.append(commaSplit[0].casefold())
            elif len(commaSplit) == 2:
                wordList.append(commaSplit[0].casefold())
        if regionName in wordList:
            regionNameList.append(nameNode)
            regionIdx = ccfList.index(nameNode)
            ccfIdList.append(ccfInfo['id'][regionIdx])

    ccfIdList = np.array(ccfIdList, dtype='int64')
    sample = {}
    sample['regionNameList'] = regionNameList
    print(f'Importing data from {h5adLocation}')
    f = h5py.File(h5adLocation)
    x_attrs = dict(f['X'].attrs)
    # slice_codes = (f['obs']['brain_section_barcode']['codes'][:])
    # slice_codes = f['obs']['brain_section_label']['codes'][:]
    # sliceIdx = np.where(slice_codes == sliceNumber)[0]
    sample['original_filename'] = h5adLocation.split(os.sep)[-1]
    tempCCFID = f['obs']['ccf_id'][:]
    ccfIDMask = np.empty([0,1])
    for i in ccfIdList:
        ccfIDMask = np.array(np.append(ccfIDMask,np.where(tempCCFID == i)[0]), dtype='int32') 
    ccfIDMask = np.sort(ccfIDMask)

    sample['z_coord'] = f['obs']['z_coord'][ccfIDMask]
    # sex and age are stored as strings in the categories section
    sample['sex'] = f['obs']['sex']['categories'][0].decode("utf-8")
    sample['age'] = f['obs']['age']['categories'][0].decode("utf-8")
    sample['polygon_center'] = np.array([f['obs']['polygon_center_x'][ccfIDMask], f['obs']['polygon_center_y'][ccfIDMask], f['obs']['polygon_center_z'][ccfIDMask]]).T

    sample['cell_label'] = f['obs']['cell_label'][ccfIDMask]
    for i in range(len(sample['cell_label'])):
        sample['cell_label'][i] = sample['cell_label'][i].decode("utf-8")

    sample['geneList'] = list(np.array(f['var']['_index']))
    # convert byte strings to strings
    for i in range(len(sample['geneList'])):
        sample['geneList'][i] = sample['geneList'][i].decode("utf-8")
    # ccfIDMask masks gene matrix to only include cells contained in desired ccf region
    try:
        gene_matrix = csr_matrix((f['X']['data'],f['X']['indices'],f['X']['indptr']),x_attrs['shape'])
        gene_matrix = gene_matrix[ccfIDMask,:]
    except ValueError:
        gene_matrix = csr_matrix(f['X'][ccfIDMask,:])
    # sample['geneMatrix'] = np.squeeze(np.array(gene_matrix[:,geneIdx].todense()))
    sample['tissuePositionCoordinates'] = f['obsm']['spatial'][ccfIDMask,0:2]
    sample['ccfCoordinates'] = np.array([f['obs']['x_CCF'][ccfIDMask], f['obs']['y_CCF'][ccfIDMask], f['obs']['z_CCF'][ccfIDMask]]).T
    # cell type class for each cell is stored in the 'codes' section of the 'hrc_mmc_class_name' key
    sample['cellTypeClassInt'] = f['obs']['hrc_mmc_class_name']['codes'][ccfIDMask]
    # the mapping of cell type class integer to a string is stored in the 'categories' section of 'hrc_mmc_class_name'
    sample['cellTypeClassNames'] = f['obs']['hrc_mmc_class_name']['categories'][:]
    for i in range(len(sample['cellTypeClassNames'])):
        sample['cellTypeClassNames'][i] = sample['cellTypeClassNames'][i].decode("utf-8")
    sample['cellTypeClassNames'] = list(sample['cellTypeClassNames'])
    sample['cellTypeClassColorHex'] = np.zeros_like(f['obs']['hrc_mmc_class_color']['codes'][ccfIDMask], dtype=np.dtype('U7'))
    for i in range(len(f['obs']['hrc_mmc_class_color']['categories'][:])):
        hexIdx = np.where(f['obs']['hrc_mmc_class_color']['codes'][ccfIDMask] == i)[0]
        sample['cellTypeClassColorHex'][hexIdx] = f['obs']['hrc_mmc_class_color']['categories'][i].decode("utf-8")
    sample['ccfID'] = f['obs']['ccf_id'][ccfIDMask]
    f.close()
    return sample

h5adLocation = os.path.join('/', 'media','zjpeters', 'Expansion', 'merscopeDataFromAllenInstitute', 'sourcedata', 'mouse_609882_registered_082725.h5ad')
amygSample = loadRegionalCellsFromSample(h5adLocation, 'amygdalar')

#%% create mask of genes present in overlappingList

geneMask = [gene in overlappingGenes for gene in amygSample['geneList']]
maskedGeneList = np.array(amygSample['geneList'])[geneMask].tolist()
#%% plot data in 3d
plt.close('all')
fig = plt.figure(figsize=(7.2,7.2))
ax = fig.add_subplot(projection='3d')
ax.scatter(amygSample['ccfCoordinates'][:,0], amygSample['ccfCoordinates'][:,1], amygSample['ccfCoordinates'][:,2], linewidth=0, s=2, marker='.', alpha=0.7)
ax.axis('equal')
ax.set_xticks([])
ax.set_yticks([])
ax.set_zticks([])
plt.show()
#%% function to create individual samples

def saveRegionalCellsAsSample(h5adLocation, regionName, save_path, gene_list):
    
    sample = {}
    sample['original_filename'] = h5adLocation.split(os.sep)[-1]
    outputBasename = f"{sample['original_filename'].split('.')[0]}_{regionName}"
    sample['sampleID'] = outputBasename
    outputPath = os.path.join(save_path, outputBasename)
    if os.path.exists(os.path.join(outputPath)):
        print(f"Sample {outputBasename} has already been created.")
        print("Loading data now.")
    else:
        os.mkdir(outputPath)
        # prepare CCF info to look for specific region
        ccfInfo = pd.read_csv(os.path.join(ccfLocation, 'allen_ccf_annotation.csv'))
        # overlappingGenes = pd.read_csv(os.path.join(currentLocation, 'data','list_of_overlapping_genes.csv'))
        # regionIdx = np.where(ccfInfo['name'] == regionName)[0]
        # structureID = int(ccfInfo['id'][regionIdx].iloc[0])
        # combStructID = []
        
        # for StructureIDpath in enumerate(ccfInfo['structure_id_path']):
        #     splitStructureID = StructureIDpath[1].split('/')
    
        #     if np.any(np.array(splitStructureID) == f'{structureID}'):
        #         combStructID.append(ccfInfo['id'][StructureIDpath[0]])
        
        # combStructID = np.array(combStructID)
        ccfList = list(ccfInfo["name"])
        ccfIdList = []
        regionNameList = []
        for nameNode in ccfList:
            nameNodeSplit = nameNode.split(" ")
            # create an empty list to store individual words
            wordList = []
            for splitName in nameNodeSplit: 
                commaSplit = splitName.split(",")
                if len(commaSplit) == 1:
                    wordList.append(commaSplit[0].casefold())
                elif len(commaSplit) == 2:
                    wordList.append(commaSplit[0].casefold())
            if regionName in wordList:
                regionNameList.append(nameNode)
                regionIdx = ccfList.index(nameNode)
                ccfIdList.append(ccfInfo['id'][regionIdx])
    
        ccfIdList = np.array(ccfIdList, dtype='int64')
        
        sample['regionNameList'] = regionNameList
        print(f'Importing data from {h5adLocation}')
        f = h5py.File(h5adLocation)
        x_attrs = dict(f['X'].attrs)
        # slice_codes = (f['obs']['brain_section_barcode']['codes'][:])
        # slice_codes = f['obs']['brain_section_label']['codes'][:]
        # sliceIdx = np.where(slice_codes == sliceNumber)[0]
        
        tempCCFID = f['obs']['ccf_id'][:]
        ccfIDMask = np.empty([0,1])
        for i in ccfIdList:
            ccfIDMask = np.array(np.append(ccfIDMask,np.where(tempCCFID == i)[0]), dtype='int32') 
        ccfIDMask = np.sort(ccfIDMask)
    
        sample['z_coord'] = f['obs']['z_coord'][ccfIDMask]
        # sex and age are stored as strings in the categories section
        sample['sex'] = f['obs']['sex']['categories'][0].decode("utf-8")
        sample['age'] = f['obs']['age']['categories'][0].decode("utf-8")
        sample['polygon_center'] = np.array([f['obs']['polygon_center_x'][ccfIDMask], f['obs']['polygon_center_y'][ccfIDMask], f['obs']['polygon_center_z'][ccfIDMask]]).T
    
        sample['cell_label'] = f['obs']['cell_label'][ccfIDMask]
        for i in range(len(sample['cell_label'])):
            sample['cell_label'][i] = sample['cell_label'][i].decode("utf-8")
    
        sample['geneList'] = list(np.array(f['var']['_index']))
        # convert byte strings to strings
        geneIdxMask = []
        for i in range(len(sample['geneList'])):
            sample['geneList'][i] = sample['geneList'][i].decode("utf-8")
            if sample['geneList'][i] in gene_list:
                geneIdxMask.append(i)
        geneIdxMask = np.array(geneIdxMask, dtype='int64')
        sample['geneList'] = np.array(sample['geneList'])[geneIdxMask].tolist()
        # geneMask = [gene in overlappingGenes for gene in sample['geneList']]
        
        # ccfIDMask masks gene matrix to only include cells contained in desired ccf region
        try:
            gene_matrix = csr_matrix((f['X']['data'],f['X']['indices'],f['X']['indptr']),x_attrs['shape'])
            gene_matrix = gene_matrix[ccfIDMask,:]
        except ValueError:
            gene_matrix = csr_matrix(f['X'][ccfIDMask,:])
        sample['geneMatrix'] = gene_matrix[:,geneIdxMask]

        sample['tissuePositionCoordinates'] = f['obsm']['spatial'][ccfIDMask,0:2]
        sample['ccfCoordinates'] = np.array([f['obs']['x_CCF'][ccfIDMask], f['obs']['y_CCF'][ccfIDMask], f['obs']['z_CCF'][ccfIDMask]]).T
        # cell type class for each cell is stored in the 'codes' section of the 'hrc_mmc_class_name' key
        sample['cellTypeClassInt'] = f['obs']['hrc_mmc_class_name']['codes'][ccfIDMask]
        # the mapping of cell type class integer to a string is stored in the 'categories' section of 'hrc_mmc_class_name'
        sample['cellTypeClassNames'] = f['obs']['hrc_mmc_class_name']['categories'][:]
        for i in range(len(sample['cellTypeClassNames'])):
            sample['cellTypeClassNames'][i] = sample['cellTypeClassNames'][i].decode("utf-8")
        sample['cellTypeClassNames'] = list(sample['cellTypeClassNames'])
        sample['cellTypeClassColorHex'] = np.zeros_like(f['obs']['hrc_mmc_class_color']['codes'][ccfIDMask], dtype=np.dtype('U7'))
        for i in range(len(f['obs']['hrc_mmc_class_color']['categories'][:])):
            hexIdx = np.where(f['obs']['hrc_mmc_class_color']['codes'][ccfIDMask] == i)[0]
            sample['cellTypeClassColorHex'][hexIdx] = f['obs']['hrc_mmc_class_color']['categories'][i].decode("utf-8")
        sample['ccfID'] = f['obs']['ccf_id'][ccfIDMask]
        f.close()
        # save json file
        sampleDataDict = {
            "sampleID": sample['sampleID'],
            "original_filename": sample['original_filename'],
            "age": sample['age'],
            "sex": sample['sex'],
            "geneList": sample['geneList'],
            "cellTypeClassNames": sample['cellTypeClassNames'],
            "ccfRegionsIncluded": sample['regionNameList']
        }
        # Serializing json
        json_object = json.dumps(sampleDataDict, indent=4)
         
        # Writing to processing_information.json
        with open(os.path.join(outputPath, f"{outputBasename}_processing_information.json"), "w") as outfile:
            outfile.write(json_object)
            
        pts = {'x_CCF': sample['ccfCoordinates'][:,0], 'y_CCF': sample['ccfCoordinates'][:,1], 
               'z_CCF': sample['ccfCoordinates'][:,2], "cellTypeClassInt": sample['cellTypeClassInt'], 
               "cellTypeClassColorHex": sample['cellTypeClassColorHex'], "ccfID": sample['ccfID'], 
               "cell_label": sample['cell_label']}
        pts = pd.DataFrame(pts)
        pts.to_csv(os.path.join(outputPath, f"{outputBasename}_ccf_coordinates.csv"), index=False)
        sp_sparse.save_npz(os.path.join(outputPath, f"{outputBasename}_gene_matrix"), sample['geneMatrix'])
    return sample

amygSample = saveRegionalCellsAsSample(h5adLocation, 'amygdalar', os.path.join('/','media','zjpeters','Expansion','eklein','rawdata'), overlappingGenes)

#%% loop to generate data for all 4 samples

filenames = ['mouse_609882_registered_082725.h5ad', 'mouse_609889_registered_082725.h5ad',
             'mouse_638850_registered_082725.h5ad', 'mouse_687997_registered_082725.h5ad']

for filename in filenames:
    fileLocation = h5adLocation = os.path.join('/', 'media','zjpeters', 'Expansion', 'merscopeDataFromAllenInstitute', 'sourcedata', filename)
    amygSample = saveRegionalCellsAsSample(h5adLocation, 'amygdalar', os.path.join('/','media','zjpeters','Expansion','eklein','rawdata'), overlappingGenes)

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
    geneMatrixLoc = glob(os.path.join(save_path, "*_gene_matrix.npz"))[0]
    sample['geneMatrix'] = sp_sparse.load_npz(geneMatrixLoc)
    return sample
        
amygSample = load3dRegionalSample(os.path.join('/','media','zjpeters','Expansion','eklein','rawdata','mouse_609882_registered_082725_amygdalar'))







