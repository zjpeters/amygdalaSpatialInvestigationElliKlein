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
import umap 
from matplotlib import pyplot as plt
import matplotlib.patches as mpatches 
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
#%% write a function to find all cells of a cell type
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



#%% function to count cells of different cell types

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
totalNumberOfCells = countAllCellsOfCellType(allSamples['mouse_609882_registered_082725_amygdalar'], "Glut")

# 1. recreate a for Loop for each sample to count Glut cells
# 2. rename x = totalNumberOfCells and append it to a list 
# 3. Write a code that converts raw numbers to ratio %


totalNumberOfCellsList = []
cellsRatioList = []
for sampleID in foldernames:
    sample = allSamples[sampleID]
   # allSamples[sampleID]
    totalNumberOfCells = countAllCellsOfCellType(sample, "Glut")
    totalNumberOfCellsList.append(totalNumberOfCells)
    cellsRatio = (totalNumberOfCells/len(sample["cellTypeClassInt"]))* 100
    cellsRatioList.append(cellsRatio)

# raw numbers to ratio with Glut cells

#cellsRatioList = []
#len(sample["cellTypeClassInt"])
#ratio = ("totalNumberOfCells"/len(sample["cellTypeClassInt"]))* 100
#%%
def createCellTypeMask(sample, cellType):
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
           
    allCellTypeIdx = np.empty(0, dtype="int64")
    for cellTypeIdx in cellTypeIdxList:
        allCellTypeIdx = np.append(allCellTypeIdx, np.where(sample['cellTypeClassInt'] == cellTypeIdx)[0])
        print(sample["cellTypeClassNames"][cellTypeIdx])
    allCellTypeIdx = np.sort(allCellTypeIdx)
    return allCellTypeIdx
cellIdxArray = createCellTypeMask(allSamples['mouse_609882_registered_082725_amygdalar'], "Glut")
#%% set up umap and save it

umapSaveLoc = os.path.join(derivatives, f"{allSamples['mouse_609882_registered_082725_amygdalar']['sampleID']}_umap_embedding.csv")
reducer = umap.UMAP(random_state=9)
if os.path.exists(umapSaveLoc):
    embedding = pd.read_csv(umapSaveLoc)
    embedding = np.array(embedding, dtype= "float32")
else: 
    embedding = reducer.fit_transform(allSamples['mouse_609882_registered_082725_amygdalar']["geneMatrix"][cellIdxArray, :])
    embeddingdf = pd.DataFrame(embedding)
    embeddingdf.to_csv(umapSaveLoc, index = False)
#%% plotting Umap

plt.close("all")
plt.scatter(embedding[:, 0], embedding[:, 1], s = 1, linewidth= 0)
plt.show()
#%% plot with colors based on cellType

hexdf = pd.read_csv(os.path.join(currentLocation, "data", "ccf_color_hex_codes.csv"))
plt.close("all")
fig, ax = plt.subplots(1,1)
ax.scatter(embedding[:, 0], embedding[:, 1], s = 1, linewidth= 0, c = np.array(allSamples['mouse_609882_registered_082725_amygdalar']["cellTypeClassColorHex"])[cellIdxArray])
uniqueCellColors = np.unique(np.array(allSamples['mouse_609882_registered_082725_amygdalar']["cellTypeClassColorHex"])[cellIdxArray])
uniqueCellIdx = np.unique(np.array(allSamples['mouse_609882_registered_082725_amygdalar']["cellTypeClassInt"])[cellIdxArray])

"""
in order to display proper colors on the UMAP scatter plots, need to add the 
scatters group by group, rather than all at once. To do this, the ax.scatter 
function above needs to be put inside the loop, and cells only added when they
are from a particular mask. 

1. Loop over the uniqueCellIdx and use the np.where(...)[0] to create a mask to 
    use to mask the cells in your scatter plot, such as:
        ax.scatter(embedding[newMask, 0], embedding[newMask,1]...)
2. In the same loop, use the uniqueCellIdx to select the correct color from the 
    hexdf dataframe and use that as a single input for c= in your ax.scatter
    as well as for the color in mpatches.Patch
"""

#%%
umapFigLoc = os.path.join(derivatives, f"{allSamples['mouse_609882_registered_082725_amygdalar']['sampleID']}_Glut_umap_figure.pdf")
plt.close("all")
fig, ax = plt.subplots(1,1, figsize = [7.2, 7.2])
handles, labels= ax.get_legend_handles_labels()
for i, currentCell in enumerate(uniqueCellIdx):
    color = uniqueCellColors[i]
    newMask = np.where(np.array(allSamples['mouse_609882_registered_082725_amygdalar']["cellTypeClassInt"])[cellIdxArray] == currentCell)[0]
    print(f"n of {hexdf['cell_type'][currentCell]}: ",len(newMask))
    ax.scatter(embedding[newMask, 0], embedding[newMask, 1], c = color, s = 1, linewidth = 0)
    #hexdf[hexdf["cell_type"] == currentCell]
    # h = hexdf["color_hex"][i].lstrip("#")
    # rgbCode = np.array(tuple(int(h[i:i+2], 16) for i in (0,2,4)))/255
    patch = mpatches.Patch(color = color, label = hexdf["cell_type"][currentCell])
    print(patch)
    handles.append(patch)
    print(handles)
ax.legend(handles = handles, ncols = 3, fontsize = 7)
plt.title(f"UMAP of all glutaminergic cells in amygdala of: \n {allSamples['mouse_609882_registered_082725_amygdalar']['sampleID']}, {allSamples['mouse_609882_registered_082725_amygdalar']['age']}, {allSamples['mouse_609882_registered_082725_amygdalar']['sex']} ")
plt.show()
plt.savefig(umapFigLoc, dpi = 300, bbox_inches = "tight")

#%% create empty array for all samples then for loop to add gene matrices to the array
# with each sample identifiy only Glut cells, append only Glut cells of Gene Matrix, use same mask 
#-> cellIdxArray = createCellTypeMask(allSamples['mouse_609882_registered_082725_amygdalar'], "Glut")

allSubGeneMatrix = np.empty([0, 279])
allSubSexArray = np.empty(0, dtype= "str")
allSubColorArray = np.empty(0, dtype= "str")
for sampleID in foldernames:
    cellIdxArray = createCellTypeMask(allSamples[sampleID], "Glut")
    allSubGeneMatrix = np.append(allSubGeneMatrix, allSamples[sampleID]["geneMatrix"][cellIdxArray, :].todense(), axis= 0)
    subSexArray = np.full(len(cellIdxArray), dtype= "str", fill_value = allSamples[sampleID]["sex"])
    allSubSexArray = np.append(allSubSexArray,np.array(subSexArray), axis =0)
    allSubColorArray = np.append(allSubColorArray, allSamples[sampleID]["cellTypeClassColorHex"][cellIdxArray], axis= 0)
allSubGeneMatrix = sp_sparse.csr_matrix(allSubGeneMatrix)

#%% create umap with allSubGeneMatrix 
# 1. with cellTypeColors, one umap for all samples for just one global cell type, Glut or GABA
# 2. plot male and female, either different colors or two different umaps side by side
# 3. save umap and images as pdf
    
    
        
    












