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
# import matplotlib.colors as pltcolors
import matplotlib.cm as cm
currentLocation = os.path.dirname(os.path.abspath(__file__))
# added a section that checks whether code is running on your computer or mine
if os.path.exists(os.path.join('C:',os.sep,'Users','Rena')):
    rawdata = os.path.join('C:',os.sep,'Users','Rena','Downloads','Uni','Spatial Transcriptomics',"amygdalaSpatialInvestigationElliKlein","rawdata")
    derivatives = os.path.join('C:',os.sep,'Users','Rena','Downloads','Uni','Spatial Transcriptomics',"amygdalaSpatialInvestigationElliKlein","derivatives")
elif os.path.exists(os.path.join(os.sep,'media','zjpeters','Expansion')):
    rawdata = os.path.join(os.sep,'media','zjpeters','Expansion', 'eklein','rawdata')
    derivatives = os.path.join(os.sep,'media','zjpeters','Expansion', 'eklein','derivatives')
#%% function to load data for regional samples

def load3dRegionalSample(save_path, remove_blanks = True):
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
    if remove_blanks == True:
        geneIdxArray = []
        geneList = []
        for i, gene in enumerate(sample["geneList"]):
            if "Blank-" not in gene:
                geneList.append(gene)
                geneIdxArray.append(i)
        geneIdxArray = np.array(geneIdxArray, dtype = "int32")
        sample["geneList"] = geneList
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
    if remove_blanks == True:
        sample["geneMatrix"]= sample["geneMatrix"][:, geneIdxArray]
    return sample

#%% load samples

foldernames = ['mouse_609882_registered_082725_amygdalar', 'mouse_609889_registered_082725_amygdalar',
             'mouse_638850_registered_082725_amygdalar', 'mouse_687997_registered_082725_amygdalar']

allSamples = {}
for i, foldername in enumerate(foldernames):
    print(i, foldername)
    allSamples[foldername] = load3dRegionalSample(os.path.join(rawdata, foldername))
    
#%%
geneIdxArray = []
for i, gene in enumerate(allSamples['mouse_609882_registered_082725_amygdalar']["geneList"]):
    if "Blank-" not in gene:
    
        print(i, gene)
        geneIdxArray.append(i)
geneIdxArray = np.array(geneIdxArray, dtype = "int32")

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

allSubIDArray = np.empty(0, dtype = "object")
allSubGeneMatrix = np.empty([0, len(allSamples['mouse_609882_registered_082725_amygdalar']["geneList"])])
allSubSexArray = np.empty(0, dtype= "str")
allSubColorArray = np.empty(0, dtype= "str")
for sampleID in foldernames:
    cellIdxArray = createCellTypeMask(allSamples[sampleID], "Glut")
    allSubGeneMatrix = np.append(allSubGeneMatrix, allSamples[sampleID]["geneMatrix"][cellIdxArray, :].todense(), axis= 0)
    subSexArray = np.full(len(cellIdxArray), dtype= "str", fill_value = allSamples[sampleID]["sex"])
    allSubSexArray = np.append(allSubSexArray,np.array(subSexArray), axis =0)
    allSubColorArray = np.append(allSubColorArray, allSamples[sampleID]["cellTypeClassColorHex"][cellIdxArray], axis= 0)
    subIDArray = np.full(len(cellIdxArray), dtype= "object", fill_value = allSamples[sampleID]["sampleID"])
    allSubIDArray = np.append(allSubIDArray, np.array(subIDArray), axis = 0)
allSubGeneMatrix = sp_sparse.csr_matrix(allSubGeneMatrix)
del(subSexArray, subIDArray, cellIdxArray)

#%% create umap with allSubGeneMatrix 
# 1. with cellTypeColors, create one umap for all samples for just one global cell type -> Glut or GABA
# 2. plot male and female, either different colors or two different umaps side by side
# 3. save umap and save images as pdf
   
# creating umap for glut cells in all samples
globalUmapCsvLoc = os.path.join(derivatives, "global_Glut_umap_embedding.csv") 
if os.path.exists(globalUmapCsvLoc):
    embedding = pd.read_csv(globalUmapCsvLoc)
    embedding = np.array(embedding, dtype= "float32")
else: 
    embedding = reducer.fit_transform(allSubGeneMatrix)
    embeddingdf = pd.DataFrame(embedding)
    embeddingdf.to_csv(globalUmapCsvLoc, index = False)
globalUmapFigLoc = os.path.join(derivatives, "global_Glut_umap_figure.pdf")
plt.close("all")
plt.scatter(embedding[:, 0], embedding[:, 1], s = 1, linewidth= 0, c = allSubColorArray)
plt.title("global UMAP of Glut cells")
plt.savefig(globalUmapFigLoc)
plt.show()

#%%
# creating one umap with all Glut cells in male and female, using different colors for M and F
maleIdxArray = allSubSexArray == "M"
femaleIdxArray = allSubSexArray == "F"

#embedding[maleIdxArray, :]
#embedding[femaleIdxArray, :]
umapMaleFemaleFigLoc = os.path.join(derivatives, "male_female_Glut_umap_figure.pdf")
plt.close("all")
plt.scatter(embedding[maleIdxArray][:,0], embedding[maleIdxArray][:,1], s = 1, linewidth= 0, c = "blue", label = "Male" )
plt.scatter(embedding[femaleIdxArray][:,0], embedding[femaleIdxArray][:,1], s = 1, linewidth= 0, c = "red", label = "Female" )
plt.title("UMAP of Glut cells in Male and Female")
# adding markerscale argument allows us to change the size of the markers in the legend
plt.legend(markerscale=5)
plt.savefig(umapMaleFemaleFigLoc)
plt.show()

#%%
# creating two umaps side by side of all Glut cells in male and female
umapMaleFemaleSepFigLoc = os.path.join(derivatives, "male_female_seperated_Glut_umap_figure.pdf")
plt.close("all")
fig, ax = plt.subplots(1,2, figsize = [7.2, 7.2])
ax[0].scatter(embedding[maleIdxArray][:,0], embedding[maleIdxArray][:,1], s = 1, linewidth= 0, c = "blue", label = "Male")
ax[1].scatter(embedding[femaleIdxArray][:,0], embedding[femaleIdxArray][:,1], s = 1, linewidth= 0, c = "red", label = "female")
ax[0].set_title("Male")
ax[1].set_title("Female")
fig.suptitle("UMAPS of Glut cells in Male and Female")
plt.savefig(umapMaleFemaleSepFigLoc)
plt.show()

#%% plotting each subject with assigned color
umapAllSamples = os.path.join(derivatives, "all_samples_umap_figure.pdf")
plt.close("all")
for sampleID in foldernames:
    sampleIDIdxArray = allSubIDArray == sampleID
    #embedding[maleIdxArray, :]
    #embedding[femaleIdxArray, :]
    plt.scatter(embedding[sampleIDIdxArray][:,0], embedding[sampleIDIdxArray][:,1], s = 1, linewidth= 0, label = sampleID )
    
    plt.title("UMAP of all samples")
plt.legend(markerscale=5)
plt.savefig(umapAllSamples)
plt.show()

#%% one sample slices
plt.close('all')
uniqueSampleSlices = np.unique(allSamples['mouse_609882_registered_082725_amygdalar']["z_coord"])
"""
realized the previous code I'd written using this approach used matplotlib.cm
instead of matplotlib.colors, which the stack overflow we saw used. this works
"""
### can uncomment the line below and comment out the three lines after "create colormap..."
### to use a colormap that is spread across a single colormap chosen from:
### https://matplotlib.org/stable/users/explain/colors/colormaps.html
# colors = cm.Spectral(np.arange(len(uniqueSampleSlices)) / len(uniqueSampleSlices))
# create colormap that combines tab20b and tab20c to get a total of 40 (mostly) unique colors
colors = np.empty([0,4])
colors = np.append(colors, cm.tab20b(np.arange(20)/20), axis=0)
colors = np.append(colors, cm.tab20c(np.arange(20)/20), axis=0)
cellIdxArray = createCellTypeMask(allSamples['mouse_609882_registered_082725_amygdalar'], "Glut")
sampleIDIdxArray = allSubIDArray == 'mouse_609882_registered_082725_amygdalar'
sampleEmbedding = embedding[sampleIDIdxArray, :]
fig, ax = plt.subplots(1,1, figsize=[7.2, 5.2])
for i, z_cord in enumerate(uniqueSampleSlices):
    sliceIdxArray = allSamples['mouse_609882_registered_082725_amygdalar']["z_coord"][cellIdxArray] == z_cord
    ax.scatter(sampleEmbedding[sliceIdxArray][:,0], sampleEmbedding[sliceIdxArray][:,1], s = 1, linewidth= 0, label = z_cord, color=colors[i])

plt.title("UMAP of all samples")
plt.legend(markerscale=3, ncols = 3, fontsize=6)
#plt.savefig(umapAllSamples)
plt.show()
    










