#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 27 07:35:40 2026

@author: zjpeters
"""
import os
# checks whether using windows or linux
if os.path.exists(os.path.join('/','media','zjpeters','Expansion')):
    externalDrivePath = os.path.join('/','media','zjpeters','Expansion')
else:
    externalDrivePath = os.path.join('D:',os.sep)
import h5py
import numpy as np
from matplotlib import pyplot as plt
import pandas as pd
import sys
sys.path.insert(0, os.path.join(externalDrivePath, 'merscopeDataFromAllenInstitute','code'))
from scipy.sparse import csr_matrix
import imageio
"""
since the allen data is all contained as a single dataset for all volumes, could
create the FOVs using the CCF coordinates, therefore allowing for 3d analysis 
"""
rawdata = os.path.join(externalDrivePath, 'asymmetryAllenData','rawdata')
derivatives = os.path.join(externalDrivePath, 'asymmetryAllenData','derivatives')

if os.path.exists(os.path.join(os.sep, 'home','zjpeters', 'Documents', 'stanly','code', 'data', 'ccf')):
    stanlyLocation = os.path.join(os.sep, 'home','zjpeters', 'Documents', 'stanly','code', 'data', 'ccf') 
else:
    stanlyLocation = os.path.join('C:', os.sep, 'Users','onyh19ug', 'Documents', 'STANLY','code', 'data', 'ccf') 
    
#%% load gene list and cell type lists

h5adLocation = os.path.join(externalDrivePath, 'merscopeDataFromAllenInstitute', 'sourcedata', 'mouse_609882_registered_082725.h5ad')
f = h5py.File(h5adLocation)
x_attrs = dict(f['X'].attrs)
geneList = list(np.array(f['var']['_index']))
# convert byte strings to strings
for i in range(len(geneList)):
    geneList[i] = geneList[i].decode("utf-8")

cellTypeClassNames  = f['obs']['hrc_mmc_class_name']['categories'][:]
for i in range(len(cellTypeClassNames)):
    cellTypeClassNames[i] = cellTypeClassNames[i].decode("utf-8")
#%% create necessary functions
def createRegionalCCFMask(sample, regionName):
    ccfInfo = pd.read_csv(os.path.join(stanlyLocation, 'allen_ccf_annotation.csv'))
    
    regionIdx = np.where(ccfInfo['name'] == regionName)[0]
    structureID = int(ccfInfo['id'][regionIdx].iloc[0])
    combStructID = []
    
    for StructureIDpath in enumerate(ccfInfo['structure_id_path']):
        splitStructureID = StructureIDpath[1].split('/')

        if np.any(np.array(splitStructureID) == f'{structureID}'):
            combStructID.append(ccfInfo['id'][StructureIDpath[0]])
    
    combStructID = np.array(combStructID)
    
    ccfIDMask = np.empty([0,1])
    for i in combStructID:
        ccfIDMask = np.array(np.append(ccfIDMask,np.where(sample['ccfID'] == i)[0]), dtype='int32')        

    return ccfIDMask

def loadRegionalCellsFromSample(h5adLocation, regionName, geneName):
    
    # prepare CCF info to look for specific region
    ccfInfo = pd.read_csv(os.path.join(stanlyLocation, 'allen_ccf_annotation.csv'))
    
    regionIdx = np.where(ccfInfo['name'] == regionName)[0]
    structureID = int(ccfInfo['id'][regionIdx].iloc[0])
    combStructID = []
    
    for StructureIDpath in enumerate(ccfInfo['structure_id_path']):
        splitStructureID = StructureIDpath[1].split('/')

        if np.any(np.array(splitStructureID) == f'{structureID}'):
            combStructID.append(ccfInfo['id'][StructureIDpath[0]])
    
    combStructID = np.array(combStructID)
    
    sample = {}
    print(f'Importing data from {h5adLocation}')
    f = h5py.File(h5adLocation)
    x_attrs = dict(f['X'].attrs)
    # slice_codes = (f['obs']['brain_section_barcode']['codes'][:])
    # slice_codes = f['obs']['brain_section_label']['codes'][:]
    # sliceIdx = np.where(slice_codes == sliceNumber)[0]
    sample['original_filename'] = h5adLocation.split(os.sep)[-1]
    tempCCFID = f['obs']['ccf_id'][:]
    ccfIDMask = np.empty([0,1])
    for i in combStructID:
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
    geneIdx = sample['geneList'].index(geneName)
    try:
        gene_matrix = csr_matrix((f['X']['data'],f['X']['indices'],f['X']['indptr']),x_attrs['shape'])
        gene_matrix = gene_matrix[ccfIDMask,:]
    except ValueError:
        gene_matrix = csr_matrix(f['X'][ccfIDMask,:])
    sample['geneMatrix'] = np.squeeze(np.array(gene_matrix[:,geneIdx].todense()))
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

#%% load single samples data
hippSample = loadRegionalCellsFromSample(h5adLocation, 'Cortical amygdalar area', 'Pvalb')

#%%

listOfFiles = ['mouse_609882_registered_082725.h5ad', 'mouse_609889_registered_082725.h5ad', 
               'mouse_638850_registered_082725.h5ad', 'mouse_687997_registered_082725.h5ad']

geneChoice = 'Maob'
allHippCoords = np.empty((0,3))
allHippGene = np.empty((0,1), dtype='uint16')
allHippColors = np.empty((0,1), dtype='str')
for h5adFilename in listOfFiles:
    h5adLocation = os.path.join(externalDrivePath, 'merscopeDataFromAllenInstitute', 'sourcedata', h5adFilename)
    hippSample = loadRegionalCellsFromSample(h5adLocation, 'Cortical amygdalar area', geneChoice)
    posMask = np.where(hippSample['geneMatrix'] > 0)[0]
    allHippCoords = np.append(allHippCoords, hippSample['ccfCoordinates'][posMask,:], axis=0)
    minMax = (hippSample['geneMatrix'][posMask] - np.min(hippSample['geneMatrix'][posMask]))/(np.max(hippSample['geneMatrix'][posMask])-np.min(hippSample['geneMatrix'][posMask]))
    allHippGene = np.append(allHippGene, minMax)
    # allHippCoords = np.append(allHippCoords, hippSample['ccfCoordinates'], axis=0)
    # allHippColors = np.append(allHippColors, hippSample['cellTypeClassColorHex'])


# #%% plot coordinates
# allHippSortIdx = np.sort(allHippColors)
# plt.close('all')
# fig = plt.figure(figsize=(7.2,7.2))
# ax = fig.add_subplot(projection='3d')
# ax.view_init(elev=-45, azim=-80, roll=15)
# ax.scatter(allHippCoords[:,0], allHippCoords[:,1], allHippCoords[:,2], linewidth=0, s=1, c=allHippGene, cmap='Reds', vmin=0, vmax=15, marker='.')
# # ax.scatter(allHippCoords[:,0], allHippCoords[:,1], allHippCoords[:,2], linewidth=0, s=1, c=allHippColors, marker='.')
# ax.axis('equal')
# plt.show()

#%% try loop

plt.close('all')
n = 0
for i in range(360, 0, -10):
    fig = plt.figure(figsize=(7.2,7.2))
    ax = fig.add_subplot(projection='3d')
    ax.view_init(elev=i, azim=0, roll=-90)
    ax.scatter(allHippCoords[:,0], allHippCoords[:,1], allHippCoords[:,2], linewidth=0, s=2, c=allHippGene, cmap='Reds', vmin=0, marker='.', alpha=0.7)
    # ax.scatter(allHippCoords[:,0], allHippCoords[:,1], allHippCoords[:,2], linewidth=0, s=1, c=allHippColors, marker='.')
    ax.axis('equal')
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_zticks([])
    plt.show()
    plt.savefig(os.path.join(derivatives, 'gifImages', f'image_{n}_{geneChoice}.png'), bbox_inches='tight', dpi=200)
    n += 1
    plt.close()
    
allImages = []
for n in range(0, 36):
    allImages.append(imageio.v3.imread(os.path.join(derivatives, 'gifImages', f'image_{n}_{geneChoice}.png')))

imageio.v3.imwrite(os.path.join(derivatives, 'gifImages', f'{geneChoice}_gif.gif'), allImages, duration = 500, loop = 1)



#%% create pngs for cell  type images
# allHippCoords = np.empty((0,3))
# allHippGene = np.empty((0,1), dtype='uint16')
# allHippColors = np.empty((0,1), dtype='str')
# for h5adFilename in listOfFiles:
#     h5adLocation = os.path.join(externalDrivePath, 'merscopeDataFromAllenInstitute', 'sourcedata', h5adFilename)
#     hippSample = loadRegionalCellsFromSample(h5adLocation, 'Hippocampal formation', 'Gfap')
#     allHippCoords = np.append(allHippCoords, hippSample['ccfCoordinates'], axis=0)
#     allHippColors = np.append(allHippColors, hippSample['cellTypeClassColorHex'])

# create gif using cell types

# plt.close('all')
# n = 0
# for i in range(360, 0, -10):
#     fig = plt.figure(figsize=(7.2,7.2))
#     ax = fig.add_subplot(projection='3d')
#     ax.view_init(elev=i, azim=0, roll=-90)
#     # ax.scatter(allHippCoords[:,0], allHippCoords[:,1], allHippCoords[:,2], linewidth=0, s=1, c=allHippGene, cmap='Reds', vmin=0, vmax=15, marker='.')
#     ax.scatter(allHippCoords[:,0], allHippCoords[:,1], allHippCoords[:,2], linewidth=0, s=1, c=allHippColors, marker='.')
#     ax.axis('equal')
#     ax.set_xticks([])
#     ax.set_yticks([])
#     ax.set_zticks([])
#     plt.show()
#     plt.savefig(os.path.join(derivatives, 'gifImages', f'image_{n}.png'), bbox_inches='tight', dpi=200)
#     n += 1
#     plt.close()

# # write gif
# allImages = []
# for n in range(0, 36):
#     allImages.append(imageio.v3.imread(os.path.join(derivatives, 'gifImages', f'image_{n}.png')))

# imageio.v3.imwrite(os.path.join(derivatives, 'gifImages', 'cell_type_gif.gif'), allImages, duration = 500, loop = 1)

