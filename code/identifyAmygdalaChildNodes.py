#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 27 07:35:40 2026

@author: zjpeters
"""
import os
import pandas as pd

# next line identifies the path for the directory where this script is located,
# since it will be different between computers, i.e. my computer vs yours
currentLocation = os.path.dirname(os.path.abspath(__file__))
# this line combines the above location with where directory where the 
# "allen_ccf_annotation.csv" file is located
ccfLocation = os.path.join(currentLocation,'data','ccf')
    
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

for nameNode in ccfList:
    print(nameNode)
