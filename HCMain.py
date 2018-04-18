# -*- coding: utf-8 -*-
"""
Created on Wed Jan 17 23:58:05 2018

Principal script to run the unidimensional heat transfer test.

@author: Leandro
"""


"""
Import statements:
""" 

import HCClass
import HCConductionLoop
import HCLibraryMaterial


"""
Setting constants:
"""

# Bar constants.
lengthBar = 1.0
alphaHeat = 1.0
nodes = 5
initTemp = 1.0

# Time constants.
timeEnd = 0.5
r = 0.05

# Heat conduction process.
T0 = 1.0
qIn = 0.0
model = "kCN" # kSimple, kProm, kCN

# Main function to run the scheme.
def main():

    """
    Creating the objects.
    """
    
    # Material.
    libraryWater = HCLibraryMaterial.water
    water = HCClass.Material(libraryWater)
    water.plotK()
    
    # Bar.
    bar = HCClass.Bar(lengthBar, alphaHeat, nodes, initTemp, timeEnd, r)
    bar.setInitialDist([1.0 + bar.nodes], [0])
    
    
    """
    Starting the main loop.
    """
    
    # model, material, bar, time, T0, qIn
    lastTemperature = HCConductionLoop.heatConductionLoop(model, water, bar, T0, qIn)
    
main()