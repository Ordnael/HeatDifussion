# -*- coding: utf-8 -*-
"""
Created on Thu Jan 18 00:19:38 2018

Script with constructor classes for time, bar and material objects.

@author: Leandro
"""


"""
Import statements:
"""        

import matplotlib.pyplot as plt
import numpy as np


"""
Class definitions:
"""        

# Bar:
class Bar:
    
    # Inputs: Lenght, thermal expansion coefficient, number of nodes,
    # initial temperature of bar, end time, r value.
    def __init__(self, lengthBar, alphaHeat, nodes, temperature, timeEnd, rVal):
        self.lengthBar = lengthBar
        self.alphaHeat = alphaHeat
        self.nodes = nodes
        self.timeEnd = timeEnd
        self.rVal = rVal
        # Constructor operations for length discretization.
        self.lenVector = np.linspace(0, lengthBar, nodes)
        self.increment = lengthBar/(nodes - 1.0)
        # Constructor operations for time discretization.
        self.delta = 2*rVal*np.power(self.increment,2.0)/alphaHeat
        self.timeVector = np.linspace(0, timeEnd, ((timeEnd/self.delta) + 1))
        # Constructor operations for temperature matrix.
        iCoordTemp = len(self.timeVector)
        jCoordTemp = len(self.lenVector)
        self.temperature = temperature*np.ones([iCoordTemp, jCoordTemp])
        
        
    # Method to input a temperature distribution.
    # Inputs: Temperature values, nodes to apply.     
    def setInitialDist(self, newTemperature, nodesDist):
        for i in range(len(newTemperature)):
            self.temperature[0][nodesDist[i]] = newTemperature[i]
                     
# Material:
class Material:
    
    # Input: Material values for T, K and degree from library.    
    def __init__(self, matLibrary):
        self.temperatures = matLibrary['temp']
        self.condCoef = matLibrary['coef']
        self.degree = matLibrary['degr']
        self.poly = np.polyfit(self.temperatures, self.condCoef, self.degree)
        self.eqn = np.poly1d(self.poly)
        
    # Method to plot the regression line and the input points.
    def plotK(self):
        xCoord = np.linspace(self.temperatures[0], self.temperatures[-1], 100)
        yCoord = self.eqn(xCoord)
        
        plt.figure()
        plt.plot(xCoord, yCoord, label="Regression curve")
        plt.plot(self.temperatures, self.condCoef, 'b.', label="Input points")
        plt.xlabel('Temperature'); plt.ylabel('Condution Coefficient')
        plt.legend()
        plt.axis("tight")
        plt.title("Fitted curve against table points for K")