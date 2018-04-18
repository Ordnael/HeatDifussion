# -*- coding: utf-8 -*-
"""
Created on Thu Jan 18 00:45:05 2018

Script with principal loop to iterate temperature values over time.

@author: Leandro
"""


"""
Import statements:
""" 

import matplotlib.pyplot as plt
import numpy as np
import HCClass
import HCConductionFunction as HCC


"""
Main loop of heat conduction:
"""

# Inputs: Model of conduction, material model, bar and time objects,
# reference temperature, input flash heat.
def heatConductionLoop(model, material, bar, T0, qIn):
    
    print(bar.lenVector)
    print(bar.temperature)
    plt.figure()
    plt.plot(bar.lenVector, bar.temperature[0], label=str(0))
    # Outer loop through all the time series.
    for i in range(len(bar.temperature)-1):
        q = 0
        if i == 1:
            q = qIn
            
        # Setting the temperature on the extremes.
        bar.temperature[i+1][0] = (1-(2*bar.rVal))*bar.temperature[i][0] + 2*bar.rVal*(bar.temperature[i][1] + q)
        bar.temperature[i+1][-1] = (1-(2*bar.rVal))*bar.temperature[i][-1] + 2*bar.rVal*(bar.temperature[i][-2] + q)
        
        # Calculating the internal products.
        HCC.condEquation(model, material, bar, i, T0)
        
        plotVec = [1, int(len(bar.timeVector)/4), int(len(bar.timeVector)/2), int(3*len(bar.timeVector)/4), len(bar.timeVector)-1]
        if i in plotVec:
            plt.plot(bar.lenVector, bar.temperature[i], label=str(bar.delta*i))


    plt.xlabel('Length'); plt.ylabel('Temperature')
    plt.legend()
    plt.axis("tight")
    plt.title("K regression - Second approach on finite differences")
    plt.savefig("conduction1DKreg.pdf")
    
    print(bar.temperature)
    
    return bar.temperature