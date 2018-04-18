# -*- coding: utf-8 -*-
"""
Created on Thu Jan 18 01:15:20 2018

Script with numerical schemes for heat conduction at inner nodes.

@author: Leandro
"""


"""
Import statements:
""" 

import matplotlib.pyplot as plt
import numpy as np
import HCClass


"""
Models for heat conduction.
"""

# Inputs: model of heat conduction, material model, bar and time objects,
# matrix position (i,j), reference temperature.
def condEquation(model, material, bar, i, T0):
    
    # Simple heat conduction.
    if model == "kSimple":
        for j in range(1, len(bar.temperature[0])-1):
            bar.temperature[i+1][j] = bar.rVal*(bar.temperature[i][j-1] + bar.temperature[i][j+1]) + (1-(2*bar.rVal))*bar.temperature[i][j]
        return 0
        
    # Finite difference for d(KpromT').
    elif model == "kProm":
        for j in range(1, len(bar.temperature[0])-1):
            kprom_r = (material.eqn(bar.temperature[i][j+1]-T0) + material.eqn(bar.temperature[i][j]-T0))/2.0
            kprom_l = (material.eqn(bar.temperature[i][j-1]-T0) + material.eqn(bar.temperature[i][j]-T0))/2.0
            temp_der_r = (bar.temperature[i][j+1] - bar.temperature[i][j])/bar.increment
            temp_der_l = (bar.temperature[i][j] - bar.temperature[i][j-1])/bar.increment
            inner = kprom_r*temp_der_r - kprom_l*temp_der_l
            bar.temperature[i+1][j] = bar.temperature[i][j] + (bar.delta/(bar.increment))*(inner)
        return 0
    
    # Crank-Nicolson method.
    elif model == "kCN":
        
        innerTemp = np.zeros([1, len(bar.temperature[0])])
        innerTemp[0] = np.copy(bar.temperature[i])
        innerTemp[0][0] = bar.temperature[i+1][0]
        innerTemp[0][-1] = bar.temperature[i+1][-1]
        
        errorL = 10
        # Counter of implicit loop.
        l = 0
        while (errorL > 0.0000001 and l < 20):
            # Counter increased for iteration.
            l += 1
            matrixK = np.zeros([len(bar.temperature[0]), len(bar.temperature[0])])
            matrixC = np.zeros(len(bar.temperature[0]))
            # Loop over length.
            for j in range(0, len(bar.temperature[0])):
                # Calculate k values for each bar point.
                
                # Default cases for first and last element.
                if j == 0:
                    matrixK[j][j] = 1
                    matrixC[j] = innerTemp[0][0]
                elif j == len(bar.temperature[0]) - 1:
                    matrixK[j][j] = 1
                    matrixC[j] = innerTemp[0][-1]
                else:
                    # Weighted left and right temperatures.
                    tempLeft = 0.25*(innerTemp[l-1][j-1] + innerTemp[l-1][j] + innerTemp[0][j-1] + innerTemp[0][j])
                    tempRight = 0.25*(innerTemp[l-1][j] + innerTemp[l-1][j+1] + innerTemp[0][j] + innerTemp[0][j+1])
                    kValueL = material.eqn(tempLeft)
                    kValueR = material.eqn(tempRight)
                    
                    # Matrix values.
                    matrixK[j][j-1] = -kValueL/(2.0*bar.increment)
                    matrixK[j][j] = 1.0/bar.delta + (kValueL + kValueR)/(2.0*bar.increment)
                    matrixK[j][j+1] = -kValueR/(2.0*bar.increment)
                    
                    # Constant vector values.
                    valueCLeft = innerTemp[0][j-1]*kValueL/(2.0*bar.increment)
                    valueCRight = innerTemp[0][j+1]*kValueR/(2.0*bar.increment)
                    valueCMid = innerTemp[0][j]*(1/bar.delta - (kValueL+kValueR)/(2.0*bar.increment))
                    
                    matrixC[j] = valueCMid + valueCLeft + valueCRight
                    
            # Solving the linear system KT + C = 0 for T vector.
            innerTempNew = np.linalg.solve(matrixK, matrixC)
            innerTemp = np.vstack([innerTemp, innerTempNew])
            
            # Find the total error
            errorL = 0.0
            errorSum = 0.0
            p = 2.0
            for lIter in range(0, len(innerTempNew)):
                errorSum += np.power((innerTemp[l][lIter]-innerTemp[l-1][lIter]),p)
            errorL = np.power(errorSum/len(innerTempNew),(1.0/p))
            
        bar.temperature[i+1] = innerTemp[-1]
        print('Iterations = ' + str(l))
        return 0

