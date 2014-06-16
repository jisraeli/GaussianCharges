import numpy as np
from math import*


def getProperty(datalines, lineIndex):
    properties = ['Density', 'Enthalpy', 'Isothermal', 'Isobaric', 'Thermal', 'Dielectric']
    for i in range(4):
        for elem in properties:
            if elem in datalines[lineIndex-i]:
                return elem



#read data into array
lines = open('simple_liquidTest.out').readlines()
datalines = []
for line in lines:
    dataline = line.split()
    datalines.append(dataline)

#find values and put in np arrays

initTemp = '249.15'
finalTemp = '298.15'
reachedInit = False
reachedFinal = False
passedFinal = False
data = []
count = 0
for i in range(len(datalines)):
    if not reachedInit and initTemp in datalines[i]:
        reachedInit = True
        passedFinal = False
    elif reachedInit and finalTemp in datalines[i]:
        reachedFinal = True
        passedFinal = False
    elif reachedFinal and not finalTemp in datalines[i]:
        passedFinal = True
        reachedInit = False
        reachedFinal = False

    if reachedInit and not passedFinal:
        datarow = []
        for elem in dataline:
            datarow.append(float(elem))
        datablock.


