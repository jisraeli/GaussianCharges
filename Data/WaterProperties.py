import numpy as np
from math import*
from matplotlib import pyplot as plt


def getProperty(datalines, lineIndex):
    properties = ['Density', 'Enthalpy', 'Isothermal', 'Isobaric', 'Thermal', 'Dielectric']
    for i in range(4):
        for elem in properties:
            if elem in datalines[lineIndex-i]:
                print elem
                return elem



#read data into array

#lines = open('simple_liquidTest.out').readlines()
lines = open('tip3g_liquidTest.out').readlines()
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
data = {}
datablock = []
count = 0
for i in range(len(datalines)):
    if not reachedInit and initTemp in datalines[i]:
        reachedInit = True
        passedFinal = False
        prop = getProperty(datalines, i)
    elif reachedInit and finalTemp in datalines[i]:
        reachedFinal = True
        passedFinal = False
    elif reachedFinal and not finalTemp in datalines[i]:
        passedFinal = True
        reachedInit = False
        reachedFinal = False
        data[prop] = datablock
        datablock = []

    if reachedInit and not passedFinal:
        datarow = []
        for elem in datalines[i]:
            if not elem in ['atm', '+-']:
                datarow.append(float(elem))
        datablock.append(datarow)


#for prop in data.keys()[:1]:
prop = 'Density'
Data = np.asarray(data[prop])
tempVals = Data[:, 0]
experimentalVals = Data[:, 2]
calcVals = Data[:, 3]
calcDevs = Data[:, 4]
plt.plot(tempVals, experimentalVals, color='red')
plt.plot(tempVals, calcVals, color='blue')
plt.show()
