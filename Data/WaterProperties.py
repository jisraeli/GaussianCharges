import numpy as np
from math import*
from matplotlib import pyplot as plt


def getProperty(datalines, lineIndex):
    properties = [['Liquid', 'Density'],
                ['Liquid', 'Enthalpy', 'of', 'Vaporization'],
                ['Liquid', 'Isothermal', 'Compressibility'],
                ['Liquid', 'Isobaric', 'Heat', 'Capacity'],
                ['Liquid', 'Thermal', 'Expansion', 'Coefficient'],
                ['Liquid', 'Dielectric', 'Constant']]
    for i in range(4):
        for elem in properties:
            N = len(elem)
            if [True]*N == [name in datalines[lineIndex-i] for name in elem]:
                elemName = ' '.join(elem)
                print elemName
                return elemName

def extractData(InFile):
    lines = open(InFile).readlines()
    datalines = []
    for line in lines:
        dataline = line.split()
        datalines.append(dataline)

    initTemp = '249.15'
    finalTemp = '298.15'
    reachedInit = False
    reachedFinal = False
    passedFinal = False
    data = {}
    datablock = []
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
    return data

#for prop in data.keys()[:1]:
data = extractData('simple_liquidTest.out')
dataFull = extractData('tip3g_liquidTest.out')
data['legend'] = '3Params'
dataFull['legend'] = '6Params'
figures = []
for i in range(len(data.keys()[:6])):
    prop = data.keys()[:6][i]
    Data = np.asarray(data[prop])
    DataFull = np.asarray(dataFull[prop])
    tempVals = Data[:, 0]
    experimentalVals = Data[:, 2]
    calcVals = Data[:, 3]
    calcDevs = Data[:, 4]
    calcFullVals = DataFull[:, 3]
    calcFullDevs = DataFull[:, 4]
    #fig = plt.figure()
    figure = plt.subplot(2,3,i+1)
    plt.title(prop)
    lineExp = plt.plot(tempVals, experimentalVals, color='red')
    line = plt.errorbar(tempVals, calcVals, color='green', yerr=calcDevs)
    lineFull = plt.errorbar(tempVals, calcFullVals, color='blue', yerr=calcFullDevs)
    #plt.legend([lineExp, line, lineFull], ['Experiment', data['legend'], dataFull['legend']])
    if i>2:
        plt.xlabel('Temperaure')
    if i==3:
        plt.legend([lineExp, line, lineFull], ['Experiment', data['legend'], dataFull['legend']], loc='best')
plt.show()
