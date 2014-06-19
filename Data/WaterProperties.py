import numpy as np
from math import*
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

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

def extractTestData(lines):
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


def SplitIterations(InFile):
    lines = open(InFile).readlines()
    badIter = ['Objective', 'function', 'rises!']
    iterLines = []
    iterList = []
    i = 1
    for line in lines:
        dataline = line.split()
        if 'Iteration' in dataline and str(i)+':' in dataline:
            print "saved iteration ", i
            iterList.append(iterLines)
            iterLines = []
            i += 1
        elif [True]*len(badIter) == [word in dataline for word in badIter]:
            print "skipping bad iteration!"
            iterLines = []
            i += 1
        iterLines.append(line)
    return iterList


def extractData(InFile, Test=True, Optimize=False, OnlyLast=False):
    if Test and not Optimize:
        lines = open(InFile).readlines()
        return extractTestData(lines)
    if not Test and Optimize:
        iterLines = SplitIterations(InFile)
        if OnlyLast:
            return extractTestData(iterLines.pop())
        else:
            data = []
            for lines in iterLines:
                data.append(extractTestData(lines))
            return data


def plotProperties(dataList, fileName, Test=True, Optimize=False, OnlyLast=False):
    if Test and not Optimize:
        data = dataList
        pp = PdfPages(fileName)
        for i in range(len(data.keys()[:6])):
            prop = data.keys()[:6][i]
            Data = np.asarray(data[prop])
            tempVals = Data[:, 0]
            experimentalVals = Data[:, 2]
            calcVals = Data[:, 3]
            calcDevs = Data[:, 4]
            figure = plt.subplot(1, 2, (i % 2) + 1)
            plt.xlabel('Temperaure')
            if i%2==0:
                plt.ylabel(prop)
                plt.title('Test Results:')
            elif i%2==1:
                plt.title(prop)
            lineExp = plt.errorbar(tempVals, experimentalVals, color='red')
            line = plt.errorbar(tempVals, calcVals, color='blue', yerr=calcDevs)
            if i%2==0:
                plt.legend([lineExp, line], ['Experiment', 'Calculated'], loc='best')
                pp.savefig()
                plt.clf()
        pp.close()
    elif not Test and Optimize:
        pp = PdfPages(fileName)
        for j in range(len(dataList)):
            data = dataList[j]
            for i in range(len(data.keys()[:6])):
                prop = data.keys()[:6][i]
                Data = np.asarray(data[prop])
                tempVals = Data[:, 0]
                experimentalVals = Data[:, 2]
                calcVals = Data[:, 3]
                calcDevs = Data[:, 4]
                plt.subplot(1, 2, (i % 2) + 1)
                plt.xlabel('Temperaure')
                if i % 2 == 0:
                    plt.ylabel(prop)
                    plt.title('Iteration '+str(j)+' Results:')
                elif i % 2 == 1:
                    plt.title(prop)
                lineExp = plt.errorbar(tempVals, experimentalVals, color='red')
                line = plt.errorbar(tempVals, calcVals, color='blue', yerr=calcDevs)
                if i % 2 == 0:
                    plt.legend([lineExp, line], ['Experiment', 'Calculated'], loc='best')
                    pp.savefig()
                    plt.clf()
        pp.close()


def generateAnalysis(InFile):
    if 'Test' in InFile and not 'Optimize' in InFile:
        test = True
        data = extractData(InFile, Test)
        splitName = InFile.split('Test')
        fileName = splitName[0]+'_Report.pdf'
        plotProperties(data, fileName, Test=test)
    if not 'Test' in InFile and 'Optimize' in InFile:
        test = False
        optimize = True
        data = extractData(InFile, Test=test, Optimize=True)
        splitName = InFile.split('Optimize')
        fileName = splitName[0]+'_Report.pdf'
        plotProperties(data, fileName, Test=test, Optimize=optimize)

InFile = '6ParamOptimize.out'
generateAnalysis(InFile)

