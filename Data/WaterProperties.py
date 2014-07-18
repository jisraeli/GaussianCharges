#!/Users/jisraeli/Library/Enthought/Canopy_64bit/User/bin/python

import numpy as np
from math import*
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import sys, getopt

def getProperty(datalines, lineIndex):
    properties = [['Density'],
                 ['Enthalpy', 'of', 'Vaporization'],
                 ['Isothermal', 'Compressibility'],
                 ['Isobaric', 'Heat', 'Capacity'],
                 ['Thermal', 'Expansion', 'Coefficient'],
                 ['Dielectric', 'Constant']]
    for i in range(4):
        for elem in properties:
            N = len(elem)
            if [True]*N == [name in datalines[lineIndex-i] for name in elem]:
                elemName = ' '.join(elem)
                return elemName


def extractTestData(lines):
    datalines = []
    for line in lines:
        dataline = line.split()
        datalines.append(dataline)

    initTemp = '249.15'
    finalTemp = '373.15'
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
    '''
    return a list called iterList where iterList[i] has the file lines
    for iteration i in the FB optimize output
    '''

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
    '''
    extractData from FB output
    for Test output:
        return data dictionary where data["property name"] contains
        parsed FB output for that property

    for Optimize output:
        return a dictionary where data["property name"] contains
        a propertData dictionary where propertyDictionary["k"]
        contains parsed FB output for that property on iteration k
    '''

    if Test and not Optimize:
        lines = open(InFile).readlines()
        return extractTestData(lines)
    if not Test and Optimize:
        iterLines = SplitIterations(InFile)
        if OnlyLast:
            return extractTestData(iterLines.pop())
        else:
            dataByIteration = []
            for lines in iterLines:
                dataByIteration.append(extractTestData(lines))
            dataByProperty = {}
            for prop in dataByIteration[0].keys():
                tempDict = {}
                for i in range(len(dataByIteration)):
                    tempDict[str(i)] = dataByIteration[i][prop]
                dataByProperty[prop] = tempDict
            return dataByProperty


def plotFigure(pp, data, Test=True, Optimize=False, Prop=None,
               CompareData=None):
    if Test and not Optimize:
        for i in range(len(data.keys()[:6])):
                prop = data.keys()[:6][i]
                Data = np.asarray(data[prop])
                tempVals = Data[:, 0]
                experimentalVals = Data[:, 2]
                calcVals = Data[:, 3]
                calcDevs = Data[:, 4]
                plt.subplot(1, 1, 1)
                plt.xlabel('Temperaure')
                if i % 2 == 0:
                    plt.ylabel(prop)
                    if Optimize:
                        plt.title('Iteration '+str(j)+' Results:')
                    else:
                        plt.title('Test Results:')
                elif i % 2 == 1:
                    plt.title(prop)
                lineExp = plt.errorbar(tempVals, experimentalVals, color='red')
                line = plt.errorbar(tempVals, calcVals, color='blue', yerr=calcDevs)
                if i % 1 == 0:
                    plt.legend([lineExp, line], ['Experiment', 'Calculated'], loc='best')
                    pp.savefig()
                    plt.clf()
    if Optimize and not Test:
        print "plotting ", Prop, "..."
        plt.subplot(1, 1, 1)
        N = len(data.keys())
        colorInds = np.linspace(0.0, 1.0, num=N+1)[::-1]
        i = N - 3
        if i<0:
            i = 0
        while i < N-1:
            Data = np.asarray(data[str(i)])
            tempVals = Data[:, 0]
            calcVals = Data[:, 3]
            col = colorInds[i+1]
            plt.errorbar(tempVals, calcVals, color=str(col))
            i += 1
        Data = np.asarray(data[str(i)])
        tempVals = Data[:, 0]
        experimentalVals = Data[:, 2]
        calcVals = Data[:, 3]
        calcDevs = Data[:, 4]
        lineExp = plt.errorbar(tempVals, experimentalVals, color='red', marker='D',lw='2')
        line_tip3g = plt.errorbar(tempVals, calcVals, color='0.0', lw='2',
                                  yerr=calcDevs, marker='s')
            
        NameList = ['Model', 'Experiment']
        LineList = [line_tip3g, lineExp]
        prop = Prop
        colorList = ['green', 'blue', 'orange']
        markerList = ['v', '*', '^']
        i = 0
        for model in CompareData.keys():
            modelData = CompareData[model]
            propData = np.asarray(modelData[prop])
            propVals = propData[:, 3]
            propDevs = propData[:, 4]
            tempLine = plt.errorbar(tempVals, propVals, color=colorList[i], lw='5',
                                  yerr=propDevs, ls=':', marker=markerList[i])
            NameList.append(model)
            LineList.append(tempLine)
            i += 1
        print "NameList: ", NameList
        plt.xlabel('Temperaure')
        plt.title(prop)
        plt.legend(LineList, NameList, loc='best')



        
        '''
        tip3pData = np.asarray(tip3p[prop])
        tip3pVals = tip3pData[:, 3]
        tip3pDevs = tip3pData[:, 4]
        tip3pfbData = np.asarray(tip3pfb[prop])
        tip3pfbVals = tip3pfbData[:, 3]
        tip3pfbDevs = tip3pfbData[:, 4]
        
        line_tip3p = plt.errorbar(tempVals, tip3pVals, color='green', lw='5',
                                  yerr=tip3pDevs, ls=':', marker='v')
        line_tip3pfb = plt.errorbar(tempVals, tip3pfbVals, color='blue', marker='*',
                                    lw='2', yerr=tip3pfbDevs, ls=':')
        plt.xlabel('Temperaure')
        plt.title(prop)
        plt.legend([lineExp, line_tip3g, line_tip3p, line_tip3pfb],
                   ['Experiment', 'TIP3G', 'TIP3P', 'TIP3P-FB'], loc='best')
        '''
        if 'Density' in prop:
            plt.ylim(960, 1010)
            if 'TIP4P-FB' in CompareData.keys():
                plt.ylim(960, 1002)
        elif 'Isothermal' in prop:
            plt.ylim(30, 70)
        elif 'Dielectric' in prop:
            ymin, ymax = plt.ylim()
            plt.ylim(ymin, 120)
        pp.savefig()
        plt.clf()


def plotProperties(dataList, fileName, Test=True, Optimize=False,
                   OnlyLast=False, compare=None):
    if Test and not Optimize:
        data = dataList
        pp = PdfPages(fileName)
        plotFigure(pp, data)
        pp.close()
    elif not Test and Optimize:
        print("creating pdf...")
        pp = PdfPages(fileName)
        for prop in dataList.keys():
            data = dataList[prop]
            plotFigure(pp, data, Test, Optimize, prop, CompareData=compare)
        print ("Saving pdf...")
        pp.close()


def generateAnalysis(InFile):
    if 'Test' in InFile and not 'Optimize' in InFile:
        test = True
        data = extractData(InFile, Test=test)
        splitName = InFile.split('Test')
        fileName = splitName[0]+'_TestReport.pdf'
        plotProperties(data, fileName, Test=test)
    if not 'Test' in InFile and 'Optimize' in InFile:
        test = False
        optimize = True
        data = extractData(InFile, Test=test, Optimize=True)
        splitName = InFile.split('Optimize')
        fileName = splitName[0]+'_OptimizationReport.pdf'
        CompareData = {}
        IncludeTIP4Pfb = raw_input("Include TIP4P-FB in the report (y/n)? ")
        if IncludeTIP4Pfb == 'y':
            tip4pfbData = extractData('tip4p-fb.out')
            CompareData['TIP4P-FB'] = tip4pfbData
        IncludeTIP3P = raw_input("Include TIP3P in the report (y/n)? ")
        if IncludeTIP3P == 'y':
            tip3pData = extractData('tip3p.txt')
            CompareData['TIP3P'] = tip3pData
        IncludeTIP3Pfb = raw_input("Include TIP3P-FB in the report (y/n)? ")
        if IncludeTIP3Pfb == 'y':
            tip3pfbData = extractData('tip3p-fb.txt')
            CompareData['TIP3P-FB'] = tip3pfbData
        IncludeiAmoeba = raw_input("Include iAmoeba in the report (y/n)? ")
        if IncludeiAmoeba == 'y':
            iAmoebaData = extractData('iamoeba.out')
            CompareData['iAmoeba'] = iAmoebaData
        plotProperties(data, fileName, Test=test, Optimize=optimize,
                       compare=CompareData)


def main(argv):
    inputfile = ''
    try:
        opts, args = getopt.getopt(argv, "hi:o:", ["ifile="])
    except getopt.GetoptError:
        print 'test.py -i <inputfile>'
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print 'test.py -i <inputfile>'
            sys.exit()
        elif opt in ("-i", "--ifile"):
            inputfile = arg
    print 'Input file is "', inputfile
    generateAnalysis(inputfile)
print("Report is Ready!")


if __name__ == "__main__":
    main(sys.argv[1:])

#InFile = '6Param_LiquidOptimize.out'
#generateAnalysis(InFile)
#print("Report is Ready!")
