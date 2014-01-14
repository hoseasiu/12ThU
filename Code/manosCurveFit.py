''' manosCurveFit.py - Hosea Siu 2014
    Takes light curve data as input and fits it to a function accoring to Harris et al. 1989
'''

from lmfit import minimize, Parameters, Parameter, report_fit
from operator import itemgetter
import matplotlib.pyplot as plt, numpy as np, os.path, sys

''' Class used to store and manipulate MANOS light curve data
'''
class lightCurveData:
    ''' Gets the txt data files and converts them into a lightCurveData object
    '''
    def getData(self, fileName, formatSpec, offsetsList, passNumber):
        basepath = os.path.abspath(os.path.dirname(sys.argv[0]))
        filepath = os.path.abspath(os.path.join(basepath, "..", 'Data', fileName))

        textFile = np.loadtxt(filepath,'a12')          # extract data from text file

        # separate data into appropriately-named arrays
        # TODO- currently assumes that everything is a float- generalize later (needed?)
        tempData = {}

        for i in range(len(formatSpec)):
            tempData[formatSpec[i][0]] = textFile[:,formatSpec[i][1]].astype(float)

        if offsetsList is not None:     # if offset values exist, then perform an offset
            tempData = self.offsetMags(tempData, offsetsList[passNumber])
            
        for i in range(len(formatSpec)):
            
            if passNumber == 0:
                self.data[formatSpec[i][0]] = tempData[formatSpec[i][0]]
            else:
                self.data[formatSpec[i][0]] = np.append(self.data[formatSpec[i][0]],tempData[formatSpec[i][0]])

    ''' Given a data type (i.e. 'jd', 'diffMag', etc), sorts the light curve data by that data type
        For fitting plots, we always sort by jd, though it may be switched to check other patterns
    '''
    def sortByDataType(self, key):
        if key is None:
            return
        
        dict = self.data
        keyIndex = dict.keys().index(key)
        numKeys = len(dict.keys())
        numEl = -1                  # checking to make sure that all dict entries have the same number of elements
        for i in range(numKeys):
            if numEl == -1:
                numEl = len(dict[dict.keys()[i]])
            else:
                if numEl != len(dict[dict.keys()[i]]):
                    print 'Error: invalid dictionary- nonrectangular'
                    return

        keyList = dict.keys()
        dictArray = np.zeros((numKeys,numEl))          # form a np array to use the transpose function
        for i in range(numKeys):
            dictArray[i] = dict[dict.keys()[i]]

        sortedList = sorted(dictArray.T.tolist(), key = itemgetter(keyIndex))

        temp = np.array(sortedList).T   # temporary variable for transpose
        for i in range(numKeys):
            dict[keyList[i]] = temp[i]
        self.data = dict

    ''' Given a dict of nights and corresponding magnitude offsets, offset the magnitude data by that much
        Also center the data by a weighted average
        "dataFile" refers to a data set from a specific file
    ''' # TODOLATER - figure out a convention for offseting "night," since observatories getting data at the same night will still have offsets
    def offsetMags(self, dataFile, offsets = None):
        
        if offsets is None:
            return
        if len(np.unique(dataFile['night'])) == len(offsets):  # check if each night has an associated offset value
            for i in range(len(dataFile['diffMag'])):
                finishedOffset = False
                for j in range(len(offsets)):
                    if finishedOffset:
                        break
                    elif dataFile['night'][i] == offsets.keys()[j]:
                        dataFile['diffMag'][i] += offsets[offsets.keys()[j]]
            return dataFile
        else:
            print 'Error: number of nights does not match the number of offsets'
            return
            
    ''' lightCurveData initializer
        requires an object name, file names; can optionally take offsets and sorting
        (sorting should always be by jd for light curve plots)
    '''
    def __init__(self, objectName, fileNamesAndFormat, offsetsList = None):
        if fileNamesAndFormat is None:
            print 'Error: file name(s) and format(s) not specified'
            return

        # check for errors in the format specification
        for i in range(len(fileNamesAndFormat)):
            fileName = fileNamesAndFormat[i].keys()[0]
            formatSpec = fileNamesAndFormat[i][fileName]

            if len(formatSpec) <= 0:
                print 'Error: formatSpec input incorrect'
                return
            for j in range(len(formatSpec)):
                if len(formatSpec[j]) !=2:
                    print 'Error: formatSpec input incorrect'
                    return
        
        self.name = objectName
        self.data = {}      # initialize an empty dictionary with the keys in formatSpec
        for i in range(len(formatSpec)):
            self.data[formatSpec[i][0]] = {}

        # check for missing (required)specifications
        if 'jd' not in self.data.keys() or 'diffMag' not in self.data.keys() or 'magErr' not in self.data.keys():
            print 'Error: Insufficient data for light curve'
            print '       Did you include \'jd\', \'diffMag\', and \'magErr\'?'
            print '       Right now, I have ' + str(self.data.keys())
        else:
            for j in range(len(fileNamesAndFormat)):
                # TODO - this structure is mess... find a better way to do it
                fileName = fileNamesAndFormat[j].keys()[0]
                formatSpec = fileNamesAndFormat[j][fileName]
                self.getData(fileName, formatSpec, offsetsList, j)
            self.sortByDataType('jd')
            self.data['diffMag'] -= np.average(self.data['diffMag'], weights = 1.0/self.data['magErr'])        # center the data around zero with a weighted average

''' Creates range lists of floats
    Has some error handling, but still has machine precision issues (shouldn't be an actual problem)
'''
def floatRange(start, stop, stepSize=0):    
    if start == stop:
        return []
    elif stepSize == 0:
        return [start, stop]
    elif (start < stop and stepSize > 0) or (start > stop and stepSize < 0):
        numbers = []
        start = float(start)
        stop = float(stop)
        current = start
        while True:
            numbers.append(current)
            if (start < stop):                  # allows for forward and backwards listing
                current = current + stepSize
                if current > stop:
                    break
            else:
                current = current - stepSize
                if current < stop:
                    break
        return numbers
    else:
        print 'Error: floatRange input error'
        return []

''' Generate the light curve model using Parameters
'''
def makeModel(params, t, mag=None, err=None):
    P = params['P'].value
    H = np.zeros(len(t))  
    for L in range(1,len(params.keys())/2+1):
        trigTerm = 2.0*np.pi*float(L)/float(P)
        offset = t-t[0]
        A = params['A'+str(L)].value
        B = params['B'+str(L)].value
        H = H + A*np.sin(trigTerm*offset)+B*np.cos(trigTerm*offset)
    if mag is None and err is None:     # no data given- the model values are returned
        return H
    elif err is None:                   # magnitude data is given, but no errors- the chi squared is returned
        return H-mag
    else:                               # both magnitude and error are given- the reduced chi squared is returned
        return (H-mag)/err

''' Fit the data to a model using scipy.optimize.leastsq.
    Takes the dataset(a lightCurveData object), the range of orders to check (a two-element list)
    Three period guess cases:
        - no guess (Nyquist sampling criterion used)        # TODO - figure out step size
        - range of periods provided (three-element list of [minGuess, maxGuess, step])
        - single initial guess provided (three-element list of [min, max, guess]) where minPeriod and maxPeriod are hard limits
    The range of orders to check over may also be provided
'''
# TODO - change initPeriod to take in a specific number, a specific range, or a default range
def fitData(lightCurveData, method = None, periodGuess = None, minPeriod = None, maxPeriod = None, orderRange = [2,6]):
    if method is None:         # Nyquist sampling criterion is used
         pass
        # TODO - Nyquist sampling
    if method == 'single':
        if type(periodGuess) is not float and type(periodGuess) is not int:
            print 'Error: \'single\' period guess method requires a guess input of type int or float'
            return
        else:
            initPeriodList = [float(periodGuess)]
    elif method == 'range':
        if type(periodGuess) is not list:
            print 'Error: \'range\' period guess method requires a guess input of type [min, max, stepSize] in floats'
            return
        elif len(periodGuess) != 3:
            print 'Error: \'range\' period guess method requires a guess input of type [min, max, stepSize] in floats'
            return            
        else:
            initPeriodList = floatRange(periodGuess[0], periodGuess[1], periodGuess[2])
    else:
        print 'Error: invalid period guess method'
        return

    initPeriodList = np.array(initPeriodList)/24.0      # convert periods from hours to JD
    
    # information from light curve for plotting
    time = lightCurveData.data['jd']-lightCurveData.data['jd'][0]       # use the minimum JD as a reference point of 0
    mag = lightCurveData.data['diffMag']
    err = lightCurveData.data['magErr']

    # storing the best fitting result
    [bestResult, bestParams, bestVar, bestOrder, bestPeriod] = [None, None, None, None, None]

    # store the periods and errors
    periodsTested = []
    periodErrors = []
    
    for p in range(len(initPeriodList)):
        initPeriod = initPeriodList[p]
        print 'attempting P = ' + str(initPeriod*24.0) + ' hours'
        for m in range(orderRange[0],orderRange[1]+1):

            params = Parameters()   # create the Parameters object and add P, A, and B values to it
            params.add('P', value = initPeriod, min = minPeriod, max = maxPeriod)
            for i in range(m):
                params.add('A' + str(i+1), value = 0)
                params.add('B' + str(i+1), value = 0)
            
            fitResult = minimize(makeModel, params, args = (time, mag, err))
            periodsTested.append(params['P'].value*24.0)       # add the period to a list for the error plot
            periodErrors.append(np.mean(np.sqrt(fitResult.residual**2)))

            if fitResult.success:
                n = float(len(time))       # number of observations
                k = float(2*m+1)           # total free parameters
                var = 1/(n-k)*sum(makeModel(params, time, mag, err)**2)
                
                if bestResult is None:
                    [bestResult, bestParams, bestVar, bestOrder, bestPeriod] = [fitResult, params, var, m, initPeriod]
                elif bestVar > var:
                    [bestResult, bestParams, bestVar, bestOrder, bestPeriod] = [fitResult, params, var, m, initPeriod]
            else:
                print 'optimization failed for P = ' + str(initPeriod) + ', m = ' + str(m)
        
    return (bestResult, bestOrder, periodsTested, periodErrors)

''' Phase folds the magnitude data given a set that has the phase offset at the first data point
'''
def phaseFold(time, period):
    for i in range(len(time)):
        periodOffset = int(time[i]/period)
        if periodOffset >= 1:
            time[i] -= period*periodOffset

    print 'fitted period: ' + str(period*24.0) + ' h'
    return time


''' Takes the best fit output and the data and plots them
    plotting error bars is optional (default is True)
    plotting the full phase is optional (default is True)
'''
def plotAndPrintResults(fit, m, lightCurveData, printReport = True, plotFullPeriod = True, plotErrors = True, phaseFoldData = True, plotResiduals = False, periodErrors = None):
    if printReport:
        report_fit(fit.params)
    time = lightCurveData.data['jd']-lightCurveData.data['jd'][0]
    mag = lightCurveData.data['diffMag']
    err = lightCurveData.data['magErr']

    fittedPeriod = fit.params['P'].value

    if phaseFoldData:
        time = phaseFold(time, fittedPeriod)
    
    if plotFullPeriod and time.min()+fittedPeriod < time.max():
        modelTime = np.linspace(time.min(), time.max(),100)
    elif plotFullPeriod:
        modelTime = np.linspace(time.min(), time.min()+fittedPeriod,100)
    else:
        modelTime = np.linspace(time.min(), time.max(),100)

    plt.figure(1)
    if plotResiduals:
        plt.subplot(211)
    modelMags = makeModel(bestFit.params, modelTime)
    plt.plot(modelTime, modelMags, 'b-')     # plot the model fit
    amp = np.ptp(modelMags)
    print 'amplitude = ' + str(amp)
    if (plotErrors):                                                    # plot the data with error bars (default)
        plt.errorbar(time, mag, err, fmt = 'rx')
    else:                                                               # plot the data without error plots
        plt.plot(time, mag, 'rx')

    if phaseFoldData:
        plt.xlabel('\delta JD')
    else:
        plt.xlabel('JD + ' + str(int(min(lightCurveData.data['jd']))))    # the subtraction simplifies display of the JD
    plt.ylabel('Differential Magnitude')
    plt.gca().invert_yaxis()        # flip the y axis
    plt.title(objectName + ' Light Curve, ' + ' P = ' + str(fit.params['P'].value*24.0) + 'h, ' + 'a = ' + str(amp) + ', m = ' + str(m))

    lightCurveAxis = plt.axis()     # used to make sure that the residuals plots the x limits the same way
    
    if plotResiduals:
        plt.subplot(212)
        plt.plot(time, fit.residual, 'rx')
        plt.title('Residuals for ' + objectName + ' Fit')
        plt.ylabel('Residual Magnitude')
        plt.xlim(lightCurveAxis[0:2])
        if phaseFoldData:
            plt.xlabel('\delta JD')
        else:
            plt.xlabel('JD + ' + str(int(min(lightCurveData.data['jd']))))

    plt.subplots_adjust(hspace = 0.5)
    print 'mean RMS of residuals = ' + str(np.mean(np.sqrt(fit.residual**2)))

    if periodErrors is not None:
        if type(periodErrors) is not list:
            print 'Error: period errors must be given as a list'
        else:
            if len(periodErrors) != 2:
                print 'Error: invalid period errors format'
            else:
                plt.figure(2)
                plt.plot(periodErrors[0], periodErrors[1], 'x')
                plt.title('Mean RMS of Residuals')
                plt.ylabel('mean RMS')
                plt.xlabel('period (h)')
        
    plt.show()

'''         start of the run routine        '''

##objectName = 'Spartacus20090130'
##fileName = 'Spartacus20090130_MANOS.txt'
##formatSpec = [['jd',3],['diffMag',6],['magErr',7]]
##fileNamesAndFormat = [{fileName:formatSpec}]
##offsets = None
##T0 = [3, 4, 0.1]                                # initial guess at period, in hours

##objectName = 'da14'
##fileName = 'da14.120301.568.mag.txt'
##formatSpec = [['jd',0],['diffMag',1],['magErr',2]]
##fileNamesAndFormat = [{fileName:formatSpec}]
##offsets = None
##T0 = [4, 6, 0.25]

objectName = 'Elisa'
mine = {'elisa_mine_mod.txt':[['night',0], ['jd',2],['diffMag',5],['magErr',6]]}
his = {'elisa_his.txt':[['night',0], ['jd',2],['diffMag',4],['magErr',5]]}
fileNamesAndFormat = [mine,his]
offsets = [{1:0.0,2:-0.04,3:0.464},{1:-0.324,2:-0.257,3:-0.237,4:-0.194,5:-0.223,6:-0.321,7:-0.246,8:-0.372,9:-0.15}]
T0 = [14, 18, 0.25]

# note that these functions are all overloaded (there are extra options you can set- see the function definitions for examples
lcd = lightCurveData(objectName, fileNamesAndFormat, offsets)         # read in data from the text file and create a lightCurveData object
(bestFit, m, periodsTested, periodErrors) = fitData(lcd, method = 'range', periodGuess = T0)       # fit the data to a model (can also add min and max periods in JD)
plotAndPrintResults(bestFit, m, lcd, plotResiduals = True, periodErrors = [periodsTested, periodErrors])             # plot and print the results
