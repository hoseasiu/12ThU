''' manosCurveFit.py - Hosea Siu 2014
    Takes light curve data as input and fits it to a function accoring to Harris et al. 1989
'''

from lmfit import minimize, Parameters, Parameter, report_fit
from operator import itemgetter
import matplotlib.pyplot as plt, os.path, sys, numpy as np

''' Class used to store and manipulate MANOS light curve data
'''
class lightCurveData:
    ''' Gets the txt data files and converts them into a lightCurveData object
    '''
    def getData(self, filename, formatSpec):
        basepath = os.path.abspath(os.path.dirname(sys.argv[0]))
        filepath = os.path.abspath(os.path.join(basepath, "..", 'Data', fileName))

        textFile = np.loadtxt(filepath,'a12')          # extract data from text file

        # separate data into appropriately-named arrays
        # TODO- currently assumes that everything is a float- generalize later (needed?)
        for i in range(len(formatSpec)):
            self.data[formatSpec[i][0]] = textFile[:,formatSpec[i][1]].astype(float)

        self.data['diffMag'] -= np.average(self.data['diffMag'], weights = 1.0/self.data['magErr'])        # center the data around zero with a weighted average

    ''' Given a data type (i.e. 'jd', 'diffMag', etc), sorts the light curve data by that data type
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
    ''' # TODO - figure out a convention for offseting "night," since observatories getting data at the same night will still have offsets
    def offsetMags(self, offsets = None):
        if offsets is None:
            return
        if len(np.unique(self.data['night'])) == len(offsets):  # check if each night has an associated offset value
            for i in range(len(self.data['diffMag'])):
                finishedOffset = False
                for j in range(len(offsets)):
                    if finishedOffset:
                        break
                    elif self.data['night'][i] == offsets.keys()[j]:
                        self.data['diffMag'][i] += offsets[offsets.keys()[j]]
        else:
            print 'Error: number of nights does not match the number of offsets'
            return
            
    ''' must pass this object a filename, and an n x 2 list of data types and associated columns
    '''
    def __init__(self, objectName, fileName, formatSpec, offsets = None, sortby = 'jd'):     # TODO - fix the naming of offsets here
        # check for errors in the format specification
        if len(formatSpec) <= 0:
            print 'Error: formatSpec input incorrect'
            return
        for i in range(len(formatSpec)):
            if len(formatSpec[i]) !=2:
                print 'Error: formatSpec input incorrect'
                return
        
        self.name = objectName
        self.data = {}      # initialize an empty dictionary with the keys in formatSpec
        for i in range(len(formatSpec)):
            self.data[formatSpec[i][0]] = {}
        if 'jd' not in self.data.keys() or 'diffMag' not in self.data.keys() or 'magErr' not in self.data.keys():
            print 'Error: Insufficient data for light curve'
            print '       Did you include \'jd\', \'diffMag\', and \'magErr\'?'
            print '       Right now, I have ' + str(self.data.keys())
        else:
            self.getData(fileName, formatSpec)
            self.sortByDataType(sortby)
            self.offsetMags(offsets)

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
        print 'floatRange input error'
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
    Takes the dataset(a lightCurveData object), the range of orders to check (a two-element list), and an intial guess at the period.
    Minimum and maximum values for the period may also be inputs
'''
def fitData(lightCurveData, initPeriod, minPeriod = None, maxPeriod = None, orderRange = [2,6]):
    # information from light curve for plotting
    time = lightCurveData.data['jd']-lightCurveData.data['jd'][0]
    mag = lightCurveData.data['diffMag']
    err = lightCurveData.data['magErr']

    bestResult = None       # storing the best fitting result
    bestParams = None
    bestVar = None
    bestOrder = None

    for m in range(orderRange[0],orderRange[1]+1):

        params = Parameters()   # create the Parameters object and add P, A, and B values to it
        params.add('P', value = initPeriod, min = minPeriod, max = maxPeriod)
        for i in range(m):
            params.add('A' + str(i+1), value = 0)
            params.add('B' + str(i+1), value = 0)
        
        fitResult = minimize(makeModel, params, args = (time, mag, err))

        if fitResult.success:
            n = float(len(time))       # number of observations
            k = float(2*m+1)           # total free parameters
            var = 1/(n-k)*sum(makeModel(params, time, mag, err)**2)
            if bestResult is None:
                bestResult = fitResult
                bestParams = params
                bestVar = var
                bestOrder = m
            elif bestVar > var:
                bestResult = fitResult
                bestParams = params        
                bestVar = var
                bestOrder = m
            print 'optimization succeeded for m = ' + str(m)
        else:
            print 'optimization failed for m = ' + str(m)
        
    return (bestResult, bestOrder)

''' Phase folds the magnitude data given a set that has the phase offset at the first data point
'''
def phaseFold(time, period):
    print 'fitted period: ' + str(period)
    for i in range(len(time)):
        periodOffset = int(time[i]/period)
        if periodOffset >= 1:
            time[i] -= period*periodOffset

    return time


''' Takes the best fit output and the data and plots them
    plotting error bars is optional (default is True)
    plotting the full phase is optional (default is True)
'''
def plotAndPrintResults(fit, m, lightCurveData, printReport = True, plotFullPeriod = True, plotErrors = True, phaseFoldData = True):
    if printReport:
        report_fit(fit.params)
    time = lightCurveData.data['jd']-lightCurveData.data['jd'][0]
    mag = lightCurveData.data['diffMag']
    err = lightCurveData.data['magErr']

    fittedPeriod = fit.params['P'].value

    if phaseFoldData:
        time = phaseFold(time, fittedPeriod)
    
    if plotFullPeriod:
        modelTime = np.linspace(time.min(), time.min()+fittedPeriod,100)
    else:
        modelTime = np.linspace(time.min(), time.max(),100)

    plt.plot(modelTime, makeModel(bestFit.params, modelTime), 'b-')    # plot the model fit
    if (plotErrors):                                        # plot the data with error bars (default)
        plt.errorbar(time, mag, err, fmt = 'rx')
    else:
        plt.plot(time, mag, 'rx')

    if phaseFoldData:
        plt.xlabel('\delta JD')
    else:
        plt.xlabel('JD + ' + str(int(min(lightCurveData.data['jd']))))    # the subtraction simplifies display of the JD
    plt.ylabel('Differential Magnitude')
    ax = plt.gca()
    ax.invert_yaxis()                                       # flip the y axis
    plt.title(objectName + ' Light Curve, ' + ' P = ' + str(fit.params['P'].value*24.0) + 'h, m = ' + str(m))

    plt.show()

###### start of the run routine ######

# the only lines that should need editing from object to object
##objectName = 'Spartacus20090130'
##fileName = 'Spartacus20090130_MANOS.txt'
##formatSpec = [['jd',3],['diffMag',6],['magErr',7]]
##offsets = None
##T0 = 0.25                                # initial guess at period

objectName = 'da14'
fileName = 'da14.120301.568.mag.txt'
formatSpec = [['jd',0],['diffMag',1],['magErr',2]]
offsets = None
T0 = 0.25

##objectName = 'Elisa'
##fileName = 'elisa_mine_mod.txt'
##formatSpec = [['night',0], ['jd',2],['diffMag',5],['magErr',6]]
##offsets = {1:0.0,2:-0.04,3:0.464}
##T0 = 0.25

# note that these functions are all overloaded (there are extra options you can set- see the function definitions for examples
lcd = lightCurveData(objectName, fileName, formatSpec, offsets)         # read in data from the text file and create a lightCurveData object
(bestFit, m) = fitData(lcd, T0)                                         # fit the data to a model (can also add min and max periods in JD)
plotAndPrintResults(bestFit, m, lcd)                                    # plot and print the results
plt.plot(bestFit.residual, 'rx')
plt.show()
