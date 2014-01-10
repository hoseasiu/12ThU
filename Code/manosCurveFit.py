''' manosCurveFit.py - Hosea Siu 2014
    Takes light curve data as input and fits it to a function accoring to Harris et al. 1989
'''

from lmfit import minimize, Parameters, Parameter, report_fit
import matplotlib.pyplot as plt, os.path, sys, numpy as np

''' Class used to store and manipulate MANOS light curve data
'''
class lightCurveData:
            
    def getData(self, filename, formatSpec, magOffset=0):
        basepath = os.path.abspath(os.path.dirname(sys.argv[0]))
        filepath = os.path.abspath(os.path.join(basepath, "..", 'Data', fileName))

        textFile = np.loadtxt(filepath,'a12')          # extract data from text file

        # separate data into appropriately-named arrays
        # TODO- currently assumes that everything is a float- generalize later (needed?)
        for i in range(len(formatSpec)):
            if formatSpec[i][0] == 'diffMag':
                self.data[formatSpec[i][0]] = textFile[:,formatSpec[i][1]].astype(float) + magOffset
            else:
                self.data[formatSpec[i][0]] = textFile[:,formatSpec[i][1]].astype(float)

    # must pass this object a filename, and an n x 2 list of data types and associated columns
    def __init__(self, objectName, fileName, formatSpec, magOffset = 0.0):
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
            self.getData(fileName, formatSpec, magOffset)


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
def makeModel2(params, t, mag=None, err=None):
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
def fitData(lightCurveData, orderRange, initPeriod, minPeriod = None, maxPeriod = None):
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
        
        fitResult = minimize(makeModel2, params, args = (time, mag, err))

        if fitResult.success:
            n = float(len(time))       # number of observations
            k = float(2*m+1)           # total free parameters
            var = 1/(n-k)*sum(makeModel2(params, time, mag, err)**2)
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


''' Takes the best fit output and the data and plots them
    plotting error bars is optional (default is True)
    plotting the full phase is optional (default is True)
'''
def plotAndPrintResults(fit, m, lightCurveData, plotFullPeriod = True, plotErrors = True):
    report_fit(fit.params)
    time = lightCurveData.data['jd']-lightCurveData.data['jd'][0]
    mag = lightCurveData.data['diffMag']
    err = lightCurveData.data['magErr']

    fittedPeriod = fit.params['P'].value

    if plotFullPeriod:
        modelTime = np.linspace(time.min(), time.min()+fittedPeriod,100)
    else:
        modelTime = np.linspace(time.min(), time.max(),100)

    plt.plot(modelTime, makeModel2(bestFit.params, modelTime), 'b-')    # plot the model fit
    if (plotErrors):                                        # plot the data with error bars (default)
        plt.errorbar(time, mag, err, fmt = 'rx')
    else:
        plt.plot(time, mag, 'rx')
    plt.xlabel('Julian Date +' + str(int(min(lightCurveData.data['jd']))))    # the subtraction simplifies display of the JD
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
##magOffset = 0
##T0 = 0.25                                # initial guess at period

objectName = 'da14'
fileName = 'da14.120301.568.mag.txt'
formatSpec = [['jd',0],['diffMag',1],['magErr',2]]
magOffset = -18.5
T0 = 0.25

lcd = lightCurveData(objectName, fileName, formatSpec, magOffset)       # read in data from the text file
(bestFit, m) = fitData(lcd, [2,6], T0)                                  # fit the data to a model (can also add min and max periods in JD)
plotAndPrintResults(bestFit, m, lcd)                                    # plot and print the results
