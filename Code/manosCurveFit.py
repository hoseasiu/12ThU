''' manosCurveFit.py - Hosea Siu 2014
    Takes light curve data as input and fits it to a function accoring to Harris et al. 1989
'''

from numpy import *
from scipy.optimize import leastsq
import matplotlib.pyplot as plt, os.path, sys

'''class used to store and manipulate MANOS light curve data
'''
class lightCurveData:
            
    def getData(self,filename, formatSpec):
        textFile = loadtxt(filename,'a12')          # extract data from text file

        # separate data into appropriately-named arrays
        # TODO- currently assumes that everything is a float- generalize later (needed?)
        for i in range(len(formatSpec)):
            self.data[formatSpec[i][0]] = textFile[:,formatSpec[i][1]].astype(float)

    # must pass this object a filename, and an n x 2 list of data types and associated columns
    def __init__(self, objectName, fileName, formatSpec):
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

''' Reads in the data file (.txt format) and returns it as a lightCurveData object
'''
def readData (objectName, fileName, formatSpec):
    basepath = os.path.abspath(os.path.dirname(sys.argv[0]))
    filepath = os.path.abspath(os.path.join(basepath, "..", 'Data', fileName))
    return lightCurveData(objectName, filepath, formatSpec)


''' Generates the light curve model (called by a lambda function later)
'''
def model(params,t):
    if len(params)<5:
        print 'Error: not enough parameters- model must be at least second order'
        return
    else:
        params = params.tolist()
        P = params.pop(0)
        A = params[0:len(params)/2]
        B = params[len(params)/2:len(params)+1]

        H = zeros(len(t))
        for L in range(1,len(A)+1):
            trigTerm = 2*pi*L/float(P)
            offset = t-t[0]
            H = H + A[L-1]*sin(trigTerm*offset)+B[L-1]*cos(trigTerm*offset)
        return H

'''lambda functions used for the least squares fit
'''
fitfunc = lambda params, t: model(params,t)                         # target function from the model
errfunc = lambda params, t, mag, err: (fitfunc(params,t)-mag)/err   # distance to the target function

''' Fit the data to a model using scipy.optimize.leastsq.
    Takes the dataset(a lightCurveData object), the range of orders to check (a two-element list), and an intial guess at the period.
'''
def fitData(lightCurveData, orderRange, initPeriod):
    # information from light curve for plotting
    time = lightCurveData.data['jd']-lightCurveData.data['jd'][0]
    mag = lightCurveData.data['diffMag']
    err = lightCurveData.data['magErr']

    bestResult = (-1, -1, array([]))       # storing the best order, error, array

    for m in range(orderRange[0],orderRange[1]):
        params0 = [initPeriod] + [0,0]*m    # format the initial parameter guesses
        params1, success = leastsq(errfunc, params0[:], args=(time,mag,err))    # least squares fit

        if success:
            n = float(len(time))       # number of observations
            k = float(2*m+1)           # total free parameters
            var = 1/(n-k)*sum(errfunc(params1, time, mag, err)**2)
            print 'optimization succeeded for m = ' + str(m) + ', bias-corrected variance = ' + str(var)
            if bestResult[1] == -1 or bestResult[1] > var:          # check to see if this is a better result
                bestResult = (m,var,params1)
        else:
            print 'optimization failed for m = ' + str(m)
    return bestResult

''' Takes the best fit output and the data and plots them
    plotting error bars is optional (default is True)
    plotting the full phase is optional (default is True)
'''
def plotAndPrintResults(fit, lightCurveData, plotFullPeriod = True, plotErrors = True):
    time = lightCurveData.data['jd']-lightCurveData.data['jd'][0]
    mag = lightCurveData.data['diffMag']
    err = lightCurveData.data['magErr']

    if plotFullPeriod:
        modelTime = linspace(time.min(), time.min()+fit[2][0],100)
    else:
        modelTime = linspace(time.min(), time.max(),100)
    plt.plot(modelTime, fitfunc(fit[2],modelTime), 'r-')    # plot the model
    if (plotErrors):                                        # plot the data with error bars (default)
        plt.errorbar(time, mag, err, fmt = 'rx')
    else:
        plt.plot(time, mag, 'rx')
    plt.xlabel('Julian Date +' + str(int(min(lightCurveData.data['jd']))))    # the subtraction simplifies display of the JD
    plt.ylabel('Differential Magnitude')
    ax = plt.gca()
    ax.invert_yaxis()                                       # flip the y axis
    plt.title(objectName + ' Light Curve, ' + ' T = ' + str(fit[2][0]*24.0) + 'h, m = ' + str(fit[0]))

    print '\n' + 'FITTING RESULTS:'                         # display results on console
    print 'order = ' + str(fit[0])
    print 'bias-corrected variance = ' + str(fit[1])
    parameters = fit[2].tolist()
    print 'period = ' + str(parameters.pop(0)*24.0) + ' hours'
    print 'A values: ' + str(parameters[0:len(parameters)/2])
    print 'B values: ' + str(parameters[len(parameters)/2:len(parameters)+1])
    plt.show()

###### start of the run routine ######

# the only lines that should need editing from object to object
objectName = 'da14' #'Spartacus20090130'
fileName = 'da14.120301.568.mag.txt' #'Spartacus20090130_MANOS.txt'
formatSpec = [['jd',0],['diffMag',1],['magErr',2]] #[['jd',3],['diffMag',6],['magErr',7]]
T0 = 0.25                                # initial guess at period

lcd = readData(objectName, fileName, formatSpec)        # read in data from the text file
bestFit = fitData(lcd, [2,7], T0)                       # fit the data to a model
plotAndPrintResults(bestFit, lcd, False)                       # plot and print the results
