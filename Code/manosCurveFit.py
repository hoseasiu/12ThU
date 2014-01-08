''' manosCurveFit.py - Hosea Siu 2014
    Takes light curve data as input and fits it to a function
'''

from numpy import *
from scipy.optimize import leastsq
import matplotlib.pyplot as plt, os.path, sys

# class used to store and manipulate MANOS light curve data
class lightCurveData:
            
    def getData(self,filename):

# columns are as follows:
# Image number, Flag, Filter, JD, Airmass, Aperture size, Differential mag, mag error

        data = loadtxt(filename,'a12')          # extract data from text file

        # separate data into appropriately-named arrays
        self.flag = data[:,1].astype(int)
        self.filter = data[:,2]
        self.jd = data[:,3].astype(float)
        self.airmass = data[:,4].astype(float)
        self.apSize = data[:,5].astype(float)
        self.diffMag = data[:,6].astype(float)
        self.magErr = data[:,7].astype(float)

        #da14
##        self.jd = data[:,0].astype(float)
##        self.diffMag = data[:,1].astype(float)
##        self.magErr = data[:,2].astype(float)
        
    def __init__(self, filename):
        self.imageNumber = []
        self.flag = []
        self.filter = []
        self.jd = []
        self.airmass = []
        self.apSize = []
        self.diffMag = []
        self.magErr = []
        self.getData(filename)        

'''function that allows for range lists of floats
has some error handling, but still has machine precision issues (shouldn't be an actual problem)
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
            if (start < stop):
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

'''generates the light curve model (called by a lambda function later)
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

fitfunc = lambda params, t: model(params,t)             # target function from the model
errfunc = lambda params, t, mag, err: (fitfunc(params,t)-mag)/err  # distance to the target function


''' Takes the best fit output and the data and plots them
'''
def plotAndPrintResults(fit, data):
    time = data.jd-data.jd[0]
    mag = data.diffMag
    err = data.magErr

    modelTime = linspace(time.min(), time.max(),100)
    plt.plot(modelTime, fitfunc(fit[2],modelTime), 'r-')
    plt.errorbar(time, mag, err, fmt = 'rx')
    plt.xlabel('Julian Date +' + str(int(min(data.jd)))) # the subtraction simplifies the display of the JD
    plt.ylabel('Differential Magnitude')
    ax = plt.gca()
    ax.invert_yaxis()   # flip the y axis
    plt.title(objectName + ' Light Curve, ' + ' T = ' + str(fit[2][0]*24.0) + 'h, m = ' + str(fit[0]))

    print '\n' + 'FITTING RESULTS:'
    print 'order = ' + str(fit[0])
    print 'bias-corrected variance = ' + str(fit[1])
    parameters = fit[2].tolist()
    print 'period = ' + str(parameters.pop(0)*24.0) + ' hours'
    A = parameters[0:len(parameters)/2]
    B = parameters[len(parameters)/2:len(parameters)+1]
    print 'A values: ' + str(A)
    print 'B values: ' + str(B)
    plt.show()


###### start of the run routine ######

# the only two lines that should need editing from object to object
objectName = 'Spartacus20090130'
##objectName = 'da14.120301.568.mag'
filename = 'Spartacus20090130_MANOS.txt'
# initial guess at period
T0 = 0.2

# read in and format the light curve data from a text file
basepath = os.path.abspath(os.path.dirname(sys.argv[0]))
filepath = os.path.abspath(os.path.join(basepath, "..", 'Data', filename))
lcd = lightCurveData(filepath)

# information from light curve for plotting
time = lcd.jd-lcd.jd[0]
mag = lcd.diffMag
err = lcd.magErr

bestResult = (-1, -1, array([]))       # storing the best order, error, array

for m in range(2,7):        # try 2nd to 6th order fits
        
    # format the initial parameter guesses
    params0 = [T0] + [0,0]*m

    params1, success = leastsq(errfunc, params0[:], args=(time,mag,err))

    if success:
        n = float(len(time))       # number of observations
        k = float(2*m+1)           # total free parameters (number of days not counted because we're dong relative photometry only)
        var = 1/(n-k)*sum(errfunc(params1, time, mag, err)**2)
        print 'optimization succeeded for m = ' + str(m) + ', bias-corrected variance = ' + str(var)
        # check to see if this is a better result
        if bestResult[1] == -1 or bestResult[1] > var:
            bestResult = (m,var,params1)
    else:
        print 'optimization failed for m = ' + str(m)

plotAndPrintResults(bestResult, lcd)    # plot and print the results
