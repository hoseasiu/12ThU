# various MANOS curvefitting routines
from numpy import *
import matplotlib.pyplot as plt
from scipy import stats

# Python implementation of Matlab function realhist.m, by Eran O. Ofek, 1995
# http://wise-obs.tau.ac.il/~eran/matlab.html
def realhist(V, Range = [0], Option = 'n', Norm = 'a'):
    # V should be a numpy column vector
    # Range, Option, and Norm are optional arguments, and already have defaults
    # Range is calculated from V if not specified
    Range = array(Range)
    Eps = 1e-8
    NBinDef = 10
    OneSig = 0.6827

    if Range.shape == (1,):
        MinVal = min(V)
        MaxVal = max(V)
        NormEps = (MaxVal-MinVal)*Eps
        Range = array([MinVal-NormEps, MaxVal+NormEps, NBinDef])

    if len(Range.shape) == 1:
        # Range case
        if len(Range) == 2:
            Range[2] = NBinDef
        Step = (float(Range[1])-float(Range[0]))/float(Range[2])

        # TODO - not sure if the -1 below is necessary
        temp = array(range(0,int(Range[2])-1,1), 'float64')
        X = Range[0] + Step*0.5 + Step*temp
        N = plt.hist(V,X)
        plt.draw()
        Edges = linspace(Range[0], Range[1], Step)
        EL = len(Edges)
    else:
        # Edges case
        # bin centers
        EL = max(Range.shape)

        # TODO - not sure if the -1 below is necessary
        X  = 0.5*(Range[0:EL-1] + Range[1:EL]);

        # TODO - not sure what the histc_debug thing does

        Edges = Range

    # options
    if Option == 'n':
        pass # do nothing, default
    elif Option == 'c+':
        # cumulative (increasing)
        N = cumsum(N)
    elif Option == 'c-':
        CN = cumsum(rot90(N,2))
        N = rot90(CN,2)
    else:
        print 'Unknown option'

    # always calculate the Poisson errors, since nargout isn't an option in Python
    Bounds = zeros([len(N), 2])
##    for i in range(len(N)):
##        B = stats.poisson.interval(1-OneSig, N[i])
        # TODO - Poisson fitting function in python? - [M, B] = poissfit(N(I),1-OneSig);
        # not sure if this works the same way

        # TODO- check array shape
##        Bounds[i,0:1] = B

    # TODO - normalizations
    if Norm == 'a':
        pass # do nothing, default
    elif Norm == 'p':
        SumN = sum(N)
        N = N/SumN
        Bounds = Bounds/SumN
    elif Norm == 's':
        # TODO not sure if the -1 is needed here
        N = N/abs(sin(Edges[0:EL-1]) - sin(Edges[1:EL]))
        Bounds = Bounds/abs(sin(Edges[0:EL-1]) - sin(Edges[1:EL]));
    else:
        print 'Unknown normalization'
    
    return [X,N,Edges,Bounds]

## realhist test case
R = [random.random() for _ in xrange(1000)]
realhist(R, [0,1,50])
plt.show()
