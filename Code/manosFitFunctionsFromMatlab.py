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

##------------------------------------------------------##

### Python implementation of Matlab function periodis.m, by Eran O. Ofek, March 1994
### http://wise-obs.tau.ac.il/~eran/matlab.html
##def periodis(x,l_f,h_f,df,pr=0):
### input  : - Timeseries x, two-column array of [Time, Mag], in which the first column
###            is the time and the second column is the magnitude.
###          - Lowest frequency (h_l).
###          - Highest frequency (h_f).
###          - Frequency interval (df).
###          - Optional probability cutoff (pr), default no cutoff (pr = 0).
###            Power spectra peaks with probability smaller than this
###            cutoff are not listed.
##    
##    c_x = 1
##    c_y = 2
##
#### START MATLAB
##
##    noise = std(x(:,c_y)).*std(x(:,c_y))
##    N0 = length(x(:,c_x));
##    Ni = -6.362 + 1.193.*N0 + 0.00098.*N0.*N0;
##    tem0 = x(:,c_y) - mean(x(:,c_y));
##    f_ind = l_f;
##    k = 1;
##    Pxw = zeros((h_f-l_f)./df,2);
##    while f_ind<h_f,
##       Pxw(k,1) = f_ind; 
##       om = 2.*pi.*f_ind;
##       tau = sum(sin(2.*om.*x(:,c_x)))./sum(cos(2.*om.*x(:,c_x)));
##       tau = atan(tau)./(2.*om);
##       Axc = cos(om.*(x(:,c_x) - tau));
##       Axs = sin(om.*(x(:,c_x) - tau));
##       Ax1 = sum(tem0.*Axc);
##       Ax1 = Ax1.*Ax1;
##       Ax2 = sum(tem0.*Axs);
##       Ax2 = Ax2.*Ax2;
##       temp = Ax1./sum(Axc.*Axc) + Ax2./sum(Axs.*Axs);
##       Pxw(k,2) = 0.5.*temp;
##       f_ind = f_ind + df;
##       k = k + 1;
##    end
##    % normalization of the periodogram
##    Pxw(:,2) = Pxw(:,2)./noise;
##
##
##
##    if (nargout==2),
##       DiffVec = diff(sign(diff([0;Pxw(:,2);0])));
##       IdV     = find(DiffVec==-2);
##
##       Fm      = [Pxw(IdV,1), Pxw(IdV,2), 1./Pxw(IdV,1), (1-exp(-Pxw(IdV,2))).^Ni];
##       Fm      = sortrows(Fm,2);
##    end
##
##
#### END MATLAB
##
##
##    return [Pxw, Fm]


def fitharmo(X,Y,DelY,Har,Deg,PlotPar,File):
    return [Par,Par_Err,Cov,Chi2,Freedom,Par1,Resid]

def fitFonlyByHarmo(Frange, Har, Deg, time, Mag, MagErr, cosmicErr, figureNum):
    return [bestF, dFreq, Chi2vec]

def fitFGbyHarmo(Grange, phi1, phi2, Frange, Har, Deg, time, Mag, MagErr, cosmicErr, plotNum):
    return [bestG,dG,bestF,dFreq,Chi2mat]

def run_calSpin4(Field, Chip, Date, AstName, SuperArcs, FigureNum):
    return SuperArcs


## realhist test case
R = [random.random() for _ in xrange(1000)]
realhist(R, [0,1,50])
plt.show()




