from numpy import *

def fitharmo(X,Y,DelY,Har,Deg=1,PlotPar,File):
##
## fitharmo function      LSQ harmonies fitting.
##                      fit harmonies of the form:
##                      Y= a_1*sin(w1*t)     + b_1*cos(w1*t)   +
##                         a_2*sin(2*w1*t)   + b_2*cos(2*w1*t) + ...
##                         a_n*sin(n_1*w1*t) + b_n*cos(n_1*w1*t) + ...
##                         c_1*sin(w2*t)     + d_1*cos(w2*t) + ...
##                         s_0 + s_1*t + ... + s_n.*t.^n_s
##                         (note that w is angular frequncy, w=2*pi*f,
##                          the program is working with frequncy "f").
##                      to a set of N data points. Return the parameters,
##                      the errors on the parameters,
##                      the Chi squars, and the covariance matrix.
## input  : - Column vector of the independent variable.
##          - Column Vector of the dependent variable.
##          - Vector of the std error in the dependent variable.
##            If only one value is given, the points
##            are taken to be with equal weight. and Std error
##            equal to the value given.
##          - matrix of harmonies to fit.
##            N*2 matrix, where N is the number of different frequncies.
##            Each row should contain two numbers, the first is the
##            frequency to fit and the second is the number of harmonies
##            of that frequncy to fit. If there is more then one row
##            then all the frequncies and their harmonics will be fitted
##            simoltanusly.
##          - Degree of polynomials to fit. (Default is 0).
##          - Vector of plots control characters.
##            If argument is given then X vs. Y graph is plotted.
##            If equal to empty string (e.g. '') then plot X vs. Y
##            with red fitted function line and yellow circs for
##            the observations.
##            If one or two character are given then the first character
##            is for the observations sign and the second for the fitted
##            function line.
##            If third character is given then histogram of resdiual
##            is plotted. when the third character should contain the
##            number of bins.
##          - File name in which summary table of the fit will be written.
##            The summary table includes all information regarding the
##            fit parameters and Chi2 test.
## output : - Fitted parameters [a_1,b_1,...,a_n,b_n,c_1,d_1,...,s_0,...]
##            The order of the parameters is like the order of the
##            freqencies matrix, and then the constant + linear terms.
##          - Fitted errors in the parameters [Da_1,Db_1,...]
##          - The covariance matrix.
##          - Chi2 of the fit.
##          - Degrees of freedom.
##          - sine/cosine parameters in form od Amp. and phase (in fraction),
##            pairs of lines for [Amp, Amp_Err; Phase, Phase_Err]...
##            phase are given in the range [-0.5,0.5].
##          - The Y axis residuals vector.
## See also : fitharmonw.m
## Tested : Matlab 5.1
##     By : Eran O. Ofek                    May 1994
##                             Last Update  Mar 1999
##    URL : http://wise-obs.tau.ac.il/~eran/matlab.html

    N_X = len(X)
    N_Y = len(Y)
    N_DY = len(DelY)
    if(N_X==N_Y):
        print 'Error: X and Y must have the same length'
    if(N_X==N_Y):
        if(N_DY==1):
            # take equal weights
            if (DelY<=0):
                print 'Error: DelY must be positive'
            else:
                DelY = DelY*ones((N_X,1))
        else:
            print 'Error: Y and DelY must have the same length'

    # number of parameters
    N_Pars = Deg+1+2.*sum(Har[:,2])

    # degree of freedom
    Freedom = N_X - N_Pars

    # the size of the harmonies matrix
    Srow_Har = size(Har,0)
    Scol_Har = size(Har,1)
    if Scol_Har == 2:
        print 'Error: Number of columns in the harmonic freq. should be two'
        
    

    return [Par,Par_Err,Cov,Chi2,Freedom,Par1,Resid]
