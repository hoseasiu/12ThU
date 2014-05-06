function [ newPoints ] = chauvenets3to10points( points )
%chauvenets3to10points Performs Chauvenet's Criterion rejection for 3 to 10
%points

% populate the max ratios list up to 10
maxRatio = [0 0 1.38 1.54 1.65 1.75 1.8 1.87 1.91 1.96];

% get rid of NaNs
points(isnan(points)) = [];
N = length(points);
if N < 3
    newPoints = points;
else
    
    stdPoints = std(points);
    meanPoints = mean(points);
    
    newPoints = [];
    
    for i = 1:N
        tau = abs(points(i)-meanPoints)/stdPoints;
        if tau < maxRatio(N)% && points(i)>(-2*stdPoints)
            newPoints(length(newPoints)+1) = points(i);
        end
    end
end

end

