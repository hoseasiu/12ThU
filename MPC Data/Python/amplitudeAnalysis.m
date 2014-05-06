% MPC light curve analysis
clear all; close all; clc;
file = 'parsedMPCdata.xlsx';

amplitude = xlsread(file,1,'K2:K275');
moid = xlsread(file,1,'L2:L275');
H = xlsread(file,1,'M2:M275');

% % minH = min(h);
% minH = 22;
% % maxH = max(h);
% maxH = 27;

for i = fliplr(1:length(moid))
    if isnan(amplitude(i)) || isnan(moid(i)) || isnan(H(i))
        amplitude(i) = [];
        moid(i) = [];
        H(i) = [];
    end
end

moidLD = moid/0.002655;     % convert MOID from AU to LD

combinedData = sortrows([amplitude moidLD H],2);         % amplitude, MOID, H
movingAvg = zeros(length(combinedData)-29,6);       % average, std period, midpoint, std frequency
cumMovAvg = zeros(length(combinedData)-29,5);
for edge = 1:(length(combinedData)-29)
    movingAvg(edge,1) = mean(combinedData(edge:edge+29,1));     % average amp of the set
    movingAvg(edge,2) = std(combinedData(edge:edge+29,1));      % std of amp of the set
    movingAvg(edge,3) = mean(combinedData(edge:edge+29,2));     % midpoint of MOID
    movingAvg(edge,6) = median(combinedData(edge:edge+29,1));   % median
    % cumulative moving average
    cumMovAvg(edge,1) = mean(combinedData(1:edge+29,1));
    cumMovAvg(edge,2) = std(combinedData(1:edge+29,1));
    cumMovAvg(edge,3) = mean(combinedData(1:edge+29,2));
end

% figure;
% scatter(combinedData(:,2),combinedData(:,1));

% new way of finding MOID breaks
LDbreaks = zeros(round(max(combinedData(:,2))),1);        % contains the breakpoints for MOID in the dataset
breaksCounter = 1;
while breaksCounter < length(LDbreaks)
    for j = 2:length(combinedData)
        if combinedData(j-1,2) <= breaksCounter && combinedData(j,2) > breaksCounter
            LDbreaks(breaksCounter) = j-1;
            breaksCounter = breaksCounter + 1;
        end
    end
end
LDbreaks(end) = [];     % TODO - kinda sketchy, fix this later

maxLDdivide = 10;

tP = zeros(maxLDdivide,1);
fP = zeros(maxLDdivide,1);
n1std = zeros(maxLDdivide,1);
n2std = zeros(maxLDdivide,1);
n1mean = zeros(maxLDdivide,1);
n2mean = zeros(maxLDdivide,1);

for i = 1:maxLDdivide
    % counting sample sizes
    n1(i) = length(combinedData(1:LDbreaks(i),1));
    n2(i) = length(combinedData(LDbreaks(i):end,1));
    
    n1std(i) = std(combinedData(1:LDbreaks(i),1));
    n2std(i) = std(combinedData(LDbreaks(i):end,1));
    n1mean(i) = mean(combinedData(1:LDbreaks(i),1));
    n2mean(i) = mean(combinedData(LDbreaks(i):end,1));

    % t-testing
    [tH, tP(i), tci(i,:), tstats(i)] = ttest2(combinedData(1:LDbreaks(i),1),combinedData(LDbreaks(i):end,1),'Tail','right'); % one-tailed t-test of amplitude
    
    % f-testing
    [fH, fP(i), fci(i,:), fstats(i)] = vartest2(combinedData(1:LDbreaks(i),1),combinedData(LDbreaks(i):end,1),'Tail','right'); % one-tailed f-test of amplitude
    
end


%%

statTest = figure;
plot(1:length(tP),tP,'d-','MarkerFaceColor','b');
hold on;
plot(1:length(fP),fP,'ro-','MarkerFaceColor','r');
legend('t-test', 'f-test','Location','Best');
ylabel('p value');
xlabel('bin edge (LD)');
plot(xlim,[0.1 0.1],'g--');
plot(xlim,[0.05 0.05],'g--');
plot(xlim,[0.01 0.01],'g--');
title('t-test and f-test p values');
hold off;
% saveas(statTest,'..\..\Documentation\Thesis\statTest_amp.png');

moidVsAmp = figure;
shadedErrorBar(movingAvg(:,3),movingAvg(:,1),movingAvg(:,2));
hold on;
plot([1 1],ylim,'r');
plot(movingAvg(:,3),movingAvg(:,1), 'kx');
xlabel('MOID (LD)');
ylabel('Light Curve Amplitude (mag)');
title('Moving Average Amplitude (30-Object Window)');
hold off;
% saveas(moidVsAmp,'..\..\Documentation\Thesis\moving_average_amp.png');


moidVsStdAmp = figure;
plot(movingAvg(:,3),movingAvg(:,2),'x-');
hold on;
plot([1 1],ylim,'r');
xlabel('MOID (LD)');
ylabel('Standard Deviation of the Mean of Light Curve Amplitude (mag)');
title('Standard Deviation of Moving Average Amplitude (30-Object Window)');
hold off;
% saveas(moidVsStdAmp,'..\..\Documentation\Thesis\std_moving_average_amp.png');

amplitudeVsMoid = figure;
semilogx(moidLD, amplitude, 'o');
title('NEO Light Curve Amplitude');
xlabel('MOID (LD)');
ylabel('Amplitude (mags)');
hold on;
plot([1 1],ylim,'r');
hold off;
% saveas(amplitudeVsMoid,'..\..\Documentation\Thesis\amplitude_vs_moid.png');

HVsAmplitude = figure;
scatter(H, amplitude, 'o');
title('Magnitude vs Amplitude');
xlabel('Absolute Magnitude');
ylabel('Amplitude (mags)');
% saveas(HVsAmplitude,'..\..\Documentation\Thesis\amplitude_vs_h.png');