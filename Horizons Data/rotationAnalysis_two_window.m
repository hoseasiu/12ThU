% JPL Horizons light curve analysis
clear all; close all; clc;
file = 'horizons_data_analysis.xlsx';

period = xlsread(file,1,'G2:G10938');
moid = xlsread(file,1,'J2:J10938');
h = xlsread(file,1,'D2:D10938');

minH = min(h);
% minH = 25;
maxH = max(h);
% maxH = 30;

% trim the longer array to be the same size as the shorter array
if length(period)~=length(moid)
    if length(period)>length(moid)
        period(length(moid)+1:end) = [];
        h(length(moid)+1:end) = [];
    else
        moid(length(period)+1:end) = [];
        h(length(period)+1:end) = [];
    end
end        

for i = fliplr(1:length(period))
    if isnan(period(i))
        period(i) = [];
        moid(i) = [];
        h(i) = [];
    end
end

% % TODO - justify this removal
maxIndex = find(period == max(period));
period(maxIndex) = [];
moid(maxIndex) = [];
h(maxIndex) = [];

maxIndex = find(period == max(period));
period(maxIndex) = [];
moid(maxIndex) = [];
h(maxIndex) = [];

moidLD = moid/0.002655;     % convert MOID from AU to LD

arr = fliplr(1:length(h));
for i = arr
    if h(i) < minH || h(i) > maxH
        period(i) = [];
        moidLD(i) = [];
        h(i) = [];
    end
end

combinedData = sortrows([period moidLD h 1./period],2);         % period, MOID, H, freqeuncy
movingAvg = zeros(length(combinedData)-29,6);       % average, std period, midpoint, std frequency
cumMovAvg = zeros(length(combinedData)-29,5);
for edge = 1:(length(combinedData)-29)
    movingAvg(edge,1) = mean(combinedData(edge:edge+29,1));     % average period of the set
    movingAvg(edge,2) = std(combinedData(edge:edge+29,1));      % std of period of the set
    movingAvg(edge,3) = mean(combinedData(edge:edge+29,2));     % midpoint of MOID
    movingAvg(edge,4) = std(combinedData(edge:edge+29,4));      % std of frequency of the set
    movingAvg(edge,5) = mean(combinedData(edge:edge+29,4));     % average frequency of the set
    movingAvg(edge,6) = median(combinedData(edge:edge+29,1));   % median
    % cumulative moving average
    cumMovAvg(edge,1) = mean(combinedData(1:edge+29,1));
    cumMovAvg(edge,2) = std(combinedData(1:edge+29,1));
    cumMovAvg(edge,3) = mean(combinedData(1:edge+29,2));
end

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

freq = 1./period;

%% draw figures

binningExample = figure;
hist(moidLD,174);
hold on;
plot([5 5],ylim,'r','LineWidth',2);
xlim([0 60]);
xlabel('MOID (LD)');
ylabel('number of objects');
title('Binning Example');


%%
load('data20to25');
load('data25to30');

freqVsMoid = figure;
loglog(moidLD,freq,'o');
title('NEO Rotation Frequency');
xlabel('MOID (LD)');
ylabel('Rotation Frequency (1/h)');
hold on;
plot([1 1],[1E-3 1E4],'r');

freqVsH = figure;
semilogy(h20to25,freq20to25,'xb');
hold on;
semilogy(h25to30,freq25to30,'og');
set(gca,'XDir','reverse')
xlabel('Absolute Magnitude (H)');
ylabel('Rotation Frequency (1/h)');
title('Magnitude vs Rotational Frequency');
legend('Magnitude 20 to 25', 'Magnitude 25 to 30','Location','Best');
hold off;

freqVsMoidComparison = figure;
loglog(moidLD20to25,freq20to25,'xb');
hold on;
loglog(moidLD25to30,freq25to30,'og');
plot([1 1], ylim,'r');
xlabel('MOID (LD)');
ylabel('Rotation Frequency (1/h)');
title('MOID vs Rotational Frequency');
legend('Magnitude 20 to 25', 'Magnitude 25 to 30');
hold off;

hVsMoid = figure;
semilogx(moidLD,h,'o');
hold on;
plot([1 1],ylim,'r');
set(gca,'YDir','reverse')
xlabel('MOID (LD)');
ylabel('Absolute Magnitude (H)');
title('Magnitude vs MOID');
hold off;

% moidVsFreqZoomed = figure;
% shadedErrorBar(movingAvg(:,3),movingAvg(:,5),movingAvg(:,4));
% hold on;
% plot([1 1],ylim,'r');
% plot(movingAvg(:,3),movingAvg(:,5), 'kx');
% xlabel('MOID (LD)');
% ylabel('Rotation Frequency (1/h)');
% title('Moving Average Rotational Frequency (Magnitude 20 to 25, 30-Object Window)');
% % xlim([0 10]);
% hold off;

moidVsFreqStdZoomed = figure;
plot(movingAvg20to25(:,3),movingAvg20to25(:,4),'x-b');
hold on;
plot(movingAvg25to30(:,3),movingAvg25to30(:,4),'o-g');
plot([1 1],ylim,'r');
xlabel('MOID (LD)');
ylabel('Standard Deviation of the Mean of Rotation Frequency (1/h)');
title('Standard Deviation of Moving Average Rotational Frequency (30-Object Window)');
legend('Magnitude 20 to 25', 'Magnitude 25 to 30');
xlim([0 10]);
hold off;

% moidVsPeriodZoomed = figure;
% % shadedErrorBar(movingAvg(:,3),movingAvg(:,1),movingAvg(:,2));
% plot(movingAvg(:,3),movingAvg(:,2),'x');
% hold on;
% plot([1 1],ylim,'r');
% % plot(movingAvg(:,3),movingAvg(:,1), 'kx');
% xlabel('MOID (LD)');
% ylabel('Rotation Period (h)');
% title('Moving Average Rotational Period (30-Object Window)');
% % xlim([0 10]);
% hold off;

%%

maxLDdivide = 10;

ttestP = zeros(maxLDdivide,1);
ftestP = zeros(maxLDdivide,1);
bin1std = zeros(maxLDdivide,1);
bin2std = zeros(maxLDdivide,1);

for i = 1:maxLDdivide
    % counting sample sizes
    n1(i) = length(combinedData(1:LDbreaks(i),1));
    n2(i) = length(combinedData(LDbreaks(i):end,1));
    
    bin1std(i) = std(combinedData(1:LDbreaks(i),4));
    bin2std(i) = std(combinedData(LDbreaks(i):end,4));

    % t-testing
    [h,ttestP(i)] = ttest2(combinedData(1:LDbreaks(i),4),combinedData(LDbreaks(i):end,4),'Tail','right');
    
    % f-testing
    [h,ftestP(i)] = vartest2(combinedData(1:LDbreaks(i),4),combinedData(LDbreaks(i):end,4),'Tail','right'); % one-tailed f-test of frequency
    
end

statTest = figure;
plot(1:length(ttestP),ttestP,'d-','MarkerFaceColor','b');
hold on;
plot(1:length(ftestP),ftestP,'ro-','MarkerFaceColor','r');
legend('t-test', 'f-test');
ylabel('p value');
xlabel('bin edge (LD)');
plot(xlim,[0.1 0.1],'g--');
plot(xlim,[0.05 0.05],'g--');
plot(xlim,[0.01 0.01],'g--');
% xlim([0 5.5]);
ylim([0 0.25]);
title(horzcat('t-test and f-test p values (absolute magnitude ',int2str(minH),'-',int2str(maxH),')'));
hold off;

varPlot = figure;
plot(1:length(bin1std),bin1std,'o-');
hold on;
plot(1:length(bin2std),bin2std, 'rx-');
hold off;
xlabel('bin edge (LD)');
ylabel('standard deviation of the mean');
legend('smaller MOID', 'larger MOID');
title(horzcat('Standard Deviation of Frequency vs Dividing MOID Value (absolute magnitude ',int2str(minH),'-',int2str(maxH),')'));       % need a better title

%%

% % continuous f-testing
% for i = 1:length(LDbreaks)-1
%     [h, p] = vartest2(combinedData(1:LDbreaks(1)),combinedData(LDbreaks(i):LDbreaks(i+1)));
%     fPercent(i) = p;     % some of these values will be NaN
%     [h, p] = ttest2(combinedData(1:LDbreaks(1)),combinedData(LDbreaks(i):LDbreaks(i+1)));
%     tPercent(i) = p;
% end
% figure;
% % semilogy(1:length(LDbreaks)-1,fPercent,'o');
% scatter(1:length(LDbreaks)-1,fPercent,'o');
% hold on;
% plot(xlim,[0.1 0.1]);
% plot(xlim,[0.05 0.05]);
% plot(xlim,[0.01 0.01]);
% title('f test results');
% % ylim([0 0.1])
% 
% figure;
% % semilogy(1:length(LDbreaks)-1,tPercent,'o');
% scatter(1:length(LDbreaks)-1,tPercent,'o');
% hold on;
% plot(xlim,[0.1 0.1]);
% plot(xlim,[0.05 0.05]);
% plot(xlim,[0.01 0.01]);
% title('t test results');