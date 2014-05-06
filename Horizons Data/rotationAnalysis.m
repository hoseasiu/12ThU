% JPL Horizons light curve analysis
clear all; close all; clc;
% file = 'horizons_data_analysis.xlsx';
% 
% period = xlsread(file,1,'G2:G10938');
% moid = xlsread(file,1,'J2:J10938');
% H = xlsread(file,1,'D2:D10938');
% 
% % trim the longer array to be the same size as the shorter array
% if length(period)~=length(moid)
%     if length(period)>length(moid)
%         period(length(moid)+1:end) = [];
%         H(length(moid)+1:end) = [];
%     else
%         moid(length(period)+1:end) = [];
%         H(length(period)+1:end) = [];
%     end
% end        
% 
% for i = fliplr(1:length(period))
%     if isnan(period(i)) || isnan(H(i)) || isnan(moid(i))
%         period(i) = [];
%         moid(i) = [];
%         H(i) = [];
%     end
% end
% 
% moidLD = moid/0.002655;     % convert MOID from AU to LD

load('rotationAnalysisData');       % reload the full dataset
freq = 1./period;
% plot the full range
freqVsMoid = figure;
loglog(moidLD,freq,'o');
hold on;
plot([1 1], ylim,'r');
title('NEO Rotational Frequency');
xlabel('MOID (LD)');
ylabel('Rotational Frequency (1/h)');
% saveas(freqVsMoid,'..\Documentation\Thesis\frequency_vs_moid.png');

magVsDia = figure;
mags = linspace(10,30,100);
dias = 10.^((15.618-mags-2.5*log10(0.14))/5);
plot(mags,dias);
grid on;
set(gca,'XDir','reverse');
xlabel('Absolute Magnitude (H)')
ylabel('Diameter (km)');
title('Conversion from Absolute Magnitude to Diameter for Mean NEO Albedo');
% saveas(magVsDia,'..\Documentation\Thesis\mag_vs_dia.png');

magVsDiaSubset = figure;
mags2 = linspace(24.5,29.5,100);
dias2 = 10.^((15.618-mags2-2.5*log10(0.14))/5);
plot(mags2,dias2);
grid on;
set(gca,'XDir','reverse');
xlabel('Absolute Magnitude (H)')
ylabel('Diameter (km)');
title('Conversion from Absolute Magnitude to Diameter for Mean NEO Albedo (H = 24.5 to 29.5)');
% saveas(magVsDiaSubset,'..\Documentation\Thesis\mag_vs_dia_subset.png');

magVsMoid = figure;
% hLim = [5 35];
% diaLim = 10.^((15.618-hLim-2.5*log10(0.14))/5);
semilogx(moidLD,H,'o');
hold on;
plot([1 1],ylim,'r');
% [ax, H1, H2] = plotyy(moidLD,H,[1,1],diaLim,'semilogx');
% set(H1,'LineStyle','o','MarkerEdgeColor','b');
% set(H2,'Color','r');
% set(ax,'YColor','k')
% set(ax(1),'YLim',hLim);
title('Magnitude vs MOID');
xlabel('MOID (LD)');
ylabel('Absolute Magnitude (H)');
set(gca,'YDir','reverse');
hold off;
saveas(magVsMoid,'..\Documentation\Thesis\magnitude_vs_moid.png');


% find the biggest magnitude window from magnitude 20 to 30 that doesn't
% have a statistically significant correlation
preMinH = 20;
preMaxH = 30;

maxWindowSize = 0;
bestMin = 0;
bestMax = 0;
for minH = linspace(preMinH,preMaxH-0.5,(preMaxH-preMinH)*2)
    for maxH = linspace(minH+0.5,30,(30-minH)*2)
        load('rotationAnalysisData');       % reload the full dataset
        % select subsection of H based on minH and maxH
        arr = fliplr(1:length(H));
        for i = arr
            if H(i) < minH || H(i) > maxH
                period(i) = [];
                moidLD(i) = [];
                H(i) = [];
            end
        end
        freq = 1./period;           % get freqeuncy from period
        [R,P]=corrcoef(H,freq);
        if maxWindowSize < maxH-minH && P(1,2) > 0.05        % find a bigger window of non-significant correlations
            maxWindowSize = maxH-minH;
            bestMin = minH;
            bestMax = maxH;
        end
    end
end

% load the best set
load('rotationAnalysisData');       % reload the full dataset
freq = 1./period;

% plot the full range
magVsFreq = figure;
semilogy(H,freq,'o');
title('Magnitude vs Rotational Frequency');
xlabel('Absolute Magnitude (H)');
ylabel('Rotational Frequency (1/h)');
set(gca,'XDir','reverse');
% saveas(magVsFreq,'..\Documentation\Thesis\mag_vs_frequency.png');

% select subsection of H based on pre-selection (H>20)
arr = fliplr(1:length(H));
for i = arr
    if H(i) < 20
        period(i) = [];
        moidLD(i) = [];
        H(i) = [];
    end
end
freq = 1./period;           % get freqeuncy from period

% plot H > 20
magVsFreqPreselect = figure;
semilogy(H,freq,'o');
hold on;
plot([bestMin bestMin], ylim, 'k');
plot([bestMax bestMax], ylim, 'k');
title('Magnitude vs Rotational Frequency (H > 20)');
xlabel('Absolute Magnitude (H)');
ylabel('Rotational Frequency (1/h)');
set(gca,'XDir','reverse');
% saveas(magVsFreqPreselect,'..\Documentation\Thesis\mag_vs_frequency_preselect.png');

% select subsection of H based on minH and maxH
arr = fliplr(1:length(H));
for i = arr
    if H(i) < bestMin || H(i) > bestMax
        period(i) = [];
        moidLD(i) = [];
        H(i) = [];
    end
end
freq = 1./period;           % get freqeuncy from period

combinedData = sortrows([moidLD H freq],1);         % MOID, H, freqeuncy
moidMidpoints = zeros(length(combinedData)-29,1);
stdFreq = zeros(length(combinedData)-29,1);
meanFreq = zeros(length(combinedData)-29,1);

for edge = 1:(length(combinedData)-29)
    moidMidpoints(edge) = mean(combinedData(edge:edge+29,1));       % midpoint of MOID
    stdFreq(edge) = std(combinedData(edge:edge+29,3));              % std of the frequency
    meanFreq(edge) = mean(combinedData(edge:edge+29,3));            % mean of the frequency
end

%% draw figures

moidVsFreqStd = figure;
plot(moidMidpoints,stdFreq,'x-b');
hold on;
plot([1 1],ylim,'r');
xlabel('MOID (LD)');
ylabel('Standard Deviation of the Mean of Rotation Frequency (1/h)');
title('Standard Deviation of Moving Average Rotational Frequency (30-Object Window)');
hold off;
% saveas(moidVsFreqStd,'..\Documentation\Thesis\std_moving_average_rot_freq.png');

freqVsMoidComparison = figure;
loglog(moidLD,freq,'o');
hold on;
plot([1 1], ylim,'r');
xlabel('MOID (LD)');
ylabel('Rotation Frequency (1/h)');
title('MOID vs Rotational Frequency');
hold off;
% saveas(freqVsMoidComparison,'..\Documentation\Thesis\frequency_vs_moid_comparison.png');

maxLDdivide = 10;

ttestP = zeros(maxLDdivide,1);
ftestP = zeros(maxLDdivide,1);
th = zeros(maxLDdivide,1);
fh = zeros(maxLDdivide,1);
% mean and std refer to frequency
n1mean = zeros(maxLDdivide,1);
n2mean = zeros(maxLDdivide,1);
n1std = zeros(maxLDdivide,1);
n2std = zeros(maxLDdivide,1);
n1 = zeros(maxLDdivide,1);
n2 = zeros(maxLDdivide,1);

for i = 1:maxLDdivide
    
    n1(i) = find(combinedData(:,1)<i,1,'last');
    n2(i) = length(combinedData)-n1(i);
    
    % frequencies
    sample1 = combinedData(1:n1(i),3);
    sample2 = combinedData((n1(i)+1):end,3);
    
    if i == 1 || i == 2
        hCheckMeanN1(i) = mean(combinedData(1:n1(i),2));
        hCheckStdN1(i) = std(combinedData(1:n1(i),2));
        hCheckMeanN2(i) = mean(combinedData((n1(i)+1):end,2));
        hCheckStdN2(i) = std(combinedData((n1(i)+1):end,2));
    end
    
    n1mean(i) = mean(sample1);
    n2mean(i) = mean(sample2);
    
    n1std(i) = std(sample1);
    n2std(i) = std(sample2);
    
    [th(i),ttestP(i),tci(i,:),tstats(i)] = ttest2(sample1,sample2,'Tail','right');
    [fh(i),ftestP(i),fci(i,:),fstats(i)] = vartest2(sample1,sample2,'Tail','right');
    
end

hCheck = figure;
errorbar([1 4],hCheckMeanN1,hCheckStdN1,'x');
hold on;
errorbar([2 5],hCheckMeanN2,hCheckStdN2,'xr');
title('Mean and Standard Deviation of Mean Magnitude');
ylabel('Absolute Magnitude (H)');
xlabel('LD of Samples');
hold off;
saveas(hCheck,'..\Documentation\Thesis\hCheckRaw.png');

statTest = figure;
plot(1:maxLDdivide,ttestP,'d-','MarkerFaceColor','b');
hold on;
plot(1:length(ftestP),ftestP,'ro-','MarkerFaceColor','r');
legend('t-test', 'f-test');
ylabel('p value');
xlabel('bin edge (LD)');
plot(xlim,[0.1 0.1],'g--');
plot(xlim,[0.05 0.05],'g--');
plot(xlim,[0.01 0.01],'g--');
title('t-test and f-test p values');
hold off;
% saveas(statTest,'..\Documentation\Thesis\statTest_freq.png');

% statTest = figure;
% plot(1:length(ttestP),ttestP,'d-','MarkerFaceColor','b');
% hold on;
% plot(1:length(ftestP),ftestP,'ro-','MarkerFaceColor','r');
% legend('t-test', 'f-test');
% ylabel('p value');
% xlabel('bin edge (LD)');
% plot(xlim,[0.1 0.1],'g--');
% plot(xlim,[0.05 0.05],'g--');
% plot(xlim,[0.01 0.01],'g--');
% % xlim([0 5.5]);
% ylim([0 0.25]);
% title(horzcat('t-test and f-test p values (absolute magnitude ',int2str(minH),'-',int2str(maxH),')'));
% hold off;
% 
% %%
% 
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