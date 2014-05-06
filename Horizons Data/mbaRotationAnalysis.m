% MBA analysis
clear all; close all; clc;
% mbaFile = 'horizons_main_belt_with_rot.csv';
% neoFile = 'horizons_data_analysis.xlsx';
% 
% mbaPeriod = xlsread(mbaFile,1,'G2:G4032');
% mbaFreq = 1./mbaPeriod;
% mbaMoidLD = xlsread(mbaFile,1,'K2:K4032');
% mbaH = xlsread(mbaFile,1,'D2:D4032');
% 
% neoPeriod = xlsread(neoFile,1,'G2:G10938');
% neoMoid = xlsread(neoFile,1,'J2:J10938');
% neoH = xlsread(neoFile,1,'D2:D10938');
% 
% % trim the longer array to be the same size as the shorter array
% if length(neoPeriod)~=length(neoMoid)
%     if length(neoPeriod)>length(neoMoid)
%         neoPeriod(length(neoMoid)+1:end) = [];
%         neoH(length(neoMoid)+1:end) = [];
%     else
%         neoMoid(length(neoPeriod)+1:end) = [];
%         neoH(length(neoPeriod)+1:end) = [];
%     end
% end        
% 
% for i = fliplr(1:length(neoPeriod))
%     if isnan(neoPeriod(i)) || neoPeriod(i) == 0
%         neoPeriod(i) = [];
%         neoMoid(i) = [];
%         neoH(i) = [];
%     end
% end
% 
% neoFreq = 1./neoPeriod;

load('mbaAndNeoData');

% count the number of objects from 5th to 30th magnitude
for mag = 1:30
    count = 1;
    for i = 1:length(mbaH)      % counting MBAs
        if mag <= mbaH(i) && mbaH(i) <= mag+1
            mbaMagSubset(mag,count,:) = [mbaH(i) mbaFreq(i)];       % mag subset dimensions represent (magnitude bin, object number (for reference), and [H, freq])
            count = count + 1;
        end
    end
end

for mag = 1:30
    count = 1;
    for i = 1:length(neoH)      % counting NEOs
        if mag <= neoH(i) && neoH(i) <= mag+1
            neoMagSubset(mag,count,:) = [neoH(i) neoFreq(i) neoMoid(i)/0.002655];
            count = count + 1;
        end
    end
end

% count the number of elements in each magnitude bin
for i = 1:(length(mbaMagSubset(:,1,1)))
    mbaNumVec(i) = nnz(mbaMagSubset(i,:,1));
end
mbaNumVec = [mbaNumVec zeros(1,30-length(mbaNumVec))]; % fill in extra zeros

for i = 1:(length(neoMagSubset(:,1,1)))
    neoNumVec(i) = nnz(neoMagSubset(i,:,1));
end
neoNumVec = [neoNumVec zeros(1,30-length(neoNumVec))]; % fill in extra zeros

% rebin for analysis
neoHVec = reshape(neoMagSubset(14:17,:,1),1,[]);
neoHVec(neoHVec==0) = [];
neoFreqVec = reshape(neoMagSubset(14:17,:,2),1,[]);
neoFreqVec(neoFreqVec==0) = [];
neoMoidVec = reshape(neoMagSubset(14:17,:,3),1,[]);
neoMoidVec(neoMoidVec==0) = [];

mbaFreqVec = reshape(mbaMagSubset(14:17,:,2),1,[]);
mbaFreqVec(mbaFreqVec==0) = [];

nNeo = zeros(1,10);
meanFreqNeo = zeros(1,10);
stdFreqNeo = zeros(1,10);
% remove non-encounter NEOs
for j = fliplr(1:10)
    for i = fliplr(1:length(neoMoidVec))
        if neoMoidVec(i) > j
            neoHVec(i) = [];
            neoFreqVec(i) = [];
            neoMoidVec(i) = [];
        end
    end
    nNeo(j) = length(neoHVec);
    meanFreqNeo(j) = mean(neoHVec);
    stdFreqNeo(j) = std(neoHVec);

    % f and t test
    [fH, fP(j), fci(j,:), fstats(j)] = vartest2(neoFreqVec,mbaFreqVec,'Tail','right');
    [tH, tP(j), tci(j,:), tstats(j)] = ttest2(neoFreqVec,mbaFreqVec,'Tail','right');
end

%%
% plot the comparison of the two counts
countComparison = figure;
width = 0.8;
bar(1:30,mbaNumVec,width);
hold on;
bar(1:30,neoNumVec,width/3,'FaceColor','r','EdgeColor','r');
plot(xlim,[30, 30],'g--');
title('Number of Objects in 1-Magnitude-Wide Bins');
xlabel('Absolute Magnitude (H)');
ylabel('Number of Objects');
legend('MBAs', 'NEOs', '30-object reference line');
ylim([0 100]);
hold off;
% saveas(countComparison,'..\Documentation\Thesis\mbaNeoCountComparison.png');

magVsFreq = figure;
semilogy(mbaH,1./mbaPeriod,'o');
ylim([1E-3 1E3]);
xlim([5 35]);
hold on;
semilogy(neoH,1./neoPeriod,'or');
xlabel('Absolute Magnitude (H)');
ylabel('Rotation Frequency (1/h)');
legend('Main Belt Objects','Near-Earth Objects','Location','Best');
title('Main Belt and Near-Earth Objects Magnitude vs Rotational Frequency');
hold off;
% saveas(magVsFreq,'..\Documentation\Thesis\mag_vs_frequency_mba_neo.png');

% plot MOID vs absolute magnitude
moidCheck = figure;
cstring = 'gbcm';
mag = 14:17;
for i = 1:length(mag)
    neoMoidVec = reshape(neoMagSubset(mag(i),:,3),1,[]);
    neoMoidVec(neoMoidVec==0) = [];
    neoHVec = reshape(neoMagSubset(mag(i),:,1),1,[]);
    neoHVec(neoHVec==0) = [];
    semilogx(neoMoidVec,neoHVec,horzcat('o',cstring(i)));
    hold on;
end
semilogx([1 1], ylim, 'r');
hold off;
set(gca,'YDir','reverse');
xlabel('MOID (LD)');
ylabel('Absolute Magnitude (H)');
legend('H = 14-15','H = 15-16','H = 16-17','H = 17-18','Location','Best');
title('MOID vs Magnitude for 14th to 18th Magnitude NEOs');
% saveas(moidCheck,'..\Documentation\Thesis\neoMoidMag14to18.png');

statTest = figure;
plot(1:length(tP),tP,'d-','MarkerFaceColor','b');
hold on;
plot(1:length(fP),fP,'ro-','MarkerFaceColor','r');
legend('t-test', 'f-test');
ylabel('p value');
xlabel('bin edge (LD)');
plot(xlim,[0.1 0.1],'g--');
plot(xlim,[0.05 0.05],'g--');
plot(xlim,[0.01 0.01],'g--');
title('t-test and f-test p values');
hold off;
% saveas(statTest,'..\Documentation\Thesis\statTestMba.png');