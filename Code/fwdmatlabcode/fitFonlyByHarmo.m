function [bestF, dFreq, Chi2vec] = fitFonlyByHarmo(Frange, Har, Deg, time, Mag, MagErr, cosmicErr, figureNum);

tic
t0 = getT0;

% Find the frequency
Chi2vec = zeros(length(Frange), 1);
for j=1:1:length(Frange),
   [Par,Par_Err,Cov,Chi2,Freedom,Par1,Resid]=fitharmo(time, Mag, MagErr.*cosmicErr, [Frange(j), Har], Deg);
   if(isreal(Chi2) == 1),
      Chi2vec(j) = Chi2;
   else
      Chi2vec(j) = nan;
   end
end
%Chi2vec
[MinChi2, minInx] = min(Chi2vec);%
bestF = Frange(minInx);%

% Old condition:
%[j]=find(Chi2vec(:,:)<MinChi2+2.2958);

NumberOfDataPoints = length(Mag)
%NumberOfDataPoints = 1;

% Using one-sigma as a limit for good matches.
DeltaChi2_1sig = chi2inv(0.68268, 1 + 2.*Har);
DeltaChi2_2sig = chi2inv(0.95449, 1 + 2.*Har);
DeltaChi2_3sig = chi2inv(0.99730, 1 + 2.*Har);
DeltaChi2_5sig = chi2inv(0.99999942, 1 + 2.*Har);
fprintf('DeltaChi2 1-sigma: %f\n', DeltaChi2_1sig);
fprintf('DeltaChi2 2-sigma: %f\n', DeltaChi2_2sig);
fprintf('DeltaChi2 3-sigma: %f\n', DeltaChi2_3sig);
fprintf('DeltaChi2 5-sigma: %f\n', DeltaChi2_5sig);

if (figureNum > 0),
    figure(figureNum)
    if (length(Frange) > 1),
       plot(Frange, Chi2vec./NumberOfDataPoints, 'b-', 'LineWidth', 6); hold on;
    elseif (length(Frange) == 1),
      plot(bestF, MinChi2, 'ko');
      dFreq = 0;
   end;
    ylabel('chi square values'); %%hold on; temp 20/02/2006
    xlabel('Frequencies f [rotation per day]'); %% hold off; temp 20/02/2006
    title ('chi square fittness for frequency'); %% hold on; temp 20/02/2006
end
jindex5 = find(Chi2vec < MinChi2 + DeltaChi2_5sig);
jindex3 = find(Chi2vec < MinChi2 + DeltaChi2_3sig);
jindex2 = find(Chi2vec < MinChi2 + DeltaChi2_2sig);
jindex1 = find(Chi2vec < MinChi2 + DeltaChi2_1sig);
if (figureNum > 0),
   plot(Frange(:), ones(size(Frange)).*(MinChi2 + DeltaChi2_1sig)./NumberOfDataPoints, '-k', 'LineWidth', 2);
   plot(Frange(:), ones(size(Frange)).*(MinChi2 + DeltaChi2_2sig)./NumberOfDataPoints, '--k', 'LineWidth', 1);
   plot(Frange(:), ones(size(Frange)).*(MinChi2 + DeltaChi2_3sig)./NumberOfDataPoints, '-.k', 'LineWidth', 1);
   plot(Frange(:), ones(size(Frange)).*(MinChi2 + DeltaChi2_5sig)./NumberOfDataPoints, ':k', 'LineWidth', 1);
   legend('Chi2', '1 \sigma', '2 \sigma', '3 \sigma', '5 \sigma', 3);
end

fprintf('Using 3-sigma to determine error\n');
jindex = jindex3;
if (isempty(jindex) == 0),
   dFreq = max(Frange(jindex)) - min(Frange(jindex));
   if (dFreq == 0),
      dFreq = abs(2.*(max(Frange(jindex+1))-max(Frange(jindex))));
   end
else
   dFreq = Frange(end)-Frange(end-1);
end


% Old:
% Using one-sigma as a limit for good matches.
%DeltaChi2 = chi2inv(0.6804, 1 + 2.*Har);
% Using five-sigma as a limit for good matches.
%DeltaChi2 = chi2inv(0.94, 1 + 2.*Har);
%
%fprintf('DeltaChi2: %f\n', DeltaChi2);
%[j]=find(Chi2vec(:,:)<minChi2 + DeltaChi2);

%if (isempty(j) == 0),
%    dFreq = max(Frange(j)) - min(Frange(j))
%%%     if (dFreq == 0),
%%%        if (length(Frange) > 1),
%%%            dFreq = 2.*(max(Frange(j+1))-max(Frange(j)));
%%%         else,
%%%            dFreq = 0;
%%%         end;
%%%     end;
%else
%   if (length(Frange) > 1),
%      dFreq = Frange(end)-Frange(end-1);
%   else
%      dFreq = 0;
%   end;
%end;

[Freq_, dFreq_] = writePrintableErr(bestF, dFreq);
fprintf('Freq is %s\n', num2str(Freq_));
fprintf('dFreq is %s\n', num2str(dFreq_));
Period = 24./bestF;
dPeriod = abs(-24.*dFreq./bestF.^2);
[Period_, dPeriod_] = writePrintableErr(Period, dPeriod);
fprintf('Period is %s\n', num2str(Period_));
fprintf('dPeriod is %s\n', num2str(dPeriod_));

% Plot Chi2 result
if (figureNum > 0),
%   figure(figureNum+10)
%   if (length(Frange) > 1),
%      plot(Frange, Chi2vec, 'b-'); hold on;
%      ylabel('chi square values'); %%hold on; temp 20/02/2006
%      xlabel('Frequencies f [rotation per day]'); %% hold off; temp 20/02/2006
%      title ('chi square fittness for frequency'); %% hold on; temp 20/02/2006
%   elseif (length(Frange) == 1),
%      plot(bestF, minChi2, 'ko');
%      dFreq = 0;
%      ylabel('min Chi square'); %% hold on; temp 20/02/2006
%      xlabel('Frequencies f [rotation per day]'); %% hold off; temp 20/02/2006
%      title ('Frequency and min Chi square'); %% hold on; temp 20/02/2006
%   end;

   % Plot Folded lightcurve by the frequency with legend
   % Translate JD to date
   Inx = 1;
   [date] = jd2date(time(1));
   nights(Inx) = round(date(3)*10000+date(2)*100+date(1)+date(4))-1;
   beginInx(Inx) = Inx;
   for i=2:1:length(time),
       [date] = jd2date(time(i));
       day = round(date(3)*10000+date(2)*100+date(1)+date(4))-1;
       if (nights(Inx) ~= day),
           endInx(Inx) = i-1;
           Inx = Inx +1;
           nights(Inx) = day;
           beginInx(Inx) = i;
       end;
   end;
   endInx(Inx) = i;
   numNights = length(nights);
   colorMat=jet;
   idx = 1 : floor( size( colorMat, 1 ) / numNights ) : size( colorMat, 1 );
   plotMark = ['+','p','*','d','x','o','s','v','^','<','>','h','.','+','p','*','d','x','o','s','v','^','<','>','h','.','+','p','*','d','x','o','s','v','^','<','>','h','.'];
   plotColor = ['k','r','b','g','m','y','k','b','r','g','m','y','k','b','r','g','m','y','k','b','r','g','m','y','k','b','r','g','m','y','k','b','r','g','m','y','k','b','r','g','m'];
   figure(figureNum+1);
   % Plot the graph with plot just for having a nice legend
   for i=1:1:numNights,
       F = folding([time(beginInx(i):endInx(i))-t0,Mag(beginInx(i):endInx(i)),MagErr(beginInx(i):endInx(i))],1./bestF);
   %    ColorMark = strcat(plotColor(i), plotMark(i));
   %    plot(F(:,1),F(:,2), ColorMark); hold on;

       h = plot(F(:,1).*24./bestF,F(:,2), plotMark(i)); hold on;
%       set( h, 'color', colorMat(idx(i),:));
       set( h, 'color', plotColor(i));
       nightsNames(i,:) = num2str(nights(i));
   end;
   legend(nightsNames, 3);
   axis ij;
%   xlabelStr = strcat('synodic phase  [epoch = ',num2str(t0), '].  Period = ',num2str(24./bestF),' hours');
   xlabelStr = strcat('Time [hours].  epoch = ',num2str(t0), '.  Period = ',num2str(24./bestF),' hours');
   xlabel(xlabelStr); hold on;
   ylabel('arbitrary mag [1,0]'); hold on;

   % Plot Folded lightcurve by the frequency
   F = folding([time-t0,Mag,MagErr],1./bestF);
   figure(figureNum+2); hold on; xlim([0 1.2]);
   errorxy([F(:,1),F(:,2),MagErr], 'EdgeColor', 'r', 'FaceColor','r','Marker', '.'); hold on;
   errorxy([F(:,1)+1,F(:,2),MagErr], 'EdgeColor', 'r', 'FaceColor','r','Marker', '.'); hold on;
   axis ij;
   xlabelStr = strcat('synodic phase  [epoch = ',num2str(t0), '].  Period = ',num2str(24./bestF),' hours');
   xlabel(xlabelStr); hold on;
   ylabel('arbitrary mag [1,0]'); hold on;

end;

toc