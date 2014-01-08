function [bestG,dG,bestF,dFreq,Chi2mat] = fitFGbyHarmo(Grange, phi1, phi2, Frange, Har, Deg, time, Mag, MagErr, cosmicErr, plotNum);

Chi2mat = zeros(length(Grange),length(Frange));
minFindex = zeros(length(Grange));

for i=1:1:length(Grange)
   MagG = Mag + 2.5.*log10((1-Grange(i)).*phi1 + Grange(i).*phi2);
   for j=1:1:length(Frange)
      [Par,Par_Err,Cov,Chi2,Freedom,Par1,Resid]=fitharmo(time, MagG, MagErr.*cosmicErr, [Frange(j), Har], Deg);
      if(isreal(Chi2) == 1),
         Chi2mat(i,j) = Chi2;
      else
         Chi2mat(i,j) = nan;
      end
   end
   [mintemp, FreqInx] = min(Chi2mat(i,:));
   % Take here the minimal frquency and print it at the next line
   fprintf('mintemp is %6.6f of freq: %6.5f. Now at checking %d out of %d\n', mintemp, Frange(FreqInx), (i-1).*j+j, length(Grange).*length(Frange));
   [minChi2F(i), minFindex(i)] = min(Chi2mat(i,:));
end
[minChi2G, minGindex] = min(minChi2F);
fprintf('minChi2k is %f at minkIndex num. %d\n', minChi2G, minGindex);
Chi2matMin = Chi2mat(minGindex,minFindex(minGindex));
bestG = Grange(minGindex);
bestF = Frange(minFindex(minGindex));

MinChi2=min(min(Chi2mat));
fprintf('MinChi2 is %f\n', MinChi2);
if (plotNum > 0),
    figure(plotNum);
end;

% Using one-sigma as a limit for good matches.
DeltaChi2_1sig = chi2inv(0.68268, 1 + 2.*Har);
DeltaChi2_2sig = chi2inv(0.95449, 1 + 2.*Har);
DeltaChi2_3sig = chi2inv(0.99730, 1 + 2.*Har);
DeltaChi2_5sig = chi2inv(0.99999942, 1 + 2.*Har);
fprintf('DeltaChi2 1-sigma: %f\n', DeltaChi2_1sig);
fprintf('DeltaChi2 2-sigma: %f\n', DeltaChi2_2sig);
fprintf('DeltaChi2 3-sigma: %f\n', DeltaChi2_3sig);
fprintf('DeltaChi2 5-sigma: %f\n', DeltaChi2_5sig);

fprintf('Using 1-sigma to determine error\n');

if (length(Grange) > 1 & length(Frange) > 1),
%   contour(Frange,Grange,Chi2mat,MinChi2+[0;2.2958;6.1801;11.829]) % original line - temp 13.3.2006
   if (plotNum > 0),
% OLD VALUES
%       contour(Frange,Grange,Chi2mat,MinChi2+[0], 'k-'); hold on;
%       contour(Frange,Grange,Chi2mat,MinChi2+[2.2958], 'k-'); hold on;
%       contour(Frange,Grange,Chi2mat,MinChi2+[6.1801], 'k-.'); hold on;
%       contour(Frange,Grange,Chi2mat,MinChi2+[11.829], 'k:'); hold on;
       contour(Frange,Grange,Chi2mat,MinChi2+[0], 'k.'); hold on;
       contour(Frange,Grange,Chi2mat,MinChi2+[DeltaChi2_1sig], 'k-'); hold on;
       contour(Frange,Grange,Chi2mat,MinChi2+[DeltaChi2_2sig], 'k-.'); hold on;
       contour(Frange,Grange,Chi2mat,MinChi2+[DeltaChi2_3sig], 'k:'); hold on;
       contour(Frange,Grange,Chi2mat,MinChi2+[DeltaChi2_5sig], 'k:'); hold on;
       hold off;
       colorbar;
       ylabel('Slope Parameter G'); hold on;
       xlabel('Frequencies f [rotation per day]'); %% hold off; temp 20/02/2006
       title ('chi square surface for frequency and G parameter'); %% hold on; temp 20/02/2006
   end;
   % Old condition:
   %[i,j]=find(MinChi2<Chi2mat(:,:) & Chi2mat(:,:)<MinChi2+2.2958);
   fprintf('Using 1 sigma to determine error.\n');
   [i,j]=find(Chi2mat(:,:)<MinChi2 + DeltaChi2_1sig);
   if (isempty(i) == 0),
      dG = max(Grange(i)) - min(Grange(i))
%      if (dG == 0),
%          if (Grange(i) ~= Grange(end)),
%              dG = 2.*(max(Grange(i+1))-max(Grange(i)));
%          else
%              dG = 10.*(max(Grange) - min(Grange));
%          end;
%      end
   else
       dG = Grange(end)-Grange(end-1);
   end;
   if (isempty(j) == 0),
       dFreq = max(Frange(j)) - min(Frange(j));
       if (dFreq == 0),
           dFreq = 2.*(max(Frange(j+1))-max(Frange(j)));
       end;
   else
       dFreq = Frange(end)-Frange(end-1);
   end;
elseif (length(Grange) > 1 & length(Frange) == 1),
   index = find(MinChi2<Chi2mat & Chi2mat<MinChi2+1);
   if (isempty(index) == 0),
      dG = max(Grange(index)) - min(Grange(index));
      if (dG == 0),
         dG = abs(2.*(max(Grange(index+1))-max(Grange(index))));
      end
   else
      dG = Grange(end)-Grange(end-1);
   end
   dFreq = 0;
   if (plotNum > 0),
       plot(Grange, Chi2mat(:,1), 'ro');
       xlabel('Slope Parameter G'); %% hold on; temp 20/02/2006
       ylabel('chi square values'); %% hold on; temp 20/02/2006
       title ('chi square fittness for G parameter'); %% hold on; temp 20/02/2006
   end
elseif (length(Grange) == 1 & length(Frange) > 1),
   if (plotNum > 0),
       plot(Frange, Chi2mat(1,:), 'b-', 'LineWidth', 6); hold on;
       ylabel('chi square values'); %%hold on; temp 20/02/2006
       xlabel('Frequencies f [rotation per day]'); %% hold off; temp 20/02/2006
       title ('chi square fittness for frequency'); %% hold on; temp 20/02/2006
   end
   % Old version
%   jindex = find(MinChi2<Chi2mat & Chi2mat<MinChi2+1);
%   size(jindex)
   jindex5 = find(Chi2mat(:,:) < MinChi2 + DeltaChi2_5sig);
   jindex3 = find(Chi2mat(:,:) < MinChi2 + DeltaChi2_3sig);
   jindex2 = find(Chi2mat(:,:) < MinChi2 + DeltaChi2_2sig);
   jindex1 = find(Chi2mat(:,:) < MinChi2 + DeltaChi2_1sig);
   if (plotNum > 0),
%      plot(Frange(jindex5), Chi2mat(1,jindex5), 'LineWidth', 4, 'Color', 'k'); hold on;
%      plot(Frange(jindex3), Chi2mat(1,jindex3), 'LineWidth', 3, 'Color', 'g'); hold on;
%      plot(Frange(jindex2), Chi2mat(1,jindex2), 'LineWidth', 2, 'Color', 'r'); hold on;
%      plot(Frange(jindex1), Chi2mat(1,jindex1), 'LineWidth', 8, 'Color', 'c'); hold on;

      plot(Frange(:), ones(size(Frange)).*(MinChi2 + DeltaChi2_1sig), '-k', 'LineWidth', 2);
      plot(Frange(:), ones(size(Frange)).*(MinChi2 + DeltaChi2_2sig), '--k', 'LineWidth', 1);
      plot(Frange(:), ones(size(Frange)).*(MinChi2 + DeltaChi2_3sig), '-.k', 'LineWidth', 1);
      plot(Frange(:), ones(size(Frange)).*(MinChi2 + DeltaChi2_5sig), ':k', 'LineWidth', 1);
%      legend('Chi2', '5 \sigma', '3 \sigma', '2 \sigma', '1 \sigma', 3);
      legend('Chi2', '1 \sigma', '2 \sigma', '3 \sigma', '5 \sigma', 3);
   end
%   size(jindex)
   jindex = jindex1;
   if (isempty(jindex) == 0),
      dFreq = max(Frange(jindex)) - min(Frange(jindex));
      if (dFreq == 0),
         dFreq = abs(2.*(max(Frange(jindex+1))-max(Frange(jindex))));
      end
   else
      dFreq = Frange(end)-Frange(end-1);
   end
   dG = 0;
elseif (length(Grange) == 1 & length(Frange) == 1),
    dG = 0;
    dFreq = 0;
   if (plotNum > 0),
       plot(bestG,bestF, 'ko');
       ylabel('Slope Parameter G'); %% hold on; temp 20/02/2006
       xlabel('Frequencies f [rotation per day]'); %% hold off; temp 20/02/2006
       title ('Frequency and G parameter'); %% hold on; temp 20/02/2006
   end
end

