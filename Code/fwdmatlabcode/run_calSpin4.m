function [SuperArcs] = run_calSpin4(Field, Chip, Date, AstName, SuperArcs, FigureNum)

Def.Chip = '';
Def.Date = '';
Def.AstName = '';
Def.SuperArcs = [];
Def.FigureNum = 0;

if (nargin==1),
   Chip = Def.Chip;
   Date = Def.Date;
   AstName = Def.AstName;
   SuperArcs = Def.SuperArcs;
   FigureNum = Def.FigureNum;
elseif (nargin==2),
   Date = Def.Date;
   AstName = Def.AstName;
   SuperArcs = Def.SuperArcs;
   FigureNum = Def.FigureNum;
elseif (nargin==3),
   AstName = Def.AstName;
   SuperArcs = Def.SuperArcs;
   FigureNum = Def.FigureNum;
elseif (nargin==4),
   SuperArcs = Def.SuperArcs;
   FigureNum = Def.FigureNum;
elseif (nargin==5),
   FigureNum = Def.FigureNum;
elseif (nargin==6),
   % do nothing
else
   error('Illegal number of input arguments');
end

% Select relevant files
if (isempty(SuperArcs)),
   [SuperArcs] = CreateSuperArcs(Field, Chip, Date, FigureNum, []);
end;
% Stick with a single specific object if needed
SuperArcsOrig = SuperArcs;
IsAstFound = 0;
if (strcmp(AstName, '') == 0),
   for Isarc=1:1:length(SuperArcs),
      if (isfield(SuperArcs{Isarc},'Properties')),
         if(~isempty(findstr(SuperArcs{Isarc}.Properties.Name,AstName))),
            SuperArcsName{1} = SuperArcs{Isarc};
            SuperArcs = SuperArcsName;
            fprintf('Arc location in cell array: %d\n', Isarc);
            IsAstFound = 1;
            break;
         end;
      end;
   end;
   if (~IsAstFound),
      % Try and search for the PTFdesig
      for Isarc=1:1:length(SuperArcs),
         if (IsAstFound),
            break;
         end;
         for Iarc=1:1:length(SuperArcs{Isarc}.Arc),
            if (isfield(SuperArcs{Isarc}.Arc{Iarc},'PTFdesig')),
               if(~isempty(findstr(SuperArcs{Isarc}.Arc{Iarc}.PTFdesig, AstName))),
                  SuperArcsName{1} = SuperArcs{Isarc};
                  SuperArcs = SuperArcsName;
                  fprintf('Arc location in cell array: %d\n', Isarc);
                  IsAstFound = 1;
                  break;
               end;
            end;
         end;
      end;
   end;
   if (~IsAstFound),
      fprintf ('%s does not exist within SuperArcs\n', AstName);
      return;
   end;
end;

% FrangeMax = 24./1.5;                % Maximal frequency is 1.5 hours - due to the rubble pile spin barrier
% FrangeInterval= 0.001;
Har = 2; % How to determine the number of harmonies?
Deg = 0; % How to determine the order?
MaxRadNeighbor = 5; % arcsec
MinImageNum = 8;

%MarkerVec = ['o', 's', '^', 'p', 'd', 'h', 'v', '<', '>'];

for Isarc=1:1:length(SuperArcs),
   
   fprintf('Now at arc number: %d\n', Isarc);
   SuperArcs{Isarc}.Spin.Status = 'running spin analysis';
   % Determine the name to be shown in the legend
   legendVec{Isarc} = getArcTitleForLegend(SuperArcs{Isarc}.Arc, 1);

   AstNameStrForFiles = strrep(legendVec{Isarc}, ' ', '_');
   AstNameStrForFiles = strrep(AstNameStrForFiles, '(', '');
   AstNameStrForFiles = strrep(AstNameStrForFiles, ')', '');
   Fid = fopen(strcat('SpinAnalysisLog_', AstNameStrForFiles),'w');
   fprintf(Fid, '-------------------------------------------------------------------------------------\n');
   fprintf(Fid, '%s\n', legendVec{Isarc});
   fprintf('-------------------------------------------------------------------------------------\n');
   fprintf('%s\n', legendVec{Isarc});

   % Prepare the measurements for spin analysis
   Time = []; Mag = []; MagErr = []; ChipVec = [];
   MagRange = []; MagErrMed = []; MagRange = 0.9;
   
%   length(SuperArcs{Isarc}.Arc)
%   for Iarc=1:1:length(SuperArcs{Isarc}.Arc),%
%      SuperArcs{Isarc}.Arc{Iarc}%
%   end;%
   ChipInx = 0;
   for Iarc=1:1:length(SuperArcs{Isarc}.Arc),
      if (isfield(SuperArcs{Isarc}.Arc{Iarc},'ZP')),
         % DON'T USE MEASUREMENTS WITH NEIGHBORING STARS !!! SP 20110505 TEMP
         InxNoNeigh = find(SuperArcs{Isarc}.Arc{Iarc}.NeighDist.*180./pi*3600 > MaxRadNeighbor);

         if (isfield(SuperArcs{Isarc}.Arc{Iarc}, 'JD_LightTimeCorrected')),
            VecLen = length(SuperArcs{Isarc}.Arc{Iarc}.JD_LightTimeCorrected(InxNoNeigh));
            Time(end+1:end+VecLen) = (SuperArcs{Isarc}.Arc{Iarc}.JD_LightTimeCorrected(InxNoNeigh)');
         else
            VecLen = length(SuperArcs{Isarc}.Arc{Iarc}.JD(InxNoNeigh));
            Time(end+1:end+VecLen) = (SuperArcs{Isarc}.Arc{Iarc}.JD(InxNoNeigh)');
         end;
         
         if (isfield(SuperArcs{Isarc}.Arc{Iarc}.Mag, 'MagCorrectionTo1AU')),
            Mag(end+1:end+VecLen) = SuperArcs{Isarc}.Arc{Iarc}.Mag.Auto.Mag(InxNoNeigh)' ...
                                  - SuperArcs{Isarc}.Arc{Iarc}.ZP(InxNoNeigh)' ...
                                  + SuperArcs{Isarc}.Arc{Iarc}.Mag.MagCorrectionTo1AU(InxNoNeigh)';
                              
            MagVec                = SuperArcs{Isarc}.Arc{Iarc}.Mag.Auto.Mag(InxNoNeigh)' ...
                                  - SuperArcs{Isarc}.Arc{Iarc}.ZP(InxNoNeigh)' ...
                                  + SuperArcs{Isarc}.Arc{Iarc}.Mag.MagCorrectionTo1AU(InxNoNeigh)';
         else
            Mag(end+1:end+VecLen) = SuperArcs{Isarc}.Arc{Iarc}.Mag.Auto.Mag(InxNoNeigh)' ...
                                  - SuperArcs{Isarc}.Arc{Iarc}.ZP(InxNoNeigh)';

            MagVec                = SuperArcs{Isarc}.Arc{Iarc}.Mag.Auto.Mag(InxNoNeigh)' ...
                                  - SuperArcs{Isarc}.Arc{Iarc}.ZP(InxNoNeigh)';
         end;

         Range = err_cl(MagVec, MagRange);
         MagRange(Iarc) = Range(2)-Range(1);


         MagErr(end+1:end+VecLen) = sqrt((SuperArcs{Isarc}.Arc{Iarc}.Mag.Auto.Err(InxNoNeigh).^2)' + (SuperArcs{Isarc}.Arc{Iarc}.ErrZP(InxNoNeigh).^2)');
         
         MagErrVec = sqrt((SuperArcs{Isarc}.Arc{Iarc}.Mag.Auto.Err(InxNoNeigh).^2)' + (SuperArcs{Isarc}.Arc{Iarc}.ErrZP(InxNoNeigh).^2)');
         MagErrVec = MagErrVec(find(~isnan(MagErrVec)));
         MagErrMed(Iarc) = median(MagErrVec);


         ChipInx = ChipInx + 1;
         ChipVec(end+1:end+VecLen) = ChipInx;

%          VecLen = length(SuperArcs{Isarc}.Arc{Iarc}.JD_LightTimeCorrected)
%          Time(end+1:end+VecLen) = (SuperArcs{Isarc}.Arc{Iarc}.JD_LightTimeCorrected');
%          Mag(end+1:end+VecLen) = SuperArcs{Isarc}.Arc{Iarc}.Mag.Auto.Mag' ...
%                                - SuperArcs{Isarc}.Arc{Iarc}.ZP' ...
%                                + SuperArcs{Isarc}.Arc{Iarc}.Mag.MagCorrectionTo1AU';
%          MagErr(end+1:end+VecLen) = sqrt((SuperArcs{Isarc}.Arc{Iarc}.Mag.Auto.Err.^2)' + (SuperArcs{Isarc}.Arc{Iarc}.ErrZP.^2)');

      end
   end
   
   % Remove mag/magerr with NaN
   Time = Time(find(~isnan(Mag(:))))';
   MagErr = MagErr(find(~isnan(Mag(:))))';
   Mag = Mag(find(~isnan(Mag(:))))';
   ChipVec = ChipVec(find(~isnan(Mag(:))))';

   if (isempty(Time) | isempty(Mag) | isempty(MagErr) | isempty(ChipVec)),
      fprintf(Fid, '%s has no valid measurements\n', legendVec{Isarc});
      fprintf('%s has no valid measurements\n', legendVec{Isarc});
      SuperArcs{Isarc}.Spin.Status = 'failed spin analysis';
      continue;
   end

   if (length(Time) <= MinImageNum),
      fprintf(Fid, '%s has too few measurements: %d\n', legendVec{Isarc}, length(Time));
      fprintf('%s has too few measurements: %d\n', legendVec{Isarc}, length(Time));
      SuperArcs{Isarc}.Spin.Status = 'failed spin analysis';
      continue;
   end


   MagVal = median(MagRange);
   MagErrMed = median(MagErrMed);
   if (MagErrMed < MagVal),
      AmpMin = sqrt(MagVal.^2 - MagErrMed.^2);
   else
      AmpMin = 0;
   end
   SuperArcs{Isarc}.Spin.AmpMin = round(AmpMin.*100)./100;

   
   
   % The minimal frequency is the overall observed time of the asteroid.
   FrangeMin = floor(1./(max(Time) - min(Time)).*10000)./10000;
   T0 = floor(min(Time));
   % The maximal frequency is the time difference between images
   TimeSorted = sort(Time);
   TimeMin = TimeSorted(2)-TimeSorted(1);
   for Iimage=3:1:length(TimeSorted),
      if (TimeSorted(Iimage)-TimeSorted(Iimage-1) < TimeMin & TimeSorted(Iimage)-TimeSorted(Iimage-1) ~= 0),
         TimeMin = TimeSorted(Iimage)-TimeSorted(Iimage-1);
      end;
   end;
   FrangeMax = ceil(1./TimeMin./2./4);                     % divided by 2 to get bi-model lightcurve
                                                          % divided by 4 to get 4 points in the lightcurve.

   % Frequncy interval is same as Frange Max divided by 4 to smooth the power spectrum
   FrangeInterval = floor(TimeMin./4.*10000)./10000;

   if (FrangeMax <= FrangeMin),
      fprintf(Fid, '%s has too few measurements: %d\n', legendVec{Isarc}, length(Time));
      fprintf('%s has too few measurements: %d\n', legendVec{Isarc}, length(Time));
      SuperArcs{Isarc}.Spin.Status = 'failed spin analysis';
      continue;
   end;


   % Search for the spin by power spectrum (periodis). This is a fast
   % search over a wide range of frequencies, with one order of magnitude.
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   [FreqPS, PeaksPS] = periodis([Time-T0, Mag, MagErr], FrangeMin, FrangeMax, FrangeInterval);
   % Check and plot
   FavFreqs = [];
   if (isempty(PeaksPS) ~= 1),
      fprintf(Fid, '~~~~~~~~~~~~~~~ periodis ~~~~~~~~~~~~~~~\n');
      fprintf('~~~~~~~~~~~~~~~ periodis ~~~~~~~~~~~~~~~\n');

      ProbabilityLimit = 0.8;
      FavFreqs = PeaksPS(find(PeaksPS(:,4) > ProbabilityLimit & PeaksPS(:,1) > 0), :);
      if (isempty(FavFreqs)),
         FavFreqs = PeaksPS(end, :); % Frequency
         % or leave this asteroid alone!
         fprintf(Fid, '%s has no solution. Highest probability: %5.3f\n', legendVec{Isarc}, PeaksPS(end,4));
         fprintf('%s has no solution. Highest probability: %5.3f\n', legendVec{Isarc}, PeaksPS(end,4));
         SuperArcs{Isarc}.Spin.Status = 'periodis failed';
         continue;
%%%%%%%%%%         continue; % but plot and check each case and see if some good lightcurves are not thrown away - DP
      end;

      [NumFavFreq, Dummy] = size(FavFreqs);
      fprintf(Fid, 'Per [hou]   Freq [1/day] power     probability\n');
      fprintf('Per [hou]   Freq [1/day] power     probability\n');
      for IfavFreq=NumFavFreq:-1:1,
         fprintf(Fid, '%f   %f     %f %f\n', 24.*FavFreqs(IfavFreq, 3), FavFreqs(IfavFreq, 1), FavFreqs(IfavFreq, 2), FavFreqs(IfavFreq, 4));
         fprintf('%f   %f     %f %f\n', 24.*FavFreqs(IfavFreq, 3), FavFreqs(IfavFreq, 1), FavFreqs(IfavFreq, 2), FavFreqs(IfavFreq, 4));
      end;

%       Period = 24./PeaksPS(end,1);
%       probability = PeaksPS(end,4);
%       [Period_, probability_] = writePrintableErr(Period, probability);

%       if (FigureNum),
%          % Plot the power spectrum
%          figure;
%          plot (24./FreqPS(:,1),FreqPS(:,2),'k'); hold on; grid on;
%          graphTitle = ['Power Spectrum of ',legendVec{Isarc}, '. Best Period = ',num2str(24./PeaksPS(end,1)),' hours']; title (graphTitle); hold on;
%          ylabelText = ['power']; ylabel(ylabelText); hold on;
%          xlabelStr = ['period [hours]']; xlabel(xlabelStr); hold off;
%          set(gcf,'position',[1120 510 560 180]);
% 
%          % Plot the folded lightcurve
%          F = folding([Time-T0,Mag,MagErr],1./PeaksPS(end,1));
%          figure;
%          errorxy([(24./PeaksPS(end,1)).*F(:,1), F(:,2), F(:,3)],'EdgeColor','k','FaceColor','k'); hold on;
% 
%          % Plot the model on the folded lightcurve
%          FreqPhase = 1;
%          [Par,Par_Err,Cov,Chi2,Freedom,Par1,Resid]=fitharmo(F(:,1), F(:,2), MagErr, [FreqPhase, Har],Deg);
%          X = [0:0.01:F(end,1)];
%          Y = getModel(X, FreqPhase, Par, Har, Deg);
%          plot((24./PeaksPS(end,1)).*X,Y,'k'); hold on;
%          axis ij;
%          xTitle = ['Time [hours].  epoch = ',num2str(T0),'.  Period = ',num2str(Period_),' hours'];
%          yLabelStr = ['absolute mag'];
%          xlabel(xTitle); ylabel(yLabelStr); hold on; hold off;
%          graphTitle = ['Folded lightcurve of ', legendVec{Isarc}];   
%          title (graphTitle); hold on;
%       end;

   end;

   % Search for the spin by Furrier series (fitharmo) %
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   fprintf(Fid, '~~~~~~~~~~~~~~~ fitharmo ~~~~~~~~~~~~~~\n');
   fprintf('~~~~~~~~~~~~~~~ fitharmo ~~~~~~~~~~~~~~\n');
   
   % Calculate the delta chi2 that are coresponding to one standard
   % deviation. This should include the errors of all the free parameters
   DeltaChi2 = chi2inv(0.6804, 1 + 2.*Har + length(unique(ChipVec)))

   % Arrange the checked frequencies - SHOULD IT BE REMOVED ??? 2011/Apr/26 DP 
   Frange = [];
%    FavFreqEdge = 0.5;
%    if (~isempty(FavFreqs)),
%       for IfavFreq=1:1:NumFavFreq,
%          Frange(end+1,:) = [FavFreqs(IfavFreq)-FavFreqEdge:FrangeInterval./10:FavFreqs(IfavFreq)+FavFreqEdge];
%          % Also check the double period (half the frequency) due to the
%          % double-peak lightcurve
%          FrangeLen = length(Frange);
%          Frange(end,end+1:FrangeLen.*2) = [FavFreqs(IfavFreq)./2-FavFreqEdge:FrangeInterval./10:FavFreqs(IfavFreq)./2+FavFreqEdge];
%          Frange = sort(Frange);
%       end;
%    else
%       Frange = [FrangeMin:FrangeInterval:FrangeMax];
%    end;
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   Frange = [FrangeMin:FrangeInterval:FrangeMax];

   % Run the fit
   CosmicErr = 1; % ???
   [Chi2vec] = PTF_fitFreqByHarmo_fitMagShift(Frange, Har, Deg, Time-T0, Mag, MagErr, ChipVec, CosmicErr);
   %%[Chi2vec] = PTF_fitFreqByHarmo2(Frange, Har, Deg, Time-T0, Mag, MagErr, CosmicErr);
   [InxGoodChi]=find(Chi2vec(:,:)<min(Chi2vec)+DeltaChi2);

   % Check and plot

   % If the best freq is at the edge of the freq range than increase the
   % range and run again
   while (FrangeMin > 1./100 | FrangeMax < 24./(1/60)),
      if (length(Frange) == InxGoodChi(end)),
         fprintf('--->>> The best freq is at the MAXIMAL edge of the freq range than increase the\n');
         FrangeMax = FrangeMax.*2;
         Frange = [FrangeMin:FrangeInterval:FrangeMax];
         [Chi2vec] = PTF_fitFreqByHarmo_fitMagShift(Frange, Har, Deg, Time-T0, Mag, MagErr, ChipVec, CosmicErr);
         %%[Chi2vec] = PTF_fitFreqByHarmo2(Frange, Har, Deg, Time-T0, Mag, MagErr, CosmicErr);
         [InxGoodChi]=find(Chi2vec(:,:)<min(Chi2vec)+DeltaChi2);
      elseif (InxGoodChi(1) == 1),
%      elseif (BestFreq < min(Frange)+FrangeInterval.*10),
         fprintf('--->>> The best freq is at the MINIMAL edge of the freq range than increase the\n');
         FrangeMin = FrangeMin./2;
         Frange = [FrangeMin:FrangeInterval:FrangeMax];
         [Chi2vec] = PTF_fitFreqByHarmo_fitMagShift(Frange, Har, Deg, Time-T0, Mag, MagErr, ChipVec, CosmicErr);
         %%[Chi2vec] = PTF_fitFreqByHarmo2(Frange, Har, Deg, Time-T0, Mag, MagErr, CosmicErr);
         [InxGoodChi]=find(Chi2vec(:,:)<min(Chi2vec)+DeltaChi2);
      else
         [minChi2, InxMinChi] = min(Chi2vec);
         BestFreq = Frange(InxMinChi);
         break;
      end
   end

   if (FrangeMin < 1./100 | FrangeMax > 24./(1/60)),
      fprintf(Fid, 'The found frequency for %s is not physical. fitharmo failed\n', legendVec{Isarc});
      fprintf('The found frequency for %s is not physical. fitharmo failed\n', legendVec{Isarc});
      if (strcmp(SuperArcs{Isarc}.Spin.Status, 'periodis failed')),
         SuperArcs{Isarc}.Spin.Status = 'failed spin analysis';
         continue;
      end
      SuperArcs{Isarc}.Spin.Status = 'fitharmo failed';
      if (24.*(max(Time)-min(Time)) < 24./PeaksPS(end,1).*2),
         fprintf(Fid, 'The time range is too short to give a solid spin value for %s. Abort\n', legendVec{Isarc});
         fprintf('The time range is too short to give a solid spin value for %s. Abort\n', legendVec{Isarc});
         continue;
      end
   end

   if (~strcmp(SuperArcs{Isarc}.Spin.Status, 'fitharmo failed')),
      % Determine the frequency error
      FreqSolutions = []; dFreq = [];

      % Make the frequenct range more dense so it will be possible to get the extramum
      if (length(Chi2vec(InxGoodChi)) < 10),
         Frange = [FrangeMin:FrangeInterval./10:FrangeMax];
         [Chi2vec] = PTF_fitFreqByHarmo_fitMagShift(Frange, Har, Deg, Time-T0, Mag, MagErr, ChipVec, CosmicErr);
         %%[Chi2vec] = PTF_fitFreqByHarmo2(Frange, Har, Deg, Time-T0, Mag, MagErr, CosmicErr);
         [InxGoodChi]=find(Chi2vec(:,:)<min(Chi2vec)+DeltaChi2);
      end;
      ExtramChi = find_local_extramum(Frange(InxGoodChi)', Chi2vec(InxGoodChi));

      % The frequency range is not completely checked - this should not happen
      if (isempty(ExtramChi)),
         dFreq = 0;
      % Only one solution is better significantly than all of the others
      elseif(length(find(ExtramChi(:,3)>0)) == 1),
         dFreq = max(Frange(InxGoodChi)) - min(Frange(InxGoodChi));
      % There are more than one possible solutions
      elseif (length(find(ExtramChi(:,3)>0)) > 1),
         fprintf(Fid, 'There are more than one possible solutions\n');
         fprintf('There are more than one possible solutions\n');
         if (strcmp(SuperArcs{Isarc}.Spin.Status, 'periodis failed')),
            fprintf(Fid, 'periodis and fitharmo have failed. Abort\n');
            fprintf('periodis and fitharmo have failed. Abort\n');
            SuperArcs{Isarc}.Spin.Status = 'failed spin analysis';
            continue;
         end;
         fprintf(Fid, 'Frequency      Error      Chi2\n');
         fprintf('Frequency  Error      Chi2\n');
         FreqErrorPosMax = []; FreqErrorPosMin = [];
         FreqSolutions = ExtramChi(find(ExtramChi(:,3)>0),:);
         [NumSolutions, dummy] = size(FreqSolutions);
         % If there are too many solutions - there is no solution...
         if (NumSolutions > 5),
            fprintf(Fid, '%s: There are too many solutions -> so there is no solution. fitharmo failed\n', legendVec{Isarc});
            fprintf('%s: There are too many solutions -> so there is no solution. fitharmo failed\n', legendVec{Isarc});
            SuperArcs{Isarc}.Spin.Status = 'fitharmo failed';
         end;
         SuperArcs{Isarc}.Spin.Status = 'many solutions from fitharmo';

         for Ifreq=1:1:NumSolutions,
            % Get the index of the frequency solution within the Frange vector
            FreqErrorPosMax(Ifreq) = 0; FreqErrorPosMin(Ifreq) = 0;
            FreqSolutionInx = find(abs(Frange - FreqSolutions(Ifreq,1)) < 0.005);
            for Inx=FreqSolutionInx(1):1:length(Frange),
               if (Chi2vec(Inx) > Chi2vec(FreqSolutionInx(1))+DeltaChi2), 
                  FreqErrorPosMax(Ifreq) = Inx;
                  break;
               end;
            end;
            for Inx=FreqSolutionInx(1):-1:1,
               if (Chi2vec(Inx) > Chi2vec(FreqSolutionInx(1))+DeltaChi2),
                  FreqErrorPosMin(Ifreq) = Inx;
                  break;
               end;
            end;
            if (FreqErrorPosMax(Ifreq) == 0 | FreqErrorPosMin(Ifreq) == 0),
               dFreq(Ifreq) = 0;
            else
               dFreq(Ifreq) = Frange(FreqErrorPosMax(Ifreq)) - Frange(FreqErrorPosMin(Ifreq));
               fprintf(Fid, '%6.4f     %6.4f   %6.4f \n', FreqSolutions(Ifreq, 1), dFreq(Ifreq), FreqSolutions(Ifreq, 2));
               fprintf('%6.4f     %6.4f   %6.4f \n', FreqSolutions(Ifreq, 1), dFreq(Ifreq), FreqSolutions(Ifreq, 2));
%               foldCurvesByFreqs(legendVec{Isarc}, SuperArcs, [FreqSolutions(Ifreq, 1)]);
            end;
         end;
      elseif(length(find(ExtramChi(:,3)>0)) < 1),
         % Problem - no minima found - can it happend?
         fprintf(Fid, 'No minima found with fitharmo\n');
         fprintf('No minima found with fitharmo\n');
         dFreq = 0;
      end;

      if (dFreq == 0),
         if (strcmp(SuperArcs{Isarc}.Spin.Status, 'periodis failed')),
            SuperArcs{Isarc}.Spin.Status = 'failed spin analysis';
            continue;
         end;
         SuperArcs{Isarc}.Spin.Status = 'fitharmo failed';
         fprintf(Fid, 'There is a problem with the frequency error: is the frequency value worthless ?\n');
         fprintf(Fid, 'Best frequency is %f +- %f\n', BestFreq, dFreq);
         fprintf(Fid, 'Period is %f +- %f\n', Period_, dPeriod_);
         fprintf('There is a problem with the frequency error: is the frequency value worthless ?\n');
         fprintf('Best frequency is %f +- %f\n', BestFreq, dFreq);
         fprintf('Period is %f +- %f\n', Period_, dPeriod_);
      end;

      % If the period is not bi-model - devide the freq by two.
      F = folding([Time-T0,Mag,MagErr,ChipVec],1./BestFreq);
      BestFreqOrig = BestFreq;
      FreqPhase = 1;
      [Par,Par_Err,Cov,Chi2,Freedom,Par1,Resid]=fitharmo_magShift(F(:,1), F(:,2), F(:,3), F(:,4), [FreqPhase, Har],Deg);
      %%[Par,Par_Err,Cov,Chi2,Freedom,Par1,Resid]=fitharmo(F(:,1), F(:,2), F(:,3), [FreqPhase, Har],Deg);
      X = [0:0.001:F(end,1)];

      Par_sincos = Par(1:4);
      Par_sincos(end+1) = mean(Par(5:end));
      
      
      %Par';
      MeanConst = mean(Par(5:end));
      fprintf('Mean Y-value of the model: %f\n', MeanConst);
      
      Y = getModel(X, FreqPhase, Par_sincos, Har, Deg);
%      figure; plot(X,Y,'g-'); %
      Extram=find_local_extramum(X',Y);
      [Row,Col] = size(Extram);
      if (Row == 2),
         fprintf(Fid, 'Double the period to fit a bi-model lightcurve!\n');
         fprintf('Double the period to fit a bi-model lightcurve!\n');
         BestFreq = BestFreq./2;
         F = folding([Time-T0,Mag,MagErr],1./BestFreq);
      end;

      Period = 24./BestFreq;
      if (length(dFreq) >1),
         [MinVal, MinInx] = min(FreqSolutions(:, 2));
         dFreq = dFreq(MinInx);
      end;
      dPeriod = abs(-24.*dFreq./BestFreq.^2);
      [Period_, dPeriod_] = writePrintableErr(Period, dPeriod);
      
      % Correct the Mag of each night+field+chip by the fit result
      UniqueChip = unique(ChipVec);
      fprintf('Number of night-field-chip sets: %d\n', length(UniqueChip));
      for Ichip=1:1:length(UniqueChip),
         Mag(find(ChipVec == UniqueChip(Ichip))) = Mag(find(ChipVec == UniqueChip(Ichip))) - (Par(4+Ichip) - mean(Par(5:end)));
      end;

   % Old and not relevant condition for a result within the checked range
%   if (FrangeMin + dFreq < BestFreq & BestFreq < FrangeMax - dFreq & dFreq < BestFreq),

      % Condition for a result within the checked range.
      % In the condition use BestFreqOrig (before moving to the bi-model frequency),
      % so it will be consistent with the chi2vec.
      if(BestFreqOrig < max(Frange)-FrangeInterval.*10 & BestFreqOrig > min(Frange)+FrangeInterval.*10 & dFreq < BestFreqOrig.*0.5),
         fprintf(Fid, 'Minimal Chi2 is %f\n', min(Chi2vec));
         fprintf(Fid, 'Best frequency is %f +- %f\n', BestFreq, dFreq);
         fprintf(Fid, 'Period is %f +- %f\n', Period_, dPeriod_);
         fprintf('Minimal Chi2 is %f\n', min(Chi2vec));
         fprintf('Best frequency is %f +- %f\n', BestFreq, dFreq);
         fprintf('Period is %f +- %f\n', Period_, dPeriod_);
      else
         if (strcmp(SuperArcs{Isarc}.Spin.Status, 'periodis failed')),
            SuperArcs{Isarc}.Spin.Status = 'failed spin analysis';
            continue;
         end;
         SuperArcs{Isarc}.Spin.Status = 'fitharmo failed';
         fprintf(Fid, 'There is a problem with the frequency error: is the frequency value worthless ?\n');
         fprintf(Fid, 'Best frequency is %f +- %f\n', BestFreq, dFreq);
         fprintf(Fid, 'Period is %f +- %f\n', Period_, dPeriod_);
         fprintf('There is a problem with the frequency error: is the frequency value worthless ?\n');
         fprintf('Best frequency is %f +- %f\n', BestFreq, dFreq);
         fprintf('Period is %f +- %f\n', Period_, dPeriod_);
      end;
   end;
   
   if (strcmp(SuperArcs{Isarc}.Spin.Status, 'many solutions from fitharmo') | strcmp(SuperArcs{Isarc}.Spin.Status, 'fitharmo failed')),
      if (24.*(max(Time)-min(Time)) < Period_.*2),
         fprintf(Fid, 'The time range is too short to give a solid spin value for %s. Abort\n', legendVec{Isarc});
         fprintf('The time range is too short to give a solid spin value for %s. Abort\n', legendVec{Isarc});
         continue;
      end;
   end;
%end; %%%% REMOVE THIS LINE !!! DP 20110927
   if (FigureNum),

      % Plot the calibrated lightcurve before folding
      plotSuperAstByName(SuperArcs, legendVec{Isarc});
      %   PlotPosition = get(gcf,'position');
      %   set(gcf,'position',[1 550 560 420]);
      set(gcf,'position',[1 800 1700 200]);
      legend off;

      if (~strcmp(SuperArcs{Isarc}.Spin.Status, 'periodis failed')),
         % Plot the power spectrum
         figure;
         plot (24./FreqPS(:,1),FreqPS(:,2),'k'); hold on; grid on;
         graphTitle = ['Power Spectrum of ',legendVec{Isarc}, '. Best Period = ',num2str(24./PeaksPS(end,1)),' hours']; title (graphTitle); hold on;
         ylabelText = ['power']; ylabel(ylabelText); hold on;
         xlabelStr = ['period [hours]']; xlabel(xlabelStr); hold off;
         set(gcf,'position',[1100 310 560 180]);

         % Plot the folded lightcurve and the model
         figure;
         plotFoldedCurveLegend(legendVec{Isarc}, Time, Mag, MagErr, PeaksPS(end,1), 0, Har, Deg);
         set(gcf,'position',[1100 528 560 420]);
      end;

      if (~strcmp(SuperArcs{Isarc}.Spin.Status, 'fitharmo failed')),
         % Plot Chi2 result
         figure;
         if (length(Frange) > 1),
            plot(Frange, Chi2vec, 'b-'); hold on;
            [InxGoodChi]=find(Chi2vec(:,:)<min(Chi2vec)+DeltaChi2);
            H = plot(Frange(InxGoodChi), Chi2vec(InxGoodChi), 'g.'); hold on;
            set(H, 'LineWidth', 3);
            ylabel('chi square values');
            xlabel('Frequencies f [rotation per day]');
            graphTitle = [legendVec{Isarc}, ': Chi square fittness vs. frequency. Best period = ', num2str(Period_)];
            title (graphTitle); hold on;
         elseif (length(Frange) == 1),
            plot(BestFreq, min(Chi2vec), 'ko');
            dFreq = 0;
            ylabel('min Chi square');
            xlabel('Frequencies f [rotation per day]');
            graphTitle = [legendVec{Isarc}, ': Chi square fittness vs. frequency. Best period = ', num2str(Period_)];
            title (graphTitle); hold on;
         end;
         set(gcf,'position',[1100 50 560 180]);
         saveas(H, strcat('PTF_SpinChi2_', AstNameStrForFiles), 'fig');
      
         if (~isempty(FreqSolutions)),
            if (NumSolutions <= 10),
               for Ifreq=1:1:NumSolutions,
                  foldCurvesByFreqs(legendVec{Isarc}, SuperArcs, [FreqSolutions(Ifreq, 1)]);
               end
            else
               fprintf ('--->>> Too many possible solutions - do not print them!\n');
            end
         end

         % Plot Folded lightcurve by the frequency with legend
         figure;
         [H] = plotFoldedCurveLegend(legendVec{Isarc}, Time, Mag, MagErr, BestFreq, dFreq, Har, Deg);
         saveas(H, strcat('PTF_LC_', AstNameStrForFiles), 'fig');
         
      end;
   end;


   if (~strcmp(SuperArcs{Isarc}.Spin.Status, 'fitharmo failed')),
      SuperArcs{Isarc}.Spin.Status = 'determine manually spin analysis';
      Period = 24./BestFreq;
      dPeriod = abs(-24.*dFreq./BestFreq.^2);
      if (dFreq == 0),
         Period_ = Period;
      else
         [Period_, dPeriod_] = writePrintableErr(Period, dPeriod);
      end;

      SuperArcs{Isarc}.Spin.Period = Period_;
      SuperArcs{Isarc}.Spin.PeriodErr = dPeriod_;
      SuperArcs{Isarc}.Spin.Par = Par;
      
      checkSuperArcsPath(SuperArcs, Isarc, Isarc);
      set(gcf,'position',[850 310 250 180]);

      fprintf(Fid, 'Please determine manually spin analysis for asteroid %s\n', legendVec{Isarc});
      fprintf('Please determine manually spin analysis for asteroid %s\n', legendVec{Isarc});

      if (1),
         InputQuestion = ['Asteroid ', legendVec{Isarc}, '-> Determine U (3, 2, 1, 11=Limit, 22=no, Enter=do not decide, 99=leave now):\n'];
         Reply = input(InputQuestion);
         if (Reply == 3 | Reply == 2 | Reply == 1),
            SuperArcs{Isarc}.Spin.Status = ['U=',num2str(Reply)];
         elseif (Reply == 11)
            Spin = [];
            Spin.Status = 'Limit';
            Spin.AmpMin = SuperArcs{Isarc}.Spin.AmpMin;
            InputQuestion = ['Asteroid ', legendVec{Isarc}, '-> Determine Period limit (Enter=do not decide)\n'];
            Reply = input(InputQuestion);
            if (Reply),
               Spin.Period = Reply;
            end;
            SuperArcs{Isarc}.Spin = Spin;
         elseif (Reply == 22)
            Spin = [];
            Spin.Status = 'Failed';
            Spin.AmpMin = SuperArcs{Isarc}.Spin.AmpMin;
            SuperArcs{Isarc}.Spin = Spin;
         elseif (Reply == 99)
            Spin = [];
            Spin.Status = 'Run again';
            Spin.AmpMin = SuperArcs{Isarc}.Spin.AmpMin;
            SuperArcs{Isarc}.Spin = Spin;
            return;
         end;

         if (Reply ~= 22),
            InputQuestion = ['Asteroid ', legendVec{Isarc}, '-> Determine AmpManu (Enter=do not decide)\n'];
            Reply = input(InputQuestion);
            if (Reply),
               SuperArcs{Isarc}.Spin.AmpManu = Reply;
            end;

            InputQuestion = ['Asteroid ', legendVec{Isarc}, '-> Determine AmpManuErr (Enter=do not decide)\n'];
            Reply = input(InputQuestion);
            if (Reply),
               SuperArcs{Isarc}.Spin.AmpManuErr = Reply;
            end;
         end;
         
         % Close graphs automatically
         for I=1:1:50,
            H=gcf; close(H);
         end;
         
      end;

   end;

end; %%% RETURN THIS LINE!!! DP 20110927

% % A patch! Insert the result to the SuperArcs - should be done some place else - DP
% SuperArcs=SuperArcsOrig;
% SuperArcs{ArcLocation}.FoldedLC.F = F;
% SuperArcs{ArcLocation}.FoldedLC.Period = Period;
% SuperArcs{ArcLocation}.FreqPS = FreqPS;
% SuperArcs{ArcLocation}.PeaksPS = PeaksPS;
% ArcLocation
