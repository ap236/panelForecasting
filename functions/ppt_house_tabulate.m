function ppt_house_tabulate (dataPath, doBayes)

% Analysing housing market forecasts for panel forecasting project

% Andreas Pick - 2025

% choices that don't change any more:
noFore = 128; % number of forecasts per MSA

sumBest2One = 0;  % ensure that best forecasts methods sum to one
% (if several same, they get proportional score)
wselindc = 1; % min no of adjacent neighbours for MSA
h = 1;        % forecast horizon - set up to be direct forecasts for given
rollwin = 60;
N = 377;

if doBayes == 1      % include equal weighted forecasts
  nFreq = 9;            % no. frequentist forecasts
  nBayes = 4;           % no. Bayesian forecasts
  noMethods = nFreq + nBayes;
  reOrder = [1 2 4 10:13 5 6 9 7 8]; %
  loadMethods = 1:13;
else
  nFreq = 9;            % no. frequentist forecasts
  nBayes = 1;           % no. Bayesian forecasts
  noMethods = nFreq + nBayes;
  reOrder = [1 2 4 10 5 6 9]; % sortMethods
  loadMethods = 1:10;
end
benchmark = 3;          % individual forecasts as benchmark


% ------ Load Data ---------------------

savePath = [dataPath 'results/house/'];
load([savePath 'house_' num2str(doBayes)]);
msfe_ave_benchmark = nan(1,3);
msfe_ave = nan(noMethods,3);
%msfe_per_unit_ratio = nan(noMethods,noModels);
msfe_ave_ratio = nan(noMethods,3);
beat_benchmark = nan(noMethods,3);
fraction_worst = nan(noMethods,3);
fraction_best = nan(noMethods,3);

dmstat = nan(noMethods-1,1);
noObs = nan(noMethods,3,N);


ctr = 0;
for meanWselInd = 0:2  % 0 - all forecasts, 1 - w*theta close to mean,
  % 2 - w*theta 1 std away
  ctr = ctr+1;

  if meanWselInd == 1
    WselText = 'mean';
  elseif meanWselInd == 2
    WselText = 'std';
  elseif meanWselInd == 0
    WselText = 'all';
  else
    error('meanWselInd incorrectly set')
  end

  % select forecasts that are close to the mean or are one std away from the
  % mean. Where close is c times the standard deviation away.
  c = 0.1;
  N = length(wiTp1{1,3});

  yerr_meanW = nan(noFore,N,noMethods);
  yerr_stdW = nan(noFore,N,noMethods);
  in_meanW = zeros(noFore,N);
  in_stdWu = zeros(noFore,N);
  in_stdWl = zeros(noFore,N);
  for tt = 1:noFore
    for ii = 1:N % over 1:N
      meanW = mean(wthetaiT{1,tt}(:,ii)); % mean of estimation sample
      stdW = std(wthetaiT{1,tt}(:,ii)); % std of estimation sample
      % close to mean
      if abs((wthetaiTp1{1,tt}(ii) - meanW)/(stdW*c)) < 1
        yerr_meanW(tt,ii,:) = yerr(tt,ii,loadMethods);
        in_meanW(tt,ii) = 1;
      end
      % one std away from mean
      if (abs((wthetaiTp1{1,tt}(ii) - meanW + stdW)/(c*stdW)) < 1)
        yerr_stdW(tt,ii,:) = yerr(tt,ii,loadMethods);
        in_stdWu(tt,ii) = 1;
      end
      if (abs((wthetaiTp1{1,tt}(ii) - meanW - stdW)/(c*stdW)) < 1)
        yerr_stdW(tt,ii,:) = yerr(tt,ii,loadMethods);
        in_stdWl(tt,ii) = 1;
      end
    end
  end
  yerr_all = yerr;

  if meanWselInd == 1
    yerrThis = yerr_meanW;
  elseif meanWselInd == 2
    yerrThis = yerr_stdW;
  else
    yerrThis = yerr_all;
  end

  % ------ Output ---------------------

  yerri = yerrThis(:,:,loadMethods);
  nForecasts = size(yerrThis,1);
  msfe_per_unit = nan(Nthis,noMethods);
  %msfe_per_period = nan(nForecasts,noMethods);
  for cMethod = [benchmark 1 2 4:noMethods]
    msfe_per_unit(Wsel >= wselindc,cMethod) = nanmean(yerri(:,Wsel >= wselindc,cMethod).^2,1);
    noObs(cMethod,ctr,Wsel >= wselindc) = sum(~isnan(yerrThis(:,Wsel >= wselindc,cMethod)),1);
  end
  msfe_ave(:,ctr) = nanmean(msfe_per_unit,1);
  msfe_ave_benchmark(ctr) = nanmean(msfe_per_unit(:,benchmark));

  % ratio of average MSFEs
  msfe_ave_ratio(:,ctr) = nanmean(msfe_per_unit,1)'./(nanmean(msfe_per_unit(:,3),1)*ones(noMethods,1));

  % --- DM test statistic -------------------------------------------------
  if ctr == 1
    h = 1;
    indDMstat = nan(Nthis,noMethods);
    indDMsumStats = nan(3,noMethods);
    dctr = 0;
    for cMethod = reOrder
      if cMethod == benchmark
        continue
      end
      dctr=dctr+1;
      % panel DM test of Pesaran and Pick (2013)
      zmat = nan(size(yerrThis,1),Nthis);
      zmat(:,:) = yerrThis(:,Wsel >= wselindc,cMethod).^2 - ...
        yerrThis(:,Wsel >= wselindc,benchmark).^2;

      zmat(:,(sum(~isnan(zmat))<=3)) = [];                % need more than 3 obs per individual

      z = vec_ap(zmat);
      dmstat(dctr,ctr) = PanelDieboldMarianoPPP_matlab (z, h, size(zmat,2));
      % DM test per unit
      for ind = 1:size(zmat,2)
        indDMstat(ind,dctr) = modifiedDieboldMariano(zmat(:,ind),1);
      end
      indDMsumStats(:,dctr) = [ ...
        sum(indDMstat(:,dctr) < -1.96) ...
        sum((indDMstat(:,dctr) > -1.96).*(indDMstat(:,dctr) < 1.96)) ...
        sum(indDMstat(:,dctr) > 1.96)];
    end

  end

  % --- how often beat benchmark/how often best/how often worst -----------
  for n = 1:noMethods
    beat_benchmark(n,ctr) = nanmean(msfe_per_unit(:,n) < msfe_per_unit(:,benchmark));
  end
  fraction_worst(:,ctr) = nanmean(msfe_per_unit == max(msfe_per_unit,[],2))';
  fraction_best(:,ctr) = nanmean(msfe_per_unit == min(msfe_per_unit,[],2))';


  msfeiout = nan(noMethods,Nthis);
  msfeiout(:,:) = msfe_per_unit';
  save('-ascii',[dataPath 'results/house/msfei_' WselText],'msfeiout');

end

% printing to screen ===================================================
absMSFEmult = 1;

matstr1 = '%1.3f & %1.3f & %1.3f &  --  &  --  &  --  & %1.3f & %1.3f & %1.3f & %1.3f & %1.3f & %1.3f';
matstr = '%1.3f';
for j = 1:11
  matstr = [matstr ' & %1.3f'];
end

if doBayes == 0
  disp('--- Upper panel of Table 2: House price inflation forecasts ---')
else
  disp('--- Upper panel of Table S.6: House price inflation forecasts ---')
end

printStr = sprintf(matstr1,[msfe_ave(3,1)*absMSFEmult ...
  msfe_ave(3,2)*absMSFEmult msfe_ave(3,3)*absMSFEmult ...
  fraction_best(benchmark,:) fraction_worst(benchmark,:)]);
disp(printStr)


for jj = reOrder
  printStr = sprintf(matstr,[msfe_ave_ratio(jj,:) ...
    beat_benchmark(jj,:) ...
    fraction_best(jj,:) fraction_worst(jj,:)]);
  disp(printStr)
end

if doBayes == 0
  matstr2 = ' Panel DM & %1.2f & %1.2f & %1.2f & %1.2f & %1.2f & %1.2f & %1.2f\n';
  matstr3 = 'DM<-1.96/DM>1.96 & %d/%d & %d/%d & %d/%d & %d/%d & %d/%d & %d/%d & %d/%d';

  disp('--- Upper pannel of Table 3: --- ')
  printStr = sprintf(matstr2,dmstat(1:7));
  disp(printStr)
  printStr = sprintf(matstr3,indDMsumStats([1 3],1:7));
  disp(printStr)
else
  disp('--- Upper pannel of Table S.7: --- ')
  matstr2 = ' Panel DM & %1.2f & %1.2f & %1.2f & %1.2f & %1.2f\n';
  matstr3 = 'DM<-1.96/DM>1.96 & %d/%d & %d/%d & %d/%d & %d/%d & %d/%d';
  printStr = sprintf(matstr2,dmstat([5:7 11:12]));
  disp(printStr)
  printStr = sprintf(matstr3,indDMsumStats([1 3],[5:7 11:12]));
  disp(printStr)
  
end