function ppt_cpi_tabulate (dataPath, doBayes)

% Tabulates results for panel forecasting project for CPI application 
% for Pesaran, Pick, Timmermann (2025, QE) 

% Andreas Pick - 2025

nFore = 599;
N = 187;

% ======= LOAD and TRANSFORM DATA
load([dataPath 'results/cpi/cpi_forecasts_doBayes_' num2str(doBayes)],'-mat');% all forecasts
yerr = yerr/100;
yhat = yhat/100;
noMethods = size(yerr,3) - (doBayes==0)*2;
if doBayes == 0
  subsetMethods = [1:6 9 10];
else
  subsetMethods = 1:noMethods;
end


% select forecasts that are close to the mean or are one std away from the
% mean. Where close is c times the standard deviation away.
c = 0.1;

msfe_ave       = cell(3,1);
msfe_ave_ratio = cell(3,1);
beat_benchmark = cell(3,1);
fraction_worst = cell(3,1);
fraction_best  = cell(3,1);
noObs          = cell(3,1);

mctr = 0;
for meanWselInd = 0:2  % 0 - all, 1 - close to mean, 2 - ones stdev away
  mctr = mctr+1;

  if meanWselInd == 1
    WselText = 'mean';
  elseif meanWselInd == 2
    WselText = 'std';
  elseif meanWselInd == 0
    WselText = 'all';
  else
    error('meanWselInd incorrectly set')
  end

  yerr_meanW = nan(N,nFore,noMethods);
  yerr_stdW = nan(N,nFore,noMethods);
  in_meanW = nan(nFore,N);
  in_stdWu = nan(nFore,N);
  in_stdWl = nan(nFore,N);
  for tt = 1:nFore
    for ii = 1:N
      meanW = nanmean(wthetaiT{tt,1}(:,ii)); % mean of estimation sample
      if isnan(meanW)
        continue
      end
      stdW = std(wthetaiT{tt,1}(~isnan(wthetaiT{tt,1}(:,ii)),ii)); % std of estimation sample
      % close to mean
      if abs((wthetaiTp1{tt,1}(ii) - meanW)/(stdW*c)) < 1
        yerr_meanW(ii,tt,:) = yerr(ii,tt,subsetMethods);
        in_meanW(tt,ii) = 1;
      else
        in_meanW(tt,ii) = 0;
      end
      % one std away from mean
      if (abs((wthetaiTp1{tt,1}(ii) - meanW + stdW)/(c*stdW)) < 1)
        yerr_stdW(ii,tt,:) = yerr(ii,tt,subsetMethods);
        in_stdWu(tt,ii) = 1;
      else
        in_stdWu(tt,ii) = 0;
      end
      if (abs((wthetaiTp1{tt,1}(ii) - meanW - stdW)/(c*stdW)) < 1)
        yerr_stdW(ii,tt,:) = yerr(ii,tt,subsetMethods);
        in_stdWl(tt,ii) = 1;
      else
        in_stdWl(tt,ii) = 0;
      end
    end
  end

  if meanWselInd == 1
    yerrThis = yerr_meanW;
  elseif meanWselInd == 2
    yerrThis = yerr_stdW;
  else
    yerrThis = yerr(:,:,subsetMethods);
  end

  maxNoModels = size(yerrThis,3);
  benchmark = 3;
  if doBayes == 0
    % subsetMethods = 1:8;
    % noSubsetMethods = length(subsetMethods);
    reOrder = [1 2 4 8 5:7];
    reOrder2 = [3 reOrder];
  else
    % subsetMethods = 1:noMethods;
    reOrder = [1 2 4 10 11 12 13 5 6 9 7 8];
    reOrder2 = [3 reOrder];
  end

  % allocate memory
  msfe_ave_ratio{mctr} = nan(noMethods,1);
  beat_benchmark{mctr} = nan(noMethods,1);
  fraction_worst{mctr} = nan(noMethods,1);
  fraction_best{mctr}  = nan(noMethods,1);

  % --- MSFE---------------------------------------------------------------
  % Calculate the MSFE per unit
  msfe_per_unit = nan(N,noMethods);
  msfe_per_unit(:,:) = nanmean(yerrThis.^2,2);
  msfe_ave{mctr} = nanmean(msfe_per_unit(:,benchmark)); % average msfe of benchmark model
  for iMethod = 1:noMethods
    noObs{mctr}(iMethod,:) = sum(~isnan(yerrThis(:,:,iMethod)),2);
  end

  msfe_ave_ratio{mctr}(:) = nanmean(msfe_per_unit,1)'./(nanmean(msfe_per_unit(:,benchmark),1)*ones(noMethods,1));

  % how often beat benchmark/how often best/how often worst
  for n = 1:noMethods
    beat_benchmark{mctr}(n) = mean(msfe_per_unit(:,n) < msfe_per_unit(:,benchmark));
  end
  fraction_worst{mctr}(:) = mean(msfe_per_unit == max(msfe_per_unit,[],2))';
  fraction_best{mctr}(:) = mean(msfe_per_unit == min(msfe_per_unit,[],2))';

  % --- DM test statistic -------------------------------------------------
  if meanWselInd == 0
    dmstat    = nan(noMethods,1);
    indDMstat = nan(N,noMethods);
    indDMsumStats = nan(3,noMethods);

    h = 1;
    for cMethod = 1:noMethods
      if cMethod == benchmark
        continue
      end
      zmat = nan(N,nFore);
      zmat(:,:) = yerrThis(:,:,cMethod).^2 - yerrThis(:,:,benchmark).^2;
      z = vec_ap(zmat');
      dmstat(cMethod) = PanelDieboldMarianoPPP_matlab (z, h, N);

      % DM test per unit
      for ind = 1:N
        indDMstat(ind,cMethod) = modifiedDieboldMariano(zmat(ind,:)',1);
      end
      indDMsumStats(:,cMethod) = [ ...
        sum(indDMstat(:,cMethod) < -1.96) ...
        sum((indDMstat(:,cMethod) > -1.96).*(indDMstat(:,cMethod) < 1.96)) ...
        sum(indDMstat(:,cMethod) > 1.96)];
    end
  end
  % --- Saving msfes for graphs -------------------------------------------
  if doBayes == 0
    msfeiout = nan(noMethods,N);
    msfeiout(:,:) = msfe_per_unit'; % ratios are calculated in the plots function
    noObsi = nan(noMethods,N);
    noObsi(:,:) = noObs{mctr};
    save('-ascii',[dataPath 'results/cpi/msfei_60_' WselText '.txt'],'msfeiout');
    save('-ascii',[dataPath 'results/cpi/noObs_60_' WselText '.txt'],'noObsi');
  end
end

% === Printing stuff =======================================================
absMSFEmult = 10^5;

matstr1 = '%1.3f & %1.3f & %1.3f &  --  &  --  &  --  & %1.3f & %1.3f & %1.3f & %1.3f & %1.3f & %1.3f';
matstr = '%1.3f';
for j = 1:11
  matstr = [matstr ' & %1.3f'];
end

if doBayes == 0
  disp('--- Lower panel of Table 2: CPI inflation forecasts ---')
else
  disp('--- Lower panel of Table S.6: CPI inflation forecasts ---')
end

printStr = sprintf(matstr1,[msfe_ave{1}*absMSFEmult ...
  msfe_ave{2}*absMSFEmult msfe_ave{3}*absMSFEmult ...
  fraction_best{1}(reOrder2(1),:) fraction_best{2}(reOrder2(1),:) ...
  fraction_best{3}(reOrder2(1),:) fraction_worst{1}(reOrder2(1),:) ...
  fraction_worst{2}(reOrder2(1),:) fraction_worst{3}(reOrder2(1),:)]);
disp(printStr)


for jj = 2:noMethods
  printStr = sprintf(matstr,[msfe_ave_ratio{1}(reOrder(jj-1),:) ...
    msfe_ave_ratio{2}(reOrder(jj-1),:) msfe_ave_ratio{3}(reOrder(jj-1),:) ...
    beat_benchmark{1}(reOrder(jj-1),:) beat_benchmark{2}(reOrder(jj-1),:) ...
    beat_benchmark{3}(reOrder(jj-1),:) ...
    fraction_best{1}(reOrder2(jj),:) fraction_best{2}(reOrder2(jj),:) ...
    fraction_best{3}(reOrder2(jj),:) fraction_worst{1}(reOrder2(jj),:) ...
    fraction_worst{2}(reOrder2(jj),:) fraction_worst{3}(reOrder2(jj),:)]);
  disp(printStr)
end

if doBayes == 0
  matstr2 = ' Panel DM & %1.2f & %1.2f & %1.2f & %1.2f & %1.2f & %1.2f & %1.2f\n';
  matstr3 = 'DM<-1.96/DM>1.96 & %d/%d & %d/%d & %d/%d & %d/%d & %d/%d & %d/%d & %d/%d\n';

  disp('--- Lower pannel of Table 3: --- ')
  printStr = sprintf(matstr2,dmstat(reOrder,:));
  disp(printStr)
  printStr = sprintf(matstr3,indDMsumStats([1 3],reOrder));
  disp(printStr)
else
  disp('--- Lower pannel of Table S.7: --- ')
  matstr2 = ' Panel DM & %1.2f & %1.2f & %1.2f & %1.2f & %1.2f\n';
  matstr3 = 'DM<-1.96/DM>1.96 & %d/%d & %d/%d & %d/%d & %d/%d & %d/%d\n';
  printStr = sprintf(matstr2,dmstat([11:13 7 8],:));
  disp(printStr)
  printStr = sprintf(matstr3,indDMsumStats([1 3],[11:13 7 8]));
  disp(printStr)
  
end
